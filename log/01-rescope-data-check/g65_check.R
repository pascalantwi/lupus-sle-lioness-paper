# Data check of GSE65391 (pediatric SLE, longitudinal), outside the project. Base R.
args <- commandArgs(TRUE); D <- args[1]; O <- args[2]
con <- gzfile(file.path(D, "GSE65391_series_matrix.txt.gz")); hl <- readLines(con); close(con)
b <- grep("^!series_matrix_table_begin", hl); e <- grep("^!series_matrix_table_end", hl)
head <- hl[1:b]
ch <- grep("^!Sample_characteristics_ch1", head, value = TRUE)
geo <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", head, value = TRUE), "\t")[[1]][-1])
s <- data.frame(geo = geo)
for (row in ch) {
  v <- gsub('^"|"$', "", strsplit(row, "\t")[[1]][-1])
  k <- make.names(sub(": .*", "", v[1])); s[[k]] <- sub("^[^:]*: ", "", v)
}
num <- function(x) suppressWarnings(as.numeric(x))
for (k in c("sledai","neutrophil_percent","lymphocyte_percent","neutrophil_count","lymphocyte_count","c3","c4","ds_dna","cumulative_time","visit"))
  if (k %in% names(s)) s[[k]] <- num(s[[k]])
tab <- read.delim(text = hl[(b + 1):(e - 1)], check.names = FALSE)
x <- as.matrix(tab[, -1]); rownames(x) <- tab[[1]]
cat("expression range:", round(range(x, na.rm = TRUE), 2), "\n")
sym <- read.delim(file.path(D, "gpl10558_symbol.tsv"), header = FALSE)
g <- sym$V2[match(rownames(x), sym$V1)]
genescore <- function(gs) {
  rows <- which(g %in% gs)
  z <- t(scale(t(x[rows, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
s$ifn <- genescore(c("EPSTI1","IFI44L","LY6E","OAS3","RSAD2"))
sle <- s$disease.state == "SLE"
p <- s[sle, ]
sink(file.path(O, "g65_check.txt"))
cat("SLE samples:", nrow(p), " children:", length(unique(p$subject)), " healthy samples:", sum(!sle), "\n")
cat("\nSLEDAI distribution (SLE visits):\n"); print(summary(p$sledai))
cat("visits with SLEDAI 0:", sum(p$sledai == 0, na.rm = TRUE), "\n")
cat("\nIFN score, healthy vs SLE: mean, SD\n")
print(round(rbind(healthy = c(mean(s$ifn[!sle]), sd(s$ifn[!sle])), SLE = c(mean(p$ifn), sd(p$ifn))), 3))
v <- p[, c("sledai","ifn","neutrophil_percent","lymphocyte_percent","ds_dna","c3","c4")]
cat("\nSpearman correlations, all SLE visits (between and within children):\n")
print(round(cor(v, use = "pairwise", method = "spearman"), 2))
cw <- v; for (k in names(cw)) cw[[k]] <- cw[[k]] - ave(cw[[k]], p$subject, FUN = function(z) mean(z, na.rm = TRUE))
cat("\nPearson correlations WITHIN children (each centered on the child's own mean):\n")
print(round(cor(cw, use = "pairwise"), 2))
cat("\nActivity domain (sledai_component_class):\n"); print(table(p$sledai_component_class, useNA = "ifany"))
cat("\nRenal SLEDAI item ('renal' field):\n"); print(table(p$renal, useNA = "ifany"))
cat("\nNephritis class:\n"); print(table(p$nephritis_class, useNA = "ifany"))
grp <- function(lab) {
  f <- factor(lab)
  out <- sapply(c("ifn","neutrophil_percent","lymphocyte_percent"), function(k)
    tapply(p[[k]], f, function(z) sprintf("%.2f (SD %.2f)", mean(z, na.rm = TRUE), sd(z, na.rm = TRUE))))
  print(cbind(n = table(f), out))
}
cat("\nFactors by activity domain:\n"); grp(p$sledai_component_class)
cat("\nFactors by renal item:\n"); grp(p$renal)
p$grade <- cut(p$sledai, c(-Inf, 4, 10, Inf), labels = c("low (0-4)", "moderate (5-10)", "high (>10)"))
cat("\nFactors by SLEDAI grade:\n"); grp(p$grade)
# flare: SLEDAI rise of 4 or more at the next visit of the same child
p <- p[order(p$subject, p$cumulative_time), ]
nxt <- ave(p$sledai, p$subject, FUN = function(z) c(z[-1], NA))
p$flare_next <- (nxt - p$sledai) >= 4
cat("\nVisits followed by a flare (next-visit SLEDAI rise >= 4):", sum(p$flare_next, na.rm = TRUE),
    "of", sum(!is.na(p$flare_next)), "; children with at least one:",
    length(unique(p$subject[p$flare_next %in% TRUE])), "\n")
cat("\nFactors at visits before a flare vs not:\n"); grp(ifelse(is.na(p$flare_next), NA, ifelse(p$flare_next, "before flare", "no flare")))
sink()
saveRDS(p, file.path(O, "g65_scores.rds"))
