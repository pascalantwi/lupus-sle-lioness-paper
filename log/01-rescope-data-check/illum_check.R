# Data check of ILLUMINATE baseline (GSE88884), outside the project. Base R.
args <- commandArgs(TRUE); D <- args[1]; H <- args[2]; O <- args[3]
hl <- readLines(H, warn = FALSE)
fld <- function(tag) lapply(grep(paste0("^!Sample_", tag), hl, value = TRUE),
  function(x) gsub('^"|"$', "", strsplit(x, "\t")[[1]][-1]))
s <- data.frame(geo = fld("geo_accession")[[1]],
                cel = sub(".*_(A[0-9]+)\\.CEL\\.gz$", "\\1.CEL", fld("supplementary_file")[[1]]))
for (row in fld("characteristics_ch1")) {
  k <- sub(": .*", "", row[1]); k <- make.names(sub(" \\(.*", "", k))
  s[[k]] <- sub("^[^:]*: ", "", row)
}
num <- function(x) suppressWarnings(as.numeric(x))
s$sledai <- num(s$sledai_at_baseline); s$dsdna <- num(s$antidsdna_at_baseline)
s$c3 <- num(s$c3_at_baseline); s$c4 <- num(s$c4_at_baseline)
s$sle <- grepl("SLE", s$group)
m <- read.delim(file.path(D, "markers.tsv"), check.names = FALSE, row.names = 1)
sym <- read.delim(file.path(D, "gpl17586_symbol.tsv"), header = FALSE)
rownames(m) <- sym$V2[match(sub(":1$", "", rownames(m)), sym$V1)]
m <- as.matrix(m[, s$cel])
z <- t(scale(t(m)))
sc <- function(g) colMeans(z[g, , drop = FALSE])
s$ifn <- sc(c("EPSTI1","IFI44L","LY6E","OAS3","RSAD2"))
s$bcell <- sc(c("MS4A1","CD19","CD79A","CD79B","CD22","BANK1"))
s$neut <- sc(c("FCGR3B","CXCR2","CSF3R","FPR1"))
s$lymph <- sc(c("CD3E","CD3D","CD2","MS4A1","CD79A"))
s$mix <- s$neut - s$lymph
hc <- !s$sle
thr <- mean(s$ifn[hc]) + 2 * sd(s$ifn[hc])
s$ifn_pos <- s$ifn > thr
s$grade <- cut(s$sledai, c(-Inf, 6, 11, Inf), labels = c("mild", "moderate", "severe"))
sink(file.path(O, "illum_check.txt"))
cat("Samples:", nrow(s), " SLE:", sum(s$sle), " healthy:", sum(hc), "\n")
cat("Healthy by batch:\n"); print(table(s$batch[hc]))
cat("\nIFN-positive (above healthy mean + 2 SD), SLE:", round(mean(s$ifn_pos[s$sle]), 3), "\n")
cat("By trial:\n"); print(round(tapply(s$ifn_pos[s$sle], s$batch[s$sle], mean), 3))
summ <- function(v, g) do.call(rbind, lapply(split(v, g), function(x) c(n = length(x), mean = mean(x), sd = sd(x))))
for (f in c("ifn", "bcell", "mix")) {
  cat("\n==", f, "by group (healthy vs SLE grade): mean, SD\n")
  grp <- ifelse(hc, "healthy", as.character(s$grade))
  print(round(summ(s[[f]], factor(grp, c("healthy", "mild", "moderate", "severe"))), 3))
}
cat("\nSpearman correlations among SLE patients:\n")
v <- s[s$sle, c("sledai", "dsdna", "c3", "c4", "ifn", "bcell", "mix")]
print(round(cor(v, use = "pairwise", method = "spearman"), 2))
cat("\nMismatch counts (SLE):\n")
x <- s[s$sle & !is.na(s$grade), ]
print(table(grade = x$grade, ifn_positive = x$ifn_pos, trial = x$batch))
cat("\nSerology profile x IFN-positive:\n")
x$sero <- ifelse(x$adna_at_baseline == "Positive", "dsDNA+", "dsDNA-")
x$lowc <- ifelse(x$c3base == "L" | x$c4base == "L", "lowC", "normC")
print(table(paste(x$sero, x$lowc), x$ifn_pos))
cat("\nIFN score by ancestry (SLE): mean, SD\n"); print(round(summ(x$ifn, x$race), 3))
cat("\nIFN score by sex (SLE): mean, SD\n"); print(round(summ(x$ifn, x$Sex), 3))
sink()
saveRDS(s, file.path(O, "illum_scores.rds"))
