s <- readRDS(commandArgs(TRUE)[1])
x <- s[s$sle, ]
pos <- grepl("Positive$", x$adna_at_baseline)
lowc <- x$c3base == "L" | x$c4base == "L"
prof <- ifelse(pos & lowc, "both", ifelse(pos, "dsDNA only", ifelse(lowc, "low complement only", "neither")))
print(table(profile = prof, ifn_positive = x$ifn_pos))
cat("\nIFN mean, SD by profile\n")
print(round(do.call(rbind, lapply(split(x$ifn, prof), function(v) c(n = length(v), mean = mean(v), sd = sd(v)))), 3))
cat("\nAdditional healthy vs trial healthy, IFN and mix means:\n")
h <- s[!s$sle, ]
print(round(do.call(rbind, lapply(split(h[, c("ifn", "bcell", "mix")], h$batch == "Additional HCs"), colMeans)), 3))
