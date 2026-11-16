p <- readRDS(commandArgs(TRUE)[1])
icc <- function(v, g) {
  ok <- !is.na(v); v <- v[ok]; g <- factor(g[ok])
  k <- table(g); g <- droplevels(g[g %in% names(k)[k >= 2]]); v <- v[as.character(g) %in% levels(g) | TRUE][seq_along(g)]
  a <- anova(lm(v ~ g)); n0 <- mean(table(g))
  (a[1, 3] - a[2, 3]) / (a[1, 3] + (n0 - 1) * a[2, 3])
}
keep <- ave(rep(1, nrow(p)), p$subject, FUN = length) >= 2
q <- p[keep, ]
cat("Children with 2+ visits:", length(unique(q$subject)), " visits:", nrow(q), "\n")
cat("Visits per child: "); print(summary(as.vector(table(q$subject))))
for (k in c("ifn", "neutrophil_percent", "lymphocyte_percent", "sledai"))
  cat(sprintf("Within-child stability (ICC) of %-20s %.2f\n", k, icc(q[[k]], q$subject)))
q <- q[order(q$subject, q$cumulative_time), ]
d <- ave(q$cumulative_time, q$subject, FUN = function(z) c(diff(z), NA))
cat("\nDays between consecutive visits: "); print(summary(d))
ds <- ave(q$sledai, q$subject, FUN = function(z) c(diff(z), NA))
cat("Change in SLEDAI to next visit: "); print(summary(ds)); cat("SD of change:", round(sd(ds, na.rm = TRUE), 2), "\n")
cat("\nTraining/test split (children):\n"); print(tapply(q$subject, q$set, function(z) length(unique(z))))
