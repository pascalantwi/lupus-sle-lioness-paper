p <- readRDS(commandArgs(TRUE)[1])
grp <- function(lab) {
  f <- factor(lab)
  out <- sapply(c("ifn","neutrophil_percent","lymphocyte_percent"), function(k)
    tapply(p[[k]], f, function(z) sprintf("%.2f (SD %.2f)", mean(z, na.rm = TRUE), sd(z, na.rm = TRUE))))
  print(cbind(n = table(f), out))
}
cat("Factors by nephritis class:\n"); grp(p$nephritis_class)
for (k in c("oral_steroids_category","steroid_iv_category","mycophenolate_category","hydroxychloroquine_category","cyclophosphamide_category")) {
  cat("\nFactors by", k, ":\n"); grp(p[[k]])
}
