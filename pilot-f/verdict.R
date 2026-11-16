# pilot-f/verdict.R
# Applies the pass rule of pilot-f/preregistration.md, section 2, to
# pilot-f/results/cell_summaries.rds. Writes pilot-f/results/verdict.rds and
# pilot-f/results/metrics_stability.png.
#   Rscript pilot-f/verdict.R

devtools::load_all(here::here(), quiet = TRUE)
set.seed(21261116)  # no random step
res_dir <- here::here("pilot-f", "results")
s <- readRDS(file.path(res_dir, "cell_summaries.rds"))
design <- readRDS(file.path(res_dir, "design.rds"))
get <- function(cell, q) s[s$cell == cell & s$quantity == q, ]
Judge <- function(value, mcse, threshold = NA, type = c("ge", "abs"),
                  tol = 0.02) {
  type <- match.arg(type)
  if (type == "ge") {
    if (abs(value - threshold) <= 3 * mcse) return("close")
    return(if (value >= threshold) "pass" else "fail")
  }
  if (abs(value) <= max(tol, 3 * mcse)) return("pass")
  if (abs(value) <= tol + 3 * mcse) return("close")
  "fail"
}
rows <- list()
add <- function(question, clause, cell, quantity, g, threshold, outcome) {
  rows[[length(rows) + 1]] <<- data.frame(
    question = question, clause = clause, cell = cell, quantity = quantity,
    value = g$mean, mcse = g$mcse, threshold = threshold, outcome = outcome,
    stringsAsFactors = FALSE)
}
Combine <- function(o) {
  if (any(o == "fail")) "fail" else if (any(o == "close")) "close" else "pass"
}

o1a <- character(0)
for (cl in design$cell) {
  g <- get(cl, "closedform_cor")
  o <- Judge(g$mean, g$mcse, 0.9, "ge")
  o1a <- c(o1a, o)
  add(1, "(a)", cl, "cor(IFN connectivity, IFN score^2)", g, ">= 0.9", o)
}
o1b <- character(0)
eqs <- unique(grep("^eq_", s$quantity, value = TRUE))
for (cl in design$cell) for (q in eqs) {
  g <- get(cl, q)
  o <- Judge(g$mean, g$mcse, type = "abs", tol = 0.02)
  o1b <- c(o1b, o)
  add(1, "(b)", cl, sprintf("ICC raw - ICC F1 (%s)", sub("^eq_", "", q)), g,
      sprintf("|D| <= %.4f", max(0.02, 3 * g$mcse)), o)
}
g <- get(2, "agg_ifn")
o2a <- Judge(g$mean, g$mcse, 0.05, "ge")
add(2, "(a)", 2, "ICC IFN connectivity - ICC IFN edges (raw)", g, ">= 0.05", o2a)
g <- get(2, "expl_ifn")
o2b <- Judge(g$mean, g$mcse, type = "abs", tol = 0.05)
add(2, "(b)", 2, "ICC IFN connectivity (raw) - ICC IFN score^2", g,
    sprintf("|D| <= %.4f", max(0.05, 3 * g$mcse)), o2b)
o3a <- character(0)
for (cl in c(2, 3)) {
  g <- get(cl, "fix_ifn_ownset")
  o <- Judge(g$mean, g$mcse, type = "abs", tol = 0.05)
  o3a <- c(o3a, o)
  add(3, "(a)", cl, "ICC own-module IFN connectivity - ICC own-module cross", g,
      sprintf("|D| <= %.4f", max(0.05, 3 * g$mcse)), o)
}
g <- get(4, "fix_module_ownset")
o3b <- Judge(g$mean, g$mcse, 0.05, "ge")
add(3, "(b)", 4, "ICC own-module module connectivity - ICC own-module cross", g,
    ">= 0.05", o3b)

clauses <- do.call(rbind, rows)
q <- c(q1a = Combine(o1a), q1b = Combine(o1b), q2a = o2a, q2b = o2b,
       q3a = Combine(o3a), q3b = o3b)
status <- if (any(q[c("q1a", "q1b", "q3a")] == "fail")) {
  "piloted-fail"
} else if (all(q == "pass")) {
  "piloted-pass"
} else {
  "piloted-close-call"
}
verdict <- list(clauses = clauses, status = status, reasons = q,
                date = Sys.time(),
                rule = "pilot-f/preregistration.md, sections 2, 5, and 6")
saveRDS(verdict, file.path(res_dir, "verdict.rds"))
print(clauses[!(clauses$question == 1 & clauses$clause == "(b)"), ], digits = 3)
cat("1(b): ", sum(o1b == "pass"), "of", length(o1b), "pass; largest |D|",
    format(max(abs(clauses$value[clauses$question == 1 & clauses$clause == "(b)"])), digits = 3), "\n")
cat("status:", status, "\n")

# Figure: stability (ICC) of each metric by edge type, cell 2 and cell 4
mets <- c(conn_ifn = "IFN connectivity", conn_module = "module connectivity",
          conn_cross = "cross connectivity", strength_ifn = "strength, IFN genes",
          betweenness_ifn = "betweenness, IFN genes", eigen_ifn = "eigenvector, IFN genes",
          clustering_ifn = "clustering, IFN genes", pathlen = "path length",
          efficiency = "efficiency", density = "density")
cols <- c(raw = "#B2182B", fo = "#762A83", ownset = "#2166AC")
png(file.path(res_dir, "metrics_stability.png"), width = 2200, height = 1300, res = 220)
par(mfrow = c(1, 2), mar = c(4, 11, 2, 1))
for (cl in c(2, 4)) {
  yy <- rev(seq_along(mets))
  m <- sapply(names(cols), function(t) sapply(names(mets), function(k) {
    v <- get(cl, sprintf("icc_%s_%s", t, k)); if (nrow(v)) v$mean else NA }))
  plot(NA, xlim = range(c(0, m), na.rm = TRUE) + c(-0.02, 0.02),
       ylim = c(0.5, length(mets) + 0.5), yaxt = "n", ylab = "",
       xlab = "Stability over visits (ICC)",
       main = sprintf("Cell %d%s", cl, if (cl == 4) ": real personal module" else ": measured IFN stability"))
  axis(2, at = yy, labels = mets, las = 1, cex.axis = 0.7)
  sc <- get(cl, "icc_score2")$mean
  abline(v = sc, col = "grey40", lty = 2)
  off <- c(0.2, 0, -0.2)
  for (k in seq_along(cols)) points(m[, k], yy + off[k], pch = 19, col = cols[k], cex = 0.8)
  if (cl == 2) legend("bottomright", c("raw LIONESS", "first-order", "own-module (fix)", "IFN score squared"),
                      col = c(cols, "grey40"), pch = c(19, 19, 19, NA), lty = c(NA, NA, NA, 2),
                      cex = 0.6, bg = "white")
}
dev.off()
