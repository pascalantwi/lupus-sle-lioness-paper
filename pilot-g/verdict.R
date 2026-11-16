# pilot-g/verdict.R
# Applies the pass rule of pilot-g/preregistration.md, section 2, to
# pilot-g/results/cell_summaries.rds. Writes pilot-g/results/verdict.rds and
# pilot-g/results/pathway_stability.png.
#   Rscript pilot-g/verdict.R

devtools::load_all(here::here(), quiet = TRUE)
set.seed(21361116)  # no random step
res_dir <- here::here("pilot-g", "results")
s <- readRDS(file.path(res_dir, "cell_summaries.rds"))
design <- readRDS(file.path(res_dir, "design.rds"))
get <- function(cell, q) s[s$cell == cell & s$quantity == q, ]
Judge <- function(value, mcse, threshold = NA, type = c("ge", "abs"),
                  tol = 0.05) {
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

o1a <- o1b <- o1c <- character(0)
for (cl in design$cell) {
  g <- get(cl, "min_r2_within")
  o <- Judge(g$mean, g$mcse, 0.99, "ge")
  o1a <- c(o1a, o)
  add(1, "(a)", cl, "smallest R^2, connectivity on its terms", g, ">= 0.99", o)
  g <- get(cl, "min_r2_cross")
  o <- Judge(g$mean, g$mcse, 0.99, "ge")
  o1b <- c(o1b, o)
  add(1, "(b)", cl, "smallest R^2, crosstalk on its terms", g, ">= 0.99", o)
  g <- get(cl, "cor_ifn")
  o <- Judge(g$mean, g$mcse, 0.9, "ge")
  o1c <- c(o1c, o)
  add(1, "(c)", cl, "cor(IFN connectivity, IFN score^2)", g, ">= 0.9", o)
}
g <- get(2, "expl_ifn")
o2a <- Judge(g$mean, g$mcse, type = "abs")
add(2, "(a)", 2, "ICC IFN connectivity - ICC IFN score^2", g,
    sprintf("|D| <= %.4f", max(0.05, 3 * g$mcse)), o2a)
g <- get(2, "agg_ifn")
o2b <- Judge(g$mean, g$mcse, 0.05, "ge")
add(2, "(b)", 2, "ICC IFN connectivity - mean ICC IFN edges", g, ">= 0.05",
    o2b)
pathways <- names(PathwaySetsG())
o3a <- character(0)
for (cl in design$cell) for (pw in pathways) {
  g <- get(cl, paste0("icc_resid_", pw))
  o <- Judge(g$mean, g$mcse, type = "abs")
  o3a <- c(o3a, o)
  add(3, "(a)", cl, sprintf("ICC residual connectivity, %s", pw), g,
      sprintf("|v| <= %.4f", max(0.05, 3 * g$mcse)), o)
}
g <- get(2, "personal_vs_null")
o3b <- Judge(g$mean, g$mcse, 0.05, "ge")
add(3, "(b)", 2, "ICC personal connectivity - ICC null connectivity", g,
    ">= 0.05", o3b)

clauses <- do.call(rbind, rows)
q <- c(q1a = Combine(o1a), q1b = Combine(o1b), q1c = Combine(o1c),
       q2a = o2a, q2b = o2b, q3a = Combine(o3a), q3b = o3b)
status <- if (any(q[c("q1a", "q1b", "q3a")] == "fail")) {
  "piloted-fail"
} else if (all(q == "pass")) {
  "piloted-pass"
} else {
  "piloted-close-call"
}
verdict <- list(clauses = clauses, status = status, reasons = q,
                date = Sys.time(),
                rule = "pilot-g/preregistration.md, sections 2, 5, and 6")
saveRDS(verdict, file.path(res_dir, "verdict.rds"))
print(clauses, digits = 3)
cat("status:", status, "\n")

# Figure: per pathway, the stability of connectivity, squared score,
# residual, and single edges, in cells 2 and 3
labs <- c(ifn = "P1 interferon", neut = "P2 neutrophil-like",
          lymph = "P3 lymphocyte-like", pb = "P4 plasmablast-like",
          personal = "P5 personal", null = "P6 null")
kinds <- c(conn = "connectivity", score2 = "squared score",
           edges = "single edges (mean)", resid = "connectivity residual")
cols <- c(conn = "#B2182B", score2 = "grey40", edges = "#762A83",
          resid = "#2166AC")
png(file.path(res_dir, "pathway_stability.png"), width = 2200, height = 1200,
    res = 220)
par(mfrow = c(1, 2), mar = c(4, 9, 2, 1))
for (cl in c(2, 3)) {
  m <- sapply(names(kinds), function(k) sapply(names(labs), function(pw) {
    get(cl, sprintf("icc_%s_%s", k, pw))$mean
  }))
  yy <- rev(seq_along(labs))
  plot(NA, xlim = range(c(0, m)) + c(-0.02, 0.02),
       ylim = c(0.5, length(labs) + 0.5), yaxt = "n", ylab = "",
       xlab = "Stability over visits (ICC)",
       main = sprintf("Cell %d: IFN stability %.2f", cl,
                      design$kappa_ifn[design$cell == cl]))
  axis(2, at = yy, labels = labs, las = 1, cex.axis = 0.7)
  abline(v = 0, col = "grey80")
  off <- c(0.24, 0.08, -0.08, -0.24)
  for (k in seq_along(kinds)) {
    points(m[, k], yy + off[k], pch = 19, col = cols[k], cex = 0.8)
  }
  if (cl == 2) {
    legend("bottomright", kinds, col = cols, pch = 19, cex = 0.6,
           bg = "white")
  }
}
invisible(dev.off())
