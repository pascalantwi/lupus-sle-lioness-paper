# pilot-d/verdict.R
# Applies the pass rule of pilot-d/preregistration.md, section 2, clause by
# clause, and the combination rule, to pilot-d/results/cell_summaries.rds.
# Writes pilot-d/results/verdict.rds, stability.png, and fingerprint.png.
#   Rscript pilot-d/verdict.R

devtools::load_all(here::here(), quiet = TRUE)
set.seed(21061116)  # no random step
res_dir <- here::here("pilot-d", "results")
s <- readRDS(file.path(res_dir, "cell_summaries.rds"))
design <- readRDS(file.path(res_dir, "design.rds"))
get <- function(cell, q) s[s$cell == cell & s$quantity == q, ]

# "ge": at least the threshold; "abs": |v| within max(tol, 3 MCSE)
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

g <- get(3, "x_raw_ifn")
o1a <- Judge(g$mean, g$mcse, 0.15, "ge")
add(1, "(a)", 3, "X raw (IFN)", g, ">= 0.15", o1a)

o1b <- character(0)
for (cl in design$cell) for (grp in c("ifn", "module", "cellmix", "cross")) {
  g <- get(cl, paste0("eq_", grp))
  o <- Judge(g$mean, g$mcse, type = "abs", tol = 0.02)
  o1b <- c(o1b, o)
  add(1, "(b)", cl, sprintf("I raw - I F1 (%s)", grp), g,
      sprintf("|D| <= %.4f", max(0.02, 3 * g$mcse)), o)
}

g <- get(2, "x_raw_ifn")
o2a <- Judge(g$mean, g$mcse, 0.05, "ge")
add(2, "(a)", 2, "X raw (IFN)", g, ">= 0.05", o2a)

g <- get(2, "fpgap_ifn")
o2b <- Judge(g$mean, g$mcse, 0.05, "ge")
add(2, "(b)", 2, "P raw (IFN) - P partial (IFN)", g, ">= 0.05", o2b)

o3a <- character(0)
for (cl in c(2, 3)) {
  g <- get(cl, "x_partial_ifn")
  o <- Judge(g$mean, g$mcse, type = "abs", tol = 0.03)
  o3a <- c(o3a, o)
  add(3, "(a)", cl, "X partial (IFN)", g,
      sprintf("|X| <= %.4f", max(0.03, 3 * g$mcse)), o)
}

g <- get(4, "keep_module")
o3b <- Judge(g$mean, g$mcse, 0, "ge")
add(3, "(b)", 4, "X partial (module) - 0.75 X raw (module)", g, ">= 0", o3b)

o3c <- character(0)
g <- get(6, "d1_N1_vs_E1")
o <- Judge(g$mean, g$mcse, 0.03, "ge")
o3c <- c(o3c, o)
add(3, "(c)", 6, "C(N1) - C(E1)", g, ">= 0.03", o)
for (cl in c(1, 2)) {
  g <- get(cl, "d1_N1_vs_E1")
  o <- Judge(g$mean, g$mcse, type = "abs", tol = 0.01)
  o3c <- c(o3c, o)
  add(3, "(c)", cl, "C(N1) - C(E1)", g,
      sprintf("|D| <= %.4f", max(0.01, 3 * g$mcse)), o)
}

clauses <- do.call(rbind, rows)
q <- c(q1a = o1a, q1b = Combine(o1b), q2a = o2a, q2b = o2b,
       q3a = Combine(o3a), q3b = o3b, q3c = Combine(o3c))
status <- if (any(q[c("q1b", "q3a", "q3b")] == "fail")) {
  "piloted-fail"
} else if (all(q == "pass")) {
  "piloted-pass"
} else {
  "piloted-close-call"
}
verdict <- list(clauses = clauses, status = status, reasons = q,
                date = Sys.time(),
                rule = "pilot-d/preregistration.md, sections 2, 5, and 6")
saveRDS(verdict, file.path(res_dir, "verdict.rds"))
print(clauses, digits = 3)
cat("status:", status, "\n")

lab <- sprintf("%d: kappaF %.2f, tau %.2f, %d visits%s", design$cell,
               design$kappa_f, design$tau, design$visits,
               ifelse(design$gamma > 0, ", activity", ""))
types <- c(raw = "raw LIONESS", fo = "first-order", partial = "partial (fix)")
cols <- c(raw = "#B2182B", fo = "#762A83", partial = "#2166AC")
Plot <- function(file, qfun, xlab, groups) {
  png(file.path(res_dir, file), width = 2200, height = 1300, res = 220)
  par(mfrow = c(1, length(groups)), mar = c(4, if (length(groups) > 1) 11 else 13, 2, 1))
  for (grp in groups) {
    yy <- rev(seq_len(nrow(design)))
    m <- sapply(names(types), function(t) sapply(design$cell, function(cl) {
      v <- get(cl, qfun(t, grp)); if (nrow(v)) v$mean else NA }))
    plot(NA, xlim = range(c(0, m), na.rm = TRUE) + c(-0.02, 0.02),
         ylim = c(0.5, nrow(design) + 0.5), yaxt = "n", ylab = "", xlab = xlab,
         main = grp)
    axis(2, at = yy, labels = lab, las = 1, cex.axis = 0.65)
    abline(v = 0, col = "grey50", lty = 2)
    off <- c(0.2, 0, -0.2)
    for (k in seq_along(types)) points(m[, k], yy + off[k], pch = 19,
                                       col = cols[k], cex = 0.8)
    if (grp == groups[1]) legend("bottomright", types, col = cols, pch = 19,
                                 cex = 0.65, bg = "white")
  }
  dev.off()
}
Plot("stability.png", function(t, grp) sprintf("x_%s_%s", t, grp),
     "Excess stability over cross pairs (ICC)", c("ifn", "module", "cellmix"))
Plot("fingerprint.png", function(t, grp) sprintf("fp_%s_%s", t, grp),
     "Fingerprinting accuracy", c("ifn", "module", "all"))
