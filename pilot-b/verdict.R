# pilot-b/verdict.R
# Applies the pass rule of pilot-b/preregistration.md, section 2, clause by
# clause, and the combination rule, to pilot-b/results/cell_summaries.rds.
# Writes pilot-b/results/verdict.rds, shares.png, and model_gains.png.
#   Rscript pilot-b/verdict.R

devtools::load_all(here::here(), quiet = TRUE)
set.seed(20861113)  # no random step
res_dir <- here::here("pilot-b", "results")
s <- readRDS(file.path(res_dir, "cell_summaries.rds"))
design <- readRDS(file.path(res_dir, "design.rds"))
get <- function(cell, q) s[s$cell == cell & s$quantity == q, ]

# "ge": at least the threshold; "abs": |v| within max(0.01, 3 MCSE);
# "le05": at most 0.05 + max(0.01, 3 MCSE)
Judge <- function(value, mcse, threshold = NA, type = c("ge", "abs", "le05")) {
  type <- match.arg(type)
  if (type == "ge") {
    if (abs(value - threshold) <= 3 * mcse) return("close")
    return(if (value >= threshold) "pass" else "fail")
  }
  if (type == "abs") {
    if (abs(value) <= max(0.01, 3 * mcse)) return("pass")
    if (abs(value) <= 0.01 + 3 * mcse) return("close")
    return("fail")
  }
  if (value <= 0.05 + max(0.01, 3 * mcse)) return("pass")
  if (value <= 0.06 + 3 * mcse) return("close")
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

# 1(a): raw edges call false rewiring at the stress level, cell 4
g <- get(4, "paired_raw_factor")
o1a <- Judge(g$mean, g$mcse, 0.10, "ge")
add(1, "(a)", 4, "raw share, factor pairs (paired)", g, ">= 0.10", o1a)

# 1(b): partial edges at the nominal level, cells 2, 3, 4 (paired), 8, 10
# (between diseases)
o1b <- character(0)
for (cl in c(2, 3, 4, 8, 10)) {
  q <- if (cl %in% c(8, 10)) "between_par_factor" else "paired_par_factor"
  g <- get(cl, q)
  o <- Judge(g$mean, g$mcse, type = "le05")
  o1b <- c(o1b, o)
  add(1, "(b)", cl, sprintf("partial share, factor pairs (%s)",
                            if (cl %in% c(8, 10)) "between" else "paired"),
      g, sprintf("<= %.4f", 0.05 + max(0.01, 3 * g$mcse)), o)
}

# 1(c): partial edges keep a real change, cell 5
g <- get(5, "keep_module")
o1c <- Judge(g$mean, g$mcse, 0, "ge")
add(1, "(c)", 5, "partial minus 0.75 raw share, module pairs (paired)", g,
    ">= 0", o1c)

# 2(a): false rewiring at the study's levels, cell 2 or cell 3 (the larger)
g2 <- get(2, "paired_raw_factor")
g3 <- get(3, "paired_raw_factor")
g <- if (g2$mean >= g3$mean) g2 else g3
o2a <- Judge(g$mean, g$mcse, 0.10, "ge")
add(2, "(a)", g$cell, "raw share, factor pairs (paired), larger of cells 2, 3",
    g, ">= 0.10", o2a)

# 2(b): a false SLE-against-RA difference, cell 8
g <- get(8, "between_raw_factor")
o2b <- Judge(g$mean, g$mcse, 0.10, "ge")
add(2, "(b)", 8, "raw share, factor pairs (between)", g, ">= 0.10", o2b)

# 2(c): a real SLE-only change found by partial edges, cell 9
g <- get(9, "between_par_module")
o2c <- Judge(g$mean, g$mcse, 0.15, "ge")
add(2, "(c)", 9, "partial share, module pairs (between)", g, ">= 0.15", o2c)

# 3(a): edges equal their formula in prediction, cells 1 to 7
o3a <- character(0)
for (cl in 1:7) {
  g <- get(cl, "d2_N1_vs_F1")
  o <- Judge(g$mean, g$mcse, type = "abs")
  o3a <- c(o3a, o)
  add(3, "(a)", cl, sprintf("%s(N1) - %s(F1)",
                            if (cl <= 5) "AUC" else "C", if (cl <= 5) "AUC" else "C"),
      g, sprintf("|D| <= %.4f", max(0.01, 3 * g$mcse)), o)
}

# 3(b): the E1 base removes the activity artifact, cell 6
g <- get(6, "d1_N1_vs_E1")
o3b <- Judge(g$mean, g$mcse, type = "abs")
add(3, "(b)", 6, "C(N1) - C(E1)", g,
    sprintf("|E| <= %.4f", max(0.01, 3 * g$mcse)), o3b)

# 3(c): edges detect activity-linked co-expression, cell 7
g <- get(7, "d1_N1_vs_E1")
o3c <- Judge(g$mean, g$mcse, 0.03, "ge")
add(3, "(c)", 7, "C(N1) - C(E1)", g, ">= 0.03", o3c)

clauses <- do.call(rbind, rows)
q <- c(q1a = o1a, q1b = Combine(o1b), q1c = o1c, q2a = o2a, q2b = o2b,
       q2c = o2c, q3a = Combine(o3a), q3b = o3b, q3c = o3c)
status <- if (any(q[c("q1b", "q3a", "q3b")] == "fail")) {
  "piloted-fail"
} else if (all(q == "pass")) {
  "piloted-pass"
} else {
  "piloted-close-call"
}
verdict <- list(clauses = clauses, status = status, reasons = q,
                date = Sys.time(),
                rule = "pilot-b/preregistration.md, sections 2 and 5")
saveRDS(verdict, file.path(res_dir, "verdict.rds"))
print(clauses, digits = 3)
cat("status:", status, "\n")

# Figure 1: share of pairs called changed, raw against partial, by cell
lab <- sprintf("%d: aim %d, %s%s%s", design$cell, design$aim, design$levels,
               ifelse(design$aim == 3, paste(" vs", design$levels_b), ""),
               ifelse(design$module > 0, ", module +0.4", ""))
lab[design$aim == 2] <- sprintf("%d: aim 2 (no test)", design$cell[design$aim == 2])
sh <- function(cl, kind, grp) {
  q <- sprintf("%s_%s_%s", if (cl %in% 8:10) "between" else "paired", kind, grp)
  g <- get(cl, q)
  if (!nrow(g)) c(NA, NA) else c(g$mean, g$mcse)
}
tcells <- design$cell[design$aim != 2]
png(file.path(res_dir, "shares.png"), width = 2000, height = 1300, res = 220)
par(mfrow = c(1, 2), mar = c(4, 12, 2, 1))
for (grp in c("factor", "module")) {
  m <- t(sapply(tcells, function(cl) c(sh(cl, "raw", grp)[1],
                                       sh(cl, "par", grp)[1])))
  yy <- rev(seq_along(tcells))
  plot(NA, xlim = c(0, max(0.3, max(m, na.rm = TRUE) + 0.03)),
       ylim = c(0.5, length(tcells) + 0.5), yaxt = "n",
       ylab = "", xlab = "Share of pairs called changed",
       main = sprintf("%s pairs", tools::toTitleCase(grp)))
  axis(2, at = yy, labels = lab[match(tcells, design$cell)], las = 1,
       cex.axis = 0.7)
  abline(v = 0.05, col = "grey50", lty = 2)
  points(m[, 1], yy + 0.12, pch = 19, col = "#B2182B")
  points(m[, 2], yy - 0.12, pch = 17, col = "#2166AC")
  if (grp == "factor") legend("bottomright", c("raw LIONESS", "partial (fix)"),
                              pch = c(19, 17), col = c("#B2182B", "#2166AC"),
                              cex = 0.75, bg = "white")
}
dev.off()

# Figure 2: paired model differences, cells 1 to 7
qs <- c(artifact_N1R0_vs_R0 = "N1 with R0 vs R0 (the artifact)",
        d1_N1_vs_E1 = "N1 vs E1", d2_N1_vs_F1 = "N1 vs F1")
cols <- c("#B2182B", "#2166AC", "#762A83")
pcells <- 1:7
m <- sapply(names(qs), function(q) sapply(pcells, function(cl) get(cl, q)$mean))
e <- sapply(names(qs), function(q) sapply(pcells, function(cl) get(cl, q)$mcse))
png(file.path(res_dir, "model_gains.png"), width = 2000, height = 1200,
    res = 220)
par(mar = c(4, 13, 1, 1))
yy <- rev(seq_along(pcells))
plot(NA, xlim = range(c(m - 2 * e, m + 2 * e, 0)),
     ylim = c(0.5, length(pcells) + 0.5), yaxt = "n", ylab = "",
     xlab = "Paired difference in AUC (cells 1-5) or C-index (6-7)")
axis(2, at = yy, labels = lab[pcells], las = 1, cex.axis = 0.7)
abline(v = 0, col = "grey50")
off <- c(0.2, 0, -0.2)
for (k in seq_along(qs)) {
  segments(m[, k] - 2 * e[, k], yy + off[k], m[, k] + 2 * e[, k], yy + off[k],
           col = cols[k])
  points(m[, k], yy + off[k], pch = 19, col = cols[k], cex = 0.7)
}
legend("topright", legend = qs, col = cols, pch = 19, cex = 0.7, bg = "white")
dev.off()
