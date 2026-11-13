#' Factor levels of pilot-b
#'
#' The shift, spread ratio, and within-woman share of variance of the IFN
#' factor (`f`) and the composition factor (`comp`) in pregnancy, rounded
#' from `results/factor_shift.rds` (SLE, RA), plus no change and the stress
#' level.
#'
#' @param name `"none"`, `"SLE"`, `"RA"`, or `"stress"`.
#' @return A list with `f` and `comp`, each `c(shift, scale, kappa)`.
#' @export
FactorLevels <- function(name = c("none", "SLE", "RA", "stress")) {
  name <- match.arg(name)
  switch(name,
    none = list(f = c(0, 1, 0.5), comp = c(0, 1, 0.5)),
    SLE = list(f = c(0.2, 1.1, 0.9), comp = c(-0.1, 1.0, 0.55)),
    RA = list(f = c(-0.3, 0.65, 0.4), comp = c(0.45, 1.45, 0.05)),
    stress = list(f = c(1, 2.25, 0.5), comp = c(0, 1, 0.5)))
}

#' The design table of pilot-b
#'
#' The 10 cells of `pilot-b/preregistration.md`, section 3, with seeds
#' `20861113 + 1000c`. For aim 3 cells, `levels` and `module` describe the
#' SLE cohort and `levels_b`, `module_b` the RA cohort of `women_b` women.
#'
#' @return A data frame, one row per cell.
#' @export
PilotBDesign <- function() {
  d <- data.frame(
    cell = 1:10,
    aim = c(1, 1, 1, 1, 1, 2, 2, 3, 3, 3),
    levels = c("none", "SLE", "RA", "stress", "SLE", "RA", "RA", "SLE", "SLE",
               "stress"),
    module = c(0, 0, 0, 0, 0.4, 0, 0, 0, 0.4, 0),
    activity = c(rep("none", 5), "ifn", "module", rep("none", 3)),
    women = c(18, 18, 27, 18, 18, 36, 36, 18, 18, 18),
    levels_b = c(rep(NA, 7), "RA", "RA", "none"),
    module_b = c(rep(NA, 7), 0, 0, 0),
    women_b = c(rep(NA, 7), 27, 27, 27),
    stringsAsFactors = FALSE
  )
  d$seed <- 20861113L + 1000L * d$cell
  d
}

#' Simulate one cohort of pilot-b
#'
#' Four samples per woman, the first two pregnant. The IFN factor `f` (genes
#' 1 to 20, loading 0.8) and the composition factor `comp` (genes 31 to 50,
#' loading 0.6) are each the woman's level (variance `kappa`) plus sample
#' noise (variance `1 - kappa`); in pregnant samples each becomes
#' `shift + scale *` its draw ([FactorLevels()]). Every gene carries the
#' woman's offset (variance 0.5) and unit noise; the noise of genes 21 to 30
#' has correlation 0.3, plus `module` in pregnant samples. Activity is the
#' woman's level plus visit noise (variance 0.5 each); with `activity =
#' "ifn"` the IFN factor is multiplied by `exp(0.5 a)`, and with
#' `"module"` the module's correlation is `0.3 + 0.4 pnorm(a)`.
#'
#' @param women Number of women.
#' @param levels A list from [FactorLevels()].
#' @param module Change in the module's correlation in pregnancy.
#' @param activity `"none"`, `"ifn"`, or `"module"`.
#' @param p Number of genes (at least 50).
#' @return A list with `x`, `y` (1 pregnant), `woman`, `comp`, `f`, and
#'   `activity`.
#' @export
SimulateCohortB <- function(women, levels = FactorLevels("none"), module = 0,
                            activity = c("none", "ifn", "module"), p = 50) {
  activity <- match.arg(activity)
  stopifnot(p >= 50, 0.3 + module < 1)
  n <- 4L * women
  woman <- rep(seq_len(women), each = 4L)
  y <- rep(c(1L, 1L, 0L, 0L), times = women)
  draw <- function(lv) {
    v <- stats::rnorm(women, sd = sqrt(lv[3]))[woman] +
      stats::rnorm(n, sd = sqrt(1 - lv[3]))
    v[y == 1] <- lv[1] + lv[2] * v[y == 1]
    v
  }
  f <- draw(levels$f)
  comp <- draw(levels$comp)
  act <- stats::rnorm(women, sd = sqrt(0.5))[woman] +
    stats::rnorm(n, sd = sqrt(0.5))
  if (activity == "ifn") f <- f * exp(0.5 * act)
  rho <- rep(0.3, n)
  rho[y == 1] <- 0.3 + module
  if (activity == "module") rho <- 0.3 + 0.4 * stats::pnorm(act)
  offset <- matrix(stats::rnorm(women * p, sd = sqrt(0.5)), women, p)[woman, ]
  eps <- matrix(stats::rnorm(n * p), n, p)
  mod <- 21:30
  # a shared module draw per sample gives correlation rho between its genes
  m <- stats::rnorm(n)
  eps[, mod] <- sqrt(rho) * m + sqrt(1 - rho) * eps[, mod]
  load_f <- c(rep(0.8, 20), rep(0, p - 20))
  load_c <- c(rep(0, 30), rep(0.6, 20), rep(0, p - 50))
  x <- outer(f, load_f) + outer(comp, load_c) + offset + eps
  colnames(x) <- paste0("g", seq_len(p))
  list(x = x, y = y, woman = woman, comp = comp, f = f, activity = act)
}

#' Pair groups of the gene set 6 to 50
#'
#' @return A list of column indices into the 990 pairs of genes 6 to 50
#'   ([PairIndex()] of 45): `factor` (pairs within genes 6 to 20 or within 31
#'   to 50, 295 pairs) and `module` (pairs within genes 21 to 30, 45 pairs).
#' @export
PairGroups <- function() {
  pr <- PairIndex(45) + 5L
  within <- function(a, b) pr[, 1] >= a & pr[, 2] <= b & pr[, 1] <= b
  list(factor = which(within(6, 20) | (pr[, 1] >= 31)),
       module = which(within(21, 30) & pr[, 1] >= 21))
}

#' Each woman's change in edge between pregnancy and the rest
#'
#' @param e Edges, samples in rows.
#' @param y 1 pregnant, 0 not.
#' @param woman Woman of each sample.
#' @return A matrix, one row per woman, one column per pair: her mean edge in
#'   pregnancy minus her mean edge outside it.
#' @export
ChangeByWoman <- function(e, y, woman) {
  m1 <- rowsum(e[y == 1, , drop = FALSE], woman[y == 1])
  m1 <- m1 / as.vector(table(woman[y == 1]))
  m0 <- rowsum(e[y == 0, , drop = FALSE], woman[y == 0])
  m0 <- m0 / as.vector(table(woman[y == 0]))
  m1 - m0[rownames(m1), , drop = FALSE]
}

#' Share of pairs called changed by a paired test across women
#'
#' @param d Output of [ChangeByWoman()].
#' @param cols Columns (pairs) to judge.
#' @param level The test's level.
#' @return The share of the pairs with a one-sample t-test p-value below
#'   `level`.
#' @export
PairedShare <- function(d, cols, level = 0.05) {
  d <- d[, cols, drop = FALSE]
  w <- nrow(d)
  t <- colMeans(d) / (apply(d, 2, stats::sd) / sqrt(w))
  mean(2 * stats::pt(-abs(t), w - 1) < level)
}

#' Share of pairs whose change differs between two groups of women
#'
#' @param d1,d2 Outputs of [ChangeByWoman()] for the two diseases.
#' @param cols Columns (pairs) to judge.
#' @param level The test's level.
#' @return The share of the pairs with a Welch two-sample t-test p-value
#'   below `level`.
#' @export
BetweenShare <- function(d1, d2, cols, level = 0.05) {
  a <- d1[, cols, drop = FALSE]
  b <- d2[, cols, drop = FALSE]
  n1 <- nrow(a)
  n2 <- nrow(b)
  v1 <- apply(a, 2, stats::var) / n1
  v2 <- apply(b, 2, stats::var) / n2
  t <- (colMeans(a) - colMeans(b)) / sqrt(v1 + v2)
  df <- (v1 + v2)^2 / (v1^2 / (n1 - 1) + v2^2 / (n2 - 1))
  mean(2 * stats::pt(-abs(t), df) < level)
}

#' Raw and partial-correlation LIONESS edges of one disease's samples
#'
#' @param coh A cohort from [SimulateCohortB()].
#' @return A list with `raw` and `partial` edge matrices over genes 6 to 50;
#'   the partial edges condition on the IFN score (mean of genes 1 to 5
#'   standardized over the samples) and composition.
#' @export
CohortEdges <- function(coh) {
  g <- coh$x[, 6:50]
  score <- rowMeans(StandardizeN(coh$x[, 1:5])$z)
  list(raw = LionessEdges(g),
       partial = LionessPartialEdges(g, cbind(score, coh$comp)))
}

#' Concordance index of predictions for a continuous outcome
#'
#' Among pairs of samples with different outcomes, the share in which the
#' sample with the larger outcome has the larger prediction, ties in
#' prediction counted one half. For a binary outcome it equals the AUC.
#'
#' @param pred Predictions.
#' @param y Outcomes.
#' @return The concordance index.
#' @export
Concordance <- function(pred, y) {
  dy <- outer(y, y, "-")
  dp <- outer(pred, pred, "-")
  use <- dy > 0
  (sum(dp[use] > 0) + 0.5 * sum(dp[use] == 0)) / sum(use)
}

#' Out-of-fold predictions of the six models of pilot-b
#'
#' @param coh A cohort from [SimulateCohortB()].
#' @param outcome The outcome: `coh$y` (pregnancy) or `coh$activity`.
#' @param family `"binomial"` or `"gaussian"`.
#' @return A matrix of out-of-fold predictions, one column per model (`R0`,
#'   `R1`, `E1`, `F1`, `N1`, `N1R0`).
#' @export
PredictModelsB <- function(coh, outcome, family = c("binomial", "gaussian")) {
  family <- match.arg(family)
  models <- c("R0", "R1", "E1", "F1", "N1", "N1R0")
  folds <- GroupFolds(coh$woman, 5L)
  pred <- matrix(NA_real_, length(outcome), length(models),
                 dimnames = list(NULL, models))
  for (f in 1:5) {
    tr <- folds != f
    ft <- PregnancyFeatures(coh, tr)
    a <- ft$train
    b <- ft$test
    ytr <- outcome[tr]
    inner <- GroupFolds(coh$woman[tr], 5L)
    en <- function(rtr, ptr, rte, pte) {
      FitElasticNet(rtr, ptr, ytr, rte, pte, foldid = inner, family = family)
    }
    pred[!tr, "R0"] <- FitReference(a$r0, ytr, b$r0, family = family)
    pred[!tr, "R1"] <- FitReference(a$r1, ytr, b$r1, family = family)
    pred[!tr, "E1"] <- en(a$r1, a$genes, b$r1, b$genes)
    pred[!tr, "F1"] <- en(a$r1, cbind(a$genes, a$fo), b$r1,
                          cbind(b$genes, b$fo))
    pred[!tr, "N1"] <- en(a$r1, cbind(a$genes, a$edges), b$r1,
                          cbind(b$genes, b$edges))
    pred[!tr, "N1R0"] <- en(a$r0, a$edges, b$r0, b$edges)
  }
  pred
}

#' Run one replicate of a pilot-b cell
#'
#' Aim 1: the paired shares (raw and partial edges, factor and module pairs)
#' and the AUC of the six models for pregnancy stage. Aim 2: the concordance,
#' correlation, out-of-fold R squared, and root mean squared error of the
#' six models for activity. Aim 3: the between-disease shares, raw and
#' partial, factor and module pairs, and the SLE cohort's paired shares.
#'
#' @param design_row One row of [PilotBDesign()].
#' @param rep Replicate index; the seed is the cell seed plus `rep`.
#' @return A one-row data frame.
#' @export
RunReplicateB <- function(design_row, rep) {
  set.seed(design_row$seed + rep)
  pg <- PairGroups()
  out <- list(cell = design_row$cell, rep = rep)
  coh <- SimulateCohortB(design_row$women, FactorLevels(design_row$levels),
                         design_row$module, design_row$activity)
  shares <- function(ed, prefix) {
    d_raw <- ChangeByWoman(ed$raw, coh$y, coh$woman)
    d_par <- ChangeByWoman(ed$partial, coh$y, coh$woman)
    out[[paste0(prefix, "raw_factor")]] <<- PairedShare(d_raw, pg$factor)
    out[[paste0(prefix, "raw_module")]] <<- PairedShare(d_raw, pg$module)
    out[[paste0(prefix, "par_factor")]] <<- PairedShare(d_par, pg$factor)
    out[[paste0(prefix, "par_module")]] <<- PairedShare(d_par, pg$module)
    list(raw = d_raw, par = d_par)
  }
  if (design_row$aim == 1) {
    shares(CohortEdges(coh), "paired_")
    pred <- PredictModelsB(coh, coh$y, "binomial")
    for (m in colnames(pred)) out[[m]] <- Auc(pred[, m], coh$y)
  } else if (design_row$aim == 2) {
    yy <- coh$activity
    pred <- PredictModelsB(coh, yy, "gaussian")
    for (m in colnames(pred)) {
      out[[m]] <- Concordance(pred[, m], yy)
      out[[paste0("cor_", m)]] <- stats::cor(pred[, m], yy)
      out[[paste0("r2_", m)]] <- 1 - sum((yy - pred[, m])^2) /
        sum((yy - mean(yy))^2)
      out[[paste0("rmse_", m)]] <- sqrt(mean((yy - pred[, m])^2))
    }
  } else {
    d_a <- shares(CohortEdges(coh), "paired_")
    coh_a <- coh
    coh <- SimulateCohortB(design_row$women_b,
                           FactorLevels(design_row$levels_b),
                           design_row$module_b)
    ed_b <- CohortEdges(coh)
    d_b_raw <- ChangeByWoman(ed_b$raw, coh$y, coh$woman)
    d_b_par <- ChangeByWoman(ed_b$partial, coh$y, coh$woman)
    out$between_raw_factor <- BetweenShare(d_a$raw, d_b_raw, pg$factor)
    out$between_raw_module <- BetweenShare(d_a$raw, d_b_raw, pg$module)
    out$between_par_factor <- BetweenShare(d_a$par, d_b_par, pg$factor)
    out$between_par_module <- BetweenShare(d_a$par, d_b_par, pg$module)
  }
  as.data.frame(out, check.names = FALSE)
}

#' Summarize pilot-b replicates by cell
#'
#' Means and MCSEs of every column of the replicates, and of the paired
#' differences and combinations the pass rule uses.
#'
#' @param reps Data frame of replicates from [RunReplicateB()] (cells may
#'   have different columns; missing ones are `NA`).
#' @return A long data frame: `cell`, `quantity`, `mean`, `mcse`, `n_rep`.
#' @export
SummarizeCellsB <- function(reps) {
  col <- function(nm) {
    if (nm %in% names(reps)) reps[[nm]] else rep(NA_real_, nrow(reps))
  }
  base <- setdiff(names(reps), c("cell", "rep", "seconds_cell"))
  q <- stats::setNames(lapply(base, col), base)
  q$d1_N1_vs_E1 <- col("N1") - col("E1")
  q$d2_N1_vs_F1 <- col("N1") - col("F1")
  q$artifact_N1R0_vs_R0 <- col("N1R0") - col("R0")
  q$gain_R1_vs_R0 <- col("R1") - col("R0")
  q$gain_E1_vs_R1 <- col("E1") - col("R1")
  q$keep_module <- col("paired_par_module") - 0.75 * col("paired_raw_module")
  out <- do.call(rbind, lapply(split(seq_len(nrow(reps)), reps$cell),
                               function(rows) {
    do.call(rbind, lapply(names(q), function(nm) {
      v <- q[[nm]][rows]
      v <- v[!is.na(v)]
      if (!length(v)) return(NULL)
      data.frame(cell = reps$cell[rows[1]], quantity = nm, mean = mean(v),
                 mcse = if (length(v) > 1) stats::sd(v) / sqrt(length(v))
                 else NA_real_,
                 n_rep = length(v))
    }))
  }))
  rownames(out) <- NULL
  out
}

#' Known-answer checks of pilot-b
#'
#' The checks of [KnownAnswerChecks()], the partial-correlation edges against
#' the inverse-covariance formula and against Pearson edges, the concordance
#' index on known cases, and, with the run's summaries, the null cell.
#'
#' @param summaries Optional output of [SummarizeCellsB()].
#' @param seed Seed for the simulated inputs.
#' @return A data frame as [KnownAnswerChecks()].
#' @export
KnownAnswerChecksB <- function(summaries = NULL, seed = 1L) {
  k <- KnownAnswerChecks(seed = seed)
  set.seed(seed + 1L)
  add <- function(check, setting, known, computed, mcse = NA_real_,
                  type = "identity") {
    k <<- rbind(k, data.frame(check = check, setting = setting, known = known,
                              computed = computed, mcse = mcse, type = type,
                              stringsAsFactors = FALSE))
  }
  x <- matrix(stats::rnorm(30 * 4), 30, 4)
  cond <- cbind(stats::rnorm(30), stats::rnorm(30))
  x[, 1:2] <- x[, 1:2] + cond[, 1]
  lit <- apply(PairIndex(4), 1, function(ij) {
    pm <- solve(stats::cov(cbind(x[, ij], cond)))
    -pm[1, 2] / sqrt(pm[1, 1] * pm[2, 2])
  })
  add("Partial correlation equals the inverse-covariance formula",
      "30 samples, 4 genes, 2 conditioning variables; largest difference", 0,
      max(abs(PartialCorrelationPairs(x, cond) - lit)))
  add("Partial-correlation edges with no conditioning equal Pearson edges",
      "30 samples, 4 genes; largest difference", 0,
      max(abs(LionessPartialEdges(x) - LionessEdges(x))))
  add("Concordance of a perfect prediction", "10 distinct outcomes", 1,
      Concordance(1:10, 1:10))
  add("Concordance of a constant prediction", "10 distinct outcomes", 0.5,
      Concordance(rep(0, 10), 1:10))
  yb <- rep(0:1, each = 5)
  sc <- stats::rnorm(10)
  add("Concordance equals the AUC for a binary outcome", "5 against 5",
      0, abs(Concordance(sc, yb) - Auc(sc, yb)))
  if (!is.null(summaries)) {
    for (qq in c("paired_raw_factor", "paired_raw_module",
                 "paired_par_factor", "paired_par_module")) {
      s <- summaries[summaries$cell == 1 & summaries$quantity == qq, ]
      add(sprintf("Share called changed (%s) in the null cell",
                  gsub("_", " ", sub("paired_", "", qq))),
          sprintf("cell 1, %d replicates", s$n_rep), 0.05, s$mean, s$mcse,
          "mcse")
    }
    for (m in c("R0", "R1", "E1", "F1", "N1", "N1R0")) {
      s <- summaries[summaries$cell == 1 & summaries$quantity == m, ]
      add(sprintf("AUC of %s in the null cell", m),
          sprintf("cell 1, %d replicates", s$n_rep), 0.5, s$mean, s$mcse,
          "mcse")
    }
  }
  k
}
