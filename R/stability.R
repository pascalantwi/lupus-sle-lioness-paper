#' The design table of pilot-d
#'
#' The 6 cells of `pilot-d/preregistration.md`, section 3, with seeds
#' `21061116 + 1000c`.
#'
#' @return A data frame, one row per cell: `cell`, `kappa_f`, `tau`,
#'   `visits`, `gamma`, `predict` (whether the cell fits the prediction
#'   models), `people`, and `seed`.
#' @export
PilotDDesign <- function() {
  d <- data.frame(cell = 1:6,
                  kappa_f = c(0, 0.55, 0.9, 0.55, 0.55, 0.55),
                  tau = c(0, 0, 0, 0.25, 0.25, 0),
                  visits = c(6, 6, 6, 6, 3, 6),
                  gamma = c(0, 0, 0, 0, 0, 0.4))
  d$kappa_c <- 0.45  # measured 0.41 to 0.48 (amendment of 2026-11-16)
  d$predict <- d$cell %in% c(1, 2, 6)
  d$people <- 131L
  d$seed <- 21061116L + 1000L * d$cell
  d
}

#' Simulate repeated visits of the same people
#'
#' Samples are ordered person by person, visits in order. The IFN factor `f`
#' is the person's level (variance `kappa_f`) plus visit noise (variance
#' `1 - kappa_f`), and genes 1 to 20 load 0.8 on it; the cell-mix factor
#' `comp` is built the same way with `kappa_c` and is observed, and genes 31
#' to 50 load 0.6 on it. Every gene carries the person's own level (variance
#' 0.5, the same at every visit) and unit noise. The noise of genes 21 to 30
#' has correlation `0.3 + tau * s + gamma * (pnorm(activity) - 0.5)`, with
#' `s` uniform on (-1, 1) fixed for the person; activity is the person's
#' level (variance 0.39) plus visit noise (variance 0.61).
#'
#' @param people,visits Number of people and visits per person.
#' @param kappa_f,kappa_c Within-person intraclass correlations of the IFN
#'   and cell-mix factors (pilot-d uses 0.45 for the cell mix, amendment of
#'   2026-11-16).
#' @param tau Spread of the person-specific module correlation.
#' @param gamma Effect of activity on the module correlation.
#' @param p Number of genes (at least 50).
#' @return A list with `x`, `person`, `visit`, `f`, `comp`, `activity`, and
#'   `woman` (the person, the name the model functions use).
#' @export
SimulateRepeatedCohort <- function(people, visits, kappa_f, kappa_c = 0.45,
                                   tau = 0, gamma = 0, p = 50) {
  stopifnot(p >= 50, kappa_f >= 0, kappa_f < 1, kappa_c >= 0, kappa_c < 1)
  n <- people * visits
  person <- rep(seq_len(people), each = visits)
  visit <- rep(seq_len(visits), times = people)
  level <- function(k) {
    stats::rnorm(people, sd = sqrt(k))[person] +
      stats::rnorm(n, sd = sqrt(1 - k))
  }
  f <- level(kappa_f)
  comp <- level(kappa_c)
  act <- stats::rnorm(people, sd = sqrt(0.39))[person] +
    stats::rnorm(n, sd = sqrt(0.61))
  s <- stats::runif(people, -1, 1)[person]
  rho <- 0.3 + tau * s + gamma * (stats::pnorm(act) - 0.5)
  stopifnot(all(rho > 0), all(rho < 1))
  offset <- matrix(stats::rnorm(people * p, sd = sqrt(0.5)), people, p)[person, ]
  eps <- matrix(stats::rnorm(n * p), n, p)
  m <- stats::rnorm(n)
  mod <- 21:30
  eps[, mod] <- sqrt(rho) * m + sqrt(1 - rho) * eps[, mod]
  load_f <- c(rep(0.8, 20), rep(0, p - 20))
  load_c <- c(rep(0, 30), rep(0.6, 20), rep(0, p - 50))
  x <- outer(f, load_f) + outer(comp, load_c) + offset + eps
  colnames(x) <- paste0("g", seq_len(p))
  list(x = x, person = person, visit = visit, f = f, comp = comp,
       activity = act, woman = person)
}

#' Pair groups of pilot-d over genes 6 to 50
#'
#' @return A list of column indices into the 990 pairs of [PairIndex()] of
#'   45 (genes 6 to 50): `ifn` (both genes in 6 to 20), `module` (both in 21
#'   to 30), `cellmix` (both in 31 to 50), and `cross` (one in 6 to 20, one
#'   in 21 to 30).
#' @export
PairGroupsD <- function() {
  pr <- PairIndex(45) + 5L
  a <- pr[, 1]
  b <- pr[, 2]
  list(ifn = which(a >= 6 & b <= 20),
       module = which(a >= 21 & b <= 30),
       cellmix = which(a >= 31),
       cross = which(a >= 6 & a <= 20 & b >= 21 & b <= 30))
}

#' One-way intraclass correlation of every column
#'
#' Between-person mean square against within-person mean square, with the
#' usual average group size `n0` for unequal visits:
#' `(MSB - MSW) / (MSB + (n0 - 1) MSW)`.
#'
#' @param e A matrix (samples by columns) or a vector.
#' @param person The person of each sample.
#' @return A vector with one intraclass correlation per column.
#' @export
EdgeICC <- function(e, person) {
  e <- as.matrix(e)
  k <- as.vector(table(person))
  w <- length(k)
  n <- nrow(e)
  mu <- colMeans(e)
  means <- rowsum(e, person) / k
  ssb <- colSums(k * sweep(means, 2, mu)^2)
  ssw <- colSums(e^2) - colSums(k * means^2)
  msb <- ssb / (w - 1)
  msw <- ssw / (n - w)
  n0 <- (n - sum(k^2) / n) / (w - 1)
  (msb - msw) / (msb + (n0 - 1) * msw)
}

#' Fingerprinting accuracy
#'
#' For each person, the Pearson correlation of her edge vector at visit
#' `v2` with every person's vector at visit `v1`; the share of people whose
#' highest correlation is with themselves.
#'
#' @param e Edges, samples by pairs.
#' @param person,visit Person and visit of each sample.
#' @param cols Columns (pairs) to use.
#' @param v1,v2 The two visits compared.
#' @return The identification accuracy.
#' @export
FingerprintAccuracy <- function(e, person, visit, cols, v1 = 1, v2 = 2) {
  a <- e[visit == v1, cols, drop = FALSE][order(person[visit == v1]), ,
                                           drop = FALSE]
  b <- e[visit == v2, cols, drop = FALSE][order(person[visit == v2]), ,
                                           drop = FALSE]
  r <- stats::cor(t(b), t(a))
  mean(apply(r, 1, which.max) == seq_len(nrow(r)))
}

#' Run one replicate of a pilot-d cell
#'
#' Simulates a cohort, computes raw, first-order, and partial-correlation
#' edges on all samples, their mean intraclass correlation and fingerprinting
#' accuracy by pair group, and, in prediction cells, the concordance of the
#' next visit's activity for the models of [PredictModelsB()].
#'
#' @param design_row One row of [PilotDDesign()].
#' @param rep Replicate index; the seed is the cell seed plus `rep`.
#' @return A one-row data frame.
#' @export
StabilityReplicate <- function(design_row, rep) {
  set.seed(design_row$seed + rep)
  coh <- SimulateRepeatedCohort(design_row$people, design_row$visits,
                                design_row$kappa_f, design_row$kappa_c,
                                tau = design_row$tau,
                                gamma = design_row$gamma)
  g <- coh$x[, 6:50]
  score <- rowMeans(StandardizeN(coh$x[, 1:5])$z)
  edges <- list(raw = LionessEdges(g), fo = LionessEdgesFirstOrder(g),
                partial = LionessPartialEdges(g, cbind(score, coh$comp)))
  pg <- PairGroupsD()
  out <- list(cell = design_row$cell, rep = rep)
  for (t in names(edges)) {
    icc <- EdgeICC(edges[[t]], coh$person)
    for (grp in names(pg)) {
      out[[sprintf("icc_%s_%s", t, grp)]] <- mean(icc[pg[[grp]]])
      out[[sprintf("fp_%s_%s", t, grp)]] <-
        FingerprintAccuracy(edges[[t]], coh$person, coh$visit, pg[[grp]])
    }
    out[[sprintf("fp_%s_all", t)]] <-
      FingerprintAccuracy(edges[[t]], coh$person, coh$visit,
                          seq_len(ncol(edges[[t]])))
  }
  if (isTRUE(design_row$predict)) {
    keep <- coh$visit < design_row$visits
    nxt <- coh$activity[which(keep) + 1L]
    sub <- list(x = coh$x[keep, ], comp = coh$comp[keep],
                woman = coh$person[keep])
    pred <- PredictModelsB(sub, nxt, "gaussian")
    for (m in colnames(pred)) out[[m]] <- Concordance(pred[, m], nxt)
  }
  as.data.frame(out, check.names = FALSE)
}

#' Summarize pilot-d replicates by cell
#'
#' Means and MCSEs of every column, and of the excesses and differences the
#' pass rule uses: `x_<type>_<group>` (the group's intraclass correlation
#' minus that of cross pairs), `eq_<group>` (raw minus first-order
#' intraclass correlation), `fpgap_ifn` (raw minus partial fingerprinting on
#' IFN pairs), `keep_module` (partial minus 0.75 times raw module excess),
#' and `d1_N1_vs_E1`, `d2_N1_vs_F1`.
#'
#' @param reps Data frame of replicates from [StabilityReplicate()].
#' @return A long data frame: `cell`, `quantity`, `mean`, `mcse`, `n_rep`.
#' @export
SummarizeCellsD <- function(reps) {
  col <- function(nm) {
    if (nm %in% names(reps)) reps[[nm]] else rep(NA_real_, nrow(reps))
  }
  base <- setdiff(names(reps), c("cell", "rep", "seconds_cell"))
  q <- stats::setNames(lapply(base, col), base)
  for (t in c("raw", "fo", "partial")) {
    for (grp in c("ifn", "module", "cellmix")) {
      q[[sprintf("x_%s_%s", t, grp)]] <- col(sprintf("icc_%s_%s", t, grp)) -
        col(sprintf("icc_%s_cross", t))
    }
  }
  for (grp in c("ifn", "module", "cellmix", "cross")) {
    q[[paste0("eq_", grp)]] <- col(sprintf("icc_raw_%s", grp)) -
      col(sprintf("icc_fo_%s", grp))
  }
  q$fpgap_ifn <- col("fp_raw_ifn") - col("fp_partial_ifn")
  q$keep_module <- q$x_partial_module - 0.75 * q$x_raw_module
  q$d1_N1_vs_E1 <- col("N1") - col("E1")
  q$d2_N1_vs_F1 <- col("N1") - col("F1")
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

#' Known-answer checks of pilot-d
#'
#' The checks of [KnownAnswerChecksB()], the intraclass correlation against
#' an analysis of variance and against `kappa^2` for a squared factor,
#' fingerprinting on identical and on independent visits, and, with the
#' run's summaries, the null cell.
#'
#' @param summaries Optional output of [SummarizeCellsD()].
#' @param seed Seed for the simulated inputs.
#' @return A data frame as [KnownAnswerChecks()].
#' @export
KnownAnswerChecksD <- function(summaries = NULL, seed = 1L) {
  k <- KnownAnswerChecksB(seed = seed)
  set.seed(seed + 2L)
  add <- function(check, setting, known, computed, mcse = NA_real_,
                  type = "identity") {
    k <<- rbind(k, data.frame(check = check, setting = setting, known = known,
                              computed = computed, mcse = mcse, type = type,
                              stringsAsFactors = FALSE))
  }
  person <- rep(1:30, each = 4)
  v <- stats::rnorm(30)[person] + stats::rnorm(120)
  a <- stats::anova(stats::lm(v ~ factor(person)))
  lit <- (a[1, 3] - a[2, 3]) / (a[1, 3] + 3 * a[2, 3])
  add("Intraclass correlation equals the analysis-of-variance formula",
      "30 people, 4 visits; absolute difference", 0,
      abs(EdgeICC(v, person) - lit))
  kap <- 0.6
  pp <- rep(1:20000, each = 2)
  f <- stats::rnorm(20000, sd = sqrt(kap))[pp] +
    stats::rnorm(40000, sd = sqrt(1 - kap))
  add("Intraclass correlation of a squared factor equals kappa squared",
      "kappa 0.6, 20,000 people, 2 visits; absolute difference", 0,
      abs(EdgeICC(f^2, pp) - kap^2))
  e <- matrix(stats::rnorm(50 * 40), 50, 40)
  pe <- rep(1:50, 2)
  ve <- rep(1:2, each = 50)
  add("Fingerprinting of identical visits", "50 people, 40 edges", 1,
      FingerprintAccuracy(rbind(e, e), pe, ve, 1:40))
  ind <- matrix(stats::rnorm(400 * 100), 400, 100)
  add("Fingerprinting of independent visits is near chance",
      "200 people, 100 edges; accuracy minus 1/200", 0,
      abs(FingerprintAccuracy(ind, rep(1:200, 2), rep(1:2, each = 200),
                              1:100) - 1 / 200))
  if (!is.null(summaries)) {
    for (qq in c("x_raw_ifn", "x_partial_ifn", "x_raw_module",
                 "d1_N1_vs_E1")) {
      s <- summaries[summaries$cell == 1 & summaries$quantity == qq, ]
      if (!nrow(s)) next
      add(sprintf("%s in the null cell", qq),
          sprintf("cell 1, %d replicates", s$n_rep), 0, s$mean, s$mcse,
          "mcse")
    }
  }
  k
}
