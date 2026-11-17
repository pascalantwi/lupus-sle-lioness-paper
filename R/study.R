# The simulation study of protocol.md: S1 (describe and explain the LIONESS
# metrics over visits) and S2 (prediction of identity and activity).

#' The design table of the simulation study
#'
#' S1 varies one factor at a time around the base setting (the measured
#' GSE65391 values); S2 adds an activity effect to the base setting
#' (protocol.md, section 2).
#'
#' @param beta The activity effects of the S2 settings, the first 0.
#' @return A data frame, one row per setting: `cell`, `part`, `family`,
#'   `kappa_ifn`, `load_ifn`, `size_ifn`, `kappa_neut`, `kappa_lymph`,
#'   `kappa_pb`, `tau`, `visits`, `people`, `beta`, `seed`.
#' @export
StudyDesign <- function(beta = c(0, 0.15, 0.3, 0.45, 0.6)) {
  base <- list(kappa_ifn = 0.55, load_ifn = 0.8, size_ifn = 15L,
               kappa_neut = 0.41, kappa_lymph = 0.48, kappa_pb = 0.39,
               tau = 0.25, visits = 6L, beta = 0)
  rows <- list()
  add <- function(part, family, ...) {
    r <- utils::modifyList(base, list(...))
    rows[[length(rows) + 1]] <<- data.frame(part = part, family = family, r)
  }
  add("S1", "base")
  for (k in c(0, 0.25, 0.75, 0.9)) add("S1", "stability", kappa_ifn = k)
  for (l in c(0.4, 0.6)) add("S1", "loading", load_ifn = l)
  for (m in c(5L, 40L)) add("S1", "size", size_ifn = m)
  for (t in c(0, 0.1)) add("S1", "personal", tau = t)
  add("S1", "visits", visits = 3L)
  add("S1", "null", kappa_ifn = 0, kappa_neut = 0, kappa_lymph = 0,
      kappa_pb = 0, tau = 0)
  for (b in beta) add("S2", "activity", beta = b)
  d <- do.call(rbind, rows)
  d <- cbind(cell = seq_len(nrow(d)), d)
  d$people <- 131L
  d$seed <- 21461117L + 1000L * d$cell
  d
}

#' The six pathways of the study, with a given interferon pathway size
#'
#' @param size_ifn Number of interferon genes; the other pathways have 15.
#' @return A named list of column indices, as [PathwaySetsG()].
#' @export
StudySets <- function(size_ifn = 15L) {
  sizes <- c(ifn = size_ifn, neut = 15L, lymph = 15L, pb = 15L,
             personal = 15L, null = 15L)
  ends <- cumsum(sizes)
  sets <- lapply(seq_along(sizes), function(k) (ends[k] - sizes[k] + 1):ends[k])
  stats::setNames(sets, names(sizes))
}

#' Simulate one cohort of the study
#'
#' As [SimulatePathwayCohort()], with the interferon pathway's loading and
#' size set by the design row, and an activity score for every visit: the
#' activity effect times the mean of the interferon and plasmablast factors
#' (scaled to unit variance), plus a person level and visit noise with
#' within-person stability 0.39 (SLEDAI in GSE65391).
#'
#' @param design_row One row of [StudyDesign()].
#' @return A list with `x`, `person`, `visit`, `activity`, and `sets`.
#' @export
SimulateStudyCohort <- function(design_row) {
  d <- design_row
  people <- d$people
  visits <- d$visits
  n <- people * visits
  person <- rep(seq_len(people), each = visits)
  visit <- rep(seq_len(visits), times = people)
  sets <- StudySets(d$size_ifn)
  p <- max(unlist(sets))
  level <- function(k) {
    stats::rnorm(people, sd = sqrt(k))[person] +
      stats::rnorm(n, sd = sqrt(1 - k))
  }
  x <- matrix(stats::rnorm(people * p, sd = sqrt(0.5)), people, p)[person, ] +
    matrix(stats::rnorm(n * p), n, p)
  f <- list(ifn = level(d$kappa_ifn), neut = level(d$kappa_neut),
            lymph = level(d$kappa_lymph), pb = level(d$kappa_pb))
  load <- c(ifn = d$load_ifn, neut = 0.6, lymph = 0.6, pb = 0.6)
  for (k in names(f)) {
    x[, sets[[k]]] <- x[, sets[[k]]] + load[[k]] * f[[k]]
  }
  rho_p <- 0.3 + d$tau * stats::runif(people, -1, 1)[person]
  for (k in c("personal", "null")) {
    cols <- sets[[k]]
    rho <- if (k == "personal") rho_p else rep(0.3, n)
    # one shared draw per sample gives equicorrelation rho between the genes
    shared <- stats::rnorm(n)
    x[, cols] <- matrix(stats::rnorm(people * length(cols), sd = sqrt(0.5)),
                        people)[person, ] +
      sqrt(rho) * shared +
      sqrt(1 - rho) * matrix(stats::rnorm(n * length(cols)), n)
  }
  colnames(x) <- paste0("g", seq_len(p))
  activity <- d$beta * (f$ifn + f$pb) / sqrt(2) +
    stats::rnorm(people, sd = sqrt(0.39))[person] +
    stats::rnorm(n, sd = sqrt(0.61))
  list(x = x, person = person, visit = visit, activity = activity,
       sets = sets)
}

#' Graph metrics of one sample's thresholded LIONESS network
#'
#' Keeps the pairs whose absolute edge exceeds `threshold`, or, when `prop`
#' is given, the `prop` share of pairs with the largest absolute edge.
#'
#' @param e One sample's edges ([PairIndex()] order).
#' @param pr `PairIndex(p)`.
#' @param p Number of genes.
#' @param threshold Absolute edge above which a pair is kept.
#' @param prop Share of pairs kept instead of a fixed threshold.
#' @return A list: `n_edges`, `n_nodes` (genes with at least one kept
#'   edge), `transitivity` (global; 0 without connected triples), `degree`
#'   and `betweenness` (one per gene; distance the inverse absolute edge).
#' @export
GraphMetrics <- function(e, pr, p, threshold = 1, prop = NULL) {
  keep <- if (is.null(prop)) {
    which(abs(e) > threshold)
  } else {
    order(-abs(e))[seq_len(max(1L, round(prop * length(e))))]
  }
  g <- igraph::make_graph(t(pr[keep, , drop = FALSE]), n = p,
                          directed = FALSE)
  deg <- igraph::degree(g)
  tr <- igraph::transitivity(g, type = "global")
  btw <- if (length(keep)) {
    igraph::betweenness(g, directed = FALSE, weights = 1 / abs(e[keep]))
  } else {
    rep(0, p)
  }
  list(n_edges = length(keep), n_nodes = sum(deg > 0),
       transitivity = if (is.finite(tr)) tr else 0, degree = deg,
       betweenness = btw)
}

#' LIONESS network metrics of every sample
#'
#' Whole-network metrics (total edge weight; with a fixed threshold the
#' number of edges, number of nodes, and transitivity; with a proportional
#' threshold the number of nodes and transitivity) and pathway metrics
#' (connectivity, crosstalk, and the mean strength, degree, and betweenness
#' of each pathway's genes, the last two under both thresholds).
#'
#' @param e Edges, samples by pairs.
#' @param sets A named list of gene index sets.
#' @param threshold Fixed absolute-edge threshold.
#' @param prop Proportional threshold: the share of pairs kept.
#' @return A matrix, samples by metrics.
#' @export
LionessMetrics <- function(e, sets, threshold = 1, prop = 0.05) {
  p <- max(unlist(sets))
  pr <- PairIndex(p)
  blocks <- PathwayBlocks(sets, p)
  out <- PathwayMetrics(e, blocks)
  colnames(out) <- sub("^within_", "conn_", colnames(out))
  inc <- matrix(0, nrow(pr), p)
  inc[cbind(seq_len(nrow(pr)), pr[, 1])] <- 1
  inc[cbind(seq_len(nrow(pr)), pr[, 2])] <- 1
  strength <- e %*% inc
  pw_mean <- function(v, nm) {
    stats::setNames(vapply(sets, function(s) mean(v[s]), numeric(1)),
                    paste0(nm, "_", names(sets)))
  }
  graph <- t(vapply(seq_len(nrow(e)), function(q) {
    a <- GraphMetrics(e[q, ], pr, p, threshold = threshold)
    b <- GraphMetrics(e[q, ], pr, p, prop = prop)
    c(n_edges = a$n_edges, n_nodes = a$n_nodes,
      transitivity = a$transitivity, n_nodes_p = b$n_nodes,
      transitivity_p = b$transitivity,
      pw_mean(a$degree, "degree"), pw_mean(a$betweenness, "betweenness"),
      pw_mean(b$degree, "degree_p"), pw_mean(b$betweenness, "betweenness_p"))
  }, numeric(5 + 4 * length(sets))))
  st <- t(apply(strength, 1, pw_mean, nm = "strength"))
  cbind(total_weight = rowSums(e), out, st, graph)
}

#' Fourth-order expression terms of pathway connectivity and crosstalk
#'
#' The products of the first-order terms of [PathwayTerms()] that the
#' second-order expansion of the edge adds, plus the sum of fourth powers
#' for a pathway.
#'
#' @param x Expression matrix, samples by genes.
#' @param sets A named list of disjoint gene index sets.
#' @return A named list, one matrix per block: the first-order terms and the
#'   fourth-order ones.
#' @export
PathwayTerms4 <- function(x, sets) {
  t1 <- PathwayTerms(x, sets)
  z <- StandardizeN(x)$z
  q <- sapply(sets, function(a) rowSums(z[, a, drop = FALSE]^2))
  lapply(stats::setNames(names(t1), names(t1)), function(b) {
    m <- t1[[b]]
    if (startsWith(b, "within_")) {
      k <- sub("^within_", "", b)
      cbind(m, s4 = m[, "s2"]^2, s2q = m[, "s2"] * m[, "q"], q2 = m[, "q"]^2,
            z4 = rowSums(z[, sets[[k]], drop = FALSE]^4))
    } else {
      ab <- strsplit(sub("^cross_", "", b), "_")[[1]]
      cbind(m, ss2 = m[, "ss"]^2, ssqa = m[, "ss"] * q[, ab[1]],
            ssqb = m[, "ss"] * q[, ab[2]], qaqb = q[, ab[1]] * q[, ab[2]])
    }
  })
}

#' Concordance index, fast
#'
#' Harrell's concordance of predictions with an outcome, from Kendall's
#' tau-b: the share of pairs with different outcomes ordered the same way,
#' ties in the prediction counting one half. Equals [Concordance()].
#'
#' @param pred Predictions.
#' @param y Outcomes.
#' @return The concordance; 0.5 for a constant prediction.
#' @export
ConcordanceFast <- function(pred, y) {
  if (length(unique(pred)) < 2) return(0.5)
  n0 <- length(y) * (length(y) - 1) / 2
  ties <- function(v) {
    t <- table(v)
    sum(t * (t - 1) / 2)
  }
  tx <- ties(pred)
  ty <- ties(y)
  tau <- stats::cor(pred, y, method = "kendall")
  0.5 + tau * sqrt((n0 - tx) * (n0 - ty)) / (2 * (n0 - ty))
}

#' Out-of-fold predictions of an outcome from a feature set
#'
#' Gaussian elastic net (mixing 0.5), the penalty chosen by inner
#' cross-validation with folds by person, the `ref` columns unpenalized.
#'
#' @param feat Features, samples by columns.
#' @param y Outcome.
#' @param person Person of each sample.
#' @param outer Outer fold of each sample (folds by person).
#' @param ref Unpenalized columns, or `NULL`.
#' @return Out-of-fold predictions.
#' @export
OutOfFold <- function(feat, y, person, outer, ref = NULL) {
  pred <- numeric(length(y))
  if (is.null(ref)) ref <- matrix(0, length(y), 0)
  for (k in unique(outer)) {
    tr <- outer != k
    if (ncol(feat) == 0) {
      pred[!tr] <- FitReference(ref[tr, , drop = FALSE], y[tr],
                                ref[!tr, , drop = FALSE], family = "gaussian")
      next
    }
    inner <- GroupFolds(person[tr], 5L)
    pred[!tr] <- FitElasticNet(ref[tr, , drop = FALSE],
                               feat[tr, , drop = FALSE], y[tr],
                               ref[!tr, , drop = FALSE],
                               feat[!tr, , drop = FALSE], foldid = inner,
                               family = "gaussian")
  }
  pred
}

#' Child-level bootstrap tests of prediction
#'
#' Resamples people with replacement and recomputes the concordance of each
#' prediction. The superiority test declares a concordance above 0.5 when
#' the 5% bootstrap quantile exceeds 0.5; the equivalence test declares two
#' predictions equivalent when the 90% interval of their concordance
#' difference lies within the margin.
#'
#' @param preds A named list of predictions; the first is compared with the
#'   others.
#' @param y Outcome.
#' @param person Person of each sample.
#' @param n_boot Bootstrap resamples.
#' @param margin Equivalence margin.
#' @return A named numeric vector: concordance of each prediction
#'   (`c_<name>`), superiority of each (`sup_<name>`, 0 or 1), and for each
#'   comparison the difference (`d_<first>_<other>`) and equivalence
#'   (`eq_<first>_<other>`, 0 or 1).
#' @export
BootstrapTests <- function(preds, y, person, n_boot = 200L, margin = 0.05) {
  ids <- unique(person)
  rows <- split(seq_along(person), person)[as.character(ids)]
  est <- vapply(preds, ConcordanceFast, numeric(1), y = y)
  boot <- t(vapply(seq_len(n_boot), function(b) {
    take <- unlist(rows[sample.int(length(ids), replace = TRUE)],
                   use.names = FALSE)
    vapply(preds, function(pr) ConcordanceFast(pr[take], y[take]),
           numeric(1))
  }, numeric(length(preds))))
  nm <- names(preds)
  out <- c(stats::setNames(est, paste0("c_", nm)),
           stats::setNames(as.numeric(apply(boot, 2, stats::quantile, 0.05) >
                                        0.5), paste0("sup_", nm)))
  for (j in nm[-1]) {
    d <- boot[, nm[1]] - boot[, j]
    ci <- stats::quantile(d, c(0.05, 0.95))
    out[sprintf("d_%s_%s", nm[1], j)] <- est[[nm[1]]] - est[[j]]
    out[sprintf("eq_%s_%s", nm[1], j)] <- as.numeric(ci[1] > -margin &
                                                       ci[2] < margin)
  }
  out
}

#' Run one replicate of a study setting
#'
#' S1 quantities for every setting: the stability (intraclass correlation
#' over visits) of every LIONESS metric on raw and on first-order edges;
#' for connectivity and crosstalk, the R squared on the first-order terms
#' and with fourth-order terms added, and the stability of each residual;
#' the stability of each squared pathway score and of the single edges;
#' fingerprinting accuracy of the metrics (raw and first-order) and of the
#' pathway scores. S2 settings add prediction of same-visit and next-visit
#' activity from A (raw metrics), B (pathway scores and their squares), and
#' C (first-order metrics), with [BootstrapTests()].
#'
#' @param design_row One row of [StudyDesign()].
#' @param rep Replicate index; the seed is the row's seed plus `rep`.
#' @param n_boot Bootstrap resamples for the S2 tests.
#' @return A one-row data frame.
#' @export
StudyReplicate <- function(design_row, rep, n_boot = 200L) {
  set.seed(design_row$seed + rep)
  coh <- SimulateStudyCohort(design_row)
  sets <- coh$sets
  p <- ncol(coh$x)
  person <- coh$person
  e <- LionessEdges(coh$x)
  e1 <- LionessEdgesFirstOrder(coh$x)
  m_raw <- LionessMetrics(e, sets)
  m_fo <- LionessMetrics(e1, sets)
  out <- list(cell = design_row$cell, rep = rep)
  icc_raw <- EdgeICC(m_raw, person)
  icc_fo <- EdgeICC(m_fo, person)
  for (j in colnames(m_raw)) {
    out[[paste0("icc_raw_", j)]] <- icc_raw[[j]]
    out[[paste0("icc_fo_", j)]] <- icc_fo[[j]]
  }
  z <- StandardizeN(coh$x)$z
  score <- sapply(sets, function(a) rowMeans(z[, a]))
  terms4 <- PathwayTerms4(coh$x, sets)
  blocks <- PathwayBlocks(sets, p)
  edge_icc <- EdgeICC(e, person)
  for (b in names(blocks)) {
    nm <- sub("^within_", "conn_", b)
    y <- m_raw[, nm]
    f1 <- R2OnTerms(y, terms4[[b]][, 1:3])
    f4 <- R2OnTerms(y, terms4[[b]])
    out[[paste0("r2_1_", nm)]] <- f1$r2
    out[[paste0("r2_4_", nm)]] <- f4$r2
    out[[paste0("icc_res1_", nm)]] <- EdgeICC(f1$resid, person)
    out[[paste0("icc_res4_", nm)]] <- EdgeICC(f4$resid, person)
  }
  for (k in names(sets)) {
    out[[paste0("icc_score2_", k)]] <- EdgeICC(score[, k]^2, person)
    out[[paste0("icc_edges_", k)]] <- mean(edge_icc[blocks[[paste0("within_",
                                                                   k)]]])
    out[[paste0("cor_conn_score2_", k)]] <- stats::cor(m_raw[, paste0("conn_",
                                                                      k)],
                                                       score[, k]^2)
  }
  keep <- function(m) {
    m <- m[, apply(m, 2, stats::sd) > 0, drop = FALSE]
    scale(m)
  }
  feats <- list(A = keep(m_raw), B = keep(cbind(score, score^2)),
                C = keep(m_fo))
  for (s in names(feats)) {
    out[[paste0("finger_", s)]] <- FingerprintAccuracy(
      feats[[s]], person, coh$visit, seq_len(ncol(feats[[s]])))
  }
  if (design_row$part == "S2") {
    outer <- GroupFolds(person, 5L)
    y <- coh$activity
    same <- c(lapply(feats, OutOfFold, y = y, person = person, outer = outer))
    t_same <- BootstrapTests(same, y, person, n_boot)
    names(t_same) <- paste0("same_", names(t_same))
    # next visit: features and activity at visit v, activity at v + 1
    nxt <- match(paste(person, coh$visit + 1), paste(person, coh$visit))
    use <- which(!is.na(nxt))
    y_next <- y[nxt[use]]
    ref <- matrix(y[use], ncol = 1)
    nx <- c(lapply(feats, function(f) {
      OutOfFold(f[use, , drop = FALSE], y_next, person[use], outer[use],
                ref = ref)
    }), list(R0 = OutOfFold(matrix(0, length(use), 0), y_next, person[use],
                            outer[use], ref = ref)))
    t_next <- BootstrapTests(nx, y_next, person[use], n_boot)
    names(t_next) <- paste0("next_", names(t_next))
    out <- c(out, as.list(t_same), as.list(t_next))
  }
  as.data.frame(out, check.names = FALSE)
}

#' Summarize study replicates by setting
#'
#' @param reps Data frame of replicates from [StudyReplicate()].
#' @return A long data frame: `cell`, `quantity`, `mean`, `mcse`, `n_rep`.
#' @export
SummarizeStudy <- function(reps) {
  q <- setdiff(names(reps), c("cell", "rep"))
  out <- do.call(rbind, lapply(split(reps, reps$cell), function(r) {
    do.call(rbind, lapply(q, function(nm) {
      v <- r[[nm]]
      v <- v[is.finite(v)]
      if (!length(v)) return(NULL)
      data.frame(cell = r$cell[1], quantity = nm, mean = mean(v),
                 mcse = if (length(v) > 1) stats::sd(v) / sqrt(length(v))
                 else NA_real_, n_rep = length(v))
    }))
  }))
  rownames(out) <- NULL
  out
}
