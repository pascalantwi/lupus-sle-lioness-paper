#' The design table of pilot-g
#'
#' The 3 cells of `pilot-g/preregistration.md`, section 3, with seeds
#' `21361116 + 1000c`.
#'
#' @return A data frame: `cell`, `kappa_ifn`, `kappa_neut`, `kappa_lymph`,
#'   `kappa_pb`, `tau`, `people`, `visits`, `predict` (`FALSE`), and `seed`.
#' @export
PilotGDesign <- function() {
  d <- data.frame(cell = 1:3, kappa_ifn = c(0, 0.55, 0.9),
                  kappa_neut = c(0, 0.41, 0.41),
                  kappa_lymph = c(0, 0.48, 0.48),
                  kappa_pb = c(0, 0.39, 0.39), tau = c(0, 0.25, 0.25))
  d$people <- 131L
  d$visits <- 6L
  d$predict <- FALSE
  d$seed <- 21361116L + 1000L * d$cell
  d
}

#' The six pathways of pilot-g
#'
#' @return A named list of column indices: `ifn`, `neut`, `lymph`, `pb`,
#'   `personal`, `null`, 15 genes each, in that order.
#' @export
PathwaySetsG <- function() {
  nm <- c("ifn", "neut", "lymph", "pb", "personal", "null")
  idx <- lapply(seq_along(nm), function(k) (15 * (k - 1) + 1):(15 * k))
  stats::setNames(idx, nm)
}

#' Simulate repeated visits with six pathways
#'
#' Ninety genes in the six pathways of [PathwaySetsG()]. The first four
#' pathways load on their own factor (IFN loading 0.8, the others 0.6), each
#' the person's level plus visit noise with the given within-person
#' stability; the factors are independent. The personal pathway has noise
#' correlation `0.3 + tau * s` with `s` uniform on (-1, 1) fixed for the
#' person; the null pathway has noise correlation 0.3. Every gene carries
#' the person's own level (variance 0.5) and unit noise.
#'
#' @param people,visits Number of people and visits per person.
#' @param kappa_ifn,kappa_neut,kappa_lymph,kappa_pb Within-person stability
#'   of the four factors.
#' @param tau Spread of the personal pathway's correlation.
#' @return A list with `x`, `person`, and `visit`.
#' @export
SimulatePathwayCohort <- function(people, visits, kappa_ifn, kappa_neut,
                                  kappa_lymph, kappa_pb, tau = 0) {
  n <- people * visits
  person <- rep(seq_len(people), each = visits)
  visit <- rep(seq_len(visits), times = people)
  sets <- PathwaySetsG()
  p <- 90L
  level <- function(k) {
    stats::rnorm(people, sd = sqrt(k))[person] +
      stats::rnorm(n, sd = sqrt(1 - k))
  }
  x <- matrix(stats::rnorm(people * p, sd = sqrt(0.5)), people, p)[person, ] +
    matrix(stats::rnorm(n * p), n, p)
  kap <- c(ifn = kappa_ifn, neut = kappa_neut, lymph = kappa_lymph,
           pb = kappa_pb)
  load <- c(ifn = 0.8, neut = 0.6, lymph = 0.6, pb = 0.6)
  for (k in names(kap)) {
    x[, sets[[k]]] <- x[, sets[[k]]] + load[[k]] * level(kap[[k]])
  }
  # one shared draw per sample gives equicorrelation rho between the genes
  correlated <- function(cols, rho) {
    m <- stats::rnorm(n)
    sqrt(rho) * m + sqrt(1 - rho) * matrix(stats::rnorm(n * length(cols)), n)
  }
  rho_p <- 0.3 + tau * stats::runif(people, -1, 1)[person]
  stopifnot(all(rho_p > 0), all(rho_p < 1))
  # the noise of the personal and null pathways is replaced by correlated
  # noise; the person levels stay
  for (k in c("personal", "null")) {
    cols <- sets[[k]]
    rho <- if (k == "personal") rho_p else rep(0.3, n)
    x[, cols] <- matrix(stats::rnorm(people * length(cols), sd = sqrt(0.5)),
                        people)[person, ] + correlated(cols, rho)
  }
  colnames(x) <- paste0("g", seq_len(p))
  list(x = x, person = person, visit = visit)
}

#' Column indices of the within-pathway and between-pathway pairs
#'
#' @param sets A list of disjoint gene index sets.
#' @param p Number of genes.
#' @return A named list of pair-column indices ([PairIndex()] of `p`):
#'   `within_<A>` for each set and `cross_<A>_<B>` for each pair of sets.
#' @export
PathwayBlocks <- function(sets, p) {
  pr <- PairIndex(p)
  memb <- rep(NA_character_, p)
  for (k in names(sets)) memb[sets[[k]]] <- k
  a <- memb[pr[, 1]]
  b <- memb[pr[, 2]]
  out <- list()
  for (k in names(sets)) out[[paste0("within_", k)]] <- which(a == k & b == k)
  nm <- names(sets)
  for (i in seq_along(nm)) for (j in seq_along(nm)) if (i < j) {
    out[[sprintf("cross_%s_%s", nm[i], nm[j])]] <-
      which((a == nm[i] & b == nm[j]) | (a == nm[j] & b == nm[i]))
  }
  out
}

#' Pathway connectivity and crosstalk of every sample
#'
#' @param e Edges, samples by pairs ([PairIndex()] order).
#' @param blocks Output of [PathwayBlocks()].
#' @return A matrix, samples by blocks: the mean edge of each block.
#' @export
PathwayMetrics <- function(e, blocks) {
  m <- vapply(blocks, function(cols) rowMeans(e[, cols, drop = FALSE]),
              numeric(nrow(e)))
  # a single sample still gives a one-row matrix
  matrix(m, nrow(e), length(blocks), dimnames = list(NULL, names(blocks)))
}

#' The expression terms of pathway connectivity and crosstalk
#'
#' For a pathway `A`: the squared sum of its standardized genes, the sum of
#' their squares, and the sum of their squares weighted by each gene's summed
#' reference correlation with the other genes of `A`. For a pair `A`, `B`:
#' the product of the two sums and the squares of each set weighted by the
#' summed reference correlation with the other set. Standardization and
#' correlations use all samples (divisor `N`), as LIONESS does.
#'
#' @param x Expression matrix, samples by genes.
#' @param sets A list of disjoint gene index sets.
#' @return A named list, one three-column matrix per block of
#'   [PathwayBlocks()].
#' @export
PathwayTerms <- function(x, sets) {
  z <- StandardizeN(x)$z
  r <- stats::cor(x)
  diag(r) <- 0
  s <- sapply(sets, function(a) rowSums(z[, a, drop = FALSE]))
  out <- list()
  for (k in names(sets)) {
    a <- sets[[k]]
    out[[paste0("within_", k)]] <- cbind(
      s2 = s[, k]^2, q = rowSums(z[, a]^2),
      w = drop(z[, a]^2 %*% rowSums(r[a, a])))
  }
  nm <- names(sets)
  for (i in seq_along(nm)) for (j in seq_along(nm)) if (i < j) {
    a <- sets[[nm[i]]]
    b <- sets[[nm[j]]]
    out[[sprintf("cross_%s_%s", nm[i], nm[j])]] <- cbind(
      ss = s[, nm[i]] * s[, nm[j]],
      wa = drop(z[, a]^2 %*% rowSums(r[a, b])),
      wb = drop(z[, b]^2 %*% rowSums(r[b, a])))
  }
  out
}

#' R squared and residuals of a regression on expression terms
#'
#' @param y A metric, one value per sample.
#' @param terms A matrix of terms, samples by columns.
#' @return A list with `r2` and `resid`.
#' @export
R2OnTerms <- function(y, terms) {
  q <- qr(cbind(1, terms))
  res <- qr.resid(q, y)
  list(r2 = 1 - sum(res^2) / sum((y - mean(y))^2), resid = res)
}

#' Run one replicate of a pilot-g cell
#'
#' @param design_row One row of [PilotGDesign()].
#' @param rep Replicate index; the seed is the cell seed plus `rep`.
#' @return A one-row data frame: the `R^2` of every block on its terms
#'   (`r2_<block>`), their smallest values over within and cross blocks, the
#'   correlation of IFN connectivity with its squared score, and the
#'   intraclass correlations of each pathway's connectivity, squared score,
#'   residual, and single edges.
#' @export
PathwayReplicate <- function(design_row, rep) {
  set.seed(design_row$seed + rep)
  coh <- SimulatePathwayCohort(design_row$people, design_row$visits,
                               design_row$kappa_ifn, design_row$kappa_neut,
                               design_row$kappa_lymph, design_row$kappa_pb,
                               design_row$tau)
  sets <- PathwaySetsG()
  p <- ncol(coh$x)
  blocks <- PathwayBlocks(sets, p)
  e <- LionessEdges(coh$x)
  met <- PathwayMetrics(e, blocks)
  terms <- PathwayTerms(coh$x, sets)
  z <- StandardizeN(coh$x)$z
  out <- list(cell = design_row$cell, rep = rep)
  r2 <- numeric(0)
  edge_icc <- EdgeICC(e, coh$person)
  for (b in names(blocks)) {
    fit <- R2OnTerms(met[, b], terms[[b]])
    out[[paste0("r2_", b)]] <- fit$r2
    r2[b] <- fit$r2
    if (startsWith(b, "within_")) {
      k <- sub("within_", "", b)
      score <- rowMeans(z[, sets[[k]]])
      out[[paste0("icc_conn_", k)]] <- EdgeICC(met[, b], coh$person)
      out[[paste0("icc_score2_", k)]] <- EdgeICC(score^2, coh$person)
      out[[paste0("icc_resid_", k)]] <- EdgeICC(fit$resid, coh$person)
      out[[paste0("icc_edges_", k)]] <- mean(edge_icc[blocks[[b]]])
      if (k == "ifn") out$cor_ifn <- stats::cor(met[, b], score^2)
    }
  }
  out$min_r2_within <- min(r2[startsWith(names(r2), "within_")])
  out$min_r2_cross <- min(r2[startsWith(names(r2), "cross_")])
  # first-order edges: the closed forms are exact, so R^2 is 1 up to rounding
  met_fo <- PathwayMetrics(LionessEdgesFirstOrder(coh$x), blocks)
  out$min_r2_fo <- min(sapply(names(blocks), function(b) {
    R2OnTerms(met_fo[, b], terms[[b]])$r2
  }))
  as.data.frame(out, check.names = FALSE)
}

#' Summarize pilot-g replicates by cell
#'
#' Means and MCSEs of every column, plus `expl_ifn` (IFN connectivity minus
#' squared-score stability), `agg_ifn` (IFN connectivity minus single-edge
#' stability), and `personal_vs_null` (personal minus null pathway
#' connectivity stability).
#'
#' @param reps Data frame of replicates from [PathwayReplicate()].
#' @return A long data frame: `cell`, `quantity`, `mean`, `mcse`, `n_rep`.
#' @export
SummarizeCellsG <- function(reps) {
  col <- function(nm) {
    if (nm %in% names(reps)) reps[[nm]] else rep(NA_real_, nrow(reps))
  }
  base <- setdiff(names(reps), c("cell", "rep", "seconds_cell"))
  q <- stats::setNames(lapply(base, col), base)
  q$expl_ifn <- col("icc_conn_ifn") - col("icc_score2_ifn")
  q$agg_ifn <- col("icc_conn_ifn") - col("icc_edges_ifn")
  q$personal_vs_null <- col("icc_conn_personal") - col("icc_conn_null")
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

#' Known-answer checks of pilot-g
#'
#' The checks of [KnownAnswerChecksE()]; the first-order connectivity and
#' crosstalk against their closed forms (proposal section 14.1); and, with
#' the run's summaries, the null cell's residual stabilities.
#'
#' @param summaries Optional output of [SummarizeCellsG()].
#' @param seed Seed for the simulated inputs.
#' @return A data frame as [KnownAnswerChecks()].
#' @export
KnownAnswerChecksG <- function(summaries = NULL, seed = 1L) {
  k <- KnownAnswerChecksE(seed = seed)
  add <- function(check, setting, known, computed, mcse = NA_real_,
                  type = "identity") {
    k <<- rbind(k, data.frame(check = check, setting = setting, known = known,
                              computed = computed, mcse = mcse, type = type,
                              stringsAsFactors = FALSE))
  }
  set.seed(seed + 4L)
  x <- matrix(stats::rnorm(60 * 9), 60, 9) + stats::rnorm(60)
  sets <- list(a = 1:4, b = 5:9)
  blocks <- PathwayBlocks(sets, 9)
  met <- PathwayMetrics(LionessEdgesFirstOrder(x), blocks)
  z <- StandardizeN(x)$z
  r <- stats::cor(x)
  diag(r) <- 0
  terms <- PathwayTerms(x, sets)
  ca <- 4 * 3 / 2
  ra <- mean(r[1:4, 1:4][upper.tri(r[1:4, 1:4])])
  tw <- terms$within_a
  within_cf <- ra + (tw[, "s2"] - tw[, "q"]) / (2 * ca) - tw[, "w"] / (2 * ca)
  add("First-order connectivity equals its closed form",
      "9 genes, sets of 4 and 5; largest absolute difference", 0,
      max(abs(met[, "within_a"] - within_cf)))
  rab <- mean(r[1:4, 5:9])
  cross_cf <- rab + terms$cross_a_b[, "ss"] / 20 -
    (terms$cross_a_b[, "wa"] + terms$cross_a_b[, "wb"]) / 40
  add("First-order crosstalk equals its closed form",
      "9 genes, sets of 4 and 5; largest absolute difference", 0,
      max(abs(met[, "cross_a_b"] - cross_cf)))
  if (!is.null(summaries)) {
    for (pw in names(PathwaySetsG())) {
      s <- summaries[summaries$cell == 1 &
                       summaries$quantity == paste0("icc_resid_", pw), ]
      if (!nrow(s)) next
      add(sprintf("Residual connectivity stability, %s pathway, null cell", pw),
          sprintf("cell 1, %d replicates", s$n_rep), 0, s$mean, s$mcse, "mcse")
    }
  }
  k
}

#' Example networks of one simulated cohort, for figures
#'
#' Simulates one cohort of a design row, computes its LIONESS edges, and
#' keeps what the figures draw: every sample's pathway metrics and squared
#' pathway scores, and the full edge vectors of every visit of a few
#' people.
#'
#' @param design_row One row of a design table with the columns of
#'   [PilotGDesign()].
#' @param rep Replicate index; the seed is the row's seed plus `rep`.
#' @param people_kept How many people's edges to keep: the person with the
#'   highest and the one with the lowest mean interferon score, then the
#'   people nearest the median, in that order.
#' @return A list: `person`, `visit`, `metrics` (samples by blocks),
#'   `score2` (samples by pathways), `kept` (the people kept), `edges` (their
#'   samples by pairs), `edges_person`, `edges_visit`, and `sets`.
#' @export
ExampleNetworks <- function(design_row, rep = 1L, people_kept = 2L) {
  set.seed(design_row$seed + rep)
  coh <- SimulatePathwayCohort(design_row$people, design_row$visits,
                               design_row$kappa_ifn, design_row$kappa_neut,
                               design_row$kappa_lymph, design_row$kappa_pb,
                               design_row$tau)
  sets <- PathwaySetsG()
  e <- LionessEdges(coh$x)
  met <- PathwayMetrics(e, PathwayBlocks(sets, ncol(coh$x)))
  z <- StandardizeN(coh$x)$z
  score2 <- sapply(sets, function(a) rowMeans(z[, a])^2)
  ifn <- tapply(rowMeans(z[, sets$ifn]), coh$person, mean)
  ord <- order(abs(ifn - stats::median(ifn)))
  kept <- unique(c(as.integer(names(which.max(ifn))),
                   as.integer(names(which.min(ifn))),
                   as.integer(names(ifn)[ord])))[seq_len(people_kept)]
  rows <- which(coh$person %in% kept)
  list(person = coh$person, visit = coh$visit, metrics = met,
       score2 = score2, kept = kept, edges = e[rows, , drop = FALSE],
       edges_person = coh$person[rows], edges_visit = coh$visit[rows],
       sets = sets)
}
