#' The design table of pilot-e
#'
#' The 4 cells of `pilot-e/preregistration.md`, section 3, with seeds
#' `21161116 + 1000c`.
#'
#' @return A data frame: `cell`, `kappa_f`, `kappa_c`, `tau`, `visits`,
#'   `gamma`, `people`, `predict` (always `FALSE`), and `seed`.
#' @export
PilotEDesign <- function() {
  d <- data.frame(cell = 1:4, kappa_f = c(0, 0.55, 0.9, 0.55),
                  tau = c(0, 0, 0, 0.25))
  d$kappa_c <- 0.45
  d$visits <- 6L
  d$gamma <- 0
  d$people <- 131L
  d$predict <- FALSE
  d$seed <- 21161116L + 1000L * d$cell
  d
}

#' Gene sets of pilot-e, as node indices over genes 6 to 50
#'
#' @return A list of node indices (1 to 45): `ifn` (genes 6 to 20),
#'   `module` (21 to 30), `cellmix` (31 to 50).
#' @export
NodeSetsE <- function() {
  list(ifn = 1:15, module = 16:25, cellmix = 26:45)
}

#' Symmetric matrix of one sample's edges
#'
#' @param e A vector of edges in the order of [PairIndex()].
#' @param p Number of genes.
#' @return A `p` by `p` symmetric matrix with zero diagonal.
#' @export
EdgeMatrix <- function(e, p) {
  m <- matrix(0, p, p)
  pr <- PairIndex(p)
  m[pr] <- e
  m[pr[, 2:1]] <- e
  m
}

#' Network metrics of one sample's LIONESS network
#'
#' Connectivity of node sets (mean edge within a set, and between the first
#' two sets); node strength (signed sum of a node's edges); weighted
#' betweenness and average shortest path length with distance `1 / |e|`;
#' global efficiency (mean inverse distance); eigenvector centrality with
#' weight `|e|`; and, on the graph of edges with `|e| > threshold`, density
#' and local clustering (0 for a node with fewer than two neighbours).
#'
#' @param e A vector of edges ([PairIndex()] order).
#' @param p Number of genes (nodes).
#' @param sets A list of node index sets; connectivity is computed within
#'   each, and between the first two (named `cross`).
#' @param threshold Absolute-edge threshold of the unweighted graph.
#' @return A list: `global` (a named vector: `conn_<set>`, `conn_cross`,
#'   `pathlen`, `efficiency`, `density`) and `node` (a matrix, `p` rows,
#'   columns `strength`, `betweenness`, `eigen`, `clustering`).
#' @export
NetworkMetrics <- function(e, p, sets, threshold = 1) {
  m <- EdgeMatrix(e, p)
  conn <- vapply(sets, function(s) {
    sub <- m[s, s]
    mean(sub[upper.tri(sub)])
  }, numeric(1))
  names(conn) <- paste0("conn_", names(sets))
  cross <- mean(m[sets[[1]], sets[[2]]])
  w <- abs(m)
  dist <- 1 / pmax(w, 1e-6)
  diag(dist) <- 0
  g <- igraph::graph_from_adjacency_matrix(dist, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  d <- igraph::distances(g)
  ut <- d[upper.tri(d)]
  gw <- igraph::graph_from_adjacency_matrix(w, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  a <- (w > threshold) * 1
  diag(a) <- 0
  ga <- igraph::graph_from_adjacency_matrix(a, mode = "undirected")
  node <- cbind(
    strength = rowSums(m),
    betweenness = igraph::betweenness(g, directed = FALSE),
    eigen = igraph::eigen_centrality(gw)$vector,
    clustering = igraph::transitivity(ga, type = "local", isolates = "zero"))
  list(global = c(conn, conn_cross = cross, pathlen = mean(ut),
                  efficiency = mean(1 / ut),
                  density = sum(a) / (p * (p - 1))),
       node = node)
}

#' Jaccard overlap of the top hubs between two visits
#'
#' @param strength A matrix of node strengths, samples in rows.
#' @param person,visit Person and visit of each sample.
#' @param k Number of top nodes.
#' @return The mean Jaccard overlap of the top `k` nodes between visits 1
#'   and 2 of the same person, minus the mean over pairs of different people.
#' @export
HubOverlap <- function(strength, person, visit, k = 5L) {
  top <- function(v) apply(strength[visit == v, , drop = FALSE][
    order(person[visit == v]), , drop = FALSE], 1,
    function(r) order(-r)[seq_len(k)])
  t1 <- top(1)
  t2 <- top(2)
  w <- ncol(t1)
  jac <- matrix(0, w, w)
  for (i in seq_len(w)) for (j in seq_len(w)) {
    jac[i, j] <- length(intersect(t2[, i], t1[, j])) /
      length(union(t2[, i], t1[, j]))
  }
  mean(diag(jac)) - mean(jac[row(jac) != col(jac)])
}

#' Run one replicate of a pilot-e cell
#'
#' Simulates a cohort ([SimulateRepeatedCohort()]); computes raw,
#' first-order, and partial-correlation edges on all samples; computes
#' [NetworkMetrics()] for every sample; and returns the intraclass
#' correlation of every metric (node metrics averaged over the genes of each
#' set), the correlation of raw IFN connectivity with the squared IFN score,
#' the intraclass correlations of the squared IFN score and of raw IFN
#' edges, and the hub overlap.
#'
#' @param design_row One row of [PilotEDesign()].
#' @param rep Replicate index; the seed is the cell seed plus `rep`.
#' @return A one-row data frame.
#' @export
MetricsReplicate <- function(design_row, rep) {
  set.seed(design_row$seed + rep)
  coh <- SimulateRepeatedCohort(design_row$people, design_row$visits,
                                design_row$kappa_f, design_row$kappa_c,
                                tau = design_row$tau, gamma = design_row$gamma)
  g <- coh$x[, 6:50]
  p <- ncol(g)
  sc <- rowMeans(StandardizeN(coh$x[, 1:5])$z)
  edges <- list(raw = LionessEdges(g), fo = LionessEdgesFirstOrder(g),
                partial = LionessPartialEdges(g, cbind(sc, coh$comp)))
  sets <- NodeSetsE()
  out <- list(cell = design_row$cell, rep = rep)
  for (t in names(edges)) {
    mets <- lapply(seq_len(nrow(g)), function(q) {
      NetworkMetrics(edges[[t]][q, ], p, sets)
    })
    glob <- do.call(rbind, lapply(mets, `[[`, "global"))
    for (k in colnames(glob)) {
      out[[sprintf("icc_%s_%s", t, k)]] <- EdgeICC(glob[, k], coh$person)
    }
    for (k in colnames(mets[[1]]$node)) {
      nm <- do.call(rbind, lapply(mets, function(z) z$node[, k]))
      icc <- EdgeICC(nm, coh$person)
      for (s in names(sets)) {
        out[[sprintf("icc_%s_%s_%s", t, k, s)]] <- mean(icc[sets[[s]]],
                                                        na.rm = TRUE)
      }
      if (k == "strength") {
        out[[paste0("hub_", t)]] <- HubOverlap(nm, coh$person, coh$visit)
      }
    }
    if (t == "raw") {
      score <- rowMeans(StandardizeN(g[, sets$ifn])$z)
      out$closedform_cor <- stats::cor(glob[, "conn_ifn"], score^2)
      out$icc_score2 <- EdgeICC(score^2, coh$person)
      out$icc_score <- EdgeICC(score, coh$person)
      pg <- PairGroupsD()
      out$icc_edges_raw_ifn <- mean(EdgeICC(edges$raw, coh$person)[pg$ifn])
    }
  }
  as.data.frame(out, check.names = FALSE)
}

#' Summarize pilot-e replicates by cell
#'
#' Means and MCSEs of every column, plus `eq_<metric>` (raw minus
#' first-order intraclass correlation, for every metric), `agg_ifn` (raw IFN
#' connectivity minus raw IFN edges), `expl_ifn` (raw IFN connectivity minus
#' the squared IFN score), `fix_ifn` and `fix_module` (partial connectivity
#' of the IFN set or the module minus partial cross connectivity), and
#' `raw_ifn_excess` (the same for raw edges).
#'
#' @param reps Data frame of replicates from [MetricsReplicate()].
#' @return A long data frame: `cell`, `quantity`, `mean`, `mcse`, `n_rep`.
#' @export
SummarizeCellsE <- function(reps) {
  col <- function(nm) {
    if (nm %in% names(reps)) reps[[nm]] else rep(NA_real_, nrow(reps))
  }
  base <- setdiff(names(reps), c("cell", "rep", "seconds_cell"))
  q <- stats::setNames(lapply(base, col), base)
  raw_names <- grep("^icc_raw_", base, value = TRUE)
  for (nm in raw_names) {
    m <- sub("^icc_raw_", "", nm)
    q[[paste0("eq_", m)]] <- col(nm) - col(paste0("icc_fo_", m))
  }
  q$agg_ifn <- col("icc_raw_conn_ifn") - col("icc_edges_raw_ifn")
  q$expl_ifn <- col("icc_raw_conn_ifn") - col("icc_score2")
  q$fix_ifn <- col("icc_partial_conn_ifn") - col("icc_partial_conn_cross")
  q$fix_module <- col("icc_partial_conn_module") -
    col("icc_partial_conn_cross")
  q$raw_ifn_excess <- col("icc_raw_conn_ifn") - col("icc_raw_conn_cross")
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

#' Known-answer checks of pilot-e
#'
#' The checks of [KnownAnswerChecksD()]; network metrics against hand
#' computation on small graphs; module connectivity against its definition;
#' and, with the run's summaries, the null cell's IFN connectivity excess.
#'
#' @param summaries Optional output of [SummarizeCellsE()].
#' @param seed Seed for the simulated inputs.
#' @return A data frame as [KnownAnswerChecks()].
#' @export
KnownAnswerChecksE <- function(summaries = NULL, seed = 1L) {
  k <- KnownAnswerChecksD(seed = seed)
  add <- function(check, setting, known, computed, mcse = NA_real_,
                  type = "identity") {
    k <<- rbind(k, data.frame(check = check, setting = setting, known = known,
                              computed = computed, mcse = mcse, type = type,
                              stringsAsFactors = FALSE))
  }
  # nodes 1-2-3 with |e| = 1 on 1-2 and 2-3 and 0.4 on 1-3: distances 1,
  # 1, and min(2, 2.5) = 2, the shortest 1-3 path passing through node 2
  e <- c(1, 0.4, 1)
  nm <- NetworkMetrics(e, 3, list(a = 1:2, b = 3), threshold = 0.9)
  add("Average shortest path length", "three nodes, distances 1, 1, 2",
      4 / 3, nm$global[["pathlen"]])
  add("Global efficiency", "three nodes, distances 1, 1, 2", 5 / 6,
      nm$global[["efficiency"]])
  add("Density above the threshold", "two of three links above 0.9", 2 / 3,
      nm$global[["density"]])
  add("Betweenness of the middle node", "path graph of three nodes", 1,
      nm$node[2, "betweenness"])
  add("Strength is the signed sum of edges", "node 1: 1 + 0.4", 1.4,
      nm$node[1, "strength"])
  set.seed(seed + 3L)
  e4 <- stats::rnorm(6)
  nm4 <- NetworkMetrics(e4, 4, list(a = 1:3, b = 4))
  add("Connectivity is the mean edge within the set", "4 nodes, set 1 to 3",
      mean(e4[c(1, 2, 4)]), nm4$global[["conn_a"]])
  if (!is.null(summaries)) {
    s <- summaries[summaries$cell == 1 & summaries$quantity == "raw_ifn_excess", ]
    if (nrow(s)) {
      add("IFN connectivity stability over cross in the null cell",
          sprintf("cell 1, %d replicates", s$n_rep), 0, s$mean, s$mcse, "mcse")
    }
  }
  k
}
