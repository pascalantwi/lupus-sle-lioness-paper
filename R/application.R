# The application to GSE65391 (protocol.md, section 8, entry of 2026-11-17).

#' Read the BloodGen3 module lists from the BloodGen3Module source package
#'
#' The package keeps the module lists in `R/sysdata.rda`; this extracts that
#' file to a temporary directory and loads it as data.
#'
#' @param tarball Path to `BloodGen3Module_<version>.tar.gz`.
#' @return A list: `modules` (a data frame `module`, `gene`, `function`)
#'   and `annotation` (`module`, `function`, `cluster`).
#' @export
ReadBloodGen3 <- function(tarball) {
  dir <- tempfile("bg3")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  utils::untar(tarball, files = "BloodGen3Module/R/sysdata.rda", exdir = dir)
  e <- new.env(parent = emptyenv())
  load(file.path(dir, "BloodGen3Module", "R", "sysdata.rda"), envir = e)
  m <- e$Module_listGen3
  a <- e$Gen3_ann
  list(modules = data.frame(module = m$Module, gene = m$Gene,
                            "function" = m$Function, check.names = FALSE,
                            stringsAsFactors = FALSE),
       annotation = data.frame(module = a$Module, "function" = a$Function,
                               cluster = a$Cluster, check.names = FALSE,
                               stringsAsFactors = FALSE))
}

#' Expression by gene symbol
#'
#' Maps probes to gene symbols and keeps, for each gene, the probe with the
#' highest mean expression. Values whose 99th percentile exceeds 50 are
#' taken as linear and log2-transformed after shifting the minimum to 1.
#'
#' @param table Probes by samples.
#' @param platform Data frame `id`, `symbol` ([ReadPlatformSymbols()]).
#' @param genes Symbols wanted.
#' @return Samples by genes (the genes found, in the order of `genes`).
#' @export
ExpressionByGene <- function(table, platform, genes) {
  if (stats::quantile(table, 0.99, na.rm = TRUE) > 50) {
    table <- log2(table - min(table, na.rm = TRUE) + 1)
  }
  sym <- platform$symbol[match(rownames(table), platform$id)]
  keep <- !is.na(sym) & sym %in% genes
  tb <- table[keep, , drop = FALSE]
  sym <- sym[keep]
  mu <- rowMeans(tb, na.rm = TRUE)
  ord <- order(sym, -mu)
  first <- ord[!duplicated(sym[ord])]
  x <- t(tb[first, , drop = FALSE])
  colnames(x) <- sym[first]
  x[, intersect(genes, colnames(x)), drop = FALSE]
}

#' The four BloodGen3 pathways of the application
#'
#' @param bg3 Output of [ReadBloodGen3()].
#' @param genes Gene symbols on the array.
#' @return A named list of gene symbols: `ifn`, `neut`, `lymph`, `pb`; a
#'   gene in two pathways is kept in the first.
#' @export
BloodGen3Pathways <- function(bg3, genes) {
  labels <- list(ifn = c("Interferon", "Type 1 Interferon"),
                 neut = "Neutrophils", lymph = "Lymphocytes",
                 pb = "Plasma cells")
  m <- bg3$modules
  out <- lapply(labels, function(l) {
    unique(m$gene[m[["function"]] %in% l & m$gene %in% genes])
  })
  seen <- character(0)
  for (k in names(out)) {
    out[[k]] <- setdiff(out[[k]], seen)
    seen <- c(seen, out[[k]])
  }
  out
}

#' The edges of one sample from the exact closed form
#'
#' @param zq The sample's standardized expression (divisor `N`).
#' @param rp The cohort correlations of the pairs, in [PairIndex()] order.
#' @param pr `PairIndex(p)`.
#' @param n Number of samples.
#' @param first_order Return the first-order form instead.
#' @return The sample's edges, one per pair.
#' @export
SampleEdges <- function(zq, rp, pr, n, first_order = FALSE) {
  zi <- zq[pr[, 1]]
  zj <- zq[pr[, 2]]
  if (first_order) return(rp + zi * zj - rp / 2 * (zi^2 + zj^2))
  eps <- 1 / (n - 1)
  n * rp - (n - 1) * (rp - eps * zi * zj) /
    sqrt((1 - eps * zi^2) * (1 - eps * zj^2))
}

#' LIONESS metrics of every sample, one sample at a time
#'
#' The metrics of [LionessMetrics()] for networks too large to hold all
#' samples' edges: each sample's edges are computed from the exact closed
#' form (or its first-order form), its metrics taken, and the edges
#' dropped.
#'
#' @param x Expression, samples by genes.
#' @param sets A named list of gene index sets.
#' @param first_order Use the first-order form of the edges.
#' @param threshold,prop As [LionessMetrics()].
#' @param cl Optional cluster for [parallel::parLapply()].
#' @return A matrix, samples by metrics, with the columns of
#'   [LionessMetrics()].
#' @export
LionessMetricsBySample <- function(x, sets, first_order = FALSE,
                                   threshold = 1, prop = 0.05, cl = NULL) {
  # evaluate every argument here: the closure below is sent to the workers
  force(sets)
  force(first_order)
  force(threshold)
  force(prop)
  n <- nrow(x)
  p <- ncol(x)
  z <- StandardizeN(x)$z
  r <- stats::cor(x)
  pr <- PairIndex(p)
  rp <- r[pr]
  one <- function(q) {
    eq <- SampleEdges(z[q, ], rp, pr, n, first_order)
    LionessMetricsOne(eq, pr, p, sets, threshold, prop)
  }
  rows <- if (is.null(cl)) {
    lapply(seq_len(n), one)
  } else {
    parallel::parLapply(cl, seq_len(n), one)
  }
  do.call(rbind, rows)
}

#' LIONESS metrics of one sample's edges
#'
#' @param eq One sample's edges ([PairIndex()] order).
#' @param pr `PairIndex(p)`.
#' @param p Number of genes.
#' @param sets A named list of gene index sets.
#' @param threshold,prop As [LionessMetrics()].
#' @return A named vector, the columns of [LionessMetrics()].
#' @export
LionessMetricsOne <- function(eq, pr, p, sets, threshold = 1, prop = 0.05) {
  blocks <- PathwayBlocks(sets, p)
  bm <- vapply(blocks, function(cols) mean(eq[cols]), numeric(1))
  names(bm) <- sub("^within_", "conn_", names(bm))
  strength <- as.vector(rowsum(c(eq, eq), c(pr[, 1], pr[, 2]),
                               reorder = TRUE))
  pw_mean <- function(v, nm) {
    stats::setNames(vapply(sets, function(s) mean(v[s]), numeric(1)),
                    paste0(nm, "_", names(sets)))
  }
  a <- GraphMetrics(eq, pr, p, threshold = threshold)
  b <- GraphMetrics(eq, pr, p, prop = prop)
  c(total_weight = sum(eq), bm, pw_mean(strength, "strength"),
    n_edges = a$n_edges, n_nodes = a$n_nodes, transitivity = a$transitivity,
    n_nodes_p = b$n_nodes, transitivity_p = b$transitivity,
    pw_mean(a$degree, "degree"), pw_mean(a$betweenness, "betweenness"),
    pw_mean(b$degree, "degree_p"), pw_mean(b$betweenness, "betweenness_p"))
}

#' Diagnostics of one module's LIONESS connectivity
#'
#' @param x Expression of the module's genes, samples by genes.
#' @param person Person of each sample.
#' @return A named vector: number of genes; $R^2$ of connectivity on its
#'   first-order terms and with the fourth-order terms; stability of the
#'   connectivity, of each residual, and of the module score; the share of
#'   variance of the first principal component.
#' @export
ModuleDiagnostics <- function(x, person) {
  sets <- list(m = seq_len(ncol(x)))
  e <- LionessEdgesClosedForm(x)
  conn <- rowMeans(e)
  t4 <- PathwayTerms4(x, sets)$within_m
  f1 <- R2OnTerms(conn, t4[, 1:3])
  f4 <- R2OnTerms(conn, t4)
  z <- StandardizeN(x)$z
  ev <- eigen(stats::cor(x), symmetric = TRUE, only.values = TRUE)$values
  c(n_genes = ncol(x), r2_1 = f1$r2, r2_4 = f4$r2,
    icc_conn = EdgeICC(conn, person), icc_res1 = EdgeICC(f1$resid, person),
    icc_res4 = EdgeICC(f4$resid, person),
    icc_score = EdgeICC(rowMeans(z), person),
    icc_score2 = EdgeICC(rowMeans(z)^2, person), pc1_share = ev[1] / sum(ev))
}

#' Concordance on a fixed training and test split
#'
#' Fits each feature set on the training people (elastic net as
#' [OutOfFold()], inner folds by person) and returns the concordance of the
#' test people's predictions.
#'
#' @param feats A named list of feature matrices.
#' @param y Outcome.
#' @param person Person of each sample.
#' @param train Logical, the training samples.
#' @return Named concordances.
#' @export
SplitConcordance <- function(feats, y, person, train) {
  vapply(feats, function(f) {
    inner <- GroupFolds(person[train], 5L)
    pred <- FitElasticNet(matrix(0, sum(train), 0), f[train, , drop = FALSE],
                          y[train], matrix(0, sum(!train), 0),
                          f[!train, , drop = FALSE], foldid = inner,
                          family = "gaussian")
    ConcordanceFast(pred, y[!train])
  }, numeric(1))
}
