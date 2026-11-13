#' Area under the ROC curve
#'
#' The Mann-Whitney form: the probability that a random class-1 score exceeds
#' a random class-0 score, ties counted one half.
#'
#' @param score Numeric scores.
#' @param y Class labels, 0 or 1.
#' @return The AUC.
#' @export
Auc <- function(score, y) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  rk <- rank(score)
  (sum(rk[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Fold assignment by group
#'
#' Every sample of a group (a woman) goes to the same fold; groups are dealt
#' to folds in random order, so fold sizes differ by at most one group.
#'
#' @param group Group identifier of each sample.
#' @param k Number of folds.
#' @return An integer vector of folds, one per sample.
#' @export
GroupFolds <- function(group, k = 5L) {
  g <- unique(group)
  k <- min(k, length(g))
  fg <- sample(rep_len(seq_len(k), length(g)))
  fg[match(group, g)]
}

#' Features of the pilot's models, fitted on a training fold
#'
#' The IFN score is the mean of genes 1 to 5 standardized by the training
#' fold; its residual on composition uses the training regression. The gene
#' set is genes 6 to 50. Products, squares, and first-order edges use the
#' training fold's correlation and standardization for both parts. LIONESS
#' edges of training samples use the training fold as reference (exact
#' closed form); test samples are scored by the add-one form.
#'
#' @param cohort A list from [SimulatePregnancyCohort()].
#' @param train Logical vector, `TRUE` for training samples.
#' @param score_genes,gene_set Column indices of the score genes and the
#'   network gene set.
#' @return A list with `train` and `test`, each holding `r0`, `r1` (base
#'   covariates), `genes`, `quad` (products and squares), `edges`, and `fo`.
#' @export
PregnancyFeatures <- function(cohort, train, score_genes = 1:5,
                              gene_set = 6:50) {
  x <- cohort$x
  sc <- StandardizeN(x[train, score_genes, drop = FALSE])
  score <- rowMeans(StandardizeN(x[, score_genes, drop = FALSE], sc$center,
                                 sc$scale)$z)
  fit <- stats::lm.fit(cbind(1, cohort$comp[train]), score[train])
  resid <- score - drop(cbind(1, cohort$comp) %*% fit$coefficients)
  r0 <- cbind(comp = cohort$comp, ifn = resid)
  r1 <- cbind(r0, ifn2 = resid^2, comp2 = cohort$comp^2)
  g <- x[, gene_set, drop = FALSE]
  gtr <- g[train, , drop = FALSE]
  gs <- StandardizeN(gtr)
  z <- StandardizeN(g, gs$center, gs$scale)$z
  pairs <- PairIndex(length(gene_set))
  quad <- cbind(z[, pairs[, 1]] * z[, pairs[, 2]], z^2)
  r_tr <- stats::cor(gtr)[pairs]
  r_mat <- matrix(r_tr, nrow = nrow(z), ncol = length(r_tr), byrow = TRUE)
  fo <- r_mat + z[, pairs[, 1]] * z[, pairs[, 2]] -
    r_mat / 2 * (z[, pairs[, 1]]^2 + z[, pairs[, 2]]^2)
  part <- function(rows, edges) {
    list(r0 = r0[rows, , drop = FALSE], r1 = r1[rows, , drop = FALSE],
         genes = g[rows, , drop = FALSE], quad = quad[rows, , drop = FALSE],
         edges = edges, fo = fo[rows, , drop = FALSE])
  }
  list(train = part(train, LionessEdges(gtr)),
       test = part(!train, LionessTestEdges(g[!train, , drop = FALSE], gtr)))
}

#' Fit a logistic base model and score a test set
#'
#' @param x_train,x_test Covariate matrices.
#' @param y Training labels.
#' @return Linear predictors for the test samples.
#' @export
FitReference <- function(x_train, y, x_test) {
  fit <- suppressWarnings(stats::glm.fit(cbind(1, x_train), y,
                                         family = stats::binomial()))
  b <- fit$coefficients
  b[is.na(b)] <- 0
  drop(cbind(1, x_test) %*% b)
}

#' Fit an elastic net with unpenalized base terms and score a test set
#'
#' `cv.glmnet` with the given mixing, the given inner folds (each woman's
#' samples in one fold), binomial deviance, and `lambda.min`; the base
#' columns have penalty factor 0.
#'
#' @param ref_train,ref_test Base covariate matrices (unpenalized).
#' @param pen_train,pen_test Penalized feature matrices.
#' @param y Training labels.
#' @param foldid Inner fold of each training sample.
#' @param alpha The elastic-net mixing parameter.
#' @return Linear predictors for the test samples.
#' @export
FitElasticNet <- function(ref_train, pen_train, y, ref_test, pen_test,
                          foldid, alpha = 0.5) {
  xtr <- cbind(ref_train, pen_train)
  xte <- cbind(ref_test, pen_test)
  pf <- c(rep(0, ncol(ref_train)), rep(1, ncol(pen_train)))
  # cv.glmnet, or its prediction, can fail when a path ends after one
  # penalty value; only then the fallback runs
  out <- tryCatch({
    fit <- glmnet::cv.glmnet(xtr, y, family = "binomial", alpha = alpha,
                             foldid = foldid, penalty.factor = pf,
                             type.measure = "deviance")
    drop(stats::predict(fit, newx = xte, s = "lambda.min", type = "link"))
  }, error = function(e) NULL)
  if (is.null(out)) out <- ElasticNetFallback(xtr, y, xte, pf, foldid, alpha)
  out
}

#' Cross-validated elastic net that tolerates truncated paths
#'
#' The fallback of [FitElasticNet()]: the penalty sequence of the fit on all
#' training samples; the given inner folds; each fold's linear predictors from
#' its coefficients, the last available fit carried forward where a path
#' stops early, and an intercept-only model where a fit has no path; the
#' penalty with the smallest binomial deviance.
#'
#' @param xtr,xte Training and test feature matrices.
#' @param y Training labels.
#' @param pf Penalty factors.
#' @param foldid Inner fold of each training sample.
#' @param alpha The elastic-net mixing parameter.
#' @return Linear predictors for the test samples.
#' @export
ElasticNetFallback <- function(xtr, y, xte, pf, foldid, alpha = 0.5) {
  fit_path <- function(x, yy, lambda = NULL) {
    suppressWarnings(glmnet::glmnet(x, yy, family = "binomial", alpha = alpha,
                                    penalty.factor = pf, lambda = lambda))
  }
  coefs <- function(g, n_lambda, yy) {
    b <- tryCatch(as.matrix(stats::coef(g)), error = function(e) NULL)
    if (is.null(b) || ncol(b) == 0 || nrow(b) != ncol(xtr) + 1) {
      b <- matrix(c(stats::qlogis(mean(yy)), rep(0, ncol(xtr))), ncol = 1)
    }
    b[, pmin(seq_len(n_lambda), ncol(b)), drop = FALSE]
  }
  full <- fit_path(xtr, y)
  lam <- full$lambda
  dev <- matrix(NA_real_, length(y), length(lam))
  for (f in unique(foldid)) {
    held <- foldid == f
    g <- fit_path(xtr[!held, , drop = FALSE], y[!held], lambda = lam)
    eta <- cbind(1, xtr[held, , drop = FALSE]) %*%
      coefs(g, length(lam), y[!held])
    p <- pmin(pmax(stats::plogis(eta), 1e-12), 1 - 1e-12)
    dev[held, ] <- -2 * (y[held] * log(p) + (1 - y[held]) * log(1 - p))
  }
  best <- which.min(colMeans(dev))
  drop(cbind(1, xte) %*% coefs(full, length(lam), y)[, best])
}

#' Apparent rewiring: the share of pairs whose edge changes with pregnancy
#'
#' LIONESS edges on all samples; for each woman, her mean edge in pregnancy
#' minus her mean edge outside it; a one-sample t-test across women for each
#' pair; the share of pairs with `p < level`.
#'
#' @param x Samples by genes.
#' @param y 1 pregnant, 0 not.
#' @param woman Woman of each sample.
#' @param level The test's level.
#' @return The share of gene pairs called different.
#' @export
ApparentRewiring <- function(x, y, woman, level = 0.05) {
  e <- LionessEdges(x)
  m1 <- rowsum(e[y == 1, , drop = FALSE], woman[y == 1]) /
    as.vector(table(woman[y == 1]))
  m0 <- rowsum(e[y == 0, , drop = FALSE], woman[y == 0]) /
    as.vector(table(woman[y == 0]))
  d <- m1 - m0[rownames(m1), , drop = FALSE]
  w <- nrow(d)
  t <- colMeans(d) / (apply(d, 2, stats::sd) / sqrt(w))
  pv <- 2 * stats::pt(-abs(t), w - 1)
  mean(pv < level)
}

#' Run one replicate of a pilot cell
#'
#' Simulates a cohort, runs the outer 5-fold cross-validation by woman for the
#' seven models (inner folds also by woman), and returns the pooled
#' out-of-fold AUC of each, with the apparent rewiring among the pairs of
#' genes 6 to 20 and among the module pairs.
#'
#' @param design_row One row of [PilotDesign()].
#' @param rep Replicate index; the seed is the cell seed plus `rep`.
#' @param models Models to fit.
#' @return A one-row data frame.
#' @export
RunReplicate <- function(design_row, rep,
                         models = c("R0", "R1", "E1", "E5", "F1", "N1",
                                    "N1R0")) {
  set.seed(design_row$seed + rep)
  coh <- SimulatePregnancyCohort(design_row$women, design_row$factor,
                                 design_row$module)
  y <- coh$y
  folds <- GroupFolds(coh$woman, 5L)
  pred <- matrix(NA_real_, length(y), length(models),
                 dimnames = list(NULL, models))
  for (f in 1:5) {
    tr <- folds != f
    ft <- PregnancyFeatures(coh, tr)
    a <- ft$train
    b <- ft$test
    ytr <- y[tr]
    inner <- GroupFolds(coh$woman[tr], 5L)
    en <- function(rtr, ptr, rte, pte) {
      FitElasticNet(rtr, ptr, ytr, rte, pte, foldid = inner)
    }
    out <- list(
      R0 = function() FitReference(a$r0, ytr, b$r0),
      R1 = function() FitReference(a$r1, ytr, b$r1),
      E1 = function() en(a$r1, a$genes, b$r1, b$genes),
      E5 = function() en(a$r1, cbind(a$genes, a$quad), b$r1,
                         cbind(b$genes, b$quad)),
      F1 = function() en(a$r1, cbind(a$genes, a$fo), b$r1,
                         cbind(b$genes, b$fo)),
      N1 = function() en(a$r1, cbind(a$genes, a$edges), b$r1,
                         cbind(b$genes, b$edges)),
      N1R0 = function() en(a$r0, a$edges, b$r0, b$edges)
    )
    for (m in models) pred[!tr, m] <- out[[m]]()
  }
  auc <- apply(pred, 2, Auc, y = y)
  data.frame(cell = design_row$cell, rep = rep, t(auc),
             rewire_factor = ApparentRewiring(coh$x[, 6:20], y, coh$woman),
             rewire_module = ApparentRewiring(coh$x[, 21:30], y, coh$woman),
             check.names = FALSE)
}

#' Summarize replicates by cell
#'
#' Means and Monte Carlo standard errors of each model's AUC, the paired
#' differences the pass rule uses, and the apparent rewiring.
#'
#' @param reps Data frame of replicates from [RunReplicate()].
#' @return A long data frame: `cell`, `quantity`, `mean`, `mcse`, `n_rep`.
#' @export
SummarizeCells <- function(reps) {
  # exact column lookup, so that no name partially matches another
  col <- function(nm) {
    if (nm %in% names(reps)) reps[[nm]] else rep(NA_real_, nrow(reps))
  }
  q <- list(
    R0 = col("R0"), R1 = col("R1"), E1 = col("E1"), E5 = col("E5"),
    F1 = col("F1"), N1 = col("N1"), N1R0 = col("N1R0"),
    d1_N1_vs_E1 = col("N1") - col("E1"),
    d2_N1_vs_F1 = col("N1") - col("F1"),
    artifact_N1R0_vs_R0 = col("N1R0") - col("R0"),
    gain_R1_vs_R0 = col("R1") - col("R0"),
    gain_E5_vs_E1 = col("E5") - col("E1"),
    par_N1_vs_E5 = col("N1") - col("E5"),
    rewire_factor = col("rewire_factor"),
    rewire_module = col("rewire_module")
  )
  out <- do.call(rbind, lapply(split(seq_len(nrow(reps)), reps$cell),
                               function(rows) {
    do.call(rbind, lapply(names(q), function(nm) {
      v <- q[[nm]][rows]
      v <- v[!is.na(v)]
      data.frame(cell = reps$cell[rows[1]], quantity = nm,
                 mean = if (length(v)) mean(v) else NA_real_,
                 mcse = if (length(v) > 1) stats::sd(v) / sqrt(length(v))
                 else NA_real_,
                 n_rep = length(v))
    }))
  }))
  rownames(out) <- NULL
  out
}
