test_that("StudyDesign has the 13 S1 and 5 S2 settings", {
  d <- StudyDesign()
  expect_equal(sum(d$part == "S1"), 13)
  expect_equal(sum(d$part == "S2"), 5)
  expect_equal(d$seed, 21461117L + 1000L * d$cell)
  base <- d[d$family == "base", ]
  expect_equal(c(base$kappa_ifn, base$load_ifn, base$size_ifn, base$tau),
               c(0.55, 0.8, 15, 0.25))
  expect_true(all(d$beta[d$part == "S1"] == 0))
})

test_that("StudySets sizes the interferon pathway", {
  s <- StudySets(40L)
  expect_equal(unname(lengths(s)), c(40L, rep(15L, 5)))
  expect_equal(sort(unlist(s, use.names = FALSE)), 1:115)
})

test_that("the activity score has the measured within-person stability", {
  set.seed(81)
  row <- StudyDesign()[1, ]
  row$people <- 400L
  coh <- SimulateStudyCohort(row)
  expect_lt(abs(EdgeICC(coh$activity, coh$person) - 0.39), 0.05)
  expect_equal(ncol(coh$x), 90)
})

test_that("GraphMetrics matches a hand-computed graph", {
  # three genes: pairs (1,2), (1,3), (2,3); only (1,2) and (2,3) above 1
  pr <- PairIndex(3)
  g <- GraphMetrics(c(2, 0.5, -1.5), pr, 3, threshold = 1)
  expect_equal(g$n_edges, 2)
  expect_equal(g$n_nodes, 3)
  expect_equal(unname(g$degree), c(1, 2, 1))
  expect_equal(unname(g$betweenness), c(0, 1, 0))
  expect_equal(g$transitivity, 0)
  # all three pairs kept: a triangle
  expect_equal(GraphMetrics(c(2, 1.5, -1.5), pr, 3)$transitivity, 1)
  expect_equal(GraphMetrics(c(2, 0.5, -1.5), pr, 3, prop = 1 / 3)$n_edges, 1)
})

test_that("LionessMetrics totals and strengths add up", {
  set.seed(82)
  sets <- list(a = 1:3, b = 4:6)
  e <- matrix(rnorm(2 * 15), 2)
  m <- LionessMetrics(e, sets)
  expect_equal(unname(m[, "total_weight"]), rowSums(e))
  pr <- PairIndex(6)
  s1 <- sapply(1:6, function(i) sum(e[1, pr[, 1] == i | pr[, 2] == i]))
  expect_equal(unname(m[1, "strength_a"]), mean(s1[1:3]))
  expect_equal(unname(m[1, "conn_a"]), mean(e[1, pr[, 1] <= 3 & pr[, 2] <= 3]))
})

test_that("ConcordanceFast equals Concordance, with ties", {
  set.seed(83)
  y <- c(rnorm(30), rnorm(5))
  y <- c(y, y[1:5])
  pred <- round(rnorm(40), 1)
  expect_equal(ConcordanceFast(pred, y), Concordance(pred, y))
  expect_equal(ConcordanceFast(rep(1, 10), rnorm(10)), 0.5)
})

test_that("first-order terms are exact and fourth-order terms nest them", {
  set.seed(84)
  x <- matrix(rnorm(80 * 9), 80) + rnorm(80)
  sets <- list(a = 1:4, b = 5:9)
  t4 <- PathwayTerms4(x, sets)
  m <- PathwayMetrics(LionessEdgesFirstOrder(x), PathwayBlocks(sets, 9))
  expect_gt(R2OnTerms(m[, "within_a"], t4$within_a[, 1:3])$r2, 1 - 1e-10)
  expect_equal(ncol(t4$within_a), 7)
  expect_equal(ncol(t4$cross_a_b), 7)
})

test_that("StudyReplicate returns the S1 and S2 quantities", {
  d <- StudyDesign()
  row <- d[d$part == "S2", ][2, ]
  row$people <- 30L
  row$visits <- 3L
  r <- StudyReplicate(row, 1, n_boot = 20L)
  expect_equal(nrow(r), 1)
  expect_true(all(c("icc_raw_betweenness_ifn", "icc_fo_transitivity",
                    "r2_4_conn_ifn", "icc_res4_conn_ifn", "finger_A",
                    "same_c_A", "same_eq_A_B", "next_sup_R0",
                    "next_d_A_C", "finger_G", "lsame_d_N1_E1",
                    "lnext_eq_N1_F1", "lnext_c_R0") %in% names(r)))
  expect_gt(r$r2_1_conn_ifn, 0.99)
})

test_that("strength and total weight of first-order edges equal their closed forms", {
  set.seed(85)
  x <- matrix(rnorm(60 * 6), 60) + rnorm(60)
  e1 <- LionessEdgesFirstOrder(x)
  m <- LionessMetrics(e1, list(a = 1:3, b = 4:6))
  z <- StandardizeN(x)$z
  r <- cor(x)
  diag(r) <- 0
  s <- rowSums(z)
  rs <- rowSums(r)
  strength1 <- rs[1] + z[, 1] * (s - z[, 1]) - rs[1] * z[, 1]^2 / 2 -
    drop(z^2 %*% r[, 1]) / 2
  pr <- PairIndex(6)
  expect_equal(rowSums(e1[, pr[, 1] == 1 | pr[, 2] == 1]), unname(strength1))
  total <- sum(r) / 2 + (s^2 - rowSums(z^2)) / 2 - drop(z^2 %*% rs) / 2
  expect_equal(unname(m[, "total_weight"]), total)
})

test_that("the ladder fits every model and the tests read the contrasts", {
  set.seed(86)
  person <- rep(1:30, each = 3)
  genes <- matrix(rnorm(90 * 6), 90)
  y <- genes[, 1] + rnorm(90, sd = 0.5)
  r1 <- cbind(genes[, 2], genes[, 2]^2)
  m <- matrix(rnorm(90 * 3), 90)
  outer <- GroupFolds(person, 5L)
  lad <- LadderPredictions(genes, r1, m, m, y, person, outer)
  expect_named(lad, c("R1", "E1", "N1", "F1"))
  lad_c <- LadderPredictions(genes, r1, m, m, y, person, outer,
                             current = matrix(rnorm(90), ncol = 1))
  expect_true("R0" %in% names(lad_c))
  t <- BootstrapLadder(lad, y, person, n_boot = 50L)
  expect_true(all(c("c_E1", "sup_E1", "d_N1_E1", "supd_N1_E1", "eq_N1_F1",
                    "d_E1_R1") %in% names(t)))
  # the genes carry the outcome, so E1 beats R1 and its concordance is high
  expect_gt(t[["c_E1"]], 0.7)
  expect_equal(t[["supd_E1_R1"]], 1)
  # identical predictions are equivalent
  same <- BootstrapLadder(list(N1 = y, F1 = y), y, person,
                          contrasts = list(c("N1", "F1")), n_boot = 20L)
  expect_equal(same[["eq_N1_F1"]], 1)
})
