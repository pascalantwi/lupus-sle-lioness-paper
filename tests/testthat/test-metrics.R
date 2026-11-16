test_that("PilotEDesign has the 4 pre-registered cells", {
  d <- PilotEDesign()
  expect_equal(nrow(d), 4)
  expect_equal(d$seed, 21161116L + 1000L * (1:4))
})

test_that("EdgeMatrix places each pair symmetrically", {
  m <- EdgeMatrix(c(1, 2, 3), 3)
  expect_equal(m[1, 2], 1)
  expect_equal(m[3, 1], 2)
  expect_equal(m[3, 2], 3)
  expect_equal(diag(m), c(0, 0, 0))
})

test_that("NetworkMetrics matches hand computation on three nodes", {
  nm <- NetworkMetrics(c(1, 0.4, 1), 3, list(a = 1:2, b = 3),
                       threshold = 0.9)
  expect_equal(nm$global[["pathlen"]], 4 / 3)
  expect_equal(nm$global[["efficiency"]], 5 / 6)
  expect_equal(nm$global[["density"]], 2 / 3)
  expect_equal(unname(nm$node[2, "betweenness"]), 1)
  expect_equal(unname(nm$node[, "strength"]), c(1.4, 2, 1.4))
  expect_equal(nm$global[["conn_a"]], 1)
})

test_that("module connectivity is close to the squared module score", {
  set.seed(51)
  coh <- SimulateRepeatedCohort(100, 2, kappa_f = 0.55)
  g <- coh$x[, 6:50]
  e <- LionessEdges(g)
  conn <- apply(e, 1, function(r) NetworkMetrics(r, 45, NodeSetsE())$global[["conn_ifn"]])
  score <- rowMeans(StandardizeN(g[, 1:15])$z)
  expect_gt(cor(conn, score^2), 0.8)
})

test_that("HubOverlap is high for identical visits and near zero for
          independent ones", {
  set.seed(52)
  s <- matrix(rnorm(40 * 10), 40, 10)
  person <- rep(1:20, 2)
  visit <- rep(1:2, each = 20)
  same <- rbind(s[1:20, ], s[1:20, ])
  expect_gt(HubOverlap(same, person, visit), 0.5)
  expect_lt(abs(HubOverlap(s, person, visit)), 0.2)
})

test_that("MetricsReplicate returns stability for every metric and type", {
  d <- PilotEDesign()
  d$people <- 20L
  r <- MetricsReplicate(d[2, ], 1)
  expect_true(all(c("icc_raw_conn_ifn", "icc_fo_pathlen",
                    "icc_partial_betweenness_module", "closedform_cor",
                    "icc_score2", "hub_raw") %in% names(r)))
})

test_that("SummarizeCellsE gives the equivalence and fix quantities", {
  reps <- data.frame(cell = 1, rep = 1:2, icc_raw_pathlen = 0.3,
                     icc_fo_pathlen = 0.29, icc_partial_conn_ifn = 0.2,
                     icc_partial_conn_cross = 0.15)
  s <- SummarizeCellsE(reps)
  expect_equal(s$mean[s$quantity == "eq_pathlen"], 0.01)
  expect_equal(s$mean[s$quantity == "fix_ifn"], 0.05)
})

test_that("the known-answer checks of pilot-e hold", {
  k <- KnownAnswerChecksE()
  ids <- k[k$type == "identity", ]
  expect_true(all(abs(ids$computed - ids$known) < 0.02))
})

test_that("PilotFDesign uses the own-module correction and new seeds", {
  d <- PilotFDesign()
  expect_equal(d$fix, rep("ownset", 4))
  expect_equal(d$seed, 21261116L + 1000L * (1:4))
})

test_that("OwnSetResiduals are uncorrelated with their own regressors", {
  set.seed(61)
  g <- matrix(rnorm(200 * 6), 200, 6) + rnorm(200)
  cond <- cbind(rnorm(200))
  r <- OwnSetResiduals(g, list(1:3, 4:6), cond)
  z <- StandardizeN(g)$z
  loo1 <- (rowSums(z[, 1:3]) - z[, 1]) / 2
  expect_lt(abs(cor(r[, 1], loo1)), 1e-10)
  expect_lt(abs(cor(r[, 5], cond[, 1])), 1e-10)
})

test_that("MetricsReplicate computes the own-module edge type", {
  d <- PilotFDesign()
  d$people <- 20L
  r <- MetricsReplicate(d[2, ], 1)
  expect_true(all(c("icc_ownset_conn_ifn", "icc_fo_pathlen") %in% names(r)))
  expect_false("icc_partial_conn_ifn" %in% names(r))
})
