test_that("PilotGDesign has the 3 pre-registered cells", {
  d <- PilotGDesign()
  expect_equal(nrow(d), 3)
  expect_equal(d$seed, 21361116L + 1000L * (1:3))
  expect_equal(d$kappa_ifn, c(0, 0.55, 0.9))
  expect_equal(d$tau, c(0, 0.25, 0.25))
})

test_that("PathwaySetsG gives six disjoint sets of 15 genes", {
  s <- PathwaySetsG()
  expect_equal(lengths(s), rep(15L, 6), ignore_attr = TRUE)
  expect_equal(sort(unlist(s, use.names = FALSE)), 1:90)
})

test_that("PathwayBlocks counts the pairs of each block", {
  b <- PathwayBlocks(list(a = 1:3, b = 4:5), 5)
  expect_length(b$within_a, 3)
  expect_length(b$within_b, 1)
  expect_length(b$cross_a_b, 6)
  pr <- PairIndex(5)
  expect_true(all(pr[b$within_b, ] == c(4, 5)))
})

test_that("SimulatePathwayCohort gives the factor's within-person stability", {
  set.seed(61)
  coh <- SimulatePathwayCohort(300, 4, 0.9, 0, 0, 0, 0)
  sets <- PathwaySetsG()
  ifn <- rowMeans(coh$x[, sets$ifn])
  null <- rowMeans(coh$x[, sets$null])
  # the IFN score carries the person level and the stable factor
  expect_gt(EdgeICC(ifn, coh$person), EdgeICC(null, coh$person) + 0.2)
  r_null <- cor(coh$x[, sets$null])
  expect_lt(abs(mean(r_null[upper.tri(r_null)]) - 0.3 / 1.5), 0.05)
})

test_that("first-order connectivity and crosstalk equal their closed forms", {
  k <- KnownAnswerChecksG(seed = 3L)
  cf <- k[grepl("closed form", k$check), ]
  expect_equal(nrow(cf), 2)
  expect_true(all(cf$computed < 1e-10))
})

test_that("R2OnTerms is 1 for an exact linear function", {
  set.seed(62)
  t <- matrix(rnorm(60), 20, 3)
  y <- 1 + t %*% c(2, -1, 0.5)
  f <- R2OnTerms(drop(y), t)
  expect_equal(f$r2, 1)
  expect_lt(max(abs(f$resid)), 1e-10)
})

test_that("PathwayReplicate returns the pre-registered quantities", {
  row <- PilotGDesign()[2, ]
  row$people <- 30L
  row$visits <- 3L
  r <- PathwayReplicate(row, 1)
  expect_equal(nrow(r), 1)
  expect_true(all(c("min_r2_within", "min_r2_cross", "cor_ifn",
                    "icc_resid_personal", "icc_edges_ifn", "min_r2_fo") %in%
                    names(r)))
  expect_gt(r$min_r2_fo, 1 - 1e-8)
  s <- SummarizeCellsG(rbind(r, PathwayReplicate(row, 2)))
  expect_true(all(c("expl_ifn", "agg_ifn", "personal_vs_null") %in%
                    s$quantity))
})
