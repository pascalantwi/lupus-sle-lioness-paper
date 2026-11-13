test_that("PilotBDesign has the 10 pre-registered cells", {
  d <- PilotBDesign()
  expect_equal(nrow(d), 10)
  expect_equal(d$seed, 20861113L + 1000L * (1:10))
  expect_equal(sum(d$aim == 3), 3)
})

test_that("SimulateCohortB moves the factors as specified", {
  set.seed(21)
  coh <- SimulateCohortB(3000, FactorLevels("RA"))
  expect_equal(mean(coh$comp[coh$y == 1]), 0.45, tolerance = 0.05)
  expect_equal(sd(coh$comp[coh$y == 1]) / sd(coh$comp[coh$y == 0]), 1.45,
               tolerance = 0.05)
  expect_equal(mean(coh$f[coh$y == 1]), -0.3, tolerance = 0.05)
})

test_that("the module correlation follows the pregnancy change", {
  set.seed(22)
  coh <- SimulateCohortB(3000, module = 0.4)
  r1 <- cor(coh$x[coh$y == 1, 21:30])
  r0 <- cor(coh$x[coh$y == 0, 21:30])
  expect_equal(mean(r1[upper.tri(r1)]), 0.7 / 1.5, tolerance = 0.03)
  expect_equal(mean(r0[upper.tri(r0)]), 0.3 / 1.5, tolerance = 0.03)
})

test_that("PairGroups has 295 factor pairs and 45 module pairs", {
  pg <- PairGroups()
  expect_equal(length(pg$factor), 105 + 190)
  expect_equal(length(pg$module), 45)
  expect_length(intersect(pg$factor, pg$module), 0)
})

test_that("Concordance equals Auc for a binary outcome", {
  set.seed(23)
  s <- rnorm(20)
  y <- rep(0:1, 10)
  expect_equal(Concordance(s, y), Auc(s, y))
})

test_that("PairedShare and BetweenShare are near the level under no change",
          {
  set.seed(24)
  v <- replicate(15, {
    a <- SimulateCohortB(18)
    b <- SimulateCohortB(27)
    da <- ChangeByWoman(LionessEdges(a$x[, 6:50]), a$y, a$woman)
    db <- ChangeByWoman(LionessEdges(b$x[, 6:50]), b$y, b$woman)
    c(PairedShare(da, PairGroups()$module),
      BetweenShare(da, db, PairGroups()$module))
  })
  expect_lt(max(abs(rowMeans(v) - 0.05)), 0.04)
})

test_that("partial edges remove the false rewiring of the stress level", {
  set.seed(25)
  coh <- SimulateCohortB(30, FactorLevels("stress"))
  ed <- CohortEdges(coh)
  pg <- PairGroups()
  raw <- PairedShare(ChangeByWoman(ed$raw, coh$y, coh$woman), pg$factor)
  par <- PairedShare(ChangeByWoman(ed$partial, coh$y, coh$woman), pg$factor)
  expect_gt(raw, 0.15)
  expect_lt(par, raw)
})

test_that("RunReplicateB returns the measures of each aim", {
  d <- PilotBDesign()
  r1 <- suppressWarnings(RunReplicateB(d[1, ], 1))
  expect_true(all(c("paired_par_factor", "paired_raw_module", "N1") %in%
                    names(r1)))
  r6 <- suppressWarnings(RunReplicateB(d[6, ], 1))
  expect_true(all(c("N1", "cor_N1", "r2_N1", "rmse_N1") %in% names(r6)))
  r8 <- RunReplicateB(d[8, ], 1)
  expect_true("between_par_factor" %in% names(r8))
})

test_that("the known-answer checks of pilot-b hold", {
  k <- KnownAnswerChecksB()
  ids <- k[k$type == "identity", ]
  expect_true(all(abs(ids$computed - ids$known) < 0.02))
})

test_that("pilot-b's cohort is unchanged by the new visit arguments", {
  set.seed(31)
  a <- SimulateCohortB(5, FactorLevels("RA"))
  set.seed(31)
  b <- SimulateCohortB(5, FactorLevels("RA"), n_preg = 2L, n_non = 2L)
  expect_identical(a, b)
})

test_that("SimulateCohortB gives the requested visits per woman", {
  set.seed(32)
  coh <- SimulateCohortB(6, FactorLevels("SLE_c"), n_preg = 3L, n_non = 2L)
  expect_equal(nrow(coh$x), 30)
  expect_true(all(tapply(coh$y, coh$woman, sum) == 3))
})

test_that("CenterByGroup removes each group's mean", {
  m <- cbind(c(1, 3, 10, 14), c(2, 2, 5, 7))
  g <- c(2, 2, 1, 1)
  cm <- CenterByGroup(m, g)
  expect_equal(cm[, 1], c(-1, 1, -2, 2))
  expect_equal(CenterByGroup(c(1, 3, 10, 14), g), c(-1, 1, -2, 2))
})

test_that("PilotCDesign has five visits and the within-woman activity cells",
          {
  d <- PilotCDesign()
  expect_equal(d$n_preg + d$n_non, rep(5L, 10))
  expect_equal(d$within, d$aim == 2)
  expect_equal(d$seed, 20961113L + 1000L * (1:10))
  r <- suppressWarnings(RunReplicateB(d[6, ], 1))
  expect_true(all(c("N1", "r2_N1") %in% names(r)))
})
