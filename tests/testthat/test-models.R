test_that("the known-answer checks hold", {
  k <- KnownAnswerChecks()
  ids <- k[k$type == "identity", ]
  expect_true(all(abs(ids$computed - ids$known) < 0.02))
  exact <- grepl("equals", ids$check)
  expect_true(all(abs(ids$computed[exact] - ids$known[exact]) < 1e-10))
})

test_that("Auc matches the Mann-Whitney statistic", {
  set.seed(2)
  s <- rnorm(30)
  y <- rep(0:1, 15)
  w <- suppressWarnings(wilcox.test(s[y == 1], s[y == 0]))$statistic
  expect_equal(Auc(s, y), unname(w) / (15 * 15))
})

test_that("GroupFolds keeps groups together and uses k folds", {
  set.seed(3)
  g <- rep(1:12, each = 4)
  f <- GroupFolds(g, 5)
  expect_true(all(tapply(f, g, function(v) length(unique(v))) == 1))
  expect_setequal(unique(f), 1:5)
})

test_that("the simulated cohort has the stated structure", {
  set.seed(4)
  coh <- SimulatePregnancyCohort(10, "large", 0.4)
  expect_equal(dim(coh$x), c(40, 50))
  expect_equal(as.vector(table(coh$y)), c(20, 20))
  expect_true(all(tapply(coh$y, coh$woman, sum) == 2))
})

test_that("the factor shifts and the module rewires as specified", {
  set.seed(5)
  coh <- SimulatePregnancyCohort(4000, "large", 0.4, per_period = 1L)
  expect_equal(mean(coh$f[coh$y == 1]), 1, tolerance = 0.05)
  expect_equal(sd(coh$f[coh$y == 1]), 2.25, tolerance = 0.03)
  expect_equal(sd(coh$f[coh$y == 0]), 1, tolerance = 0.03)
  # within-period module correlation net of the woman offsets: 0.3 / 1.5 and
  # 0.7 / 1.5
  r0 <- cor(coh$x[coh$y == 0, 21:30])
  r1 <- cor(coh$x[coh$y == 1, 21:30])
  expect_equal(mean(r0[upper.tri(r0)]), 0.2, tolerance = 0.03)
  expect_equal(mean(r1[upper.tri(r1)]), 0.7 / 1.5, tolerance = 0.03)
})

test_that("PregnancyFeatures gives the stated dimensions", {
  set.seed(6)
  coh <- SimulatePregnancyCohort(10, "none", 0)
  tr <- GroupFolds(coh$woman, 5) != 1
  ft <- PregnancyFeatures(coh, tr)
  expect_equal(ncol(ft$train$edges), 990)
  expect_equal(ncol(ft$train$quad), 990 + 45)
  expect_equal(nrow(ft$test$fo), sum(!tr))
  expect_equal(colnames(ft$train$r1), c("comp", "ifn", "ifn2", "comp2"))
})

test_that("ApparentRewiring is near its level under no change", {
  set.seed(7)
  v <- replicate(20, {
    coh <- SimulatePregnancyCohort(20, "none", 0)
    ApparentRewiring(coh$x[, 21:30], coh$y, coh$woman)
  })
  expect_lt(abs(mean(v) - 0.05), 0.04)
})

test_that("ApparentRewiring detects the squared factor with no co-expression
          change", {
  set.seed(8)
  coh <- SimulatePregnancyCohort(40, "large", 0)
  expect_gt(ApparentRewiring(coh$x[, 6:20], coh$y, coh$woman), 0.5)
})

test_that("RunReplicate returns every model's AUC in [0, 1]", {
  set.seed(9)
  row <- PilotDesign()[5, ]
  r <- suppressWarnings(RunReplicate(row, 1))
  aucs <- unlist(r[c("R0", "R1", "E1", "E5", "F1", "N1", "N1R0")])
  expect_true(all(aucs >= 0 & aucs <= 1))
  expect_equal(r$cell, 5)
})

test_that("SummarizeCells gives the paired differences", {
  reps <- data.frame(cell = 1, rep = 1:3, R0 = 0.5, R1 = 0.5, E1 = 0.6,
                     E5 = 0.6, F1 = 0.6, N1 = c(0.7, 0.7, 0.7), N1R0 = 0.55,
                     rewire_factor = 0.05, rewire_module = 0.05)
  s <- SummarizeCells(reps)
  expect_equal(s$mean[s$quantity == "d1_N1_vs_E1"], 0.1)
  expect_equal(s$mean[s$quantity == "artifact_N1R0_vs_R0"], 0.05)
})

test_that("DescribePregnancyCohort counts usable pregnancies", {
  s <- data.frame(samplegroup = c("SLE", "SLE", "SLE", "RA", "RA"),
                  timepoint = c("1", "5", "2", "3", "2"),
                  pregnancyid = c("a", "a", "b", "c", "c"))
  d <- DescribePregnancyCohort(s)$by_group
  expect_equal(d$usable_pregnancies[d$group == "SLE"], 1)
  expect_equal(d$usable_pregnancies[d$group == "RA"], 0)
  expect_equal(d$pregnant_samples[d$group == "SLE"], 1)
})
