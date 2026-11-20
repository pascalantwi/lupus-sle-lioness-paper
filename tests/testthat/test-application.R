test_that("one-sample-at-a-time metrics equal the all-at-once metrics", {
  set.seed(101)
  x <- matrix(rnorm(40 * 8), 40) + rnorm(40)
  sets <- list(a = 1:4, b = 5:8)
  expect_equal(LionessMetricsBySample(x, sets),
               LionessMetrics(LionessEdges(x), sets), ignore_attr = TRUE,
               tolerance = 1e-10)
  expect_equal(LionessMetricsBySample(x, sets, first_order = TRUE),
               LionessMetrics(LionessEdgesFirstOrder(x), sets),
               ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("SampleEdges is one row of the closed form", {
  set.seed(102)
  x <- matrix(rnorm(30 * 5), 30)
  z <- StandardizeN(x)$z
  pr <- PairIndex(5)
  e <- SampleEdges(z[3, ], cor(x)[pr], pr, 30)
  expect_equal(e, unname(LionessEdgesDirect(x)[3, ]), tolerance = 1e-10)
})

test_that("ExpressionByGene keeps the highest-mean probe per gene", {
  tab <- rbind(p1 = c(1, 2, 3), p2 = c(5, 6, 7), p3 = c(0, 0, 1),
               p4 = c(9, 9, 9))
  colnames(tab) <- c("s1", "s2", "s3")
  plat <- data.frame(id = c("p1", "p2", "p3", "p4"),
                     symbol = c("A", "A", "B", "C"))
  x <- ExpressionByGene(tab, plat, c("B", "A", "D"))
  expect_equal(colnames(x), c("B", "A"))
  expect_equal(unname(x[, "A"]), c(5, 6, 7))
  expect_equal(rownames(x), c("s1", "s2", "s3"))
})

test_that("BloodGen3Pathways picks by label and keeps a gene once", {
  bg3 <- list(modules = data.frame(
    module = c("M1", "M1", "M2", "M3", "M4"),
    gene = c("G1", "G2", "G2", "G3", "G4"),
    "function" = c("Interferon", "Interferon", "Neutrophils", "Lymphocytes",
                   "TBD"), check.names = FALSE))
  pw <- BloodGen3Pathways(bg3, c("G1", "G2", "G3", "G4"))
  expect_equal(pw$ifn, c("G1", "G2"))
  expect_equal(pw$neut, character(0))
  expect_equal(pw$lymph, "G3")
})

test_that("ModuleDiagnostics finds connectivity linear in its terms", {
  set.seed(103)
  person <- rep(1:40, each = 4)
  f <- rnorm(40)[person] + rnorm(160)
  x <- 0.8 * f + matrix(rnorm(160 * 10), 160)
  d <- ModuleDiagnostics(x, person)
  expect_gt(d[["r2_1"]], 0.99)
  expect_gte(d[["r2_4"]], d[["r2_1"]])
  expect_gt(d[["pc1_share"]], 0.1)
})

test_that("SplitConcordance returns one concordance per feature set", {
  set.seed(104)
  person <- rep(1:30, each = 3)
  y <- rnorm(90)
  feats <- list(A = matrix(rnorm(270), 90), B = cbind(y + rnorm(90), 1))
  cc <- SplitConcordance(feats, y, person, person <= 20)
  expect_named(cc, c("A", "B"))
  expect_gt(cc[["B"]], 0.6)
})

test_that("ReadBloodGen3 reads the module lists", {
  tb <- test_path("..", "..", "data", "raw", "BloodGen3Module_1.20.0.tar.gz")
  skip_if_not(file.exists(tb), "BloodGen3Module not downloaded")
  bg3 <- ReadBloodGen3(tb)
  expect_equal(nrow(bg3$modules), 14168)
  expect_equal(length(unique(bg3$modules$module)), 382)
})

test_that("LionessMetricsBySample gives the same answer on a cluster", {
  skip_on_cran()
  set.seed(105)
  x <- matrix(rnorm(30 * 6), 30) + rnorm(30)
  sets <- list(a = 1:3, b = 4:6)
  cl <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, devtools::load_all(
    testthat::test_path("..", ".."), quiet = TRUE))
  expect_equal(LionessMetricsBySample(x, sets, cl = cl),
               LionessMetricsBySample(x, sets))
})

test_that("PathwayDiagnostics reads the right metrics and terms", {
  set.seed(106)
  person <- rep(1:30, each = 4)
  f <- rnorm(30)[person] + rnorm(120)
  x <- cbind(0.8 * f + matrix(rnorm(120 * 5), 120),
             matrix(rnorm(120 * 4), 120))
  sets <- list(a = 1:5, b = 6:9)
  m <- LionessMetrics(LionessEdges(x), sets)
  d <- PathwayDiagnostics(m, x, sets, person)
  expect_equal(d$block, c("within_a", "within_b", "cross_a_b"))
  expect_true(all(d$r2_1 > 0.99))
  expect_true(all(d$r2_4 >= d$r2_1 - 1e-12))
  # single-edge stability from the pathway's genes alone equals the
  # stability of the same edges in the whole network
  pr <- PairIndex(9)
  inside <- pr[, 1] %in% 1:5 & pr[, 2] %in% 1:5
  expect_equal(d$icc_edges[1],
               mean(EdgeICC(LionessEdges(x)[, inside], person)),
               tolerance = 1e-10)
  expect_true(is.na(d$icc_edges[3]))
})
