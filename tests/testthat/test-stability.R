test_that("PilotDDesign has the 6 pre-registered cells", {
  d <- PilotDDesign()
  expect_equal(nrow(d), 6)
  expect_equal(d$seed, 21061116L + 1000L * (1:6))
  expect_equal(d$cell[d$predict], c(1, 2, 6))
})

test_that("SimulateRepeatedCohort gives the factor stability asked for", {
  set.seed(41)
  coh <- SimulateRepeatedCohort(3000, 3, kappa_f = 0.55, kappa_c = 0.45)
  # absolute tolerance: the ICC from 3,000 people has a standard error
  # of about 0.015
  expect_lt(abs(EdgeICC(coh$f, coh$person) - 0.55), 0.05)
  expect_lt(abs(EdgeICC(coh$comp, coh$person) - 0.45), 0.05)
  expect_lt(abs(EdgeICC(coh$activity, coh$person) - 0.39), 0.05)
  expect_equal(dim(coh$x), c(9000, 50))
})

test_that("EdgeICC matches the analysis-of-variance formula", {
  set.seed(42)
  person <- rep(1:25, each = 3)
  v <- rnorm(25)[person] + rnorm(75)
  a <- anova(lm(v ~ factor(person)))
  expect_equal(EdgeICC(v, person),
               (a[1, 3] - a[2, 3]) / (a[1, 3] + 2 * a[2, 3]))
})

test_that("PairGroupsD has the stated group sizes", {
  pg <- PairGroupsD()
  expect_equal(lengths(pg), c(ifn = 105, module = 45, cellmix = 190,
                              cross = 150))
})

test_that("FingerprintAccuracy is 1 for identical visits", {
  set.seed(43)
  e <- matrix(rnorm(20 * 10), 20, 10)
  expect_equal(FingerprintAccuracy(rbind(e, e), rep(1:20, 2),
                                   rep(1:2, each = 20), 1:10), 1)
})

test_that("raw edges of IFN genes are more stable than cross pairs when IFN is
          stable, and partial edges remove it", {
  set.seed(44)
  coh <- SimulateRepeatedCohort(131, 6, kappa_f = 0.9)
  g <- coh$x[, 6:50]
  pg <- PairGroupsD()
  score <- rowMeans(StandardizeN(coh$x[, 1:5])$z)
  raw <- EdgeICC(LionessEdges(g), coh$person)
  par <- EdgeICC(LionessPartialEdges(g, cbind(score, coh$comp)), coh$person)
  expect_gt(mean(raw[pg$ifn]) - mean(raw[pg$cross]), 0.1)
  expect_lt(abs(mean(par[pg$ifn]) - mean(par[pg$cross])),
            mean(raw[pg$ifn]) - mean(raw[pg$cross]))
})

test_that("StabilityReplicate returns the measures, with predictions in
          prediction cells", {
  d <- PilotDDesign()
  d$people <- 40L
  r1 <- suppressWarnings(StabilityReplicate(d[1, ], 1))
  expect_true(all(c("icc_raw_ifn", "fp_partial_all", "N1", "E1") %in%
                    names(r1)))
  r4 <- StabilityReplicate(d[4, ], 1)
  expect_false("N1" %in% names(r4))
})

test_that("SummarizeCellsD computes the excess over cross pairs", {
  reps <- data.frame(cell = 1, rep = 1:2, icc_raw_ifn = c(0.3, 0.3),
                     icc_raw_cross = c(0.1, 0.1))
  s <- SummarizeCellsD(reps)
  expect_equal(s$mean[s$quantity == "x_raw_ifn"], 0.2)
})

test_that("the known-answer checks of pilot-d hold", {
  k <- KnownAnswerChecksD()
  ids <- k[k$type == "identity", ]
  expect_true(all(abs(ids$computed - ids$known) < 0.02))
})

test_that("DescribeGse65391 computes the within-child stability", {
  set.seed(45)
  n <- 40
  subj <- rep(sprintf("S%02d", 1:10), each = 4)
  lev <- rnorm(10)[as.integer(factor(subj))]
  tab <- rbind(matrix(rnorm(5 * n), 5, n) + rep(lev, each = 5),
               matrix(rnorm(2 * n), 2, n))
  geo <- sprintf("GSM%03d", 1:n)
  dimnames(tab) <- list(c(paste0("p", 1:5), "q1", "q2"), geo)
  samples <- data.frame(geo = geo, disease.state = "SLE", subject = subj,
                        sledai = rpois(n, 4), neutrophil_percent = 60,
                        lymphocyte_percent = 30,
                        cumulative_time = rep(c(0, 90, 180, 270), 10),
                        set = rep(c("Training", "Test"), each = 20))
  samples$neutrophil_percent <- samples$neutrophil_percent + rnorm(n)
  samples$lymphocyte_percent <- samples$lymphocyte_percent + rnorm(n)
  platform <- data.frame(id = c(paste0("p", 1:5), "q1", "q2"),
                         symbol = c("EPSTI1", "IFI44L", "LY6E", "OAS3",
                                    "RSAD2", "A", "B"))
  d <- DescribeGse65391(list(samples = samples, table = tab), platform)
  expect_equal(unname(d$counts["children_multi"]), 10)
  expect_gt(d$icc["ifn"], 0.5)
  expect_equal(unname(d$gap_days["Median"]), 90)
})
