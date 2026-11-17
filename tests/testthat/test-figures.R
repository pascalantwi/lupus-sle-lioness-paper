test_that("PathwayLayout places every gene and the pathway centres", {
  xy <- PathwayLayout(list(a = 1:3, b = 4:7))
  expect_equal(dim(xy), c(7L, 2L))
  expect_false(anyNA(xy))
  expect_equal(rownames(attr(xy, "centers")), c("a", "b"))
  # each gene sits at the inner radius from its pathway's centre
  d <- sqrt(rowSums((xy[1:3, ] -
                       matrix(attr(xy, "centers")["a", ], 3, 2,
                              byrow = TRUE))^2))
  expect_equal(d, rep(0.28, 3))
})

test_that("each figure function draws one plot without error", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(71)
  sets <- list(ifn = 1:4, null = 5:8)
  e <- rnorm(8 * 7 / 2)
  expect_length(PlotSampleNetwork(e, sets, top = 10), 10)
  blocks <- PathwayBlocks(sets, 8)
  m <- PathwayMetrics(matrix(e, 1), blocks)[1, ]
  expect_null(PlotPathwayNetwork(m[1:2], m[3]))
  expect_gt(PlotConnectivityScore(1:10 + rnorm(10, sd = 0.1), 1:10), 0.9)
  expect_null(PlotTrajectories(rnorm(6), rep(1:2, 3), rep(1:3, each = 2)))
  s <- expand.grid(cell = 2, kind = c("conn", "score2", "edges", "resid"),
                   pw = names(PathwayLabels()), stringsAsFactors = FALSE)
  s$quantity <- sprintf("icc_%s_%s", s$kind, s$pw)
  s$mean <- runif(nrow(s))
  s$mcse <- 0.01
  expect_equal(dim(PlotStabilityByPathway(s, 2)), c(6L, 4L))
})
