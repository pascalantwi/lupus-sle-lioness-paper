test_that("partial correlations equal the inverse-covariance formula", {
  set.seed(51)
  x <- matrix(rnorm(40 * 3), 40, 3)
  cond <- rnorm(40)
  x[, 1:2] <- x[, 1:2] + cond
  lit <- apply(PairIndex(3), 1, function(ij) {
    pm <- solve(cov(cbind(x[, ij], cond)))
    -pm[1, 2] / sqrt(pm[1, 1] * pm[2, 2])
  })
  expect_equal(PartialCorrelationPairs(x, cond), lit, tolerance = 1e-10)
  expect_equal(PartialCorrelationPairs(x), cor(x)[PairIndex(3)],
               tolerance = 1e-12)
})

test_that("partial-correlation LIONESS edges follow their definition", {
  set.seed(52)
  x <- matrix(rnorm(20 * 3), 20, 3)
  cond <- cbind(rnorm(20), rnorm(20))
  rho <- PartialCorrelationPairs(x, cond)
  lit <- t(sapply(1:20, function(q) {
    20 * rho - 19 * PartialCorrelationPairs(x[-q, ], cond[-q, ])
  }))
  expect_equal(LionessPartialEdges(x, cond), lit, tolerance = 1e-10)
  expect_equal(LionessPartialEdges(x), LionessEdges(x), tolerance = 1e-10)
  new <- matrix(rnorm(2 * 3), 2, 3)
  c_new <- cbind(rnorm(2), rnorm(2))
  lit_new <- t(sapply(1:2, function(q) {
    21 * PartialCorrelationPairs(rbind(x, new[q, ]), rbind(cond, c_new[q, ])) -
      20 * rho
  }))
  expect_equal(LionessPartialTestEdges(new, c_new, x, cond), lit_new,
               tolerance = 1e-10)
  expect_equal(LionessPartialTestEdges(new, NULL, x, NULL),
               LionessTestEdges(new, x), tolerance = 1e-10)
})
