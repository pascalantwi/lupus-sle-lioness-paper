#' The pilot's design table
#'
#' The 8 cells of `pilot/preregistration.md`, section 3, with the seed of
#' each cell (`20761113 + 1000c`).
#'
#' @return A data frame with one row per cell: `cell`, `change` (what
#'   pregnancy changes), `factor` (`"none"`, `"moderate"`, `"large"`),
#'   `module` (the change in the module's correlation), `women`, and `seed`.
#' @export
PilotDesign <- function() {
  d <- data.frame(
    cell = 1:8,
    change = c("null", "factor", "factor", "factor", "module", "module",
               "module", "factor and module"),
    factor = c("none", "large", "moderate", "large", "none", "none", "none",
               "large"),
    module = c(0, 0, 0, 0, 0.4, 0.2, 0.4, 0.4),
    women = c(20, 20, 20, 40, 20, 20, 40, 40),
    stringsAsFactors = FALSE
  )
  d$seed <- 20761113L + 1000L * d$cell
  d
}

#' Simulate one longitudinal pregnancy cohort
#'
#' Each woman has four samples, two pregnant (`y = 1`) and two not. The
#' shared factor `f` (composition and IFN changes of pregnancy) is the
#' woman's level plus sample noise, each of variance 0.5; in pregnant samples
#' it becomes `shift + scale * f` (`"moderate"`: 0.5 and 1.5; `"large"`: 1
#' and 2.25). Composition `comp` is the woman's level (variance 0.8) plus
#' sample noise (variance 0.2), the same in both periods. Genes 1 to 20 load
#' 0.8 on `f`, genes 31 to 50 load 0.4 on `comp`; every gene carries the
#' woman's offset (variance 0.5) and unit noise; the noise of genes 21 to 30
#' has correlation 0.3, plus `module` in pregnant samples.
#'
#' @param women Number of women.
#' @param factor `"none"`, `"moderate"`, or `"large"`.
#' @param module Change in the module's noise correlation in pregnancy.
#' @param per_period Samples per woman in each period.
#' @param p Number of genes (at least 50).
#' @return A list with `x` (samples by genes), `y` (1 pregnant, 0 not),
#'   `woman` (integer identifier), `comp`, and `f`.
#' @export
SimulatePregnancyCohort <- function(women, factor = c("none", "moderate",
                                                      "large"),
                                    module = 0, per_period = 2L, p = 50) {
  factor <- match.arg(factor)
  stopifnot(p >= 50, module >= 0, 0.3 + module < 1)
  k <- 2L * per_period
  n <- women * k
  woman <- rep(seq_len(women), each = k)
  y <- rep(rep(c(1L, 0L), each = per_period), times = women)
  f <- stats::rnorm(women, sd = sqrt(0.5))[woman] +
    stats::rnorm(n, sd = sqrt(0.5))
  sk <- switch(factor, none = c(0, 1), moderate = c(0.5, 1.5),
               large = c(1, 2.25))
  f[y == 1] <- sk[1] + sk[2] * f[y == 1]
  comp <- stats::rnorm(women, sd = sqrt(0.8))[woman] +
    stats::rnorm(n, sd = sqrt(0.2))
  offset <- matrix(stats::rnorm(women * p, sd = sqrt(0.5)), women, p)[woman, ]
  eps <- matrix(stats::rnorm(n * p), n, p)
  mod <- 21:30
  for (cl in c(0L, 1L)) {
    r <- 0.3 + if (cl == 1L) module else 0
    sigma <- matrix(r, length(mod), length(mod))
    diag(sigma) <- 1
    rows <- y == cl
    eps[rows, mod] <- eps[rows, mod, drop = FALSE] %*% chol(sigma)
  }
  load_f <- c(rep(0.8, 20), rep(0, p - 20))
  load_c <- c(rep(0, 30), rep(0.4, 20), rep(0, p - 50))
  x <- outer(f, load_f) + outer(comp, load_c) + offset + eps
  colnames(x) <- paste0("g", seq_len(p))
  list(x = x, y = y, woman = woman, comp = comp, f = f)
}
