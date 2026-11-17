# Acceptance test: the package's functions reproduce pilot-g's saved results
# with pilot-g's grid and seeds. Slow, so gated by an environment variable:
#   LIONESSPREG_ACCEPTANCE=quick  replicates 1 to 3 of every cell
#   LIONESSPREG_ACCEPTANCE=full   every replicate and every per-cell summary
test_that("the pathway functions reproduce pilot-g", {
  mode <- Sys.getenv("LIONESSPREG_ACCEPTANCE")
  skip_if(!mode %in% c("quick", "full"), "set LIONESSPREG_ACCEPTANCE")
  saved <- readRDS(test_path("fixtures", "pilot_g_replicates.rds"))
  design <- PilotGDesign()
  n_rep <- if (mode == "full") max(saved$rep) else 3L
  reps <- do.call(rbind, lapply(design$cell, function(cl) {
    do.call(rbind, lapply(seq_len(n_rep), function(r) {
      PathwayReplicate(design[design$cell == cl, ], r)
    }))
  }))
  keep <- saved$rep <= n_rep
  cols <- setdiff(names(reps), "seconds_cell")
  expect_equal(reps[, cols], saved[keep, cols], ignore_attr = TRUE,
               tolerance = 1e-12)
  if (mode == "full") {
    expect_equal(SummarizeCellsG(reps),
                 readRDS(test_path("fixtures", "pilot_g_cell_summaries.rds")),
                 ignore_attr = TRUE, tolerance = 1e-12)
  }
})
