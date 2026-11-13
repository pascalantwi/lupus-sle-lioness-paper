# pilot-c/pilot.R
# The pre-registered third pilot (pilot-c/preregistration.md). Run from the project
# root:
#   Rscript pilot-c/pilot.R --timing        ten replicates per cell, timing only
#   Rscript pilot-c/pilot.R --reps <n>      the run
# Tasks are the 10 cells, on three PSOCK workers with the BLAS pinned to one
# thread, scheduled by parLapplyLB(chunk.size = 1). Each replicate's seed is
# the cell seed plus the replicate index (PilotCDesign()).

devtools::load_all(here::here(), quiet = TRUE)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
timing <- "--timing" %in% args
n_rep <- if (timing) 10L else as.integer(args[match("--reps", args) + 1])
stopifnot(!is.na(n_rep), n_rep >= 1)
n_workers <- 3L
out_dir <- here::here("pilot-c", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

design <- PilotCDesign()
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
cl <- makePSOCKcluster(n_workers)
clusterExport(cl, "n_rep")
clusterEvalQ(cl, {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
  devtools::load_all(here::here(), quiet = TRUE)
  NULL
})
start <- Sys.time()
# the largest cells first, so the load balances
# activity cells fit twice the samples; aim 3 cells fit no models
order_cells <- order(-(design$women * ifelse(design$aim == 2, 2, 1) *
                       ifelse(design$aim == 3, 0.3, 1)))
res <- parLapplyLB(cl, split(design, design$cell)[order_cells], function(row) {
  t0 <- Sys.time()
  reps <- do.call(rbind, lapply(seq_len(n_rep), function(r) {
    suppressWarnings(RunReplicateB(row, r))
  }))
  reps$seconds_cell <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  reps
}, chunk.size = 1)
stopCluster(cl)
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
# cells of different aims have different columns
reps <- Reduce(function(a, b) merge(a, b, all = TRUE, sort = FALSE), res)
reps <- reps[order(reps$cell, reps$rep), ]

if (timing) {
  per_cell <- unique(reps[, c("cell", "seconds_cell")])
  per_cell$seconds_per_rep <- per_cell$seconds_cell / n_rep
  saveRDS(list(per_cell = per_cell, elapsed = elapsed, n_rep = n_rep,
               n_workers = n_workers, date = Sys.time()),
          file.path(out_dir, "timing.rds"))
  print(per_cell)
  cat(sprintf("elapsed %.0f s for %d replicates per cell on %d workers\n",
              elapsed, n_rep, n_workers))
} else {
  saveRDS(reps, file.path(out_dir, "replicates.rds"))
  summaries <- SummarizeCellsB(reps)
  saveRDS(summaries, file.path(out_dir, "cell_summaries.rds"))
  saveRDS(KnownAnswerChecksB(summaries),
          file.path(out_dir, "known_answers.rds"))
  saveRDS(design, file.path(out_dir, "design.rds"))
  saveRDS(list(run = "pilot", n_rep = n_rep, n_cells = nrow(design),
               date = Sys.time(), elapsed_seconds = elapsed,
               n_workers = n_workers,
               note = "pre-registered third pilot, pilot-c/preregistration.md",
               n_rep_production = NA_integer_),
          file.path(out_dir, "run_info.rds"))
  writeLines(capture.output(sessionInfo()),
             file.path(out_dir, "sessionInfo.txt"))
  cat(sprintf("done: %d cells x %d replicates in %.0f s\n", nrow(design),
              n_rep, elapsed))
}
