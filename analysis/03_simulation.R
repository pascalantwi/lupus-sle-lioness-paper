# analysis/03_simulation.R
# The simulation study of protocol.md: S1 (13 settings, describe and explain
# the LIONESS metrics over visits) and S2 (5 settings, prediction). Run from
# the project root:
#   Rscript analysis/03_simulation.R --calibrate          protocol.md, section 6
#   Rscript analysis/03_simulation.R --run quick --reps 10
#   Rscript analysis/03_simulation.R --run production
# (01_ and 02_ are the earlier pregnancy scripts, kept as the record.)
# Each replicate's seed is its setting's seed plus the replicate index
# (StudyDesign()); replicates are the tasks of parLapplyLB(chunk.size = 1) on
# PSOCK workers with the BLAS pinned to one thread.
# Writes results/: replicates.rds, cell_summaries.rds, known_answers.rds,
# design.rds, example_networks.rds, run_info.rds, sessionInfo.txt; with
# --calibrate, calibration.rds and calibration_run_info.rds.

devtools::load_all(here::here(), quiet = TRUE)
library(parallel)

# ---- parameter block (protocol.md, sections 2 and 6) -------------------------
beta <- c(0, 0.15, 0.3, 0.45, 0.6)  # S2 activity effects; set by calibration
n_rep_production <- 100L
n_workers <- 3L                     # the core cap
n_boot <- 200L                      # bootstrap resamples of the tests
calibration_beta <- c(0.05, 0.1, 0.15, 0.2, 0.25)
calibration_reps <- 12L
example_cell <- 1L                  # the base setting, for the figures
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
calibrate <- "--calibrate" %in% args
run <- if (calibrate) "quick" else opt("--run", "quick")
stopifnot(run %in% c("quick", "production"))
n_rep <- if (calibrate) calibration_reps else if (run == "production")
  n_rep_production else as.integer(opt("--reps", 10))
design <- if (calibrate) {
  d <- StudyDesign(beta = calibration_beta)
  d <- d[d$part == "S2", ]
  # calibration settings get their own seeds, apart from the study's
  d$seed <- d$seed + 500000L
  d
} else {
  StudyDesign(beta = beta)
}
set.seed(21461117L)  # no random step outside the replicates
out_dir <- here::here("results")
dir.create(out_dir, showWarnings = FALSE)

tasks <- expand.grid(rep = seq_len(n_rep), cell = design$cell)
# the S2 settings take about twice as long; start them first
tasks <- tasks[order(-(design$part[match(tasks$cell, design$cell)] == "S2"),
                     tasks$cell, tasks$rep), ]

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
cl <- makePSOCKcluster(n_workers)
clusterExport(cl, c("design", "n_boot"))
clusterEvalQ(cl, {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
  devtools::load_all(here::here(), quiet = TRUE)
  NULL
})
start <- Sys.time()
res <- parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), function(tk) {
  suppressWarnings(StudyReplicate(design[design$cell == tk$cell, ], tk$rep,
                                  n_boot = n_boot))
}, chunk.size = 1)
stopCluster(cl)
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

cols <- unique(unlist(lapply(res, names)))
reps <- do.call(rbind, lapply(res, function(r) {
  r[setdiff(cols, names(r))] <- NA_real_
  r[cols]
}))
reps <- reps[order(reps$cell, reps$rep), ]
rownames(reps) <- NULL
summaries <- SummarizeStudy(reps)
info <- list(run = run, n_rep = n_rep, n_cells = nrow(design),
             date = Sys.time(), elapsed_seconds = elapsed,
             n_workers = n_workers, n_rep_production = n_rep_production)

if (calibrate) {
  saveRDS(list(design = design, summaries = summaries, reps = reps),
          file.path(out_dir, "calibration.rds"))
  info$note <- "calibration of the S2 activity effects, protocol.md section 6"
  saveRDS(info, file.path(out_dir, "calibration_run_info.rds"))
  pw <- summaries[summaries$quantity %in% c("same_sup_A", "same_sup_B",
                                            "same_sup_C"), ]
  pw$beta <- design$beta[match(pw$cell, design$cell)]
  print(pw[, c("beta", "quantity", "mean", "mcse")], digits = 3)
} else {
  saveRDS(reps, file.path(out_dir, "replicates.rds"))
  saveRDS(summaries, file.path(out_dir, "cell_summaries.rds"))
  saveRDS(KnownAnswerChecksG(), file.path(out_dir, "known_answers.rds"))
  saveRDS(design, file.path(out_dir, "design.rds"))
  saveRDS(ExampleNetworks(design[design$cell == example_cell, ], 1L),
          file.path(out_dir, "example_networks.rds"))
  info$note <- "the simulation study of protocol.md"
  saveRDS(info, file.path(out_dir, "run_info.rds"))
  writeLines(capture.output(sessionInfo()),
             file.path(out_dir, "sessionInfo.txt"))
}
cat(sprintf("done (%s%s): %d settings x %d replicates in %.0f s\n", run,
            if (calibrate) ", calibration" else "", nrow(design), n_rep,
            elapsed))
