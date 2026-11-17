# analysis/03_simulation.R
# The simulation study. Until the design stop of stage 3 it carries pilot-g's
# design (PilotGDesign()); protocol.md replaces the design table. Run from the
# project root:
#   Rscript analysis/03_simulation.R --run quick --reps 20
#   Rscript analysis/03_simulation.R --run production
# (01_ and 02_ are the earlier pregnancy scripts, kept as the record.)
# Each replicate's seed is its cell's seed plus the replicate index; cells are
# split into seeded blocks of replicates, the tasks of parLapplyLB(chunk.size
# = 1) on PSOCK workers with the BLAS pinned to one thread.
# Writes results/: replicates.rds, cell_summaries.rds, known_answers.rds,
# design.rds, example_networks.rds, run_info.rds, sessionInfo.txt.

devtools::load_all(here::here(), quiet = TRUE)
library(parallel)

# ---- parameter block ---------------------------------------------------------
design <- PilotGDesign()
n_rep_production <- 200L
block <- 20L        # replicates per task
n_workers <- 3L     # the core cap
example_cell <- 2L  # the cell whose first replicate the network figures draw
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}
run <- opt("--run", "quick")
stopifnot(run %in% c("quick", "production"))
n_rep <- if (run == "production") n_rep_production else
  as.integer(opt("--reps", 20))
set.seed(design$seed[1])  # no random step outside the replicates
out_dir <- here::here("results")
dir.create(out_dir, showWarnings = FALSE)

starts <- seq(1L, n_rep, by = block)
tasks <- expand.grid(from = starts, cell = design$cell)
tasks$to <- pmin(tasks$from + block - 1L, n_rep)

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
cl <- makePSOCKcluster(n_workers)
clusterExport(cl, "design")
clusterEvalQ(cl, {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
  devtools::load_all(here::here(), quiet = TRUE)
  NULL
})
start <- Sys.time()
res <- parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), function(tk) {
  row <- design[design$cell == tk$cell, ]
  do.call(rbind, lapply(tk$from:tk$to, function(r) PathwayReplicate(row, r)))
}, chunk.size = 1)
stopCluster(cl)
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

reps <- do.call(rbind, res)
reps <- reps[order(reps$cell, reps$rep), ]
rownames(reps) <- NULL
summaries <- SummarizeCellsG(reps)
saveRDS(reps, file.path(out_dir, "replicates.rds"))
saveRDS(summaries, file.path(out_dir, "cell_summaries.rds"))
saveRDS(KnownAnswerChecksG(summaries), file.path(out_dir, "known_answers.rds"))
saveRDS(design, file.path(out_dir, "design.rds"))
saveRDS(ExampleNetworks(design[design$cell == example_cell, ], 1L),
        file.path(out_dir, "example_networks.rds"))
saveRDS(list(run = run, n_rep = n_rep, n_cells = nrow(design),
             date = Sys.time(), elapsed_seconds = elapsed,
             n_workers = n_workers,
             note = "pilot-g design until the design stop of stage 3",
             n_rep_production = n_rep_production),
        file.path(out_dir, "run_info.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat(sprintf("done (%s): %d cells x %d replicates in %.0f s\n", run,
            nrow(design), n_rep, elapsed))
