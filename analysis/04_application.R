# analysis/04_application.R
# The application of protocol.md (section 8, entry of 2026-11-17): LIONESS
# networks of the GSE65391 children with two or more visits, read through
# BloodGen3 pathways. Run from the project root after analysis/00_download.R:
#   Rscript analysis/04_application.R
# Writes results/: application_summary.rds, application_modules.rds,
# application_metrics.rds, application_prediction.rds,
# application_examples.rds, application_run_info.rds,
# application_sessionInfo.txt.

devtools::load_all(here::here(), quiet = TRUE)
library(parallel)

# ---- parameter block (protocol.md, sections 4 and 8) -------------------------
n_workers <- 3L
n_boot <- 200L
min_genes <- 10L
# ------------------------------------------------------------------------------

set.seed(21461118)
raw <- here::here("data", "raw")
out_dir <- here::here("results")
start <- Sys.time()
sm <- ReadSeriesMatrix(file.path(raw, "GSE65391_series_matrix.txt.gz"))
platform <- ReadPlatformSymbols(file.path(raw, "GPL10558_platform_table.txt"))
desc <- DescribeGse65391(sm, platform)
q <- desc$samples  # ordered by child, then cumulative time
person <- q$subject
visit <- stats::ave(seq_along(person), person, FUN = seq_along)
sledai <- as.numeric(q$sledai)
bg3 <- ReadBloodGen3(file.path(raw, "BloodGen3Module_1.20.0.tar.gz"))
x_all <- ExpressionByGene(sm$table[, q$geo], platform,
                          unique(bg3$modules$gene))
# a gene with no variation cannot be standardized
x_all <- x_all[, apply(x_all, 2, stats::sd) > 0, drop = FALSE]

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
cl <- makePSOCKcluster(n_workers)
clusterEvalQ(cl, {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
  devtools::load_all(here::here(), quiet = TRUE)
  NULL
})

# ---- module diagnostics (claim 5; gradient placement) ------------------------
mod_genes <- split(bg3$modules$gene, bg3$modules$module)
mod_genes <- lapply(mod_genes, function(g) intersect(unique(g),
                                                     colnames(x_all)))
mod_genes <- mod_genes[lengths(mod_genes) >= min_genes]
clusterExport(cl, c("x_all", "person"))
diag <- parLapplyLB(cl, mod_genes, function(g) {
  ModuleDiagnostics(x_all[, g, drop = FALSE], person)
})
modules <- data.frame(module = names(mod_genes), do.call(rbind, diag),
                      stringsAsFactors = FALSE)
ann <- bg3$annotation
modules$fun <- ann[["function"]][match(modules$module, ann$module)]
modules$cluster <- ann$cluster[match(modules$module, ann$module)]

# ---- the four-pathway LIONESS network ----------------------------------------
pw <- BloodGen3Pathways(bg3, colnames(x_all))
x_net <- x_all[, unlist(pw, use.names = FALSE)]
ends <- cumsum(lengths(pw))
sets <- stats::setNames(lapply(seq_along(pw), function(i) {
  (ends[i] - lengths(pw)[i] + 1):ends[i]
}), names(pw))
m_raw <- LionessMetricsBySample(x_net, sets, cl = cl)
m_fo <- LionessMetricsBySample(x_net, sets, first_order = TRUE, cl = cl)
stopCluster(cl)
icc <- data.frame(metric = colnames(m_raw),
                  raw = EdgeICC(m_raw, person), fo = EdgeICC(m_fo, person))
icc$d <- icc$raw - icc$fo
z_net <- StandardizeN(x_net)$z
score <- sapply(sets, function(s) rowMeans(z_net[, s]))
icc_score2 <- EdgeICC(score^2, person)

# ---- identity and prediction (claim 7) ----------------------------------------
keep_cols <- function(m) scale(m[, apply(m, 2, stats::sd) > 0, drop = FALSE])
feats <- list(A = keep_cols(m_raw), B = keep_cols(cbind(score, score^2)),
              C = keep_cols(m_fo))
finger <- vapply(feats, function(f) {
  FingerprintAccuracy(f, person, visit, seq_len(ncol(f)), 1, 2)
}, numeric(1))
outer <- GroupFolds(person, 5L)
same <- lapply(feats, OutOfFold, y = sledai, person = person, outer = outer)
t_same <- BootstrapTests(same, sledai, person, n_boot)
nxt <- match(paste(person, visit + 1), paste(person, visit))
use <- which(!is.na(nxt))
ref <- matrix(sledai[use], ncol = 1)
nx <- c(lapply(feats, function(f) {
  OutOfFold(f[use, , drop = FALSE], sledai[nxt[use]], person[use],
            outer[use], ref = ref)
}), list(R0 = OutOfFold(matrix(0, length(use), 0), sledai[nxt[use]],
                        person[use], outer[use], ref = ref)))
t_next <- BootstrapTests(nx, sledai[nxt[use]], person[use], n_boot)
split_c <- SplitConcordance(feats, sledai, person, q$set == "Training")

# ---- example networks for the figures -----------------------------------------
ifn_child <- tapply(q$ifn, person, mean)
kids <- c(names(which.max(ifn_child)),
          names(ifn_child)[which.min(abs(ifn_child - stats::median(ifn_child)))])
pr <- PairIndex(ncol(x_net))
rp <- stats::cor(x_net)[pr]
ex_rows <- which(person %in% kids & visit <= 2)
ex_edges <- t(vapply(ex_rows, function(i) {
  SampleEdges(z_net[i, ], rp, pr, nrow(x_net))
}, numeric(nrow(pr))))
examples <- list(kids = kids, person = person[ex_rows], visit = visit[ex_rows],
                 edges = ex_edges, metrics = m_raw[ex_rows, , drop = FALSE],
                 sets = sets, pathway_genes = pw,
                 all_metrics = m_raw, all_person = person, all_visit = visit,
                 score2 = score^2)

elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
saveRDS(list(n_children = length(unique(person)), n_visits = nrow(q),
             n_genes_matched = ncol(x_all), n_modules = nrow(modules),
             pathway_sizes = lengths(pw), split = table(q$set),
             icc_sledai = EdgeICC(sledai, person),
             icc_ifn5 = EdgeICC(q$ifn, person)),
        file.path(out_dir, "application_summary.rds"))
saveRDS(modules, file.path(out_dir, "application_modules.rds"))
saveRDS(list(icc = icc, icc_score2 = icc_score2),
        file.path(out_dir, "application_metrics.rds"))
saveRDS(list(finger = finger, same = t_same, next_visit = t_next,
             split = split_c, chance = 1 / length(unique(person))),
        file.path(out_dir, "application_prediction.rds"))
saveRDS(examples, file.path(out_dir, "application_examples.rds"))
saveRDS(list(run = "production", n_rep = NA_integer_, n_cells = NA_integer_,
             date = Sys.time(), elapsed_seconds = elapsed,
             note = "GSE65391 application, protocol.md section 8",
             n_rep_production = NA_integer_),
        file.path(out_dir, "application_run_info.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "application_sessionInfo.txt"))
cat(sprintf("done: %d modules, %d network genes, %.0f s\n", nrow(modules),
            ncol(x_net), elapsed))
