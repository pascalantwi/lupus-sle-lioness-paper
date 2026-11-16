# pilot-d/describe_data.R
# Describes GSE65391 for the pilot-d report: visits, children, the
# within-child stability of the IFN5 score, the measured neutrophil and
# lymphocyte percentages, and SLEDAI; days between visits; the authors'
# training and test split. Variable names are read from the GEO sample
# characteristics (series opened 2026-11-15): subject, disease state,
# sledai, neutrophil_percent, lymphocyte_percent, cumulative_time, set.
#   Rscript pilot-d/describe_data.R   (after Rscript analysis/00_download.R)

devtools::load_all(here::here(), quiet = TRUE)
set.seed(21061116)  # no random step
raw <- here::here("data", "raw")
sm <- ReadSeriesMatrix(file.path(raw, "GSE65391_series_matrix.txt.gz"))
platform <- ReadPlatformSymbols(file.path(raw, "GPL10558_platform_table.txt"))
desc <- DescribeGse65391(sm, platform)
out_dir <- here::here("pilot-d", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
desc$samples <- desc$samples[, c("geo", "subject", "visit", "cumulative_time",
                                 "set", "sledai", "ifn", "neutrophil_percent",
                                 "lymphocyte_percent")]
saveRDS(desc, file.path(out_dir, "data_summary.rds"))
print(desc[c("counts", "icc", "gap_days", "split")])
saveRDS(list(run = "pilot", n_rep = NA_integer_, n_cells = NA_integer_,
             date = Sys.time(), note = "GSE65391 description",
             n_rep_production = NA_integer_),
        file.path(out_dir, "describe_data_run_info.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "describe_data_sessionInfo.txt"))
