# pilot/describe_data.R
# Describes GSE235508 for the pilot report's application section: samples,
# pregnancies, and pregnancies with samples in both periods, by group.
# Variable codes are read from the GEO series record (opened 2026-11-13):
# samplegroup (SLE, SPRA seropositive RA, SNRA seronegative RA, HEALTHY),
# timepoint (0 pre-pregnancy; 1 to 3 trimesters; 4, 5, 6 at 6 weeks,
# 6 months, 12 months postpartum), pregnancyid.
#   Rscript pilot/describe_data.R   (after Rscript analysis/00_download.R)

devtools::load_all(here::here(), quiet = TRUE)
set.seed(20761113)  # no random step
sm <- ReadSeriesMatrix(here::here("data", "raw",
                                  "GSE235508_series_matrix.txt.gz"))
desc <- DescribePregnancyCohort(sm$samples)
# RA as one group (seropositive and seronegative together), as an extra row
ra <- sm$samples
ra <- ra[ra$samplegroup %in% c("SPRA", "SNRA"), ]
ra$samplegroup <- "RA (SPRA and SNRA)"
desc$by_group <- rbind(desc$by_group, DescribePregnancyCohort(ra)$by_group)
out_dir <- here::here("pilot", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(desc, file.path(out_dir, "data_summary.rds"))
print(desc)
saveRDS(list(run = "pilot", n_rep = NA_integer_, n_cells = NA_integer_,
             date = Sys.time(), note = "GSE235508 description",
             n_rep_production = NA_integer_),
        file.path(out_dir, "describe_data_run_info.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "describe_data_sessionInfo.txt"))
