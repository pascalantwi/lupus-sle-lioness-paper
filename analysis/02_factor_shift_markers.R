# analysis/02_factor_shift_markers.R
# The cell-mix change over pregnancy in GSE235508 from marker genes, within
# SLE and within RA (seropositive and seronegative together), never pooled
# across diseases. From log2(CPM + 1):
#   neutrophil score: mean of FCGR3B (CD16b), CXCR2, CSF3R, FPR1;
#   lymphocyte score: mean of CD3E, CD3D, CD2 (T cells), MS4A1 (CD20),
#     CD79A (B cells);
#   composition axis: neutrophil score minus lymphocyte score.
# Ensembl identifiers from the Ensembl REST lookup, 2026-11-13. Pregnant:
# time points 1 to 3; women with samples in both periods. Also the visits
# per woman in each period. Output: results/factor_shift_markers.rds. These
# numbers set the levels of pilot-c; no edge is computed here.
#   Rscript analysis/02_factor_shift_markers.R

devtools::load_all(here::here(), quiet = TRUE)
set.seed(20961113)  # no random step
raw <- here::here("data", "raw")
s <- ReadSeriesMatrix(file.path(raw, "GSE235508_series_matrix.txt.gz"))$samples
lc <- ReadLogCpm(file.path(raw, "GSE235508_mRNA_counts.txt.gz"))
neut <- c(FCGR3B = "ENSG00000162747", CXCR2 = "ENSG00000180871",
          CSF3R = "ENSG00000119535", FPR1 = "ENSG00000171051")
lymph <- c(CD3E = "ENSG00000198851", CD3D = "ENSG00000167286",
           CD2 = "ENSG00000116824", MS4A1 = "ENSG00000156738",
           CD79A = "ENSG00000105369")
stopifnot(all(c(neut, lymph) %in% rownames(lc)))
s$pregnant <- as.integer(s$timepoint) %in% 1:3
s$disease <- ifelse(s$samplegroup == "SLE", "SLE",
                    ifelse(s$samplegroup %in% c("SPRA", "SNRA"), "RA",
                           "Healthy"))
out <- list()
visits <- list()
for (d in c("SLE", "RA")) {
  sd_ <- s[s$disease == d, ]
  both <- tapply(sd_$pregnant, sd_$pregnancyid, function(v) any(v) && any(!v))
  sd_ <- sd_[sd_$pregnancyid %in% names(both)[both], ]
  x <- lc[, sd_$sampleID]
  nscore <- colMeans(x[neut, , drop = FALSE])
  lscore <- colMeans(x[lymph, , drop = FALSE])
  out[[d]] <- rbind(
    cbind(disease = d, score = "neutrophil",
          ShiftSummary(nscore, sd_$pregnant, sd_$pregnancyid)),
    cbind(disease = d, score = "lymphocyte",
          ShiftSummary(lscore, sd_$pregnant, sd_$pregnancyid)),
    cbind(disease = d, score = "composition axis",
          ShiftSummary(nscore - lscore, sd_$pregnant, sd_$pregnancyid)))
  vp <- tapply(sd_$pregnant, sd_$pregnancyid, sum)
  vn <- tapply(!sd_$pregnant, sd_$pregnancyid, sum)
  visits[[d]] <- data.frame(disease = d, women = length(vp),
                            mean_pregnant = mean(vp), mean_nonpregnant = mean(vn))
}
res <- list(shift = do.call(rbind, out), visits = do.call(rbind, visits))
rownames(res$shift) <- NULL
saveRDS(res, here::here("results", "factor_shift_markers.rds"))
print(res, digits = 3)
saveRDS(list(run = "production", n_rep = NA_integer_, n_cells = NA_integer_,
             date = Sys.time(),
             note = "marker-gene composition shift over pregnancy, GSE235508",
             n_rep_production = NA_integer_),
        here::here("results", "factor_shift_markers_run_info.rds"))
writeLines(capture.output(sessionInfo()),
           here::here("results", "factor_shift_markers_sessionInfo.txt"))
