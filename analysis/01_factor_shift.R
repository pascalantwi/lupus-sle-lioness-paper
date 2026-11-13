# analysis/01_factor_shift.R
# The size of the shared-factor change over pregnancy in GSE235508, within
# SLE and within RA (seropositive and seronegative together), never pooled
# across diseases. Two scores per sample, from log2(CPM + 1):
#   IFN5: the mean of EPSTI1, IFI44L, LY6E, OAS3, RSAD2 (Mai et al. 2021;
#     Ensembl identifiers from the Ensembl REST lookup, 2026-11-13; LY6E
#     maps to two identifiers, the primary-assembly one is used);
#   PC1: the first principal component of the 2,000 most variable genes
#     within the disease, a data-driven proxy for the main composition axis.
# Pregnant: time points 1 to 3; non-pregnant: 0, 4, 5, 6. Only women with
# samples in both periods. Output: results/factor_shift.rds. These numbers
# set the levels of pilot-b; no edge is computed here.
#   Rscript analysis/01_factor_shift.R

devtools::load_all(here::here(), quiet = TRUE)
set.seed(20761113)  # no random step
raw <- here::here("data", "raw")
sm <- ReadSeriesMatrix(file.path(raw, "GSE235508_series_matrix.txt.gz"))
s <- sm$samples
lc <- ReadLogCpm(file.path(raw, "GSE235508_mRNA_counts.txt.gz"))
stopifnot(all(s$sampleID %in% colnames(lc)))
ifn5 <- c(EPSTI1 = "ENSG00000133106", IFI44L = "ENSG00000137959",
          LY6E = "ENSG00000160932", OAS3 = "ENSG00000111331",
          RSAD2 = "ENSG00000134321")
stopifnot(all(ifn5 %in% rownames(lc)))
s$pregnant <- as.integer(s$timepoint) %in% 1:3
s$disease <- ifelse(s$samplegroup == "SLE", "SLE",
                    ifelse(s$samplegroup %in% c("SPRA", "SNRA"), "RA",
                           "Healthy"))
out <- list()
for (d in c("SLE", "RA")) {
  sd_ <- s[s$disease == d, ]
  both <- tapply(sd_$pregnant, sd_$pregnancyid, function(v) any(v) && any(!v))
  sd_ <- sd_[sd_$pregnancyid %in% names(both)[both], ]
  x <- lc[, sd_$sampleID]
  ifn <- colMeans(x[ifn5, , drop = FALSE])
  expressed <- rowMeans(x) > 1
  v <- apply(x[expressed, ], 1, stats::var)
  top <- names(sort(v, decreasing = TRUE))[1:2000]
  pc1 <- stats::prcomp(t(x[top, ]), scale. = TRUE)$x[, 1]
  ifn_genes_r <- stats::cor(t(x[ifn5, ]))
  out[[d]] <- rbind(
    cbind(disease = d, score = "IFN5",
          ShiftSummary(ifn, sd_$pregnant, sd_$pregnancyid)),
    cbind(disease = d, score = "PC1",
          ShiftSummary(pc1, sd_$pregnant, sd_$pregnancyid)))
  out[[d]]$mean_r_ifn5_genes <- mean(ifn_genes_r[upper.tri(ifn_genes_r)])
}
res <- do.call(rbind, out)
rownames(res) <- NULL
dir.create(here::here("results"), showWarnings = FALSE)
saveRDS(res, here::here("results", "factor_shift.rds"))
print(res, digits = 3)
saveRDS(list(run = "production", n_rep = NA_integer_, n_cells = NA_integer_,
             date = Sys.time(),
             note = "IFN5 and PC1 shift over pregnancy, GSE235508",
             n_rep_production = NA_integer_),
        here::here("results", "factor_shift_run_info.rds"))
writeLines(capture.output(sessionInfo()),
           here::here("results", "factor_shift_sessionInfo.txt"))
