# Stage 1 re-scope: data checks (2026-11-15 and 2026-11-16)

Descriptive checks run by the agent outside the package, in base R (`Rscript --vanilla`), to see what GSE65391 and ILLUMINATE (GSE88884) can support before the proposal was written. No edge was computed. The input files were downloaded to the session's scratch folder, not to `data/raw/`; stage 2 brings the download into `analysis/00_download.R`.

| Script | Data | Output |
|---|---|---|
| `g65_check.R` | GSE65391 series matrix (GEO, 91 MB) and the GPL10558 probe-to-symbol column | `g65_check.txt`: IFN5 score, measured cell percentages, SLEDAI, domains, kidney class, next-visit flares |
| `g65_more.R` | as above | `g65_more.txt`: factors by kidney class and by treatment |
| `g65_icc.R` | as above | `g65_icc.txt`: within-child stability (intraclass correlation), visit gaps, the authors' training and test split |
| `illum_check.R` | GSE88884 preprocessed matrix (GEO, 1.3 GB), marker rows only, and the GPL17586 probe-to-symbol column | `illum_check.txt`: IFN, B-cell, and composition scores against activity, labs, ancestry, sex |
| `illum_sero.R` | as above | `illum_sero.txt`: the lab profile and the healthy batches |

The scripts take their input folders as arguments; their paths are the session's, not the repository's, so they are a record, not a pipeline.

**Error found at stage 2 (2026-11-16).** `g65_icc.R`'s helper misaligned values and groups when a child had a missing count, so its neutrophil (0.22) and lymphocyte (0.30) stabilities are wrong; the tested `EdgeICC()` gives 0.41 and 0.48 (`pilot-d/results/data_summary.rds`). Its IFN and SLEDAI values are right (no missing values). Kept as recorded, with this note.
