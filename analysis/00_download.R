# analysis/00_download.R
# Fetches the public GSE235508 files into data/raw/, skipping any file that
# already exists; nothing here modifies a raw file afterwards. Run from the
# project root:  Rscript analysis/00_download.R
#
# GSE235508 series matrix (16 KB): the sample characteristics (group, time
#   point, pregnancy, disease activity, medication).
# GSE235508 processed counts (GSE235508_mRNA_counts.txt.gz, 25 MB), approved
#   by the user on 2026-11-13; not committed (.gitignore), fetched by this
#   script.
# GSE65391 series matrix (91 MB, normalized Illumina HT-12 data with the
#   clinical characteristics) and the GPL10558 platform table (probe to gene
#   symbol), for the re-scope (proposal of 2026-11-16, default D-R5); not
#   committed.
# BloodGen3Module 1.20.0 source package from Bioconductor (183 KB, GPL-2):
#   the BloodGen3 module lists (Altman et al. 2021) in R/sysdata.rda,
#   approved by the user on 2026-11-17; not committed.

devtools::load_all(here::here(), quiet = TRUE)
set.seed(20761113)  # no random step; the seed is set as every script does

raw_dir <- here::here("data", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
base <- "https://ftp.ncbi.nlm.nih.gov/geo/series"
files <- data.frame(
  accession = c("GSE235508", "GSE235508", "GSE65391", "GPL10558",
                "BloodGen3Module"),
  url = c(file.path(base, c(
    "GSE235nnn/GSE235508/matrix/GSE235508_series_matrix.txt.gz",
    "GSE235nnn/GSE235508/suppl/GSE235508_mRNA_counts.txt.gz",
    "GSE65nnn/GSE65391/matrix/GSE65391_series_matrix.txt.gz")),
    paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558",
           "&targ=self&form=text&view=data"),
    paste0("https://bioconductor.org/packages/3.23/bioc/src/contrib/",
           "BloodGen3Module_1.20.0.tar.gz")),
  stringsAsFactors = FALSE
)
files$file <- file.path(raw_dir, ifelse(files$accession == "GPL10558",
                                        "GPL10558_platform_table.txt",
                                        basename(files$url)))
files$fetched <- FALSE
for (i in seq_len(nrow(files))) {
  if (file.exists(files$file[i])) next
  for (attempt in 1:4) {
    ok <- tryCatch({
      utils::download.file(files$url[i], files$file[i], mode = "wb",
                           quiet = TRUE, method = "curl",
                           extra = "--retry 3 --fail")
      TRUE
    }, error = function(e) FALSE)
    if (ok) break
    Sys.sleep(2^attempt)
  }
  stopifnot(ok, file.exists(files$file[i]))
  files$fetched[i] <- TRUE
}
files$bytes <- file.size(files$file)
files$md5 <- unname(tools::md5sum(files$file))
files$file <- basename(files$file)
dir.create(here::here("results"), showWarnings = FALSE)
utils::write.csv(files, here::here("results", "00_download_checksums.csv"),
                 row.names = FALSE)
print(files[, c("file", "bytes", "md5", "fetched")])
saveRDS(list(run = "production", n_rep = NA_integer_, n_cells = NA_integer_,
             date = Sys.time(),
             note = "GSE235508, GSE65391, GPL10558 table, BloodGen3Module",
             n_rep_production = NA_integer_),
        here::here("results", "00_download_run_info.rds"))
writeLines(capture.output(sessionInfo()),
           here::here("results", "00_download_sessionInfo.txt"))
