# analysis/00_download.R
# Fetches the public GSE235508 files into data/raw/, skipping any file that
# already exists; nothing here modifies a raw file afterwards. Run from the
# project root:  Rscript analysis/00_download.R
#
# GSE235508 series matrix (16 KB): the sample characteristics (group, time
#   point, pregnancy, disease activity, medication). The processed counts
#   (GSE235508_mRNA_counts.txt.gz, 25 MB) are added here at stage 3, after the
#   user approves that download (proposal.md, section 9).

devtools::load_all(here::here(), quiet = TRUE)
set.seed(20761113)  # no random step; the seed is set as every script does

raw_dir <- here::here("data", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
base <- "https://ftp.ncbi.nlm.nih.gov/geo/series"
files <- data.frame(
  accession = "GSE235508",
  url = file.path(base,
                  "GSE235nnn/GSE235508/matrix/GSE235508_series_matrix.txt.gz"),
  stringsAsFactors = FALSE
)
files$file <- file.path(raw_dir, basename(files$url))
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
             date = Sys.time(), note = "GSE235508 series matrix download",
             n_rep_production = NA_integer_),
        here::here("results", "00_download_run_info.rds"))
writeLines(capture.output(sessionInfo()),
           here::here("results", "00_download_sessionInfo.txt"))
