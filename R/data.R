#' Read a GEO series matrix
#'
#' @param path Path to a `_series_matrix.txt.gz` file.
#' @return A list with `samples` (one row per sample: accession, title, and
#'   one column per characteristics key) and `table` (the data table as a
#'   matrix, or `NULL` when the matrix holds none).
#' @export
ReadSeriesMatrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  field <- function(tag) {
    l <- grep(paste0("^!Sample_", tag, "\t"), lines, value = TRUE)
    lapply(l, function(x) gsub('^"|"$', "", strsplit(x, "\t")[[1]][-1]))
  }
  acc <- field("geo_accession")[[1]]
  samples <- data.frame(geo = acc, title = field("title")[[1]],
                        stringsAsFactors = FALSE)
  # characteristics rows may hold different keys for different samples
  for (row in field("characteristics_ch1")) {
    for (i in seq_along(row)) {
      kv <- regmatches(row[i], regexpr(": ", row[i]), invert = TRUE)[[1]]
      if (length(kv) != 2) next
      key <- make.names(kv[1])
      if (!key %in% names(samples)) samples[[key]] <- NA_character_
      samples[[key]][i] <- kv[2]
    }
  }
  b <- grep("^!series_matrix_table_begin", lines)
  e <- grep("^!series_matrix_table_end", lines)
  table <- NULL
  if (length(b) && e > b + 2) {
    tab <- utils::read.delim(text = lines[(b + 1):(e - 1)],
                             check.names = FALSE,
                             na.strings = c("", "null", "NA"))
    table <- as.matrix(tab[, -1])
    rownames(table) <- tab[[1]]
  }
  list(samples = samples, table = table)
}

#' Describe the GSE235508 cohort for the pilot's application section
#'
#' Classifies each sample as pregnant (time points 1 to 3) or not (0, 4, 5,
#' 6), and counts, by group, the samples, the pregnancies, and the
#' pregnancies with at least one sample in each period (the ones a
#' within-woman comparison can use).
#'
#' @param samples The `samples` element of [ReadSeriesMatrix()] for
#'   GSE235508, with columns `samplegroup`, `timepoint`, and `pregnancyid`.
#' @return A list: `by_group` (a data frame with one row per group and
#'   `group`, `samples`, `pregnancies`, `usable_pregnancies`,
#'   `pregnant_samples`, `nonpregnant_samples`), and `by_time` (a table of
#'   samples by group and time point).
#' @export
DescribePregnancyCohort <- function(samples) {
  s <- samples
  tp <- as.integer(s$timepoint)
  s$pregnant <- tp %in% 1:3
  by_group <- do.call(rbind, lapply(sort(unique(s$samplegroup)), function(g) {
    d <- s[s$samplegroup == g, ]
    both <- tapply(d$pregnant, d$pregnancyid,
                   function(v) any(v) && any(!v))
    usable <- names(both)[both]
    du <- d[d$pregnancyid %in% usable, ]
    data.frame(group = g, samples = nrow(d),
               pregnancies = length(unique(d$pregnancyid)),
               usable_pregnancies = length(usable),
               pregnant_samples = sum(du$pregnant),
               nonpregnant_samples = sum(!du$pregnant),
               stringsAsFactors = FALSE)
  }))
  list(by_group = by_group,
       by_time = table(group = s$samplegroup, timepoint = tp))
}
