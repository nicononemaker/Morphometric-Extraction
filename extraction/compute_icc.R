# =============================================================================
# compute_icc.R
#
# Measurement-repeatability analysis. Takes THREE independent measurement passes
# of the SAME repeatability subsample (three measurements.csv files, one per
# pass) and computes the intraclass correlation coefficient (ICC) per trait.
#
# Reports ICC(3,k) -- the reliability of the MEAN of the three measurements
# (two-way mixed, consistency, average of k raters) -- which is the statistic
# the manuscript reports ("ICC for the mean of three measurements"). ICC(3,1)
# (single-measurement reliability) is also reported for completeness.
#
# Each input is a measure.py output (per-frond rows). Records are aggregated to
# the record level (mean across fronds) within each pass, then matched across
# passes by record key, so ICC is computed on record-level trait values --
# matching the analysis's unit of replication.
#
# Usage:
#   Rscript compute_icc.R pass1.csv pass2.csv pass3.csv [out.csv]
#   (any number of passes >= 2 works; 3 is the manuscript design)
# =============================================================================
# =============================================================================
# HOW TO RUN
#   RStudio (simplest): set PASS_FILES (and optional OUT) in the CONFIG block
#     below, then click Source. No terminal needed.
#   Command line:  Rscript compute_icc.R pass1.csv pass2.csv pass3.csv [out.csv]
#   (command-line .csv arguments, if given, override the CONFIG block.)
# =============================================================================
if (!requireNamespace("psych", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("psych")
}
suppressPackageStartupMessages(library(psych))

# ---- CONFIG: edit these for an RStudio run ----------------------------------
PASS_FILES <- c(
  "C:/path/to/project_root/measurements_pass1.csv",  # <<< pass 1
  "C:/path/to/project_root/measurements_pass2.csv",  # <<< pass 2
  "C:/path/to/project_root/measurements_pass3.csv"   # <<< pass 3
)
OUT <- ""     # output path; "" = icc_repeatability.csv next to pass 1

# command-line args override CONFIG when run via Rscript
args  <- commandArgs(trailingOnly = TRUE)
csvs  <- args[grepl("\\.csv$", args)]
if (length(csvs) >= 2) {
  # if a 4th csv is given and doesn't exist yet, treat it as the output path
  if (length(csvs) > 3 && !file.exists(csvs[length(csvs)])) {
    OUT <- csvs[length(csvs)]; csvs <- csvs[-length(csvs)]
  }
  PASS_FILES <- csvs
}
files <- PASS_FILES
if (length(files) < 2) stop("Need at least 2 measurement passes (set PASS_FILES or pass .csv args).")
if (!nzchar(OUT)) OUT <- file.path(dirname(files[1]), "icc_repeatability.csv")

KEY   <- "filename"                 # record key in measure.py output
FROND <- "frond"
# traits to assess (the manuscript's morphometric suite; extras skipped if absent)
TRAITS <- c("area_mm2","perimeter_mm","length_main_axis_mm","max_width_mm",
            "median_width_mm","mean_width_mm","circularity","solidity",
            "perim_over_sqrt_area","aspect_ratio","stipe_width_mm",
            "stipe_to_length_ratio","branched_fraction","peak_position_frac")

# ---- read each pass, aggregate to record level (mean across fronds) ----------
agg_pass <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!KEY %in% names(d)) stop(path, " has no '", KEY, "' column.")
  tr <- intersect(TRAITS, names(d))
  num <- data.frame(lapply(d[tr], function(x) suppressWarnings(as.numeric(x))))
  num[[KEY]] <- d[[KEY]]
  aggregate(num[tr], by = list(record = num[[KEY]]),
            FUN = function(x) mean(x, na.rm = TRUE))
}
passes <- lapply(files, agg_pass)
cat(sprintf("Read %d passes: %s\n", length(passes),
            paste(basename(files), collapse = ", ")))

# ---- records measured in ALL passes (complete cases for ICC) ----------------
common <- Reduce(intersect, lapply(passes, function(p) p$record))
cat(sprintf("Records present in all passes: %d\n", length(common)))
if (length(common) < 3) stop("Too few shared records for a meaningful ICC.")

TRAITS <- Reduce(intersect, c(list(TRAITS), lapply(passes, names)))
TRAITS <- setdiff(TRAITS, "record")

# ---- per-trait ICC ----------------------------------------------------------
rows <- lapply(TRAITS, function(tr) {
  # matrix: rows = records, cols = passes (raters)
  M <- sapply(passes, function(p) p[match(common, p$record), tr])
  M <- M[stats::complete.cases(M), , drop = FALSE]
  if (nrow(M) < 3) return(NULL)
  ic <- tryCatch(psych::ICC(M), error = function(e) NULL)
  if (is.null(ic)) return(NULL)
  r  <- ic$results
  get <- function(type) r$ICC[r$type == type]
  data.frame(
    trait = tr, n_records = nrow(M), n_passes = ncol(M),
    ICC_single = round(get("ICC3"), 3),      # ICC(3,1): one measurement
    ICC_mean   = round(get("ICC3k"), 3),     # ICC(3,k): mean of k measurements (reported)
    CI_low     = round(r$`lower bound`[r$type == "ICC3k"], 3),
    CI_high    = round(r$`upper bound`[r$type == "ICC3k"], 3),
    row.names = NULL, stringsAsFactors = FALSE)
})
res <- do.call(rbind, Filter(Negate(is.null), rows))
res <- res[order(-res$ICC_mean), ]
write.csv(res, OUT, row.names = FALSE)

cat("\n=== Measurement repeatability (ICC) ===\n")
cat("ICC_mean = reliability of the mean of the ", length(files),
    " passes (the reported statistic).\n", sep = "")
print(res, row.names = FALSE)
cat(sprintf("\nAll traits ICC_mean >= %.2f: %s\n", 0.9,
            if (all(res$ICC_mean >= 0.9, na.rm = TRUE)) "YES" else "NO (see table)"))
cat("Wrote", OUT, "\n")
