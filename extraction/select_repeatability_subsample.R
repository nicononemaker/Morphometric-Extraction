# =============================================================================
# select_repeatability_subsample.R
#
# Draws a stratified-random subsample of records to be RE-MEASURED (three
# passes) for a measurement-repeatability (ICC) analysis. Stratifies across
# frond size, collection year, and source institution so the subsample spans
# the range of average record size, age, and institutional preservation
# differences -- matching the manuscript's repeatability design.
#
# Reusable: run on any record table that has the stratification columns below.
# Reproducible: set SEED to regenerate the exact same draw.
#
# INPUT  : a records CSV (e.g. the cleaned metadata seaweed_data_clean.csv, or
#          an aggregated dataset) with columns:
#            Image.ID (record key), Average.Frond.Length, Year,
#            Herbarium.Data.Institution
# OUTPUT : repeatability_subsample.csv  (the N records to re-measure)
#
# Usage:
#   Rscript select_repeatability_subsample.R <input.csv> [N] [seed] [out.csv]
#   defaults: N = 24, seed = 42, out = repeatability_subsample.csv (next to input)
# =============================================================================

# =============================================================================
# HOW TO RUN
#   RStudio (simplest): set the four values in the CONFIG block just below,
#     then click Source (or source this file). No terminal needed.
#   Command line:  Rscript select_repeatability_subsample.R <input.csv> [N] [seed] [out.csv]
#   (command-line arguments, if given, override the CONFIG block below.)
# =============================================================================

# ---- CONFIG: edit these for an RStudio run ----------------------------------
IN   <- "C:/path/to/project_root/spec_data_w_prelim_length.csv"  # <<< records CSV with Image.ID, Average.Frond.Length, Year, Herbarium.Data.Institution
N    <- 24     # how many records to draw
SEED <- 42     # random seed (same seed -> same draw)
OUT  <- ""     # output path; "" = repeatability_subsample.csv next to input

# command-line args override the CONFIG block when the script is run via Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) IN   <- args[[1]]
if (length(args) >= 2) N    <- as.integer(args[[2]])
if (length(args) >= 3) SEED <- as.integer(args[[3]])
if (length(args) >= 4) OUT  <- args[[4]]
N <- as.integer(N); SEED <- as.integer(SEED)
if (!nzchar(OUT)) OUT <- file.path(dirname(IN), "repeatability_subsample.csv")

# ---- stratification columns (edit here if your column names differ) ---------
KEY_COL   <- "Image.ID"
SIZE_COL  <- "Average.Frond.Length"      # frond size
YEAR_COL  <- "Year"                      # collection year
INST_COL  <- "Herbarium.Data.Institution"# source institution

# ---- stratification toggle --------------------------------------------------
# TRUE  = stratify across size (terciles) x year (terciles) x institution
# FALSE = stratify across year (terciles) x institution only (drop frond size)
STRATIFY_SIZE <- TRUE

d <- read.csv(IN, stringsAsFactors = FALSE, check.names = FALSE)
strat_cols <- c(YEAR_COL, INST_COL)
if (STRATIFY_SIZE) strat_cols <- c(SIZE_COL, strat_cols)
need <- c(KEY_COL, strat_cols)
miss <- need[!need %in% names(d)]
if (length(miss)) stop("Input is missing required column(s): ", paste(miss, collapse = ", "))

# keep only records with all stratification values present
d <- d[stats::complete.cases(d[, strat_cols]), ]
if (nrow(d) < N) stop(sprintf("Only %d complete records; cannot draw %d.", nrow(d), N))

# ---- build strata -----------------------------------------------------------
# terciles give low/mid/high bins for continuous variables; institution is
# categorical. Frond size is included only when STRATIFY_SIZE is TRUE.
tercile <- function(x) cut(x, breaks = quantile(x, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                           include.lowest = TRUE, labels = c("low","mid","high"))
d$.year_bin <- tercile(as.numeric(d[[YEAR_COL]]))
d$.inst     <- as.character(d[[INST_COL]])
if (STRATIFY_SIZE) {
  d$.size_bin <- tercile(as.numeric(d[[SIZE_COL]]))
  d$.stratum  <- interaction(d$.size_bin, d$.year_bin, d$.inst, drop = TRUE)
  cat("Stratifying on: frond size (terciles) x year (terciles) x institution\n")
} else {
  d$.stratum  <- interaction(d$.year_bin, d$.inst, drop = TRUE)
  cat("Stratifying on: year (terciles) x institution  (frond size NOT used)\n")
}

# ---- proportional allocation across strata, then random pick within ---------
set.seed(SEED)
strata <- split(seq_len(nrow(d)), d$.stratum)
# target count per stratum, proportional to stratum size, rounded, >=0
prop   <- vapply(strata, length, integer(1)) / nrow(d)
alloc  <- floor(prop * N)
# distribute the remainder to the largest strata so totals hit N exactly
rem <- N - sum(alloc)
if (rem > 0) {
  ord <- order(prop * N - alloc, decreasing = TRUE)
  alloc[ord[seq_len(rem)]] <- alloc[ord[seq_len(rem)]] + 1
}
pick <- unlist(mapply(function(idx, k) if (k > 0) sample(idx, min(k, length(idx))) else integer(0),
                      strata, alloc, SIMPLIFY = FALSE), use.names = FALSE)
# if rounding left us short (tiny strata), top up randomly from the remainder
if (length(pick) < N) {
  extra <- sample(setdiff(seq_len(nrow(d)), pick), N - length(pick))
  pick  <- c(pick, extra)
}
pick <- pick[seq_len(N)]

out_cols <- c(KEY_COL, strat_cols)
sub <- d[pick, out_cols]
sub <- sub[order(sub[[YEAR_COL]]), ]
write.csv(sub, OUT, row.names = FALSE)

cat(sprintf("Drew %d records (seed %d) -> %s\n", nrow(sub), SEED, OUT))
yr <- suppressWarnings(as.numeric(sub[[YEAR_COL]]))
yr_lo <- if (all(is.na(yr))) "?" else as.character(min(yr, na.rm = TRUE))
yr_hi <- if (all(is.na(yr))) "?" else as.character(max(yr, na.rm = TRUE))
if (STRATIFY_SIZE) {
  sz <- suppressWarnings(as.numeric(sub[[SIZE_COL]]))
  cat(sprintf("  size range: %s - %s | years: %s - %s | institutions: %d\n",
              round(min(sz, na.rm = TRUE),1), round(max(sz, na.rm = TRUE),1),
              yr_lo, yr_hi, length(unique(sub[[INST_COL]]))))
} else {
  cat(sprintf("  years: %s - %s | institutions: %d\n",
              yr_lo, yr_hi, length(unique(sub[[INST_COL]]))))
}
cat("Re-measure these records three times, then run compute_icc.R.\n")
