# =============================================================================
# 00_build_climatologies.R  —  OPTIONAL Part 0.
#
# Turns raw environmental downloads into the single-layer climatology .nc files
# that Part 1 (01_build_dataset.R) reads. Each predictor is read from its
# PRECLIM_* path (config §6) and written to its CLIM_* path (config §5).
#
# SKIP this entirely if you use the paper's pre-built climatologies — just put
# them at the CLIM_* paths and run 01 directly.
#
# Every predictor follows the same shape: read the raw source, pick its
# variable, REDUCE its layers to one climatology layer, write. Predictors differ
# only in the reduction:
#   window_mean : mean of monthly layers in CLIM_MIN..CLIM_MAX  (SST, fe, no3, po4, pCO2, Kd490)
#   layer_mean  : mean of all layers/tiles, no time window      (salinity, waves, PAR)
# =============================================================================
suppressPackageStartupMessages({ library(terra) })

# --- locate and source config.R (works with the .Rproj open, via source(), or
# from the console) --------------------------------------------------------
.find_config <- function() {
  for (p in c("config.R", "../config.R")) if (file.exists(p)) return(p)
  sp <- if (requireNamespace("rstudioapi", quietly = TRUE))
    tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "") else ""
  if (!nzchar(sp)) sp <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
  if (nzchar(sp)) {
    cand <- file.path(dirname(dirname(sp)), "config.R")
    if (file.exists(cand)) return(cand)
  }
  stop("config.R not found. Open Morphometric-Extraction.Rproj, or set your ",
       "working directory to the repo root, then re-run.")
}
source(.find_config())

have     <- function(p) length(p) && all(nzchar(basename(p))) && all(file.exists(p))
find_var <- function(r, pats) { for (p in pats) { idx <- grep(p, names(r), ignore.case = TRUE)
if (length(idx)) return(r[[idx]]) }; NULL }

# --- the three reductions ----------------------------------------------------
window_mean <- function(r) {                 # monthly layers -> mean over CLIM window
  yr <- suppressWarnings(as.integer(format(time(r), "%Y")))
  keep <- which(!is.na(yr) & yr >= CLIM_MIN & yr <= CLIM_MAX)
  if (!length(keep)) keep <- seq_len(terra::nlyr(r))
  mean(r[[keep]], na.rm = TRUE)
}
layer_mean  <- function(r) if (terra::nlyr(r) > 1) mean(r, na.rm = TRUE) else r  # all layers

# --- one uniform builder for every predictor ---------------------------------
# src      : PRECLIM_* path (file, files, or a folder of .nc)
# var_pats : name patterns to pick the variable (NULL = use all layers)
# reduce   : window_mean | layer_mean
build <- function(src, var_pats, reduce, varname, out, label) {
  cat(label, "\n")
  if (!have(src)) { cat("  SKIP: source not found\n"); return(invisible()) }
  files <- if (length(src) == 1 && dir.exists(src))
    list.files(src, pattern = "\\.nc$", full.names = TRUE) else src
  r <- rast(files)
  NAflag(r) <- -1000                          # undeclared -1000 fill (HadISST)
  r <- classify(r, cbind(1e19, Inf, NA))      # also catch huge positive fills (1e20/1e30)
  v <- if (is.null(var_pats)) r else find_var(r, var_pats)
  if (is.null(v)) { cat("  SKIP: variable not found\n"); return(invisible()) }
  layer <- reduce(v)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  writeCDF(layer, out, varname = varname, overwrite = TRUE)
  cat(sprintf("  -> %-22s range %.4f .. %.4f\n", basename(out),
              global(layer, min, na.rm=TRUE)[1,1], global(layer, max, na.rm=TRUE)[1,1]))
}

cat("Building climatology files (window", CLIM_MIN, "-", CLIM_MAX, ")\n")

#        source            variable pattern        reduction     out-var          out-file        label
build(PRECLIM_SST,      NULL,                   window_mean, "SST_clim",       CLIM_SST,      "SST (HadISST)")
build(PRECLIM_IRON,     c("^fe","iron"),        window_mean, "iron_clim",      CLIM_IRON,     "iron (BIO 001_029)")
build(PRECLIM_NO3,      c("^no3","nitrate"),    window_mean, "nitrate_clim",   CLIM_NO3,      "nitrate (BIO 001_029)")
build(PRECLIM_PO4,      c("^po4","phosphate"),  window_mean, "phosphate_clim", CLIM_PO4,      "phosphate (BIO 001_029)")
build(PRECLIM_SPCO2,    c("^spco2"),            window_mean, "spco2_clim",     CLIM_SPCO2,    "pCO2 (Chau 2024)")
build(PRECLIM_KD490,    c("^KD490","kd490"),    window_mean, "kd490_clim",     CLIM_KD490,    "Kd490 (Ocean Colour 009_104)")
build(PRECLIM_SALINITY, c("^so"),               layer_mean,  "salinity_clim",  CLIM_SALINITY, "salinity (PHY 001_030)")
build(PRECLIM_WAVE,     c("VHM0","wave.?height","swh"), layer_mean, "wave.height", CLIM_WAVE, "waves (WAV 001_032)")
build(PRECLIM_PAR,      NULL,                   layer_mean,  "PAR_mean",       CLIM_PAR,      "PAR (MODIS-Aqua)")

cat("\nDone. Part 1 (01_build_dataset.R) now reads all CLIM_* files.\n")