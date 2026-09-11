# =============================================================================
# 01_build_dataset.R  —  traits + metadata + environment  (Part 1)
#
# Builds the record-level analysis table used by every downstream script:
#   aggregated_all_fronds_all_traits_w_envr.csv
#
# Pipeline:
#   1. compile per-frond ImageJ trait CSVs
#   2. merge with specimen metadata (coords + year + region-ready fields)
#   3. extract 9 environmental predictors from pre-built CLIMATOLOGY FILES at
#      each unique coordinate, join back
#   4. assign MEOW realm/province/ecoregion
#   5. aggregate fronds -> one row per herbarium record
#
# This reads single-layer climatology .nc files (one per predictor), NOT the raw
# multi-GB rasters. Get those files one of two ways:
#   - download the paper's pre-built climatologies (README > Data availability), OR
#   - run 00_build_climatologies.R (Part 0) to rebuild the self-computed ones,
#     and download the precomputed ones (salinity, waves, PAR, Kd490 recipe).
# All nine file paths are set in config.R (section 5).
#
# ENVIRONMENTAL PREDICTORS (the nine used in the manuscript):
#   SST_clim, salinity_clim, spco2_clim, wave.height, kd490_clim, PAR_mean,
#   phosphate_clim, nitrate_clim, iron_clim  (see config.R for source files)
# =============================================================================

required <- c("dplyr","terra","raster","ncdf4","sf")
new_pkgs <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(new_pkgs)) {
  options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages(new_pkgs)
}
suppressPackageStartupMessages(invisible(lapply(required, library, character.only = TRUE)))

# ------------------------------------------------------------------ config ----
# All paths live in config.R at the repo root — edit PROJECT_ROOT there, not here.
# Run with the repo root as your working directory (RStudio: Session > Set Working
# Directory > To Project/Source File location), or source() this script directly.
# --- locate and source config.R (works with the .Rproj open, via source(), or
# from the console) --------------------------------------------------------
.find_config <- function() {
  for (p in c("config.R", "../config.R")) if (file.exists(p)) return(p)
  # last resort: derive from this script's own path (RStudio or source())
  sp <- if (requireNamespace("rstudioapi", quietly = TRUE))
    tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "") else ""
  if (!nzchar(sp)) sp <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
  if (nzchar(sp)) {
    cand <- file.path(dirname(dirname(sp)), "config.R")   # repo root from analysis/
    if (file.exists(cand)) return(cand)
  }
  stop("config.R not found. Open Morphometric-Extraction.Rproj, or set your ",
       "working directory to the repo root, then re-run.")
}
source(.find_config())

# shim: map config's names to the locals this script uses. Part 1 reads the
# per-predictor CLIMATOLOGY FILES (from Part 0 or the paper's deposit).
proj_root   <- PROJECT_ROOT
stats_dir   <- STATS_DIR
raw_csv     <- RAW_CSV
out_dir     <- DATASETS_DIR
f_sst_clim       <- CLIM_SST
f_iron_clim      <- CLIM_IRON
f_nitrate_clim   <- CLIM_NO3
f_phosphate_clim <- CLIM_PO4
f_spco2_clim     <- CLIM_SPCO2
f_kd490_clim     <- CLIM_KD490
f_sal_clim       <- CLIM_SALINITY
f_wave           <- CLIM_WAVE
f_par            <- CLIM_PAR
f_meow           <- CLIM_MEOW

# Region columns are re-assigned fresh in Step 3i; drop any stale copies that
# might ride in on the metadata so the merge doesn't create .x/.y duplicates.
# (The cleaned metadata carries no environmental columns, so nothing else needs
# stripping here.)
DROP_META <- c("ECOREGION","REALM","TYPE","PROVINC","PROVINCE","BIOME")

# helpers ---------------------------------------------------------------------
have <- function(p) length(p) && all(nzchar(basename(p))) && all(file.exists(p))

KEY_CANDIDATES <- c("Image.ID","Image ID","Image_ID","ImageID","image.id","image_id",
                    "filename","file_name","file","Label","label","source_file",
                    "record.ID","record_id","ID","id")
canon_key <- function(df, side = "") {
  hit <- KEY_CANDIDATES[tolower(KEY_CANDIDATES) %in% tolower(names(df))]
  if (!length(hit)) stop(sprintf(
    "No specimen-id column found on the %s side. Columns present:\n  %s",
    side, paste(names(df), collapse = ", ")))
  col <- names(df)[tolower(names(df)) == tolower(hit[1])][1]
  if (col != "filename") names(df)[names(df) == col] <- "filename"
  df$filename <- trimws(basename(as.character(df$filename)))
  df$filename <- sub("\\.(tif|tiff|csv|png|jpg|jpeg)$", "", df$filename, ignore.case = TRUE)
  df
}

# =============================================================================
# STEP 1 — compile ImageJ trait CSVs
# =============================================================================
cat("[1] compiling ImageJ trait CSVs...\n")
trait_files <- list.files(stats_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(trait_files) > 0)
traits <- do.call(rbind, lapply(trait_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d$source_file <- basename(f); d
}))
traits <- canon_key(traits, side = "ImageJ traits")
cat(sprintf("    %d frond rows, %d unique specimens\n",
            nrow(traits), length(unique(traits$filename))))

# =============================================================================
# STEP 2 — merge with specimen metadata
# =============================================================================
cat("[2] merging specimen metadata...\n")
meta <- read.csv(raw_csv, stringsAsFactors = FALSE, check.names = FALSE)
meta <- canon_key(meta, side = "specimen metadata")
if (!all(c("Longitude","Latitude","Year") %in% names(meta)))
  stop("Metadata missing Longitude/Latitude/Year.")
meta <- meta[, !(names(meta) %in% DROP_META), drop = FALSE]
dup <- setdiff(intersect(names(traits), names(meta)), "filename")
if (length(dup)) meta <- meta[, !(names(meta) %in% dup), drop = FALSE]

fr <- merge(traits, meta, by = "filename", all.x = TRUE, sort = FALSE)
cat(sprintf("    %d frond rows | %d unmatched to metadata\n",
            nrow(fr), sum(is.na(fr$Longitude))))
write.csv(fr, file.path(out_dir, "all_specimens_all_traits.csv"), row.names = FALSE)

# =============================================================================
# STEP 3 — extract environment on UNIQUE coordinates, join back
# =============================================================================
cat("[3] extracting environment...\n")
suppressWarnings({
  fr$Longitude <- as.numeric(as.character(fr$Longitude))
  fr$Latitude  <- as.numeric(as.character(fr$Latitude))
  fr$Year      <- as.integer(as.character(fr$Year))
})
uni <- unique(fr[!is.na(fr$Longitude) & !is.na(fr$Latitude),
                 c("Longitude","Latitude","Year")])
uni$.key <- paste(uni$Longitude, uni$Latitude, uni$Year, sep = "_")
cat(sprintf("    %d unique lon/lat/year combinations\n", nrow(uni)))
coords <- uni[, c("Longitude","Latitude")]

# near-shore fallback: step the point diagonally until a non-NA ocean cell is hit
extract_ocean <- function(lon, lat, r, step = 0.05, max_iter = 20) {
  v <- raster::extract(r, cbind(lon, lat)); i <- 0
  while (is.na(v) && i < max_iter) { lon <- lon + step; lat <- lat + step
    v <- raster::extract(r, cbind(lon, lat)); i <- i + 1 }
  v
}
# extract one single-layer climatology raster at all points, with the fallback
ext_clim <- function(r) unlist(mapply(extract_ocean, uni$Longitude, uni$Latitude,
                                      MoreArgs = list(r = r)))
# read a single-layer climatology file, extract at all points, report coverage
extract_clim_file <- function(path, colname, label) {
  cat(sprintf("    %s\n", label))
  if (!have(path)) { cat("       ", basename(path), "not found; skipped\n"); return(invisible()) }
  r <- rast(path)
  if (terra::nlyr(r) > 1) r <- mean(r, na.rm = TRUE)   # safety: collapse if not 1 layer
  uni[[colname]] <<- ext_clim(r)
}

# All nine predictors are single-layer climatology files (Part 0 writes the
# self-computed ones as single-layer; salinity/wave/PAR are single-layer too).
# Extract each at the record coordinates. Ocean fields use the near-shore
# fallback; PAR uses a plain grid lookup.
extract_clim_file(f_sst_clim,       "SST_clim",       "3a SST")
extract_clim_file(f_iron_clim,      "iron_clim",      "3b iron")
extract_clim_file(f_nitrate_clim,   "nitrate_clim",   "3c nitrate")
extract_clim_file(f_phosphate_clim, "phosphate_clim", "3d phosphate")
extract_clim_file(f_spco2_clim,     "spco2_clim",     "3e pCO2")
extract_clim_file(f_kd490_clim,     "kd490_clim",     "3f Kd490")
extract_clim_file(f_sal_clim,       "salinity_clim",  "3g salinity")
extract_clim_file(f_wave,           "wave.height",    "3h wave height")

# PAR: plain grid extract (no near-shore fallback needed)
cat("    3i PAR\n")
if (have(f_par)) {
  r <- rast(f_par); if (terra::nlyr(r) > 1) r <- mean(r, na.rm = TRUE)
  uni$PAR_mean <- terra::extract(r, uni[, c("Longitude","Latitude")])[, 2]
} else cat("        PAR climatology not found; skipped\n")

# ---- MEOW ecoregions --------------------------------------------------------
# Coastal points often fall just outside the ecoregion polygons, so a strict
# point-in-polygon join leaves many NA. Restrict to MEOW ecoregions (drop PPOW
# pelagic provinces), repair invalid polygon rings, then snap each point to its
# nearest ecoregion within MAX_SNAP_KM.
cat("    3j MEOW ecoregions\n")
if (have(f_meow)) {
  meow <- st_read(f_meow, quiet = TRUE); sf::sf_use_s2(FALSE)
  type_col <- intersect(c("TYPE","Type","type"), names(meow))
  if (length(type_col)) {
    keep <- toupper(trimws(as.character(meow[[type_col[1]]]))) == "MEOW"
    if (!any(keep) && "ECOREGION" %in% names(meow))
      keep <- !is.na(meow$ECOREGION) & trimws(meow$ECOREGION) != ""
    meow <- meow[keep, ]
  } else if ("ECOREGION" %in% names(meow)) {
    meow <- meow[!is.na(meow$ECOREGION) & trimws(meow$ECOREGION) != "", ]
  }
  cat(sprintf("       %d MEOW ecoregion polygons (PPOW removed)\n", nrow(meow)))

  MAX_SNAP_KM <- 50
  meow <- st_make_valid(meow); meow <- meow[!st_is_empty(meow), ]
  meow_cols <- intersect(c("ECOREGION","PROVINCE","PROVINC","REALM","BIOME"), names(meow))
  uni_sf <- st_as_sf(uni, coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE)
  ni  <- st_nearest_feature(uni_sf, meow)
  seg <- st_nearest_points(uni_sf, meow[ni, ], pairwise = TRUE)
  ep  <- st_coordinates(seg)
  haversine_km <- function(lon1, lat1, lon2, lat2) {
    R <- 6371; d2r <- pi/180
    dlat <- (lat2-lat1)*d2r; dlon <- (lon2-lon1)*d2r
    a <- sin(dlat/2)^2 + cos(lat1*d2r)*cos(lat2*d2r)*sin(dlon/2)^2
    2*R*asin(pmin(1, sqrt(a)))
  }
  o <- ep[seq(1, nrow(ep), by = 2), ]; p <- ep[seq(2, nrow(ep), by = 2), ]
  dk <- haversine_km(o[,1], o[,2], p[,1], p[,2])
  within <- dk <= MAX_SNAP_KM
  cat(sprintf("       nearest-feature: %d of %d within %d km\n",
              sum(within), length(within), MAX_SNAP_KM))
  for (c_ in meow_cols) {
    vals <- rep(NA_character_, nrow(uni_sf))
    vals[within] <- as.character(meow[[c_]][ni[within]])
    uni[[c_]] <- vals
  }
} else cat("       f_meow not set; MEOW skipped\n")

# ---- join environment back to the frond table -------------------------------
cat("    joining environment back to fronds...\n")
fr$.key <- paste(fr$Longitude, fr$Latitude, fr$Year, sep = "_")
env_cols <- setdiff(names(uni), c("Longitude","Latitude","Year",".key"))
fr_env <- merge(fr, uni[, c(".key", env_cols)], by = ".key", all.x = TRUE, sort = FALSE)
fr_env$.key <- NULL
write.csv(fr_env, file.path(out_dir, "all_fronds_all_traits_w_envr.csv"), row.names = FALSE)

# quick sanity print for the two self-computed fields worth verifying
for (v in c("SST_clim","spco2_clim")) if (v %in% names(uni))
  cat(sprintf("    %-11s: %d/%d non-NA, mean %.3f\n", v,
              sum(!is.na(uni[[v]])), nrow(uni), mean(uni[[v]], na.rm = TRUE)))

# =============================================================================
# STEP 4 — aggregate by record (mean across fronds)
# =============================================================================
cat("[4] aggregating by record...\n")
region_cols <- intersect(c("ECOREGION","PROVINCE","PROVINC","REALM","BIOME"), names(fr_env))
const_cols  <- unique(c(env_cols, region_cols, "Longitude","Latitude","Year",
                        setdiff(names(meta), "filename")))
const_cols  <- intersect(const_cols, names(fr_env))
num_cols    <- names(fr_env)[sapply(fr_env, is.numeric)]
avg_cols    <- setdiff(num_cols, const_cols)

agg <- fr_env %>%
  group_by(filename) %>%
  summarise(n_fronds = dplyr::n(),
            across(all_of(avg_cols),   ~ mean(.x, na.rm = TRUE)),
            across(any_of(const_cols), ~ dplyr::first(.x)),
            .groups = "drop") %>%
  as.data.frame()
for (c_ in avg_cols) if (c_ %in% names(agg)) agg[[c_]][is.nan(agg[[c_]])] <- NA

# analysis scripts and the manuscript use "frond length"; the extraction schema
# calls this the main axis. Bridge the two names here, at the single assembly point.
if ("length_main_axis_mm" %in% names(agg))
  names(agg)[names(agg) == "length_main_axis_mm"] <- "length_frond_mm"

# stipe_to_length_ratio is a RATIO-OF-MEANS at the record level: recompute it from
# the aggregated stipe width and frond length, rather than averaging the per-frond
# ratios (which is fragile to a single short/damaged frond). This matches the
# manuscript's record-level values.
if (all(c("stipe_width_mm","length_frond_mm") %in% names(agg))) {
  agg$stipe_to_length_ratio <- ifelse(
    is.finite(agg$length_frond_mm) & agg$length_frond_mm > 0,
    agg$stipe_width_mm / agg$length_frond_mm, NA)
}

cat(sprintf("    %d records\n", nrow(agg)))
write.csv(agg, file.path(out_dir, "aggregated_all_fronds_all_traits_w_envr.csv"),
          row.names = FALSE)

cat("\nDone. Files written to", out_dir, ":\n")
cat("  all_specimens_all_traits.csv\n")
cat("  all_fronds_all_traits_w_envr.csv\n")
cat("  aggregated_all_fronds_all_traits_w_envr.csv\n")
