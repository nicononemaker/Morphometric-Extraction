# =============================================================================
# config.R  —  set your paths ONCE here; every analysis script sources this file.
#
# Usage: open Morphometric-Extraction.Rproj in RStudio (sets the working
# directory), edit PROJECT_ROOT below, then run the scripts in order.
#
# The ONLY line you MUST edit is PROJECT_ROOT. Everything else derives from it,
# and assumes the folder layout described in the README's "Data availability"
# and "Setup" sections. Adjust the sub-folder names in sections 5-6 only if your
# downloads use different names.
# =============================================================================

# --- 1. THE ONE PATH YOU MUST SET --------------------------------------------
# Absolute path to your PROJECT folder (holds data + outputs; NOT the code repo).
# Use forward slashes even on Windows.
PROJECT_ROOT <- "C:/path/to/project_root"        # <<< EDIT THIS ONE LINE

# --- 2. inputs (derived from PROJECT_ROOT) -----------------------------------
# ImageJ per-frond trait CSVs (measure.py outputs), and the cleaned specimen
# metadata CSV. Create these sub-folders/files under PROJECT_ROOT.
STATS_DIR <- file.path(PROJECT_ROOT, "imagej_stats")            # folder of measure.py CSVs
RAW_CSV   <- file.path(PROJECT_ROOT, "raw_spec_data.csv")       # cleaned specimen metadata

# --- 3. outputs (created automatically) --------------------------------------
DATASETS_DIR <- file.path(PROJECT_ROOT, "datasets")            # built datasets
OUT_DIR      <- file.path(PROJECT_ROOT, "analysis_out")        # figures + tables
dir.create(DATASETS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,      showWarnings = FALSE, recursive = TRUE)

# the record-level dataset 01 builds and 02-06 read (one name, one place)
AGG_CSV <- file.path(DATASETS_DIR, "aggregated_all_fronds_all_traits_w_envr.csv")

# --- 4. CLIMATOLOGY FILES  (Part 1 reads these — the only env input needed) ---
# One single-layer .nc per predictor. If you downloaded the deposited data
# (see README > Data availability), point DATASRC_DIR at that download and these
# resolve automatically. If you rebuild with Part 0 (00_build_climatologies.R),
# it writes the files to exactly these paths. Swap any single path to try a
# different data source for one predictor.
DATASRC_DIR <- file.path(PROJECT_ROOT, "environmental_data")

CLIM_SST       <- file.path(DATASRC_DIR, "climatologies", "sst_climatology.nc")        # SST_clim
CLIM_IRON      <- file.path(DATASRC_DIR, "climatologies", "iron_climatology.nc")       # iron_clim
CLIM_NO3       <- file.path(DATASRC_DIR, "climatologies", "nitrate_climatology.nc")    # nitrate_clim
CLIM_PO4       <- file.path(DATASRC_DIR, "climatologies", "phosphate_climatology.nc")  # phosphate_clim
CLIM_SPCO2     <- file.path(DATASRC_DIR, "climatologies", "spco2_climatology.nc")      # spco2_clim
CLIM_KD490     <- file.path(DATASRC_DIR, "climatologies", "kd490_climatology.nc")      # kd490_clim
CLIM_SALINITY  <- file.path(DATASRC_DIR, "climatologies", "salinity_climatology.nc")   # salinity_clim
CLIM_WAVE      <- file.path(DATASRC_DIR, "climatologies", "wave_climatology.nc")       # wave.height (VHM0)
CLIM_PAR       <- file.path(DATASRC_DIR, "climatologies", "par_climatology.nc")        # PAR_mean

# biogeographic regions (Marine Ecoregions of the World shapefile)
CLIM_MEOW      <- file.path(DATASRC_DIR, "MEOW", "WCMC-036-MEOW-PPOW-2007-2012.shp")

# --- 5. RAW SOURCES  (Part 0 only — 00_build_climatologies.R reads these) -----
# Needed ONLY to rebuild a climatology from scratch (see README > full-rebuild
# path + download/download_env_data.sh). Leave any line "" if you use the
# deposited/pre-built climatology for that predictor instead. Part 1 never
# touches these. Folder names below are the templates the README instructs you
# to create; point them at your actual downloads if they differ.
RAW_DIR <- file.path(DATASRC_DIR, "raw_products")

PRECLIM_SST       <- file.path(RAW_DIR, "hadisst",    "HadISST_sst.nc")   # monthly SST
PRECLIM_IRON      <- file.path(RAW_DIR, "cmems_bio")  # folder; var fe  picked by name
PRECLIM_NO3       <- file.path(RAW_DIR, "cmems_bio")  # folder; var no3 picked by name
PRECLIM_PO4       <- file.path(RAW_DIR, "cmems_bio")  # folder; var po4 picked by name
PRECLIM_SPCO2     <- file.path(RAW_DIR, "cmems_carbon")  # folder; monthly spco2 picked by name
PRECLIM_KD490     <- file.path(RAW_DIR, "cmems_oc")      # folder; monthly KD490 -> annual -> mean
# salinity / waves / PAR download AS climatologies but MONTHLY (12-layer) or
# multi-tile; Part 0 averages each to one layer. Point at the downloads:
PRECLIM_SALINITY  <- file.path(RAW_DIR, "cmems_phy")    # folder; var so picked by name
PRECLIM_WAVE      <- file.path(RAW_DIR, "cmems_wav", "waves.nc")   # var VHM0
PRECLIM_PAR       <- file.path(RAW_DIR, "modis_par",
                               c("par1.nc","par2.nc","par3.nc","par4.nc"))  # tiles, averaged

# climatology window for the self-computed products (SST, fe, no3, po4, pCO2, Kd490)
CLIM_MIN <- 1993
CLIM_MAX <- 2025
