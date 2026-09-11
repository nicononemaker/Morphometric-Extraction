#!/usr/bin/env bash
# =============================================================================
# download_env_data.sh — fetch the raw CMEMS products for the FULL-REBUILD path.
#
# Only needed if you are rebuilding the climatologies from scratch (Part 0 /
# 00_build_climatologies.R). Most users should instead download the pre-built
# climatologies from the data deposit — see the repo README > Data availability.
#
# This downloads the five CMEMS products. SST (HadISST), PAR (MODIS-Aqua) and
# the MEOW shapefile are NOT on CMEMS — get those from the links at the bottom.
#
# SETUP (once):
#   pip install copernicusmarine
#   copernicusmarine login          # free CMEMS account
#
# USAGE:
#   1. set DEST below to the raw_products/ folder your config.R points at
#      (PROJECT_ROOT/environmental_data/raw_products by default).
#   2. bash download_env_data.sh
#
# These are large (the Kd490 monthly series alone is hundreds of GB). Run on a
# good connection with ample disk. Folder names below match config.R's PRECLIM_*.
#
# NOTE ON DATASET IDs: verify each --dataset-id against the current CMEMS
# catalogue (`copernicusmarine describe --contains <product-number>`) before a
# long download — CMEMS occasionally revises dataset IDs when products are
# reprocessed. IDs below are the ones used for this study.
# =============================================================================
set -euo pipefail

# <<< EDIT: your raw_products folder (matches RAW_DIR in config.R)
DEST="C:/path/to/project_root/environmental_data/raw_products"

# global bounding box + common 1993-2025 climatology window (self-computed vars)
BBOX="--minimum-longitude -180 --maximum-longitude 180 --minimum-latitude -90 --maximum-latitude 90"
START="1993-01-01"
END="2025-12-31"

mkdir -p "$DEST"

# --- BIO 001_029 : iron + nitrate + phosphate (one folder: cmems_bio) --------
echo ">> BIO 001_029 (iron, nitrate, phosphate)"
mkdir -p "$DEST/cmems_bio"
copernicusmarine subset \
  --dataset-id cmems_mod_glo_bgc_my_0.25deg_P1M-m \
  --variable fe --variable no3 --variable po4 \
  $BBOX --start-datetime "$START" --end-datetime "$END" \
  --output-directory "$DEST/cmems_bio"

# --- Carbon (Chau et al. 2024) : surface ocean pCO2 (cmems_carbon) -----------
echo ">> Carbon (pCO2 / spco2)"
mkdir -p "$DEST/cmems_carbon"
copernicusmarine subset \
  --dataset-id cmems_obs-mob_glo_bgc-car_my_irr-i \
  --variable spco2 \
  $BBOX --start-datetime "$START" --end-datetime 2024-12-01 \
  --output-directory "$DEST/cmems_carbon"

# --- PHY 001_030 : salinity monthly climatology (cmems_phy) ------------------
echo ">> PHY 001_030 (salinity)"
mkdir -p "$DEST/cmems_phy"
copernicusmarine subset \
  --dataset-id cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m \
  --variable so \
  $BBOX \
  --output-directory "$DEST/cmems_phy"

# --- WAV 001_032 : significant wave height monthly climatology (cmems_wav) ----
echo ">> WAV 001_032 (wave height)"
mkdir -p "$DEST/cmems_wav"
copernicusmarine subset \
  --dataset-id cmems_mod_glo_wav_my_0.2deg-climatology_P1M-m \
  --variable VHM0 \
  $BBOX \
  --output-directory "$DEST/cmems_wav"

# --- OceanColour 009_104 : Kd490 monthly (cmems_oc) --------------------------
# 009_104 has no ready-made Kd490 climatology: download the MONTHLY 4km series;
# Part 0 (00_build_climatologies.R) reduces monthly -> climatology mean.
# WARNING: full 1997-present monthly global 4km is very large (hundreds of GB).
echo ">> OceanColour 009_104 (Kd490 monthly)"
mkdir -p "$DEST/cmems_oc"
copernicusmarine subset \
  --dataset-id cmems_obs-oc_glo_bgc-transp_my_l4-multi-4km_P1M \
  --variable KD490 \
  $BBOX --start-datetime 1997-01-01 --end-datetime 2024-12-31 \
  --output-directory "$DEST/cmems_oc"

echo ""
echo "CMEMS downloads done. Three sources are NOT on CMEMS — get these manually:"
echo "  SST  : HadISST v1.1  -> $DEST/hadisst/HadISST_sst.nc"
echo "         https://www.metoffice.gov.uk/hadobs/hadisst/"
echo "  PAR  : MODIS-Aqua L3 PAR (monthly/climatology tiles) -> $DEST/modis_par/"
echo "         https://oceancolor.gsfc.nasa.gov/"
echo "  MEOW : Marine Ecoregions of the World shapefile -> ../MEOW/"
echo "         https://data.unep-wcmc.org/datasets/38"
echo ""
echo "Then run 00_build_climatologies.R (Part 0) to build the climatology files."
