# Asparagopsis taxiformis morphometrics

Code for *Global Diversity in* Asparagopsis taxiformis *Morphology: Biogeographic Patterns and Environmental Structure*. This project analyzed the global intraspecific morphological variation in the gametophyte life stage of the red alga *Asparagopsis taxiformis* using digitized herbarium records. In this analysis we were able to identify significant biogeographic and environmental structure to morphological variation.  



The project's repository has two halves:

* **`extraction/`** - a semi-automated Fiji/ImageJ + Python pipeline that measures morphometric traits from herbarium specimen images. See [`extraction/README.md`](extraction/README.md) for full detail.
* **`analysis/`** - the R pipeline that built the morphometric dataset (trait values + environmental data + biogeographic region groupings) and runs every statistical analysis and produces each figure and table in the manuscript, in addition to exploratory results.



### Trait extraction (`extraction/`)

Three stages, run in order on a folder of specimen images:

1. **`capture.ijm`** (Fiji/ImageJ) - set pixel scale using ruler attached to image, then for each frond trace the outline, main axis (full length of stipe) and the basal stipe width. Lastly, designate the boundaries of the branched-region. Ultimately, a crop of the herbarium image is saved alongside the traced and plotted coordinates. Writes manifest.csv, \*\_crop.png, \*\_mainaxis.txt, \*\_stipe.txt, and \*\_branched.txt.

2. **`review.py`** - interactive thresholding (saturation + brightness sliders) for each cropped frond to   mask thallus tissue. Writes `thresholds.csv`.

3. **`measure.py`** — uses each frond's mask to computes the trait values for the suite. Writes `measurements.csv`.

### 

### Analysis (`analysis/`)

**Setup:** open `Morphometric-Extraction.Rproj` in RStudio setting your working directory to the repo root. Then edit `config.R` and set `PROJECT\_ROOT` to the folder containing your data and where you would like outputs to be saved to. If you would like to run the scripts outside of RStudio, set your working directory to the repo root before running any script.



**Project folder layout.** `config.R` assumes this structure under `PROJECT\_ROOT` (create the folders you need; the `datasets/` and `analysis\_out/` output folders are created automatically):

```
PROJECT\_ROOT/
├── raw\_spec\_data.csv            # specimen metadata (config: RAW\_CSV)
├── imagej\_stats/                # CSV containing per-frond trait values (config: STATS\_DIR)
├── environmental\_data/
│   ├── climatologies/           # the 9 climatology .nc files script 01 reads
│   │   ├── sst\_climatology.nc
│   │   ├── iron\_climatology.nc  ... (one per predictor)
│   ├── MEOW/                     # MEOW shapefile (WCMC-036-MEOW-PPOW-2007-2012.\*)
│   └── raw\_products/            # ONLY if rebuilding climatologies (part 0)
│       ├── hadisst/  cmems\_bio/  cmems\_carbon/  cmems\_oc/
│       └── cmems\_phy/  cmems\_wav/  modis\_par/
├── datasets/                    # built datasets (datasets generated in script 01)
└── analysis\_out/                # figures + tables (results generated in scripts 02-06; organized in 07)
```

The full analysis pipeline requires `raw\_spec\_data.csv`, `imagej\_stats/`, and `environmental\_data/climatologies/` + `MEOW/` (cleaned climatologies can be downloaded from the deposit, see *Data availability*). The `raw\_products/` tree is only for rebuilding climatologies from scratch (script 0). If your downloads use different folder names, adjust the matching paths in `config.R` sections 4–5.

The pipeline reads single layer climatology files (long term average at each coordinate pair). Pre-build climatolgies are available via download from the paper's deposit (see *Data availability* below) or can be rebuilt from multi-layer climatolgies and hindcast data products using script 00 (requires multiple multi-GB downloads). Script 00 and 01 can be skipped by downloading `aggregated\\\_all\\\_fronds\\\_all\\\_traits\\\_w\\\_envr.csv` (included in the deposit) into you datasets/ folder. 



|Script|Role|What it does|Manuscript outputs|
|-|-|-|-|
|`00\_build\_climatologies.R`|Part 0 (optional)|Rebuilds the self-computed climatology rasters (SST, iron, nitrate, phosphate, pCO₂, Kd490) and aggregates multi-layer climatolgies from raw data products and writes them to the `CLIM\_\*` paths in `config.R`. This part can be skiped if you use the deposited climatologies.|—|
|`01\_build\_dataset.R`|Part 1|Reads the nine climatology files, extracts each at every record's coordinates, assigns MEOW region, aggregates fronds to record level|builds full dataset: `aggregated\_all\_fronds\_all\_traits\_w\_envr.csv`|
|`02\_pca\_and\_pc\_regression.R`||PCA + OLS/GLS regression of PC1/PC2 on environment; biplot, scree, coefficient/Spearman plots; region \& climate-zone tests|Figs 2, 3, S1, S4, S5; Table 1; Tables S2–S4, S6|
|`03\_permanova\_analysis.R`||Region PERMANOVA/PERMDISP/pairwise, environmental PERMANOVA, variance partitioning; writes `.rds` caches|Tables S5, S7, S8|
|`04\_permanova\_figures.R`||Draws the PERMANOVA figures from the cached results|Figs 4, 5, S6, S7|
|`05\_trait\_analysis.R`||Per-trait regressions vs environment|Table S9|
|`06\_trait\_supp\_export.R`||Assembles the per-trait tables into the supplementary table|Table S9 (formatted)|

Run in order. `00` and `01` can be skipped by downloading data from the deposit.



### Environmental predictors



Nine predictors are extracted and assigned to specimens in `01\_build\_dataset.R`. Self-computed climatologies (SST, nitrate, phosphate, iron, pCO₂) are the mean of all monthly layers over an overlapping window (min 1993; max 2025). Precomputed monthly climatology products (salinity, waves, Kd490, PAR) are reduced to one layer (mean of monthly layers) and maintain their native reference periods.



|Predictor (column)|Variable|Product|Reference period|
|-|-|-|-|
|`SST\_clim`|sea surface temperature|HadISST v1.1 (Met Office Hadley Centre)|1993–2025 (self-computed)|
|`salinity\_clim`|salinity|CMEMS Global Ocean Physics Reanalysis (GLORYS12V1)|1993–2016 (precomputed)|
|`spco2\_clim`|surface ocean pCO₂|CMEMS observation-based carbon (Chau et al. 2024)|1993–2024 (self-computed)|
|`wave.height`|significant wave height|CMEMS Global Ocean Waves Reanalysis (WAVERYS)|1993–2020 (precomputed)|
|`kd490\_clim`|turbidity (Kd490)|Copernicus-GlobColour (OCEANCOLOUR\_GLO\_BGC\_L4\_MY\_009\_104)|1997–2023 (precomputed)|
|`PAR\_mean`|photosynthetically available radiation|MODIS-Aqua L3 (NASA OBPG)|2002–2025 (precomputed)|
|`phosphate\_clim`|phosphate|CMEMS Global Ocean Biogeochemistry Hindcast|1993–2025 (self-computed)|
|`nitrate\_clim`|nitrate|CMEMS Global Ocean Biogeochemistry Hindcast|1993–2025 (self-computed)|
|`iron\_clim`|iron|CMEMS Global Ocean Biogeochemistry Hindcast|1993–2025 (self-computed)|

Biogeographic regions use the Marine Ecoregions of the World scheme (Spalding et al. 2007).



## Data availability

Specimen images, coordinate files, and result CSVs are git-ignored, so the repo holds only code. Environmental data is available in two tiers.



**Compressed single-layer climatologyy (from deposited data; faster run through).** A Zenodo deposit provides the derived data needed to reproduce the analysis without the multi-GB raw environmental downloads and the time intensive compression of this data to a one layer climatolgy:

> (https://doi.org/10.5281/zenodo.22715237)



The deposit contains:

* `climatologies/` - the nine single-layer climatology rasters the analysis reads (one per environmental predictor), the direct inputs to Part 1.
* `aggregated\_all\_fronds\_all\_traits\_w\_envr.csv` - the record-level analysis dataset (traits + extracted environment + MEOW region), one row per herbarium record. Can be used to skip scripts 00 and 01, the input every analysis script (`02`–`06`) reads.
* `all\_specimens\_all\_traits.csv` - intermediate per-frond table (traits before environment/region were joined).
* `raw\_spec\_data.csv` - cleaned specimen metadata (coordinates, year, institution, catalog numbers).
* `spec\_data\_w\_prelim\_length.csv` — specimen metadata with preliminary frond-length measurements (used by the repeatability sampler).



Derived climatologies are redistributed under the source providers' open-data terms with attribution (see `download/DATA\_SOURCES\_ATTRIBUTION.md`). To run using these data products: download the deposit, point `config.R`'s `CLIM\_\*` paths at the `climatologies/` files, and either run the pipeline from `01\_build\_dataset.R`, or skip straight to `02`–`06` using the provided aggregated CSV.



**Full-rebuild (build climatolgies from raw data products; multi-hour process).** To regenerate the climatologies from the original products, download the raw sources below, run `download/download\_env\_data.sh` (CMEMS products; `pip install copernicusmarine` then `copernicusmarine login`), then run `00\_build\_climatologies.R` (Part 0) to build the self-computed ones. SST, PAR, and MEOW come from the non-CMEMS services linked below.

|Predictor(s)|Product|Variables|Reference period|Climatology|Source|
|-|-|-|-|-|-|
|iron, nitrate, phosphate|CMEMS Ocean Biogeochemistry Hindcast (BIO 001\_029)|`fe`, `no3`, `po4`|1993–2025|self-computed|Copernicus Marine|
|pCO₂|CMEMS observation-based carbon (Chau et al. 2024)|`spco2`|1993–2024|self-computed|Copernicus Marine|
|SST|HadISST v1.1|`sst`|1993–2025|self-computed|[Met Office Hadley Centre](https://www.metoffice.gov.uk/hadobs/hadisst/)|
|turbidity (Kd490)|Copernicus-GlobColour Ocean Colour (009\_104)|`KD490`|1997–2024|self-computed (monthly → annual → mean)|Copernicus Marine|
|salinity|CMEMS Ocean Physics Reanalysis (PHY 001\_030)|`so`|native|precomputed|Copernicus Marine|
|wave height|CMEMS Ocean Waves Reanalysis (WAV 001\_032)|`VHM0`|native|precomputed|Copernicus Marine|
|PAR|MODIS-Aqua L3|`par`|native|precomputed|[NASA OBPG](https://oceancolor.gsfc.nasa.gov/)|
|biogeographic region|Marine Ecoregions of the World|polygons|—|—|[UNEP-WCMC](https://data.unep-wcmc.org/datasets/38)|

Raw environmental data is licensed by its providers (Copernicus Marine, Met Office, NASA, UNEP-WCMC). The deposited climatologies are derived products redistributed under those licenses with attribution. Cite the original products (paper Table S1) when using either tier.



## Requirements

* **Fiji / ImageJ** for the capture macro.
* **Python 3.10+** for the extraction scripts: `pip install -r extraction/requirements.txt` (numpy, Pillow, scipy, scikit-image, matplotlib).
* **R 4.4+** for the analysis. Packages used across the scripts: dplyr, tidyr, terra, raster, ncdf4, sf, vegan, car, ape, nlme, ggplot2, ggrepel, scales, readr (plus optional maps, deldir, ggnewscale for some figures).



## Citation

If you use this code, please cite the manuscript (Nonemaker, Scripps Institution of Oceanography, UC San Diego). See [`LICENSE`](LICENSE) (MIT).

