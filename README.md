# Asparagopsis morphometrics pipeline

A semi-automated pipeline for extracting morphometric traits from suitable digitized macroalgae herbarium specimens via a Fiji/ImageJ macro, requiring manual tracing, and two Python scripts: threshold approval, then automated trait value computation.

The pipeline was built to quantify global morphological diversity in *Asparagopsis taxiformis* gametophytes; however, it is applicable to other taxa.

## Pipeline overview

Run each on a folder of specimen images:

**`capture.ijm`** (Fiji/ImageJ) — set scale via the attached ruler for each herbarium record. For up to three fronds per record: trace the frond outline, trace the length of the main axis, measure basal stipe width, and designate the starts and ends of branched regions. This script saves a crop of each frond and the traced coordinates. Cumulatively writes `manifest.csv`, `*_crop.png`, `*_mainaxis.txt`, `*_stipe.txt`, and `*_branched.txt`.

**`review.py`** — interactive thresholding for each frond crop, showing the cropped frond image with an initial threshold mask estimate overlaid in red and two sliders (saturation min, brightness min). Use the sliders to adjust values until only thallus is masked, then Accept (or Skip, flagging a frond as unusable). Approved cutoffs are written to `thresholds.csv`.

**`measure.py`** — fully automated portion; uses the designated masks from the previous part to compute values for the trait suite. Writes `measurements.csv`.

```
images/ ──capture.ijm──▶ manifest.csv + *_crop.png + *_mainaxis.txt
                          (+ *_stipe.txt, *_branched.txt)
        ──review.py────▶ thresholds.csv
        ──measure.py───▶ measurements.csv
```

## Traits measured

Surface area, perimeter, main-axis length, width suite (max / median / mean), circularity, solidity (area ÷ convex-hull area), perimeter-over-root-area, aspect ratio, stipe width, stipe-to-max-width ratio, stipe-to-length ratio, and branched fraction (length along the axis that is branched ÷ total length). A `qc_flags` column marks fronds whose measurements look implausible — an outlier max width (holdfast/rhizoid caught by the width march), a max width exceeding frond length, or a mask fragmented into many pieces — as a visual-check flag, never an auto-reject.

Thresholding is initialized by Otsu's method (Otsu 1979) and refined manually per frond in the review step.

## Requirements

- **Fiji / ImageJ** for the capture macro
- **Python 3.10+** with the packages in `requirements.txt`:

```
pip install -r requirements.txt
```

(numpy, Pillow, scipy, scikit-image, matplotlib)

## Usage

```bash
# 1. In Fiji: Plugins > Macros > Run… > capture.ijm
#    Choose the image folder, then the output folder when prompted.

# 2. Interactive threshold review
python review.py /path/to/output_folder

# 3. Measure
python measure.py /path/to/output_folder
```

Useful `measure.py` options:

| Option | Default | Purpose |
|---|---|---|
| `--max-gap-mm` | 1.5 | Max background gap the width march bridges |
| `--max-fill-hole-mm2` | 1.0 | Fill only interior holes smaller than this (threshold speckle) |
| `--bridge-tape-mm` | 0 (off) | Morphological close to reconnect branches severed by tape strips |
| `--overlay` | off | Save a per-frond diagnostic overlay image |

All three scripts support **resume**: re-running on the same folder skips work already completed (`capture.ijm` and `review.py` skip finished fronds; `measure.py` regenerates from what's present).

## Notes

- The `.gitignore` deliberately excludes specimen images, coordinate files, and result CSVs, so only the code is versioned. Point the scripts at a local data folder outside the repo.
- `capture.ijm` assumes a 10 mm scale line (`known_mm`) and up to 3 fronds per sheet (`expected_fronds`); both are constants at the top of the macro.

## License

MIT — see [LICENSE](LICENSE).
