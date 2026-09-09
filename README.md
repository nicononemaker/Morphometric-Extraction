# Asparagopsis morphometrics pipeline

A semi-automated pipeline for extracting morphometric traits from suitable digitized macroalgae herbarium specimens via a 
Fiji/ImageJ macro, requring manual tracing, and two Python scripts, threshold approval, then automated trait value 
computation.

The pipeline was built to quantify global morphological diversity in Asparagopsis taxiformis gametophytes however, this 
pipeline is applicable to other taxa.

## Pipeline overview

Run each on a folder of specimen images:

1. **`capture.ijm`** (Fiji/ImageJ) — set scale via attatched ruler for each herbarium record. For up to three fronds per
   record: Trace frond outline, trace the length of the main axis, measure basal stipe-width, and designate starts and
   ends of branching regions. This script will save a crop of each frond and the traced coordinates. Cumulativly writes
   `manifest.csv`, `crop.png`, `midrib.txt`, `stipe.txt`, and `bushy.txt`.

3. **`review.py`** — Interactive thresholding for each frond crop, showing the croped frond image with an initail
   threshold mask estimate overlaid in red and two sliders (saturation min, brightness min). Use sliders to adjust values
   until only thallus is masked, then Accept (or skip, flaging a frond as unusable). Approved cutoffs are written to
   `thresholds.csv`.

5. **`measure.py`** — Fully automated portion, uses deisgnated masks from previous part to computes values for the trait
   suite. Writes `measurements.csv`.

```
images/ ──capture.ijm──▶ manifest.csv + *_crop.png + *_midrib.txt
                          (+ *_stipe.txt, *_bushy.txt)
        ──review.py────▶ thresholds.csv
        ──measure.py───▶ measurements.csv
```

## Traits measured

Surface area, perimeter, length, width suite (max / median / mean), circularity, solidity (area ÷ convex-hull
area), perimeter-over-root-area, aspect ratio, stipe width, and branched fraction (length along axis branched ÷ total
length) A `qc_flag` marks fronds whose measurements look implausible

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

All three scripts support **resume**: re-running on the same folder skips work
already completed (`capture.ijm` and `review.py` skip finished fronds;
`measure.py` regenerates from what's present).

## Notes

- The `.gitignore` deliberately excludes specimen images, coordinate files, and
  result CSVs, so only the code is versioned. Point the scripts at a local data
  folder outside the repo.
- `capture.ijm` assumes a 10 mm scale line (`known_mm`) and up to 3 fronds per
  sheet (`expected_fronds`); both are constants at the top of the macro.

## License

MIT — see [LICENSE](LICENSE).
