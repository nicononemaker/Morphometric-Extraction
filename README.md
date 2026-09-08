# Asparagopsis morphometrics pipeline

A semi-automated pipeline for extracting morphometric traits from digitized
herbarium specimens of *Asparagopsis taxiformis* gametophytes. It combines a
Fiji/ImageJ macro for manual tracing with two Python scripts for interactive
threshold approval and trait computation.

The pipeline was built for a global morphological diversity study in the Smith
Lab at Scripps Institution of Oceanography, but the tools are general to any
project measuring flattened, branching thallus morphology from photographed
sheets.

## Pipeline overview

Three stages, run in order on a folder of specimen images:

1. **`capture.ijm`** (Fiji/ImageJ) — set the scale from the ruler, then per
   frond trace the outline, midrib, stipe-width line, and bushy-region
   boundaries. Saves a color crop (outside the outline blacked out) plus the
   traced coordinates. No thresholding happens here, so the thallus stays
   visible while you trace.

2. **`review.py`** — walks each crop one at a time, showing the color image with
   the current threshold mask overlaid in red and two sliders (saturation min,
   brightness min) initialized from an Otsu auto-guess. Drag until only thallus
   is masked, then Accept (or Skip to flag a frond as unusable). Approved cutoffs
   are written to `thresholds.csv`.

3. **`measure.py`** — rebuilds each frond's mask from its crop using the approved
   cutoffs and computes the trait suite in millimetres. Writes `measurements.csv`.

```
images/ ──capture.ijm──▶ manifest.csv + *_crop.png + *_midrib.txt
                          (+ *_stipe.txt, *_bushy.txt)
        ──review.py────▶ thresholds.csv
        ──measure.py───▶ measurements.csv
```

## Traits measured

Per frond, in millimetres where applicable: projected area, perimeter, midrib
length, the perpendicular width suite (max / median / mean / min and their
ratio), bounding-box dimensions, circularity, solidity (area ÷ convex-hull
area), a perimeter-over-root-area branching index, aspect ratio, stipe width,
and a manual bushy/stipe partition (bushy length and fraction, segment count,
width-profile shape descriptors). A `qc_flags` column marks fronds that need a
second look.

Thresholding is initialized by Otsu's method (Otsu 1979) and refined manually
per frond in the review step.

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
