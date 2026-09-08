#!/usr/bin/env python3
"""
Asparagopsis morphometrics — MEASUREMENT script (final stage of the pipeline).

Reads the capture + review outputs from a folder:
    manifest.csv                      (seq_index, filename, frond, px_per_mm)
    thresholds.csv                    (stem, frond, sat_min, bright_min, status)
    {stem}_frond{n}_crop.png          (color crop; outside-outline is blacked out)
    {stem}_frond{n}_midrib.txt        (midrib polyline, "x,y" per line)
    {stem}_frond{n}_stipe.txt         (optional: 2 points across the stipe)
    {stem}_frond{n}_bushy.txt         (optional: bushy-region start/end points)

For each frond it rebuilds the binary mask from the color crop using the
saturation/brightness cutoffs that review.py approved (fronds marked "skipped"
in review are ignored), then computes, in millimetres:
    area_mm2, perimeter_mm, length_midrib_mm,
    max/median/mean/min_width_mm (perpendicular-to-midrib, gap-tolerant),
    width_max_min_ratio, bbox_w_mm, bbox_h_mm,
    circularity, solidity, perim_over_sqrt_area (branching-complexity index),
    aspect_ratio, stipe_width_mm and the manual bushy-partition traits,
    n_components, qc_flags

Writes: measurements.csv

Run capture.ijm (Fiji) then review.py before this script.

Usage:
    python measure.py /path/to/output_folder [--max-gap-mm 1.5]
        [--max-fill-hole-mm2 1.0] [--bridge-tape-mm 0] [--overlay]
"""
import argparse, glob, os, sys, csv
import numpy as np
from PIL import Image
from scipy import ndimage

# ----------------------- width algorithm (validated) -----------------------
def in_mask(m, x, y):
    xi, yi = int(round(x)), int(round(y))
    if xi < 0 or yi < 0 or xi >= m.shape[1] or yi >= m.shape[0]:
        return False
    return m[yi, xi]

def crossing_width(m, px, py, dx, dy, max_reach, max_gap_px, step=1.0):
    """Full perpendicular crossing through (px,py) along unit (dx,dy).
    Bridges background gaps shorter than max_gap_px; a longer gap stops the march.
    Returns (width_px, endpoint_plus, endpoint_minus) or None if start off-tissue."""
    if not in_mask(m, px, py):
        return None
    def reach(sign):
        last, gap, d = 0.0, 0.0, step
        while d <= max_reach:
            if in_mask(m, px + sign*dx*d, py + sign*dy*d):
                last, gap = d, 0.0
            else:
                gap += step
                if gap > max_gap_px:
                    break
            d += step
        return last
    rp, rm = reach(1), reach(-1)
    return rp + rm, (px + dx*rp, py + dy*rp), (px - dx*rm, py - dy*rm)

def max_perp_width(mask, verts, max_gap_px, samples_per_seg=60):
    """Max full perpendicular crossing along the midrib polyline. Returns
    (width_px, best_segment_endpoints, profile_widths, profile_dists) where
    profile_dists is each sample's cumulative arc-distance from the base."""
    best, best_seg, profile, dists = 0.0, None, [], []
    diag = float(np.hypot(*mask.shape))
    cum = 0.0
    for i in range(len(verts) - 1):
        x0, y0 = verts[i]; x1, y1 = verts[i+1]
        sx, sy = x1 - x0, y1 - y0
        L = np.hypot(sx, sy)
        if L == 0:
            continue
        pdx, pdy = -sy / L, sx / L            # perpendicular unit vector
        for t in np.linspace(0, 1, samples_per_seg, endpoint=False):
            cx, cy = x0 + t*(x1-x0), y0 + t*(y1-y0)
            r = crossing_width(mask, cx, cy, pdx, pdy, max_reach=diag, max_gap_px=max_gap_px)
            if r:
                profile.append(r[0])
                dists.append(cum + t * L)     # arc distance from base
                if r[0] > best:
                    best, best_seg = r[0], (r[1], r[2], (cx, cy))
        cum += L
    return best, best_seg, profile, dists

# ----------------------- shape descriptors -----------------------
def polyline_length_px(verts):
    d = np.diff(verts, axis=0)
    return float(np.hypot(d[:,0], d[:,1]).sum())

def measure_frond(mask, verts, px_per_mm, max_gap_px, stipe_w_mm=-1.0,
                  bushy_pts=None, want_seg=False, keep_all=False):
    # By default keep the largest connected component (the frond body). When
    # keep_all=True (tape-bridging already cleaned the mask upstream), use the
    # whole mask so branches severed by tape and reconnected aren't discarded.
    lbl, n = ndimage.label(mask)
    if n == 0:
        return (None, None) if want_seg else None
    if keep_all:
        body = mask.astype(bool)
    else:
        sizes = ndimage.sum(np.ones_like(lbl), lbl, range(1, n+1))
        body = (lbl == (int(np.argmax(sizes)) + 1))

    area_px = float(body.sum())
    # solidity = area / convex-hull area (how filled-in vs spidery/branchy the
    # shape is; ~1.0 = solid blade, lower = more dissected/branched). Same
    # definition ImageJ uses. Convex hull of the body pixel coordinates.
    solidity = np.nan
    try:
        from scipy.spatial import ConvexHull
        ys, xs = np.nonzero(body)
        if xs.size >= 3:
            # expand each pixel to its 4 corners so the hull encloses the full
            # pixel footprint (not just centers) -> solidity stays <= 1.0
            corners = np.array([[0,0],[1,0],[0,1],[1,1]], dtype=float)
            pts = np.column_stack([xs, ys]).astype(float)
            pts = (pts[:, None, :] + corners[None, :, :]).reshape(-1, 2)
            hull_area = ConvexHull(pts).volume   # 2-D: .volume is polygon area
            if hull_area > 0:
                solidity = min(area_px / hull_area, 1.0)
    except Exception:
        solidity = np.nan
    # perimeter via marching edges (count boundary pixels' exposed sides)
    # use scipy: perimeter ~ sum over boundary; approximate with contour length
    eroded = ndimage.binary_erosion(body)
    boundary = body & ~eroded
    perim_px = float(boundary.sum())  # rough; refined below by contour if available
    try:
        from skimage import measure as skmeasure
        contours = skmeasure.find_contours(body.astype(float), 0.5)
        if contours:
            c = max(contours, key=len)
            d = np.diff(c, axis=0)
            perim_px = float(np.hypot(d[:,0], d[:,1]).sum())
    except Exception:
        pass

    ys, xs = np.nonzero(body)
    bbox_w_px = float(xs.max() - xs.min() + 1)
    bbox_h_px = float(ys.max() - ys.min() + 1)

    width_px, seg, profile, dists = max_perp_width(body, verts, max_gap_px=max_gap_px)
    length_px = polyline_length_px(verts)
    prof = np.array(profile, dtype=float) if profile else np.array([np.nan])
    min_w_px = float(np.nanmin(prof))
    mean_w_px = float(np.nanmean(prof))
    median_w_px = float(np.nanmedian(prof))

    mm = 1.0 / px_per_mm

    # ---------- STIPE WIDTH + PROFILE-SHAPE + MANUAL BUSHY PARTITION ----------
    bushy_metrics = dict(
        stipe_width_mm="NA", stipe_to_max_width_ratio="NA",
        stipe_length_mm="NA", bushy_length_mm="NA", bushy_fraction="NA",
        n_bushy_segments="NA", width_shape_factor="NA", peak_position_frac="NA")

    # stipe width as its own metric (from the manual stipe line)
    if stipe_w_mm is not None and stipe_w_mm > 0:
        bushy_metrics["stipe_width_mm"] = round(float(stipe_w_mm), 4)
        if width_px > 0:
            bushy_metrics["stipe_to_max_width_ratio"] = round(
                float(stipe_w_mm * px_per_mm) / width_px, 4)

    if len(profile) >= 3 and np.isfinite(prof).any():
        d_arr = np.array(dists, dtype=float)
        order = np.argsort(d_arr)
        d_arr = d_arr[order]; w_arr = prof[order]
        total_len_px = d_arr[-1] - d_arr[0] if d_arr[-1] > d_arr[0] else length_px

        # (3) width shape factor = mean width / max width (circular/spiky vs oval/rect)
        shape_factor = (np.nanmean(w_arr) / np.nanmax(w_arr)
                        if np.nanmax(w_arr) > 0 else np.nan)
        # (5) peak position as fraction of length from base
        peak_idx = int(np.nanargmax(w_arr))
        peak_pos = (d_arr[peak_idx] - d_arr[0]) / total_len_px if total_len_px > 0 else np.nan
        bushy_metrics["width_shape_factor"] = round(float(shape_factor), 4)
        bushy_metrics["peak_position_frac"] = round(float(peak_pos), 4)

        # ----- MANUAL bushy partition from dropped points -----
        # bushy_pts: array of (x,y) in crop coords, paired start,end,(start,end..)
        if bushy_pts is not None and len(bushy_pts) >= 2:
            # project each boundary point onto the midrib -> arc-distance from base
            def project_dist(pt):
                best_d, best_dist = None, np.inf
                cum = 0.0
                for i in range(len(verts)-1):
                    a = verts[i]; b = verts[i+1]
                    ab = b - a; L2 = float(ab @ ab)
                    if L2 == 0:
                        continue
                    t = np.clip(float((pt - a) @ ab) / L2, 0, 1)
                    proj = a + t*ab
                    dd = float(np.hypot(*(pt - proj)))
                    if dd < best_dist:
                        best_dist = dd
                        best_d = cum + t*np.sqrt(L2)
                    cum += np.sqrt(L2)
                return best_d if best_d is not None else 0.0

            bd = sorted(project_dist(np.asarray(p, float)) for p in bushy_pts)
            # pair consecutively: (0,1),(2,3),...; odd trailing start -> to tip
            bushy_len_px = 0.0
            n_seg = 0
            for i in range(0, len(bd), 2):
                start = bd[i]
                end = bd[i+1] if i+1 < len(bd) else total_len_px
                bushy_len_px += max(0.0, end - start)
                n_seg += 1
            bushy_len_px = min(bushy_len_px, total_len_px)
            stipe_len_px = max(0.0, total_len_px - bushy_len_px)
            bushy_metrics.update(
                stipe_length_mm=round(stipe_len_px*mm, 4),
                bushy_length_mm=round(bushy_len_px*mm, 4),
                bushy_fraction=round(bushy_len_px/total_len_px, 4) if total_len_px>0 else "NA",
                n_bushy_segments=n_seg,
            )
    # -------------------------------------------------------------------------
    area_mm2 = area_px * mm * mm
    perim_mm = perim_px * mm
    width_mm = width_px * mm
    length_mm = length_px * mm
    min_w_mm = min_w_px * mm
    mean_w_mm = mean_w_px * mm
    median_w_mm = median_w_px * mm

    circ = (4*np.pi*area_px) / (perim_px**2) if perim_px > 0 else np.nan

    # ---- QC FLAGS: name metrics that look implausible (semicolon-separated) ----
    flags = []
    # max width far larger than typical crossing width => likely a holdfast/
    # rhizoid crossing got sampled (the failure seen in overlays). Threshold 6x.
    if median_w_px > 0 and (width_px / median_w_px) > 6.0:
        flags.append("max_width")
    # max width larger than the frond is long is almost always wrong
    if length_px > 0 and (width_px / length_px) > 1.0:
        flags.append("max_width_gt_length")
    # fragmented mask (multiple components) can corrupt area/shape
    if n > 3:
        flags.append("fragmented")
    qc_flags = ";".join(flags) if flags else ""

    result = {
        "area_mm2": round(area_mm2, 4),
        "perimeter_mm": round(perim_mm, 4),
        "length_midrib_mm": round(length_mm, 4),
        "max_width_mm": round(width_mm, 4),
        "median_width_mm": round(median_w_mm, 4),
        "mean_width_mm": round(mean_w_mm, 4),
        "min_width_mm": round(min_w_mm, 4),
        "width_max_min_ratio": round(width_mm/min_w_mm, 4) if min_w_mm > 0 else "NA",
        "bbox_w_mm": round(bbox_w_px*mm, 4),
        "bbox_h_mm": round(bbox_h_px*mm, 4),
        "circularity": round(circ, 4) if np.isfinite(circ) else "NA",
        "solidity": round(solidity, 4) if np.isfinite(solidity) else "NA",
        "perim_over_sqrt_area": round(perim_mm/np.sqrt(area_mm2), 4) if area_mm2 > 0 else "NA",
        "aspect_ratio": round(length_mm/median_w_mm, 4) if median_w_mm > 0 else "NA",
        "stipe_width_mm": bushy_metrics["stipe_width_mm"],
        "stipe_to_max_width_ratio": bushy_metrics["stipe_to_max_width_ratio"],
        "stipe_length_mm": bushy_metrics["stipe_length_mm"],
        "bushy_length_mm": bushy_metrics["bushy_length_mm"],
        "bushy_fraction": bushy_metrics["bushy_fraction"],
        "n_bushy_segments": bushy_metrics["n_bushy_segments"],
        "width_shape_factor": bushy_metrics["width_shape_factor"],
        "peak_position_frac": bushy_metrics["peak_position_frac"],
        "n_components": int(n),
        "qc_flags": qc_flags,
    }
    if want_seg:
        return result, (body, verts, seg, profile, bushy_pts)
    return result

# ----------------------- driver -----------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="folder with manifest.csv + masks + midribs")
    ap.add_argument("--max-gap-mm", type=float, default=1.5,
                    help="bridge perpendicular gaps shorter than this (mm)")
    ap.add_argument("--overlay", action="store_true",
                    help="save a PNG per frond showing mask + midrib + widest crossing")
    ap.add_argument("--max-fill-hole-mm2", type=float, default=1.0,
                    help="fill only interior holes SMALLER than this (mm^2) — patches "
                         "threshold speckle while leaving between-branch gaps open so "
                         "area = true tissue. Set 0 to disable hole-filling entirely.")
    ap.add_argument("--bridge-tape-mm", type=float, default=0.0,
                    help="bridge thin gaps up to ~this width (mm) via morphological "
                         "close — reconnects frond pieces severed by tape strips. "
                         "0 = off (default). Try 0.5-1.5 on taped specimens.")
    args = ap.parse_args()

    man_path = os.path.join(args.folder, "manifest.csv")
    if not os.path.exists(man_path):
        sys.exit("No manifest.csv in " + args.folder)

    rows = []
    with open(man_path) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            seq, filename, frond, ppm = parts[0], parts[1], parts[2], float(parts[3])
            rows.append((seq, filename, frond, ppm))

    # load approved thresholds from review.py
    thr_path = os.path.join(args.folder, "thresholds.csv")
    if not os.path.exists(thr_path):
        sys.exit("No thresholds.csv — run review.py first to set thresholds.")
    thresholds = {}
    with open(thr_path) as tf:
        tr = csv.reader(tf); next(tr, None)
        for row in tr:
            if len(row) >= 5:
                thresholds[(row[0], row[1])] = (int(row[2]), int(row[3]), row[4])

    def hsv_channels(rgb):
        mx = rgb.max(-1); mn = rgb.min(-1)
        sat = np.where(mx == 0, 0, (mx - mn) / np.where(mx == 0, 1, mx))
        return (sat*255).astype(np.int16), (mx*255).astype(np.int16)

    out_path = os.path.join(args.folder, "measurements.csv")
    cols = ["seq_index","filename","frond","px_per_mm","area_mm2","perimeter_mm",
            "length_midrib_mm","max_width_mm","median_width_mm","mean_width_mm",
            "min_width_mm","width_max_min_ratio","bbox_w_mm","bbox_h_mm",
            "circularity","solidity","perim_over_sqrt_area","aspect_ratio",
            "stipe_width_mm","stipe_to_max_width_ratio",
            "stipe_length_mm","bushy_length_mm","bushy_fraction","n_bushy_segments",
            "width_shape_factor","peak_position_frac","n_components","qc_flags"]
    with open(out_path, "w") as out:
        out.write(",".join(cols) + "\n")
        for seq, filename, frond, ppm in rows:
            stem = os.path.splitext(filename)[0]
            crop_p = os.path.join(args.folder, f"{stem}_frond{frond}_crop.png")
            mid_p  = os.path.join(args.folder, f"{stem}_frond{frond}_midrib.txt")
            key = (stem, frond)
            if not (os.path.exists(crop_p) and os.path.exists(mid_p)):
                print(f"  MISSING files for {stem} frond {frond} — skipping")
                continue
            if key not in thresholds:
                print(f"  NO THRESHOLD for {stem} frond {frond} — run review.py; skipping")
                continue
            sat_min, bright_min, status = thresholds[key]
            if status == "skipped":
                print(f"  [{seq}] {stem} frond {frond}: flagged unusable in review — skipping")
                continue
            # build the mask from the color crop using the approved cutoffs
            rgb = np.asarray(Image.open(crop_p).convert("RGB")).astype(float)/255.0
            nonblack = rgb.max(-1) > (10/255.0)
            sat_ch, val_ch = hsv_channels(rgb)
            mask = (sat_ch >= sat_min) & nonblack
            if bright_min > 0:
                mask = mask & (val_ch >= bright_min)
            # Fill ONLY small interior holes (threshold speckle), leaving large
            # between-branch gaps open so area = true tissue, not trapped water/paper.
            # A "hole" = background enclosed by tissue; fill it only if its area is
            # below --max-fill-hole-mm2 (converted to px via this frond's ppm).
            if args.max_fill_hole_mm2 > 0:
                filled = ndimage.binary_fill_holes(mask)
                holes = filled & (~mask)              # pixels that would be filled
                if holes.any():
                    hlbl, hn = ndimage.label(holes)
                    hsizes = ndimage.sum(np.ones_like(hlbl), hlbl, range(1, hn+1))
                    max_hole_px = args.max_fill_hole_mm2 * (ppm ** 2)
                    # re-fill only holes smaller than the threshold
                    keep = np.zeros_like(mask)
                    for hi in range(1, hn+1):
                        if hsizes[hi-1] <= max_hole_px:
                            keep |= (hlbl == hi)
                    mask = mask | keep
            # Bridge thin gaps from tape strips (morphological close: dilate then
            # erode). Only runs if --bridge-tape-mm > 0, so clean specimens are
            # untouched. Reconnects pieces the tape severed without much changing
            # the overall outline. Then keep ALL frond-sized components (not just
            # the single largest) so severed branches aren't discarded.
            if args.bridge_tape_mm > 0:
                r = max(1, int(round(args.bridge_tape_mm * ppm / 2.0)))
                struct = ndimage.generate_binary_structure(2, 1)
                mask = ndimage.binary_closing(mask, structure=struct, iterations=r)
                # keep components above a small size (drop specks, keep frond parts)
                clbl, cn = ndimage.label(mask)
                if cn > 1:
                    csizes = ndimage.sum(np.ones_like(clbl), clbl, range(1, cn+1))
                    biggest = csizes.max()
                    keepmask = np.zeros_like(mask)
                    for ci in range(1, cn+1):
                        # keep any component at least 5% of the largest (real frond
                        # pieces) — discards isolated speckle/label noise
                        if csizes[ci-1] >= 0.05 * biggest:
                            keepmask |= (clbl == ci)
                    mask = keepmask
            verts = np.loadtxt(mid_p, delimiter=",", ndmin=2)
            # load manual bushy-boundary points if present
            bushy_p = os.path.join(args.folder, f"{stem}_frond{frond}_bushy.txt")
            bushy_pts = None
            if os.path.exists(bushy_p) and os.path.getsize(bushy_p) > 0:
                try:
                    bushy_pts = np.loadtxt(bushy_p, delimiter=",", ndmin=2)
                except Exception:
                    bushy_pts = None
            # load stipe line endpoints if present (for overlay + width)
            stipe_p = os.path.join(args.folder, f"{stem}_frond{frond}_stipe.txt")
            stipe_line = None
            stipe_w = -1.0
            if os.path.exists(stipe_p) and os.path.getsize(stipe_p) > 0:
                try:
                    sl = np.loadtxt(stipe_p, delimiter=",", ndmin=2)
                    if len(sl) >= 2 and np.isfinite(sl).all():
                        # use first and last point as the stipe endpoints
                        stipe_line = np.array([sl[0], sl[-1]])
                        stipe_px = float(np.hypot(stipe_line[1,0]-stipe_line[0,0],
                                                  stipe_line[1,1]-stipe_line[0,1]))
                        stipe_w = stipe_px / ppm
                except Exception:
                    stipe_line = None
                    stipe_w = -1.0
            max_gap_px = args.max_gap_mm * ppm
            keep_all = args.bridge_tape_mm > 0
            if args.overlay:
                res = measure_frond(mask, verts, ppm, max_gap_px,
                                    stipe_w_mm=stipe_w, bushy_pts=bushy_pts, want_seg=True,
                                    keep_all=keep_all)
                m, segdata = res if res else (None, None)
            else:
                m = measure_frond(mask, verts, ppm, max_gap_px,
                                  stipe_w_mm=stipe_w, bushy_pts=bushy_pts,
                                  keep_all=keep_all)
            if m is None:
                print(f"  EMPTY mask for {stem} frond {frond} — skipping")
                continue
            vals = [seq, filename, frond, round(ppm,4)] + [m[c] for c in cols[4:]]
            out.write(",".join(str(v) for v in vals) + "\n")

            # sanity check: the reported max_width should equal the profile's own
            # peak. If it's well below, something is wrong in the width reporting.
            flag = ""
            try:
                if args.overlay and segdata and segdata[3]:
                    peak_px = max(segdata[3])
                    reported_px = m["max_width_mm"] * ppm
                    if reported_px < 0.95 * peak_px:
                        flag = "  <-- WARNING: reported max_width below profile peak (reporting bug)"
                # also flag a near-flat profile (midrib may not cross the blade)
                if m["width_max_min_ratio"] != "NA" and float(m["width_max_min_ratio"]) < 1.2:
                    flag += "  <-- NOTE: width nearly constant; check midrib crosses the bushy part"
            except Exception:
                pass

            if args.overlay and segdata and segdata[2]:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                body, vv, seg, profile, bpts = segdata
                (axp,ayp),(bx,by),(cx,cy) = seg
                fig, ax = plt.subplots(1, 3, figsize=(15, 7),
                                       gridspec_kw={"width_ratios":[3,3,2]})
                # ---- panel 0: full frond (color crop) with measured AREA shaded ----
                ax[0].imshow(rgb)
                # translucent cyan over the measured region (the largest component
                # = exactly the pixels counted for area_mm2). Lets you verify the
                # threshold captured the real thallus and nothing spurious.
                area_layer = np.zeros((*body.shape, 4), dtype=float)
                area_layer[body] = [0.1, 0.8, 0.9, 0.35]   # cyan, 35% opacity
                ax[0].imshow(area_layer)
                ax[0].plot([], [], "s", color=(0.1,0.8,0.9), alpha=0.6,
                           label=f"measured area {m['area_mm2']}mm²")
                ax[0].plot(vv[:,0], vv[:,1], "c-o", ms=3, lw=1, label="midrib")
                ax[0].plot([axp,bx],[ayp,by], "r-", lw=3, label=f"max width {m['max_width_mm']}mm")
                ax[0].plot(cx,cy,"y*",ms=14)
                if bpts is not None and len(bpts) >= 1:
                    bp = np.atleast_2d(np.asarray(bpts, float))
                    for idx in range(len(bp)):
                        is_start = (idx % 2 == 0)
                        ax[0].plot(bp[idx,0], bp[idx,1], "o", ms=11,
                                   mfc=("lime" if is_start else "orange"), mec="black",
                                   label=("bushy start" if (is_start and idx==0) else
                                          ("bushy end" if (not is_start and idx==1) else None)))
                # stipe line on the full view (blue)
                if stipe_line is not None:
                    ax[0].plot(stipe_line[:,0], stipe_line[:,1], "b-", lw=2.5, label="stipe width")
                ax[0].legend(fontsize=8); ax[0].set_title(
                    f"{stem} frond {frond}  (bushy frac={m['bushy_fraction']})"
                    + (f"\nFLAGGED: {m['qc_flags']}" if m['qc_flags'] else ""),
                    color=("red" if m['qc_flags'] else "black")); ax[0].axis("off")
                # ---- panel 1: width profile ----
                if profile:
                    ax[1].plot(profile, color="purple")
                    ax[1].axhline(max(profile), ls="--", color="r")
                    ax[1].set_title("crossing-width profile (px)")
                # ---- panel 2: zoomed stipe inset ----
                if stipe_line is not None:
                    sx_c = float(stipe_line[:,0].mean()); sy_c = float(stipe_line[:,1].mean())
                    span = float(np.hypot(*(stipe_line[1]-stipe_line[0])))
                    pad = max(span * 1.5, 25)
                    ax[2].imshow(rgb)
                    ax[2].plot(stipe_line[:,0], stipe_line[:,1], "b-", lw=3)
                    ax[2].plot(stipe_line[:,0], stipe_line[:,1], "bo", ms=6)
                    ax[2].set_xlim(sx_c - pad, sx_c + pad)
                    ax[2].set_ylim(sy_c + pad, sy_c - pad)   # inverted y for image coords
                    ax[2].set_title(f"stipe width = {m['stipe_width_mm']} mm")
                    ax[2].axis("off")
                else:
                    ax[2].text(0.5, 0.5, "no stipe line", ha="center", va="center")
                    ax[2].axis("off")
                ovp = os.path.join(args.folder, f"{stem}_frond{frond}_overlay.png")
                plt.tight_layout(); plt.savefig(ovp, dpi=80, bbox_inches="tight"); plt.close()

            print(f"  [{seq}] {stem} frond {frond}: "
                  f"area={m['area_mm2']} len={m['length_midrib_mm']} "
                  f"maxW={m['max_width_mm']} medW={m['median_width_mm']} "
                  f"minW={m['min_width_mm']} (n_comp={m['n_components']}){flag}")
    print("Wrote", out_path)

if __name__ == "__main__":
    main()
