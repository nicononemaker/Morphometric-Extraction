#!/usr/bin/env python3
"""
Asparagopsis morphometrics — REVIEW tool (interactive thresholding).

Walks every *_crop.png in the capture folder ONE AT A TIME. For each crop it
shows the color image with the current threshold mask overlaid in red, plus
two sliders (Saturation min, Brightness min) starting at an Otsu auto-guess.
Drag until only thallus is red, then click Accept. Click Skip to flag a frond
as unusable.

Approved values are written to thresholds.csv:
    stem, frond, sat_min, bright_min, status   (status = ok | skipped)

measure.py reads thresholds.csv and applies these exact cutoffs.

Run:
    python review.py /path/to/capture_folder
Controls:
    drag sliders -> live mask update
    Accept (or press 'a')  -> save + next
    Skip   (or press 's')  -> flag unusable + next
    close window           -> stop (progress saved so far)
"""
import os, sys, glob, csv
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("TkAgg")   # interactive backend; falls back handled below
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

def otsu(values):
    """Otsu threshold on a 1-D array of 0..255 ints. Returns int cutoff."""
    hist, _ = np.histogram(values, bins=256, range=(0, 256))
    total = values.size
    if total == 0:
        return 128
    sum_all = np.dot(np.arange(256), hist)
    wB = 0.0; sumB = 0.0; best_var = -1.0; thresh = 0
    for t in range(256):
        wB += hist[t]
        if wB == 0:
            continue
        wF = total - wB
        if wF == 0:
            break
        sumB += t * hist[t]
        mB = sumB / wB
        mF = (sum_all - sumB) / wF
        var = wB * wF * (mB - mF) ** 2
        if var > best_var:
            best_var = var; thresh = t
    return int(thresh)

def hsv_channels(rgb):
    """Return saturation and value (brightness) as 0..255 arrays from an RGB
    image array (H,W,3) float 0..1."""
    mx = rgb.max(-1); mn = rgb.min(-1)
    sat = np.where(mx == 0, 0, (mx - mn) / np.where(mx == 0, 1, mx))
    val = mx
    return (sat * 255).astype(np.int16), (val * 255).astype(np.int16)

def load_crop(path):
    rgb = np.asarray(Image.open(path).convert("RGB")).astype(float) / 255.0
    # "outside outline" was blacked out in capture; treat near-black as outside
    nonblack = rgb.max(-1) > (10/255.0)
    sat, val = hsv_channels(rgb)
    return rgb, sat, val, nonblack

def make_mask(sat, val, nonblack, sat_min, bright_min):
    m = (sat >= sat_min) & nonblack
    if bright_min > 0:
        m = m & (val >= bright_min)
    return m

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python review.py /path/to/capture_folder")
    folder = sys.argv[1]
    crops = sorted(glob.glob(os.path.join(folder, "*_crop.png")))
    if not crops:
        sys.exit("No *_crop.png files found in " + folder)

    out_path = os.path.join(folder, "thresholds.csv")
    # resume support: load already-approved stems/fronds
    done = {}
    if os.path.exists(out_path):
        with open(out_path) as f:
            r = csv.reader(f); next(r, None)
            for row in r:
                if len(row) >= 5:
                    done[(row[0], row[1])] = row
    results = dict(done)

    # parse stem/frond from filename
    def parse(path):
        base = os.path.basename(path)[:-len("_crop.png")]
        # stem is everything before the last _frond{n}
        idx = base.rfind("_frond")
        stem = base[:idx]
        frond = base[idx + len("_frond"):]
        return stem, frond

    queue = [c for c in crops if parse(c) not in results]
    if not queue:
        print("All crops already reviewed. Delete thresholds.csv to redo.")
        return
    print(f"{len(queue)} crops to review ({len(results)} already done).")

    state = {"i": 0}

    fig, ax = plt.subplots(figsize=(8, 9))
    plt.subplots_adjust(left=0.08, bottom=0.22, right=0.97, top=0.95)
    sat_ax = plt.axes([0.15, 0.12, 0.65, 0.03])
    br_ax  = plt.axes([0.15, 0.07, 0.65, 0.03])
    acc_ax = plt.axes([0.15, 0.005, 0.3, 0.05])
    skip_ax= plt.axes([0.55, 0.005, 0.3, 0.05])

    s_sat = Slider(sat_ax, "Sat min", 0, 255, valinit=28, valstep=1)
    s_br  = Slider(br_ax,  "Bright min", 0, 255, valinit=0, valstep=1)
    b_acc = Button(acc_ax, "Accept (a)")
    b_skip= Button(skip_ax, "Skip (s)")

    cur = {}

    def load_current():
        path = queue[state["i"]]
        stem, frond = parse(path)
        rgb, sat, val, nonblack = load_crop(path)
        # Otsu init on non-black pixels only
        sat_guess = otsu(sat[nonblack]) if nonblack.any() else 28
        cur.update(dict(path=path, stem=stem, frond=frond, rgb=rgb,
                        sat=sat, val=val, nonblack=nonblack, imh=None))
        ax.clear()                       # fresh axis for the new (differently sized) crop
        # set sliders WITHOUT triggering redraw (eventson=False), then draw once
        s_sat.eventson = False; s_br.eventson = False
        s_sat.set_val(sat_guess)
        s_br.set_val(0)
        s_sat.eventson = True;  s_br.eventson = True
        redraw()
        ax.set_title(f"[{state['i']+1}/{len(queue)}]  {stem}  frond {frond}\n"
                     f"Otsu sat guess = {sat_guess}", fontsize=10)
        ax.axis("off")
        fig.canvas.draw_idle()

    def redraw(_=None):
        m = make_mask(cur["sat"], cur["val"], cur["nonblack"],
                      int(s_sat.val), int(s_br.val))
        disp = (cur["rgb"] * 255).astype(np.uint8).copy()
        disp[m] = [255, 0, 0]
        if cur.get("imh") is None:
            cur["imh"] = ax.imshow(disp)
        else:
            cur["imh"].set_data(disp)
        fig.canvas.draw_idle()

    def accept(_=None):
        results[(cur["stem"], cur["frond"])] = [
            cur["stem"], cur["frond"], int(s_sat.val), int(s_br.val), "ok"]
        save_all()
        advance()

    def skip(_=None):
        results[(cur["stem"], cur["frond"])] = [
            cur["stem"], cur["frond"], int(s_sat.val), int(s_br.val), "skipped"]
        save_all()
        advance()

    def save_all():
        with open(out_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["stem", "frond", "sat_min", "bright_min", "status"])
            for row in results.values():
                w.writerow(row)

    def advance():
        state["i"] += 1
        if state["i"] >= len(queue):
            print("Review complete. Wrote", out_path)
            plt.close(fig)
            return
        load_current()

    def on_key(event):
        if event.key == "a": accept()
        elif event.key == "s": skip()

    s_sat.on_changed(redraw)
    s_br.on_changed(redraw)
    b_acc.on_clicked(accept)
    b_skip.on_clicked(skip)
    fig.canvas.mpl_connect("key_press_event", on_key)

    load_current()
    plt.show()

if __name__ == "__main__":
    main()
