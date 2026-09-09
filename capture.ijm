// ============================================================
// Asparagopsis morphometrics — CAPTURE macro (ImageJ/Fiji)  v2
// ============================================================
// Full COLOR capture — NO thresholding here (thallus stays visible while you
// trace; fixes the red-overlay problem). Thresholding happens later in Python.
//
// Per frond it saves:
//   {stem}_frond{n}_crop.png  : color crop of the frond's bounding box, with
//                               everything OUTSIDE your freehand outline set to
//                               pure black. Python thresholds the non-black part.
//   {stem}_frond{n}_mainaxis.txt: main-axis polyline in CROP-relative coords.
// Manifest: seq_index, filename, frond, px_per_mm
//
// WORKFLOW per image:
//   1. SET SCALE: draw a line over a known distance on the ruler -> OK
//   2. Per frond: (a) freehand around the frond -> OK
//                 (b) segmented main-axis line, full blade, double-click -> OK
// ============================================================

known_mm = 10;
expected_fronds = 3;

input_dir  = getDirectory("Choose image folder");
output_dir = getDirectory("Choose OUTPUT folder");

list = getFileList(input_dir);
list = Array.sort(list);

manifest = output_dir + "manifest.csv";

// --- RESUME SUPPORT ---
// If manifest.csv already exists, read which filenames are already captured so
// we can skip them. Otherwise create it with a header. Each frond appends its
// own row immediately (below), so stopping early keeps all completed work.
done_files = "";   // newline-joined list of already-captured filenames
if (File.exists(manifest)) {
    existing = File.openAsString(manifest);
    lines = split(existing, "\n");
    for (li = 1; li < lines.length; li++) {       // skip header row
        parts = split(lines[li], ",");
        if (parts.length >= 2) done_files = done_files + parts[1] + "\n";
    }
    print("Resuming: found existing manifest with prior work.");
} else {
    File.saveString("seq_index,filename,frond,px_per_mm\n", manifest);
}

seq = 0;

for (i = 0; i < list.length; i++) {
    filename = list[i];
    if (!(endsWith(toLowerCase(filename), ".jpg") ||
          endsWith(toLowerCase(filename), ".jpeg") ||
          endsWith(toLowerCase(filename), ".tif") ||
          endsWith(toLowerCase(filename), ".tiff") ||
          endsWith(toLowerCase(filename), ".png"))) continue;

    seq++;
    // skip if this filename is already in the manifest (resume)
    if (indexOf("\n" + done_files, "\n" + filename + "\n") >= 0) {
        print("[" + seq + "] " + filename + " - already done, skipping.");
        continue;
    }
    open(input_dir + filename);
    original = getTitle();
    stem = filename;
    if (lastIndexOf(stem, ".") > 0) stem = substring(stem, 0, lastIndexOf(stem, "."));
    print("[" + seq + "] " + filename);

    // ---- SCALE ----
    setTool("line");
    waitForUser("SET SCALE", "Draw a line over exactly " + known_mm + " mm on the ruler.\nThen click OK.");
    getLine(x1, y1, x2, y2, lw);
    if (x1 == -1) {
        print("  WARNING: no scale line for " + filename + " - skipping image.");
        close(original);
        continue;
    }
    px_len = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    px_per_mm = px_len / known_mm;
    print("  px_per_mm = " + px_per_mm);

    // ---- PER-FROND CAPTURE ----
    fronds_done = 0;
    for (frond = 1; frond <= expected_fronds; frond++) {

        // For frond 2 and 3, ask whether there's another frond on this sheet.
        // (Frond 1 is assumed — a sheet with zero fronds wouldn't be in the set.)
        if (frond > 1) {
            another = getBoolean("Sheet has another frond to measure?\n\n" +
                "Yes = measure frond " + frond + "\nNo = done with this sheet, go to next image");
            if (!another) {
                print("  (only " + fronds_done + " frond(s) on this sheet)");
                break;
            }
        }

        // (a) frond outline via the POLYGON tool (click points around the frond,
        //     double-click to close). Produces a closed area selection, same as
        //     freehand, so everything downstream (crop, Clear Outside) is unchanged.
        selectWindow(original);
        run("Select None");
        setTool("polygon");
        waitForUser("FROND " + frond + " - OUTLINE",
            "Click points around frond " + frond + " to encircle it,\nthen DOUBLE-CLICK to close the loop.\nThen click OK.");
        // selectionType 2 = polygon, 3 = freehand; accept either area selection
        if (selectionType() != 2 && selectionType() != 3) {
            print("  WARNING: no outline for frond " + frond + " - skipping.");
            continue;
        }
        getSelectionBounds(bx, by, bw, bh);
        roiManager("reset");
        roiManager("add");                 // outline now safe in slot 0

        // (b) main-axis polyline (outline is stored; drawing this won't lose it)
        selectWindow(original);
        setTool("polyline");
        waitForUser("FROND " + frond + " - MAIN AXIS",
            "Draw a segmented line down frond " + frond + "'s main axis,\nthrough the full blade. Double-click to finish.\nThen click OK.");
        getSelectionCoordinates(xs, ys);
        coord_path = output_dir + stem + "_frond" + frond + "_mainaxis.txt";
        mainaxis_text = "";
        for (k = 0; k < xs.length; k++) mainaxis_text = mainaxis_text + (xs[k]-bx) + "," + (ys[k]-by) + "\n";
        File.saveString(mainaxis_text, coord_path);

        // (b2) stipe-width line. Written EXACTLY like the branched block (which
        //      works): loop the selection coords directly, no intermediate
        //      variables, no arithmetic. Width is computed later in measure.py
        //      from these two saved points. Click 2 points across the stipe,
        //      double-click to finish.
        stipe_text = "";
        stipe_ok = false;
        while (!stipe_ok) {
            selectWindow(original);
            run("Select None");
            setTool("polyline");
            waitForUser("FROND " + frond + " - STIPE WIDTH",
                "With the SEGMENTED LINE tool, click one edge of the stipe,\nthen the other edge, then DOUBLE-CLICK to finish.\n(2 clicks across the stipe width.)\nThen click OK.");
            getSelectionCoordinates(stxs, stys);
            if (stxs.length < 2) {
                proceed = getBoolean("No line detected (need 2 points).\n\nYes = try again\nNo = skip stipe for this frond");
                if (!proceed) stipe_ok = true;
            } else {
                for (sk = 0; sk < stxs.length; sk++) stipe_text = stipe_text + (stxs[sk]-bx) + "," + (stys[sk]-by) + "\n";
                stipe_ok = true;
            }
        }
        stipe_path = output_dir + stem + "_frond" + frond + "_stipe.txt";
        File.saveString(stipe_text, stipe_path);

        // (b3) branched-region boundary points, captured as START/END pairs.
        //      Each click is its own prompt (robust: a single point selection
        //      can't be lost). After each pair, asks if there's another branched
        //      region (handles branched-stipe-branched). Skip rhizoid/holdfast by not
        //      marking it. Saved in CROP-relative coords; measure.py projects
        //      them onto the main axis.
        branched_text = "";
        more_branched = true;
        region = 0;
        while (more_branched) {
            region++;
            // START point
            selectWindow(original);
            run("Select None");
            setTool("point");
            waitForUser("FROND " + frond + " - BRANCHED REGION " + region + " START",
                "Click the START of branched region " + region + "\n(where branching begins). Then click OK.");
            if (selectionType() != 10) {
                print("  frond " + frond + ": no start point for region " + region + " - stopping branched capture.");
                more_branched = false;
            } else {
                getSelectionCoordinates(sxs, sys);
                branched_text = branched_text + (sxs[0]-bx) + "," + (sys[0]-by) + "\n";

                // END point
                selectWindow(original);
                run("Select None");
                setTool("point");
                waitForUser("FROND " + frond + " - BRANCHED REGION " + region + " END",
                    "Click the END of branched region " + region + "\n(where branching stops). Then click OK.");
                if (selectionType() != 10) {
                    print("  frond " + frond + ": no end point for region " + region + " - dropping this region.");
                    // remove the orphan start we just added
                    branched_text = "";  // safest: discard partial; user can redo frond
                    more_branched = false;
                } else {
                    getSelectionCoordinates(exs, eys);
                    branched_text = branched_text + (exs[0]-bx) + "," + (eys[0]-by) + "\n";
                    // ask about another region
                    more_branched = getBoolean("Another separate branched region on this frond?\n\nYes = mark another start/end pair\nNo = done with branched regions");
                }
            }
        }
        branched_path = output_dir + stem + "_frond" + frond + "_branched.txt";
        File.saveString(branched_text, branched_path);

        // (c) build the color crop with outside-of-outline blacked out:
        //     duplicate whole image, recall outline, clear outside, crop to bbox.
        selectWindow(original);
        run("Select None");
        run("Duplicate...", "title=work");     // color copy of whole sheet
        setBackgroundColor(0, 0, 0);           // ensure cleared region is BLACK
        roiManager("select", 0);               // recall freehand outline
        run("Clear Outside");                  // everything outside ROI -> black
        run("Crop");                           // crop to the ROI bounding box
        run("Select None");
        crop_path = output_dir + stem + "_frond" + frond + "_crop.png";
        saveAs("PNG", crop_path);
        close();                               // close the work/crop image

        File.append(seq + "," + filename + "," + frond + "," + px_per_mm, manifest);
        fronds_done++;
        print("  frond " + frond + ": crop + " + xs.length + "-pt main axis + stipe + branched saved");
    }

    close(original);
}

print("CAPTURE DONE (or stopped). Completed work is saved in: " + output_dir);
print("Re-running on the same folder will resume where you left off.");
print("Next: run review.py, then measure.py.");
