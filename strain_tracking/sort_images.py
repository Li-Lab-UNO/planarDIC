"""
Sort CellScale images into one folder per deformation cycle.

Cycle boundaries come from Logs/*_PhaseTiming.csv. Each cycle runs from its
Stretch start to the next Stretch start (the last one ends at the End time).
Folders are named DefNN_<SetName>[_cNN] inside Images/. Any jpgs already in
Images/ subfolders are found and moved too, so re-running is safe.
A manifest (old path -> new path) is written to Images/sort_manifest.csv.

    ../.venv/bin/python sort_images.py --test N1739/test2 --dry-run
    ../.venv/bin/python sort_images.py --test N1739/test2
"""
import argparse
import csv
import glob
import os
import re
import shutil
from collections import Counter

import testcfg

FNAME_RE = re.compile(r"\.(\d{6})\.(\d{3})\.jpg$")


def clean(name):
    name = name.replace(" save cyle", "").replace(" save cycle", "").replace(" Save Cycle", " Save")
    return re.sub(r"[^A-Za-z0-9.\-]+", "_", name.strip()).strip("_")


def cycles(timing_csv):
    """List of (folder_name, t_start, t_end) from the PhaseTiming log."""
    rows = list(csv.DictReader(open(timing_csv)))
    stretch = [r for r in rows if r["Phase"] == "Stretch"]
    t_end = float([r for r in rows if r["Phase"] == "End"][0]["StartTime(s)"])
    n_per_set = Counter(r["Set"] for r in stretch)
    out = []
    for i, r in enumerate(stretch):
        t0 = float(r["StartTime(s)"])
        t1 = float(stretch[i + 1]["StartTime(s)"]) if i + 1 < len(stretch) else t_end + 1
        name = f"Def{i+1:02d}_{clean(r['SetName'])}"
        if n_per_set[r["Set"]] > 1:
            name += f"_c{int(r['Cycle']):02d}"
        out.append((name, t0, t1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    a.root = testcfg.get().root
    img_dir = os.path.join(a.root, "Images")
    timing = glob.glob(os.path.join(a.root, "Logs", "*_PhaseTiming.csv"))[0]
    cyc = cycles(timing)

    files = [p for p in glob.glob(os.path.join(img_dir, "**", "*.jpg"), recursive=True)
             if FNAME_RE.search(p)]
    moves, unassigned = [], []
    for p in files:
        m = FNAME_RE.search(p)
        t = int(m.group(1)) + int(m.group(2)) / 1000
        dest = next((n for n, t0, t1 in cyc if t0 - 1e-6 <= t < t1 - 1e-6), None)
        if dest is None:
            unassigned.append(p)
            continue
        new = os.path.join(img_dir, dest, os.path.basename(p))
        if os.path.abspath(p) != new:
            moves.append((p, new))

    counts = Counter(os.path.basename(os.path.dirname(n)) for _, n in moves)
    for n, t0, t1 in cyc:
        print(f"{n:32s} {t0:7.1f} - {min(t1, cyc[-1][2]):7.1f} s   {counts.get(n, 0):5d} to move")
    print(f"{len(files)} images found, {len(moves)} to move, {len(unassigned)} unassigned")
    if unassigned:
        print("  unassigned examples:", unassigned[:3])
    if a.dry_run:
        return

    if os.path.exists(os.path.join(img_dir, "sort_manifest.csv")):
        man = open(os.path.join(img_dir, "sort_manifest.csv"), "a", newline="")
    else:
        man = open(os.path.join(img_dir, "sort_manifest.csv"), "w", newline="")
        man.write("old_path,new_path\n")
    w = csv.writer(man)
    for old, new in moves:
        if os.path.exists(new):
            raise FileExistsError(new)
        os.makedirs(os.path.dirname(new), exist_ok=True)
        shutil.move(old, new)
        w.writerow([os.path.relpath(old, a.root), os.path.relpath(new, a.root)])
    man.close()
    # remove now-empty old subfolders (ignoring Finder .DS_Store files)
    for d in sorted({os.path.dirname(o) for o, _ in moves}, key=len, reverse=True):
        if d != img_dir and os.path.isdir(d):
            left = [f for f in os.listdir(d) if f != ".DS_Store"]
            if not left:
                for f in os.listdir(d):
                    os.remove(os.path.join(d, f))
                os.rmdir(d)
                print("removed empty folder", os.path.relpath(d, a.root))
    print("done; manifest at Images/sort_manifest.csv")


if __name__ == "__main__":
    main()
