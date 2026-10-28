"""
Sort CellScale images into one folder per deformation cycle.

Cycle boundaries come from Logs/*_PhaseTiming.csv. Each cycle runs from its
Stretch start to the next Stretch start (the last one ends at the End time).
Folders are named DefNN_<SetName>[_cNN] inside Images/. Any jpgs already in
Images/ subfolders are found and moved too, so re-running is safe.
A manifest (old path -> new path) is written to Images/sort_manifest.csv.

Afterwards every cycle is checked against the tester's own files: images, Data.csv rows and Coords1.csv
rows missing (gaps in their times), and images without a Data.csv or Coords1.csv row at their time. The
table goes to Images/sort_summary.csv. With --review (runDIC.sh: review = yes) it is shown in a window
that asks whether to proceed; "No, stop this test" exits with an error, so runDIC.sh skips the test's
later steps.

    .venv/bin/python strain_tracking/sort_images.py --test N1739/test2 --dry-run
    .venv/bin/python strain_tracking/sort_images.py --test N1739/test2 --review
"""
import argparse
import csv
import glob
import os
import re
import shutil
import sys
import textwrap
from collections import Counter

import numpy as np
import pandas as pd

import review
import testcfg

FNAME_RE = re.compile(r"\.(\d{6})\.(\d{3})\.jpg$")


def clean(name):
    """Folder-safe protocol name: 'Equibiaxial Save Cycle' -> 'Equibiaxial', '1:0.75 ratio' -> '1-0.75_ratio'
    style. 'save cycle' (also 'save cyle') is a label of the tester's only."""
    name = re.sub(r"\s+save\s+cyc?le$", "", name, flags=re.I)
    return re.sub(r"[^A-Za-z0-9.\-]+", "_", name.strip()).strip("_")


def kept(img_dir, n):
    """The Images/DefNN_* folder of deformation n if the test was sorted before: it keeps its name, so its
    tracks, named after it, still match (e.g. Def19_Equibiaxial_Save from earlier versions)."""
    if not os.path.isdir(img_dir):
        return None
    m = [d for d in os.listdir(img_dir) if d.startswith(f"Def{n:02d}_") and os.path.isdir(os.path.join(img_dir, d))]
    return m[0] if len(m) == 1 else None


def cycles(timing_csv):
    """List of (folder_name, t_start, t_end) from the PhaseTiming log."""
    return [c[:3] for c in phases(timing_csv)]


def phases(timing_csv):
    """List of (folder_name, t_start, t_end, t_recover) from the PhaseTiming log: t_recover is where
    the cycle's loading ends and unloading (its Recover phase) starts."""
    rows = list(csv.DictReader(open(timing_csv)))
    stretch = [r for r in rows if r["Phase"] == "Stretch"]
    recover = {(r["Set"], r["Cycle"]): float(r["StartTime(s)"]) for r in reversed(rows) if r["Phase"] == "Recover"}
    t_end = float([r for r in rows if r["Phase"] == "End"][0]["StartTime(s)"])
    n_per_set = Counter(r["Set"] for r in stretch)
    out = []
    for i, r in enumerate(stretch):
        t0 = float(r["StartTime(s)"])
        t1 = float(stretch[i + 1]["StartTime(s)"]) if i + 1 < len(stretch) else t_end + 1
        name = f"Def{i+1:02d}_{clean(r['SetName'])}"
        if n_per_set[r["Set"]] > 1:
            name += f"_c{int(r['Cycle']):02d}"
        out.append((name, t0, t1, recover.get((r["Set"], r["Cycle"]), t1)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--review", action="store_true", help="show the table and ask whether to proceed")
    a = ap.parse_args()
    if a.review:
        review.enable()
    a.root = testcfg.get().root
    img_dir = os.path.join(a.root, "Images")
    timing = glob.glob(os.path.join(a.root, "Logs", "*_PhaseTiming.csv"))[0]
    cyc = [(kept(img_dir, i + 1) or n, t0, t1) for i, (n, t0, t1) in enumerate(cycles(timing))]
    t_rec = [c[3] for c in phases(timing)]

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
    if moves:
        move(moves, a.root, img_dir)
    check(a.root, img_dir, cyc, t_rec, len(unassigned))


def move(moves, root, img_dir):
    """Move the images, recording each move in Images/sort_manifest.csv."""
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
        w.writerow([os.path.relpath(old, root), os.path.relpath(new, root)])
    man.close()
    # remove now-empty old subfolders (ignoring Finder .DS_Store files)
    for d in sorted({os.path.dirname(o) for o, _ in moves}, key=len, reverse=True):
        if d != img_dir and os.path.isdir(d):
            left = [f for f in os.listdir(d) if f != ".DS_Store"]
            if not left:
                for f in os.listdir(d):
                    os.remove(os.path.join(d, f))
                os.rmdir(d)
                print("removed empty folder", os.path.relpath(d, root))
    print("moved; manifest at Images/sort_manifest.csv")


def gaps(t, t_start, t_end, dt):
    """Times (s) of the frames missing from times t between t_start and t_end, at one every dt."""
    u = np.r_[t_start - dt, t, t_end]
    k = np.maximum(0, np.round(np.diff(u) / dt) - 1).astype(int)
    return np.array([u[i] + dt * j for i in np.nonzero(k)[0] for j in range(1, k[i] + 1)])


def where(miss, dt, t_rec):
    """Where the missing frames are: their times, as runs, and whether they fall in loading or unloading."""
    runs, start = [], 0
    for i in range(1, len(miss) + 1):
        if i == len(miss) or miss[i] - miss[i - 1] > 1.5 * dt:
            a, b = miss[start], miss[i - 1]
            runs.append(f"{a:.1f} s" if i - 1 == start else f"{a:.1f}-{b:.1f} s")
            start = i
    if len(runs) > 3:
        runs = runs[:3] + ["..."]
    n = int(np.sum(miss < t_rec - 1e-6))
    phase = ("all loading" if n == len(miss) else "all unloading" if n == 0
             else f"{n} loading, {len(miss) - n} unloading")
    return f": {', '.join(runs)} ({phase})"


def check(root, img_dir, cyc, t_recs, n_unassigned):
    """Completeness of every cycle from the tester's own files (a frame needs an image, a Data.csv row and
    a Coords1.csv row); writes Images/sort_summary.csv and, with review on, asks whether to proceed."""
    times = [frame_times(os.path.join(img_dir, n)) for n, _, _ in cyc]
    dt = float(np.median(np.concatenate([np.diff(t) for t in times])))           # image interval (s)
    tD = np.unique(pd.read_csv(glob.glob(os.path.join(root, "*Data.csv"))[0], encoding="latin1",
                               usecols=["Time_S"])["Time_S"].to_numpy(float))
    dtD = float(np.median(np.diff(tD)))
    coords = sorted(glob.glob(os.path.join(root, "Coords*.csv")))
    tC = np.unique(np.round(pd.read_csv(coords[0], skiprows=1)["Time(S)"].to_numpy(float), 1)) if coords else None
    dtC = float(np.median(np.diff(tC))) if coords else np.nan
    rows = []
    for (name, t0, t1), t_rec, t in zip(cyc, t_recs, times):
        r = dict(folder=name, images=len(t), first_image_s=np.nan, last_image_s=np.nan, cycle_length_s=round(t1 - t0, 1),
                 images_missing=0, data_rows_missing=0, coords_rows_missing=np.nan, images_unmatched=0, note="")
        rows.append(r)
        if not len(t):
            r["note"] = "NO IMAGES"
            continue
        r["first_image_s"], r["last_image_s"] = round(t[0], 1), round(t[-1], 1)
        t_end = min(t1, t[-1] + dt)
        within = lambda u: u[(u >= t0 - 1e-6) & (u < t_end - 1e-6)]  # noqa: E731
        miss = gaps(t, t0, t_end, dt)
        missD = gaps(within(tD), t0, t_end, dtD)
        missC = gaps(within(tC), t0, t_end, dtC) if coords else np.array([])
        r["images_missing"], r["data_rows_missing"] = len(miss), len(missD)
        tr = np.round(t, 1)
        unmatched = ~np.isin(tr, np.round(tD, 1))
        if coords:
            r["coords_rows_missing"] = len(missC)
            unmatched |= ~np.isin(tr, tC)
        r["images_unmatched"] = int(unmatched.sum())
        same = lambda m: len(m) == len(miss) and np.allclose(np.round(m, 1), np.round(miss, 1))  # noqa: E731
        p = []
        if len(miss):
            p.append(f"{len(miss)} image(s) missing{where(miss, dt, t_rec)}")
        for m, f, step in ((missD, "Data.csv", dtD), (missC, "Coords1.csv", dtC)):
            if len(m):
                p.append(f"{len(m)} {f} row(s) missing" + (" (the same frames)" if len(miss) and same(m)
                                                           else where(m, step, t_rec)))
        if r["images_unmatched"]:
            p.append(f"{r['images_unmatched']} image(s) without a Data/Coords row")
        r["note"] = "; ".join(p)
    rows[-1]["cycle_length_s"] = np.nan                       # the last cycle runs to the End time
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(img_dir, "sort_summary.csv"), index=False)
    flag = S["note"] != ""
    print(f"{len(S)} folders, {S['images'].sum()} images (one every {dt:.2f} s); {n_unassigned} outside every cycle; "
          f"{'gaps below' if flag.any() else 'images, Data.csv and Coords1.csv complete'} -> Images/sort_summary.csv")
    for _, r in S[flag].iterrows():
        print(f"  {r['folder']}: {r['note']}")
    if review.enabled() and show(S, flag, root, dt, n_unassigned) == "No, stop this test":
        sys.exit("stopped after sorting (see Images/sort_summary.csv)")


def frame_times(d):
    """Sorted image times (s) in folder d (none if it doesn't exist)."""
    if not os.path.isdir(d):
        return np.array([])
    return np.sort([int(m[1]) + int(m[2]) / 1000 for m in map(FNAME_RE.search, os.listdir(d)) if m])


def show(S, flag, root, dt, n_unassigned):
    """The table in a window; returns the answer."""
    fig = review.figure("runDIC: sorted images", size=(13, 8), toolbar=False)
    ax = fig.add_axes([0.02, 0.23, 0.96, 0.74])
    ax.axis("off")
    num = lambda v: "" if pd.isna(v) else f"{int(v)}"  # noqa: E731
    cells = [[r["folder"], r["images"], "" if pd.isna(r["first_image_s"]) else f"{r['first_image_s']:.1f}",
              "" if pd.isna(r["last_image_s"]) else f"{r['last_image_s']:.1f}", num(r["images_missing"]),
              num(r["data_rows_missing"]), num(r["coords_rows_missing"]),
              "\n".join(textwrap.fill(q, 62, subsequent_indent="    ") for q in r["note"].split("; ") if q)]
             for _, r in S.iterrows()]
    tab = ax.table(cellText=cells, loc="upper center", cellLoc="center",
                   colLabels=["Folder", "Images", "First\nimage (s)", "Last\nimage (s)", "Images\nmissing",
                              "Data.csv\nrows missing", "Coords1.csv\nrows missing", "Note"],
                   colWidths=[0.165, 0.055, 0.065, 0.065, 0.065, 0.085, 0.09, 0.41])
    tab.auto_set_font_size(False)
    tab.set_fontsize(10)
    tab.scale(1, 1.45)
    for (i, j), c in tab.get_celld().items():
        if i == 0:
            c.set_text_props(weight="bold")
            c.set_height(c.get_height() * 1.8)
        elif j in (0, 7):
            c.set_text_props(ha="left")
            c.PAD = 0.02
            if j == 7:
                c.set_fontsize(9.5)
        if i > 0 and flag.iloc[i - 1]:
            c.set_text_props(color="red", weight="bold")
            c.set_height(c.get_height() * (cells[i - 1][7].count("\n") + 1))   # one line per problem
    fig.text(0.02, 0.19, f"{os.path.basename(root).replace('_', ' ')}: {S['images'].sum()} images, one every "
             f"{dt:.2f} s; {n_unassigned} images outside every cycle (before the first or after the last).", fontsize=review.TEXT_FONT)
    if flag.any():
        q = ("Some folders are missing frames (red). A few missing frames are fine (the camera drops frames now and "
             "then),\nand frames missing during unloading don't affect the loading curves; a folder with no images "
             "or many missing\nmeans the test or the folders need a look. Proceed with this test?")
    else:
        q = "Every folder has all its images, Data.csv rows and Coords1.csv rows. Proceed with this test?"
    return review.ask(fig, q, ["Proceed", "No, stop this test"])


if __name__ == "__main__":
    main()
