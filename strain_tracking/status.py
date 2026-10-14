"""
Status of every test folder: what has been set up, tracked and exported, read from the files on disk.

    .venv/bin/python strain_tracking/status.py                 every sample
    .venv/bin/python strain_tracking/status.py N1739 N1757     some samples (folder prefixes)
    .venv/bin/python strain_tracking/status.py --notes         also print each test's notes from test_config.json

Columns: sorted (DefNN folders), config (test_config.json), scale (px/mm; excluded tine pairs),
widths (w_theta / w_z tine spans, mm), H (thickness, mm, from the photos or the workbook), tracked (deformations with a full-domain
track among Def11-Def19; "loading only" marks track_protocols.py --loading-only tracks), J/mean (peak_J_over_mean figure), exports (deformations with
*_fields.npz).
"""
import argparse
import glob
import json
import os
import re
import sys
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import testcfg  # noqa: E402

DEFS = [f"Def{i}" for i in range(11, 20)]


def branch(npz):
    """'loading' for a track_protocols.py --loading-only track, else 'full'."""
    try:
        with zipfile.ZipFile(npz) as z:
            if "branch.npy" in z.namelist():
                return str(np.load(z.open("branch.npy")))
    except (OSError, zipfile.BadZipFile, ValueError):
        pass
    return "full"


def def_range(found):
    """'Def11-19' for a full set, else the list."""
    if not found:
        return "none"
    if found == DEFS:
        return "Def11-19"
    return ", ".join(found)


def row(td):
    T = testcfg.Test(td)
    img = T.images
    folders = sorted(d for d in os.listdir(img) if re.match(r"Def\d\d_", d)) if os.path.isdir(img) else []
    res = T.results
    tracks = {d: [p for f in folders if f.startswith(d + "_") for p in [os.path.join(res, f + ".npz")]
                  if os.path.exists(p)] for d in DEFS}
    tracked = [d for d in DEFS if tracks[d]]
    loading = [d for d in tracked if branch(tracks[d][0]) == "loading"]
    exported = [d for d in DEFS if glob.glob(os.path.join(res, d + "_*_fields.npz"))]
    cfg = T.config
    sc = json.load(open(T.scale_file)) if os.path.exists(T.scale_file) else None
    scale = "-"
    if sc:
        scale = f"{sc['px_per_mm']:.1f}"
        if cfg.get("calibration_exclude"):
            scale += f" (excl. {', '.join(cfg['calibration_exclude'])})"
    w = cfg.get("widths_mm", {}).get("tine_span")
    th = cfg.get("thickness_mm", {})
    H = th.get("H")
    Hsrc = "" if not H else " photos" if str(th.get("source", "")).startswith("measure_thickness") \
        else " workbook" if ".xls" in str(th.get("source", "")) else " other"
    return dict(test=f"{T.sample.split('_')[0]}/{T.name.split('_')[0]}",
                sorted=f"{len(folders)} folders" if folders else "no",
                config="yes" if cfg else "no", scale=scale,
                widths=f"{w[0]:.2f} / {w[1]:.2f}" if w else "-", H=f"{H:.2f}{Hsrc}" if H else "-",
                tracked=def_range(tracked) + (" (loading only)" if loading and loading == tracked
                                              else f" ({def_range(loading)} loading only)" if loading else ""),
                jmean="yes" if os.path.exists(os.path.join(res, "peak_J_over_mean_Def11-19.png")) else "no",
                exports=def_range(exported), notes=cfg.get("notes", []), label=cfg.get("sample", ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", help="sample folder prefixes (default: all)")
    ap.add_argument("--notes", action="store_true", help="print each test's label and notes")
    a = ap.parse_args()
    samples = testcfg.samples()
    if a.samples:
        samples = [s for s in samples if any(s.startswith(p) for p in a.samples)]
    rows = [row(td) for s in samples
            for td in sorted(glob.glob(os.path.join(testcfg.ROOT, s, "test[0-9]*_*"))) if os.path.isdir(td)]
    if not rows:
        sys.exit("no test folders found")
    cols = [("test", "test"), ("sorted", "sorted"), ("config", "config"), ("scale", "scale px/mm"),
            ("widths", "w_theta / w_z mm"), ("H", "H mm"), ("tracked", "tracked"), ("jmean", "J/mean"),
            ("exports", "exports")]
    width = {k: max(len(h), *(len(str(r[k])) for r in rows)) for k, h in cols}
    print("  ".join(h.ljust(width[k]) for k, h in cols))
    for r in rows:
        print("  ".join(str(r[k]).ljust(width[k]) for k, _ in cols))
    if a.notes:
        for r in rows:
            print(f"\n{r['test']}: {r['label']}")
            for n in r["notes"]:
                print(f"  - {n}")


if __name__ == "__main__":
    main()
