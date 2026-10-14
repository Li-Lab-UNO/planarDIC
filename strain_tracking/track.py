"""
Track a grid of material points through a CellScale image sequence, referenced to its own first
frame. (For Def11-19 against the test's common reference, use track_protocols.py.)
The grid defaults to the full domain between the tine tips (--grid tine_region, output
<name>[_tag].npz); --grid standard uses roi_px.standard (output <name>_standard[_tag].npz).
An existing common-reference track (from track_protocols.py) is never overwritten.

Track a whole deformation folder; results go to <test>/DIC_results/:
    ../.venv/bin/python track.py --test N1739/test1 --def Def01
    ../.venv/bin/python track.py --test N1739/test1 --def Def13 --radius 40 --step 20 --tag s81_g20

Or an arbitrary time window across all of Images/ (saved in DIC_results/ as t<t0>-<t1>...npz):
    ../.venv/bin/python track.py --test N1739/test1 --t0 893 --t1 913
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

T = testcfg.get()
IMAGES = T.images

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--def", dest="defn", help="deformation folder or its prefix, e.g. Def13")
ap.add_argument("--t0", type=float, help="reference frame time (s); default: first frame")
ap.add_argument("--t1", type=float, help="last frame time (s); default: last frame")
ap.add_argument("--grid", choices=tuple(D.GRID_SUFFIX), default="tine_region",
                help="grid: the full domain between the tine tips (default) or the standard region")
ap.add_argument("--roi", nargs=4, type=float, default=None,
                metavar=("X0", "Y0", "X1", "Y1"),
                help="grid extent in the reference frame (px); default: the test's roi_px for --grid")
ap.add_argument("--step", type=float, default=D.DEFAULT_STEP, help="grid spacing (px)")
ap.add_argument("--radius", type=int, default=D.TrackParams.radius,
                help="subset half-width (px); subset = 2*radius+1")
ap.add_argument("--zncc-min", type=float, default=D.TrackParams.zncc_min)
ap.add_argument("--mask-saturated", action=argparse.BooleanOptionalAction,
                default=T.config.get("dic", {}).get("mask_saturated", False),
                help="leave saturated glint pixels out of the matching (default: the test's dic.mask_saturated)")
ap.add_argument("--tag", default="", help="suffix for the output name, e.g. r30")
ap.add_argument("--out", help="output .npz (default: <test>/DIC_results/<name><grid suffix>[_tag].npz)")
a = ap.parse_args()
a.roi = a.roi or T.setting("roi_px", a.grid)

t0 = -np.inf if a.t0 is None else a.t0 - 1e-6
t1 = np.inf if a.t1 is None else a.t1 + 1e-6
if a.defn:
    name = T.folder(a.defn)                                       # Def1 == Def01
    times, paths = D.find_frames([os.path.join(IMAGES, name)], t0, t1)
else:
    if a.t0 is None or a.t1 is None:
        sys.exit("give --def, or both --t0 and --t1")
    times, paths = D.find_frames([IMAGES], t0, t1, recursive=True)
    name = f"t{a.t0:g}-{a.t1:g}"
if not paths:
    sys.exit("no frames found")
if a.out is None:
    stem = name + D.GRID_SUFFIX[a.grid] + (f"_{a.tag}" if a.tag else "")
    a.out = os.path.join(T.results, stem + ".npz")
if os.path.exists(a.out) and "X_start" in np.load(a.out).files:
    sys.exit(f"{a.out} is a common-reference track (track_protocols.py); choose --tag or --out")
print(f"{len(paths)} frames, t = {times[0]:.1f} .. {times[-1]:.1f} s")
gaps = np.nonzero(np.diff(times) > 0.21)[0]
if gaps.size:
    print("  WARNING: missing frames after t =", times[gaps])

X0, shape = D.make_grid(a.roi, a.step)
print(f"grid {shape[1]} x {shape[0]} = {len(X0)} points, subset {2*a.radius+1} px")
prm = D.TrackParams(radius=a.radius, zncc_min=a.zncc_min, mask_saturated=a.mask_saturated)
t = time.time()
res = D.track_sequence(paths, X0, prm)
print(f"done in {time.time()-t:.0f} s; keyframes at t =", times[res["keyframes"]])

os.makedirs(os.path.dirname(a.out), exist_ok=True)
np.savez_compressed(a.out, X0=X0, shape=np.array(shape), times=times,
                    paths=np.array([os.path.basename(p) for p in paths]),
                    radius=a.radius, step=a.step, roi=np.array(a.roi), grid=a.grid,
                    mask_saturated=prm.mask_saturated, sat_level=prm.sat_level, sat_dilate=prm.sat_dilate, **res)
print("saved", a.out)
