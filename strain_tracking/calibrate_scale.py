"""
Infer the image scale (px per mm) from the rake motion.

The tine tips are tracked in every frame of a deformation (sub-pixel, along each
tine's axis). The change in tip-to-tip separation (px) is regressed against the
change in rake separation reported by the tester (X/Y_Size_um in *Data.csv).
X (left/right tines) and Y (top/bottom tines) give independent estimates.

The tines are pulled along their own axis, which they barely stretch, so their
tips follow the rakes. Checks:
  * per tine pair (outer tines splay sideways, shortening their footprint slightly)
  * low-force vs full range (tine bending under load would bend the curve)
  * loading vs unloading (backlash / lag)

    .venv/bin/python strain_tracking/calibrate_scale.py --test N1739/test1 --def Def01
    .venv/bin/python strain_tracking/calibrate_scale.py --test N1739/test1 --def Def01 Def19

Tine pairs that disagree with the rest are left out by an outlier rule (outliers(): more than 6% from the
median of all pairs, or Def01 and Def19 more than 5% apart, or not tracked throughout), and listed in
the output. calibration_exclude in test_config.json (e.g. ["Y3", "Y5"]) or --exclude replaces the
rule; --no-auto turns it off.
Approximate rest-state tine tips come from the test's test_config.json
(calibration_tines_px). Writes <test folder>/scale_calibration.json (used by later
scripts) and a diagnostic plot, <test folder>/DIC_results/scale_calibration.png.
"""
import argparse
import glob
import json
import os
import sys

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

TC = testcfg.get()
TEST_ROOT, IMAGES = TC.root, TC.images

def tines():
    """Approximate tine tips in the rest-state frame (px): (axis position, tip coordinate), per rake, from
    calibration_tines_px. The rod extends from the tip towards the image edge in direction `out`. Read when
    needed, not at import, so make_config.py can import find_tip before the test has a config."""
    return {side: dict(horizontal=v["horizontal"], out=v["out"], tips=[tuple(t) for t in v["tips"]])
            for side, v in TC.setting("calibration_tines_px").items()}


def find_tip(g, c, t, horizontal, out):
    """Refine centerline c and tip coordinate t of one tine in image g."""
    # centerline: darkest line across the rod, 15-45 px behind the tip
    across = np.arange(-12, 13)
    along = t + out * np.arange(15, 46)
    if horizontal:
        X, Y = np.meshgrid(along, c + across)
        cp = cv2.remap(g, X.astype(np.float32), Y.astype(np.float32), cv2.INTER_LINEAR).mean(1)
    else:
        Y, X = np.meshgrid(along, c + across)
        cp = cv2.remap(g, X.astype(np.float32), Y.astype(np.float32), cv2.INTER_LINEAR).mean(1)
    i = int(np.argmin(cp))
    if 0 < i < len(cp) - 1:
        den = cp[i - 1] - 2 * cp[i] + cp[i + 1]
        i = i + (0.5 * (cp[i - 1] - cp[i + 1]) / den if den > 0 else 0)
    c = c + across[0] + i
    # axial profile along the centerline (+-2 px), from rod into tissue
    s = np.arange(-40, 41)
    a = t - out * s
    w = np.arange(-2, 3)
    if horizontal:
        X, Y = np.meshgrid(a, c + w)
        prof = cv2.remap(g, X.astype(np.float32), Y.astype(np.float32), cv2.INTER_LINEAR).mean(0)
    else:
        Y, X = np.meshgrid(a, c + w)
        prof = cv2.remap(g, X.astype(np.float32), Y.astype(np.float32), cv2.INTER_LINEAR).mean(0)
    rod = np.median(prof[(s >= -40) & (s <= -15)])
    tis = np.median(prof[(s >= 12) & (s <= 40)])
    mid = 0.5 * (rod + tis)
    above = prof > mid if tis > rod else prof < mid
    j = np.nonzero(above & (s > -15))[0]
    if j.size == 0 or j[0] == 0:
        return c, np.nan, tis - rod
    j = j[0]
    f = (mid - prof[j - 1]) / (prof[j] - prof[j - 1])
    s_tip = s[j - 1] + f
    return c, t - out * s_tip, tis - rod


def track_tips(paths):
    """Tip coordinate of every tine in every frame -> dict side -> (T, 5) array."""
    TINES = tines()
    state = {k: [list(p) for p in v["tips"]] for k, v in TINES.items()}
    tips = {k: np.full((len(paths), len(v["tips"])), np.nan) for k, v in TINES.items()}
    for k, p in enumerate(paths):
        g = cv2.GaussianBlur(cv2.imread(p, cv2.IMREAD_GRAYSCALE).astype(np.float32), (0, 0), 0.8)
        for side, spec in TINES.items():
            for i, (c, t) in enumerate(state[side]):
                for _ in range(2):                       # refine twice (first frame: coarse start)
                    c, t_new, contrast = find_tip(g, c, t, spec["horizontal"], spec["out"])
                    if np.isfinite(t_new):
                        t = t_new
                state[side][i] = [c, t]
                tips[side][k, i] = t if contrast > 30 else np.nan
    return tips


OUTLIER_DEV = 0.06     # a tine pair more than this far from the median of all pairs is left out ...
OUTLIER_JUMP = 0.05    # ... and so is one whose deformations disagree by more than this (fractions of the median)


def outliers(T):
    """Tine pairs to leave out: a pair whose tips were not tracked in every deformation, whose scale (mean over
    the deformations) is more than OUTLIER_DEV from the median of all pairs, or whose deformations disagree by
    more than OUTLIER_JUMP of that median. The camera's pixels are square, so every pair should give the same
    scale; a pair that doesn't has tip tracking errors, or tines that bend or slip. Returns (list, reasons)."""
    ndef = T.deformation.nunique()
    ref = float(np.nanmedian(T.px_per_mm))
    out = {}
    for (ax, pr), g in T.groupby(["axis", "pair"]):
        v = g.px_per_mm.to_numpy(float)
        e = f"{ax}{pr}"
        if len(v) < ndef or not np.all(np.isfinite(v)):
            out[e] = "not tracked in every deformation"
        elif abs(v.mean() / ref - 1) > OUTLIER_DEV:
            out[e] = f"{v.mean():.2f} px/mm, {100 * (v.mean() / ref - 1):+.1f}% from the median {ref:.2f}"
        elif len(v) > 1 and np.ptp(v) / ref > OUTLIER_JUMP:
            out[e] = f"deformations disagree: {', '.join(f'{x:.2f}' for x in v)} px/mm"
    return sorted(out), out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--def", dest="defs", nargs="+", default=["Def01", "Def19"])
    ap.add_argument("--exclude", nargs="+", default=None, metavar="AXIS+PAIR",
                    help="tine pairs to leave out of the result, e.g. Y3 Y5; default: calibration_exclude of "
                         "test_config.json, else the outlier rule")
    ap.add_argument("--no-auto", action="store_true", help="don't apply the outlier rule (exclude nothing by default)")
    ap.add_argument("--low-force", type=float, default=150.0,
                    help="mN; fits restricted to frames below this rake force are reported too")
    a = ap.parse_args()
    source = "command line"
    if a.exclude is None and "calibration_exclude" in TC.config:
        a.exclude, source = list(TC.config["calibration_exclude"]), "calibration_exclude in test_config.json"

    Dd = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "*Data.csv"))[0], encoding="latin1")
    Dd = Dd.drop_duplicates("Time_S").set_index("Time_S")
    rows, plot = [], []
    for dn in a.defs:
        folder = TC.folder(dn)
        times, paths = D.find_frames([os.path.join(IMAGES, folder)])
        print(f"{folder}: tracking 20 tine tips in {len(paths)} frames")
        tips = track_tips(paths)
        dd = Dd.reindex(np.round(times, 1))
        for axis, (lo, hi), size in (("X", ("left", "right"), "X_Size_um"), ("Y", ("top", "bottom"), "Y_Size_um")):
            sep = tips[hi] - tips[lo]                              # (T, 5) pair separations
            rake = (dd[size] - dd[size].iloc[0]).to_numpy(float)   # change in rake separation (um)
            force = dd[f"{axis}_Force_mN"].to_numpy(float)
            dsep = sep - sep[0]
            load = np.arange(len(times)) <= int(np.nanargmax(rake))
            for pair in range(sep.shape[1]):
                y = dsep[:, pair]
                for label, m in (("all", np.isfinite(y)),
                                 ("loading", np.isfinite(y) & load),
                                 ("unloading", np.isfinite(y) & ~load),
                                 (f"F<{a.low_force:g}mN", np.isfinite(y) & load & (force < a.low_force))):
                    if m.sum() < 10:
                        continue
                    k, b = np.polyfit(rake[m], y[m], 1)
                    res = y[m] - (k * rake[m] + b)
                    rows.append(dict(deformation=folder, axis=axis, pair=pair + 1, subset=label,
                                     px_per_mm=1000 * k, rms_px=np.sqrt(np.mean(res ** 2)), n=int(m.sum())))
            plot.append((folder, axis, rake, dsep, force))

    T = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    piv = T[T.subset == "all"].pivot_table(index=["deformation", "axis"], columns="pair", values="px_per_mm")
    print("\npx per mm, each tine pair (1 = top/left ... 5 = bottom/right), all frames:")
    print(piv.round(2).to_string())
    summ = T.groupby(["axis", "subset"]).px_per_mm.agg(["mean", "std", "count"])
    print("\nsummary over pairs and deformations:")
    print(summ.round(3).to_string())
    print("\nfit residual RMS (px), all frames: median %.2f, max %.2f" %
          (T[T.subset == "all"].rms_px.median(), T[T.subset == "all"].rms_px.max()))

    if a.exclude is None:
        if a.no_auto:
            a.exclude, source = [], "none"
        else:
            a.exclude, reasons = outliers(T[T.subset == "all"])
            source = "outlier rule"
            for e, why in reasons.items():
                print(f"  outlier rule: {e} left out ({why})")
    print(f"excluded tine pairs ({source}): {', '.join(sorted(a.exclude)) or 'none'}")

    # use the middle three pairs (least lateral splay) over all frames, both axes; if fewer than two of
    # them are left, every pair that is left
    excl = {(e[0].upper(), int(e[1:])) for e in a.exclude}
    keep = ~np.array([(ax_, int(pr)) in excl for ax_, pr in zip(T.axis, T.pair)], bool)
    mid = T[(T.subset == "all") & T.pair.isin([2, 3, 4]) & keep]
    if mid[["axis", "pair"]].drop_duplicates().shape[0] < 2:
        mid = T[(T.subset == "all") & keep]
        print("  fewer than two middle pairs left: using every pair left")
    if mid.empty:
        sys.exit("every tine pair is excluded: check the tine tips (calibration_tines_px) and scale_calibration.png")
    s = mid.px_per_mm
    px_per_mm = float(s.mean())
    result = dict(px_per_mm=px_per_mm, mm_per_px=1 / px_per_mm, um_per_px=1000 / px_per_mm,
                  std_px_per_mm=float(s.std()), px_per_mm_X=float(mid[mid.axis == "X"].px_per_mm.mean()),
                  px_per_mm_Y=float(mid[mid.axis == "Y"].px_per_mm.mean()),
                  method="tine-tip separation vs reported rake separation (middle 3 tine pairs, all frames"
                         + (f"; excluded {', '.join(sorted(a.exclude))}" if a.exclude else "") + ")",
                  excluded=sorted(a.exclude), excluded_by=source,
                  deformations=list(T.deformation.unique()))
    print("\nRESULT: %.3f px/mm  (%.2f um/px);  X %.3f, Y %.3f px/mm;  sd across pairs/defs %.3f px/mm (%.2f%%)" %
          (px_per_mm, 1000 / px_per_mm, result["px_per_mm_X"], result["px_per_mm_Y"], s.std(),
           100 * s.std() / px_per_mm))
    with open(TC.scale_file, "w") as f:
        json.dump(result, f, indent=2)
    print("wrote", TC.scale_file)

    # diagnostic plot
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    for folder, axis, rake, dsep, force in plot:
        ls = "-" if axis == "X" else "--"
        for pair in range(dsep.shape[1]):
            ax[0].plot(rake / 1000, dsep[:, pair], ls, lw=1, color=f"C{pair}",
                       label=f"{axis} pair {pair+1}" if folder == plot[0][0] else None)
            r = dsep[:, pair] - px_per_mm * rake / 1000
            ax[1].plot(force, r - np.nanmedian(r[:5]), ls, lw=.8, color=f"C{pair}")
    ax[0].set_xlabel("reported change in rake separation (mm)"); ax[0].set_ylabel("change in tine-tip separation (px)")
    ax[0].set_title("Tine tips vs reported rake motion (solid X, dashed Y)"); ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=.3)
    ax[1].set_xlabel("rake force (mN)"); ax[1].set_ylabel("residual from single scale (px)")
    ax[1].set_title("Residual vs force (tine bending would trend here)"); ax[1].grid(alpha=.3)
    sub = T[T.subset == "all"]
    for i, axis in enumerate(("X", "Y")):
        q = sub[sub.axis == axis]
        ax[2].scatter(q.pair + (i - .5) * .15, q.px_per_mm, label=f"{axis} axis", s=30)
    ax[2].axhline(px_per_mm, color="k", lw=.8)
    ax[2].set_xlabel("tine pair (1 = top/left … 5 = bottom/right)"); ax[2].set_ylabel("px per mm")
    ax[2].set_title(f"Scale per tine pair; line = result {px_per_mm:.2f} px/mm"); ax[2].legend(); ax[2].grid(alpha=.3)
    out = os.path.join(TC.results, "scale_calibration.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=95)
    print("wrote", out)


if __name__ == "__main__":
    main()
