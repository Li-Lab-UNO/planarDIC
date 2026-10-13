"""
Animate the scale calibration for one tine pair: the tracked tine tips in every frame of a
deformation, next to the calibration plots (as in calibrate_scale.py's scale_calibration.png)
drawn up to the current frame.

Left: the specimen with every tracked tip; the chosen pair is highlighted with the line whose
length is the tip-to-tip separation. Zoomed insets of the two tips show their motion against
their rest positions. Right: change in tip separation vs the tester's reported change in rake
separation, with the calibration result as a line (its slope is px per mm), and the residual
from that line vs rake force.

    ../.venv/bin/python animate_calibration.py --test N1690/test2 --def Def01 --pair X3
Output: <test>/Images/<Def>/strain_results/scale_calibration_<pair>.mp4
"""
import argparse
import glob
import json
import os
import sys

import cv2
import imageio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import calibrate_scale as CS  # noqa: E402  (selects the test from --test)
import dic2d as D  # noqa: E402

TC = CS.TC
TINES = CS.tines()
PAIR_COLORS = [f"C{i}" for i in range(5)]          # same colors per pair as scale_calibration.png
INK, MUTED, GRID = "#1f1f1e", "#6b6b67", "#d9d9d6"


def track(paths):
    """Tip (axis position c, tip coordinate t) of every tine in every frame: side -> (T, 5, 2)."""
    state = {k: [list(p) for p in v["tips"]] for k, v in TINES.items()}
    out = {k: np.full((len(paths), len(v["tips"]), 2), np.nan) for k, v in TINES.items()}
    for k, p in enumerate(paths):
        g = cv2.GaussianBlur(cv2.imread(p, cv2.IMREAD_GRAYSCALE).astype(np.float32), (0, 0), 0.8)
        for side, spec in TINES.items():
            for i, (c, t) in enumerate(state[side]):
                for _ in range(2):
                    c, t_new, contrast = CS.find_tip(g, c, t, spec["horizontal"], spec["out"])
                    if np.isfinite(t_new):
                        t = t_new
                state[side][i] = [c, t]
                if contrast > 30:
                    out[side][k, i] = (c, t)
    return out


def xy(side, ct):
    """Image (x, y) of a tip given (c, t)."""
    return (ct[..., 1], ct[..., 0]) if TINES[side]["horizontal"] else (ct[..., 0], ct[..., 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--def", dest="defn", default="Def01")
    ap.add_argument("--pair", default="X3", help="axis and pair number, e.g. X3 (left/right tine 3) or Y2")
    ap.add_argument("--fps", type=float, default=12)
    ap.add_argument("--zoom-margin", type=int, default=45,
                    help="image px shown around each tip's path in the insets")
    a = ap.parse_args()
    axis, pair = a.pair[0].upper(), int(a.pair[1:]) - 1
    lo, hi = ("left", "right") if axis == "X" else ("top", "bottom")

    folder = TC.folder(a.defn)
    times, paths = D.find_frames([os.path.join(TC.images, folder)])
    print(f"{folder}: tracking tine tips in {len(paths)} frames")
    tips = track(paths)
    Dd = pd.read_csv(TC.first("*Data.csv"), encoding="latin1").drop_duplicates("Time_S").set_index("Time_S")
    dd = Dd.reindex(np.round(times, 1))
    size = "X_Size_um" if axis == "X" else "Y_Size_um"
    rake = (dd[size] - dd[size].iloc[0]).to_numpy(float) / 1000          # mm
    force = dd[f"{axis}_Force_mN"].to_numpy(float)
    sep = tips[hi][..., 1] - tips[lo][..., 1]                             # (T, 5) tip-to-tip separation, px
    dsep = sep - sep[0]
    scale = json.load(open(TC.scale_file))["px_per_mm"]
    res = dsep - scale * rake[:, None]
    res -= np.nanmedian(res[:5], axis=0)
    m = np.isfinite(dsep[:, pair])
    k_pair = np.polyfit(rake[m], dsep[m, pair], 1)[0]
    color = PAIR_COLORS[pair]

    # image crop around all tips, in rest-frame coordinates
    allx = np.concatenate([xy(s, tips[s][0])[0] for s in tips])
    ally = np.concatenate([xy(s, tips[s][0])[1] for s in tips])
    excursion = max(np.nanmax(np.abs(tips[s][:, :, 1] - tips[s][0, :, 1])) for s in tips)
    pad = 60 + int(excursion)
    x0, x1 = int(np.nanmin(allx) - pad), int(np.nanmax(allx) + pad)
    y0, y1 = int(np.nanmin(ally) - pad), int(np.nanmax(ally) + pad)
    rest_lo, rest_hi = [xy(s, tips[s][0, pair]) for s in (lo, hi)]
    # insets: a square around each tip's whole path (rest to farthest excursion)
    insets = []
    for s in (lo, hi):
        px, py = xy(s, tips[s][:, pair])
        cx, cy = np.nanmean([np.nanmin(px), np.nanmax(px)]), np.nanmean([np.nanmin(py), np.nanmax(py)])
        h = int(max(np.nanmax(px) - np.nanmin(px), np.nanmax(py) - np.nanmin(py)) / 2 + a.zoom_margin)
        insets.append((int(cx), int(cy), h))

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})
    fig = plt.figure(figsize=(16, 9), dpi=100, facecolor="white")
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 0.55, 1.35], height_ratios=[1, 1], left=0.03, right=0.97,
                          top=0.9, bottom=0.08, wspace=0.18, hspace=0.32)
    ax_img = fig.add_subplot(gs[:, 0])
    ax_z = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])]
    ax_cal, ax_res = fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 2])
    rmax, fmax = np.nanmax(rake), np.nanmax(force)
    others = [p for p in range(dsep.shape[1]) if p != pair]

    def render(k):
        for ax_ in (ax_img, *ax_z, ax_cal, ax_res):
            ax_.cla()
        img = cv2.cvtColor(cv2.imread(paths[k]), cv2.COLOR_BGR2RGB)
        ax_img.imshow(img[y0:y1, x0:x1], extent=(x0, x1, y1, y0))
        for s in tips:
            px, py = xy(s, tips[s][k])
            ax_img.plot(px, py, "o", ms=4, mfc="none", mec="white", mew=1)
        (ax_lo, ay_lo), (ax_hi, ay_hi) = xy(lo, tips[lo][k, pair]), xy(hi, tips[hi][k, pair])
        ax_img.plot([ax_lo, ax_hi], [ay_lo, ay_hi], "-", color=color, lw=2)
        ax_img.plot([ax_lo, ax_hi], [ay_lo, ay_hi], "o", ms=8, color=color, mec="white", mew=1.5)
        ax_img.set_xlim(x0, x1); ax_img.set_ylim(y1, y0); ax_img.axis("off")
        ax_img.set_title(f"{a.pair}: {lo} and {hi} tine {pair + 1}   separation {sep[k, pair]:.1f} px "
                         f"(change {dsep[k, pair]:+.1f} px)", fontsize=12, loc="left")
        for ax_, side, (rx, ry), (cx, cy, h) in zip(ax_z, (lo, hi), (rest_lo, rest_hi), insets):
            ax_.imshow(img[cy - h:cy + h, cx - h:cx + h], extent=(cx - h, cx + h, cy + h, cy - h))
            px, py = xy(side, tips[side][k, pair])
            if TINES[side]["horizontal"]:
                ax_.axvline(rx, color="white", ls="--", lw=1.2)
            else:
                ax_.axhline(ry, color="white", ls="--", lw=1.2)
            ax_.plot(px, py, "o", ms=10, color=color, mec="white", mew=1.5)
            ax_.set_xticks([]); ax_.set_yticks([])
            mv = (px - rx) if TINES[side]["horizontal"] else (py - ry)
            ax_.set_title(f"{side} tip: {mv:+.1f} px", fontsize=11, loc="left")
        ax_z[1].text(0, -0.12, "dashed: rest position", transform=ax_z[1].transAxes, fontsize=9, color=MUTED,
                     va="top")
        # calibration plot
        for p in others:
            ax_cal.plot(rake[:k + 1], dsep[:k + 1, p], color=GRID, lw=1)
            ax_res.plot(force[:k + 1], res[:k + 1, p], color=GRID, lw=1)
        xs = np.array([0, rmax])
        ax_cal.plot(xs, scale * xs, color=INK, lw=1, ls="--")
        ax_cal.text(0.03, 0.95, f"dashed: calibration result, {scale:.2f} px/mm\ngray: the other {axis} pairs",
                    transform=ax_cal.transAxes, ha="left", va="top", fontsize=10, color=MUTED)
        ax_cal.plot(rake[:k + 1], dsep[:k + 1, pair], color=color, lw=2)
        ax_cal.plot(rake[k], dsep[k, pair], "o", ms=9, color=color, mec="white", mew=1.5)
        ax_cal.set_xlim(-0.03 * rmax, rmax * 1.05)
        ax_cal.set_ylim(min(-5, np.nanmin(dsep) - 5), np.nanmax(dsep) * 1.08)
        ax_cal.set_xlabel("reported change in rake separation (mm)")
        ax_cal.set_ylabel("change in tine-tip separation (px)")
        ax_cal.set_title(f"Tine tips vs reported rake motion   (this pair: slope {k_pair:.1f} px/mm)",
                         fontsize=11, loc="left")
        ax_cal.grid(alpha=.3)
        ax_res.axhline(0, color=INK, lw=0.8)
        ax_res.plot(force[:k + 1], res[:k + 1, pair], color=color, lw=1.5)
        ax_res.plot(force[k], res[k, pair], "o", ms=9, color=color, mec="white", mew=1.5)
        ax_res.set_xlim(-0.03 * fmax, fmax * 1.05)
        lim = np.nanpercentile(np.abs(res[np.isfinite(res)]), 99) + 3
        ax_res.set_ylim(-lim, lim)
        ax_res.set_xlabel("rake force (mN)"); ax_res.set_ylabel("residual from single scale (px)")
        ax_res.set_title("Residual vs force (tine bending would trend here)", fontsize=11, loc="left")
        ax_res.grid(alpha=.3)
        fig.suptitle(f"{TC.sample} {TC.name.split('_')[0]}, {folder}: scale calibration from tine-tip motion   "
                     f"t = {times[k] - times[0]:5.1f} s   rake change {rake[k]:.3f} mm   force {force[k]:5.0f} mN",
                     x=0.03, ha="left", fontsize=14, weight="bold")
        fig.canvas.draw()
        return np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()

    out = os.path.join(TC.images, folder, "strain_results", f"scale_calibration_{a.pair}.mp4")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    writer = imageio.get_writer(out, fps=a.fps, codec="libx264", quality=8, macro_block_size=4)
    for k in range(len(paths)):
        writer.append_data(render(k))
    writer.close()
    print(f"pair {a.pair}: slope {k_pair:.2f} px/mm (result {scale:.2f}); wrote {out}")


if __name__ == "__main__":
    main()
