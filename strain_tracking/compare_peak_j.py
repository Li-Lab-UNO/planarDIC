"""
Cross-specimen comparison of local areal stretch at peak: J/mean maps at the reference protocol
(Def11) for each test, one row per sample, and the per-protocol summary numbers of every test side by side.

Uses the same computation as peak_j_figure.py (frame of maximum mean axial stretch, J divided by the
interior mean, points valid in fewer than half the frames gray, lowest 3 x 2 mm interior box). Use a
physical window (--window-mm) so specimens imaged at different zoom are smoothed alike.
Tests without a track for every one of Def11-Def19 are left out, with a note.

    .venv/bin/python strain_tracking/compare_peak_j.py --samples N1739 N1690      (window from analysis_settings.json)
    .venv/bin/python strain_tracking/compare_peak_j.py --samples N1739            (one sample; up to four)
    .venv/bin/python strain_tracking/compare_peak_j.py --samples N1739 N1690 --window-mm 1.4
    .venv/bin/python strain_tracking/compare_peak_j.py --samples N1739 N1690 --edge-mm 1.6 2.0   (another interior)
Output: comparisons/peak_J_<samples>_w<L>mm[_edge<Z>x<TH>mm].png / .csv
"""
import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import peak_j_figure as PJ  # noqa: E402
import testcfg  # noqa: E402

INK, MUTED, SURF, GRID = PJ.INK, PJ.MUTED, PJ.SURF, PJ.GRID
SAMPLE_COLOR = ["#2a78d6", "#eb6834", "#1baf7a", "#8a56c2"]   # categorical slots 1-4, in --samples order
DEFS = [d for _, ds in PJ.ROWS for d in ds]


def tests_of(sample):
    folder = [d for d in glob.glob(os.path.join(testcfg.ROOT, sample + "*")) if os.path.isdir(d)]
    if len(folder) != 1:
        sys.exit(f"sample {sample!r}: expected one folder, found {folder}")
    tests = []
    for t in sorted(glob.glob(os.path.join(folder[0], "test*_*"))):
        if not os.path.exists(os.path.join(t, "test_config.json")):
            continue
        folders = {f.split("_")[0]: f for f in os.listdir(os.path.join(t, "Images")) if re.match(r"Def\d+_", f)} \
            if os.path.isdir(os.path.join(t, "Images")) else {}
        missing = [d for d in DEFS if d not in folders
                   or not os.path.exists(os.path.join(t, "DIC_results", folders[d] + ".npz"))]
        if missing:
            print(f"{sample} {os.path.basename(t)}: left out, no track for {', '.join(missing)}")
            continue
        tests.append(os.path.basename(t).split("_")[0])
    if not tests:
        sys.exit(f"sample {sample!r}: no test has all of {DEFS[0]}-{DEFS[-1]} tracked")
    return tests


def specimen(sample, test, ref, window_mm, box_mm, edge):
    testcfg.select(f"{sample}/{test}")
    T = testcfg.get()
    mpp = D.load_scale()["mm_per_px"]
    tips = T.setting("tine_tips_px")
    P = {d: PJ.peak_J(T, d, window_mm, None, edge) for d in DEFS}
    z = (P[ref]["X"] - tips["z"][0]) * mpp
    th = (tips["theta"][1] - P[ref]["Y"]) * mpp
    box, outline, _ = PJ.low_spot(P[ref]["Jn"], z, th, box_mm, P[ref]["inner"])
    rows = []
    for d in DEFS:
        Jn = P[d]["Jn"]
        inn = P[d]["inner"] & np.isfinite(Jn)
        rows.append(dict(sample=sample, test=test, protocol=d, low_spot_J_over_mean=float(np.nanmean(Jn[box])),
                         interior_sd_J_over_mean=float(np.nanstd(Jn[inn])),
                         interior_pct_below_0p9=float(np.mean(Jn[inn] < 0.9) * 100),
                         interior_edge_z_mm=edge["z"], interior_edge_theta_mm=edge["theta"]))
    return dict(Jn=P[ref]["Jn"], z=z, th=th, outline=outline, lz=P[ref]["lz"], lt=P[ref]["lt"],
                prog=D.programmed_stretch(T.root, T.folder(ref)),
                spans=((tips["z"][1] - tips["z"][0]) * mpp, (tips["theta"][1] - tips["theta"][0]) * mpp),
                rows=rows, px_per_mm=1 / mpp)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", nargs="+", default=["N1739", "N1690"], help="one to four samples, one map row each")
    ap.add_argument("--ref", default="Def11")
    ap.add_argument("--box", nargs=2, type=float, default=None, metavar=("Z_MM", "THETA_MM"),
                    help="low-spot box (default: low_spot_box_mm of analysis_settings.json)")
    ap.add_argument("--window-mm", type=float, default=None, help="default: window_mm of analysis_settings.json")
    ap.add_argument("--edge-mm", nargs=2, type=float, default=None, metavar=("Z_MM", "THETA_MM"),
                    help="interior edge distances; default: edge_mm of analysis_settings.json")
    ap.add_argument("--vmin", type=float, default=0.80)
    ap.add_argument("--vmax", type=float, default=1.20)
    a = ap.parse_args()
    if not 1 <= len(a.samples) <= len(SAMPLE_COLOR):
        sys.exit(f"--samples: give 1 to {len(SAMPLE_COLOR)} samples")
    a.window_mm = a.window_mm or testcfg.analysis()["window_mm"]
    a.box = a.box or testcfg.analysis()["low_spot_box_mm"]
    edge = testcfg.edge_mm(a.edge_mm)
    S = {s: {t: specimen(s, t, a.ref, a.window_mm, a.box, edge) for t in tests_of(s)} for s in a.samples}
    df = pd.DataFrame([r for s in S.values() for sp in s.values() for r in sp["rows"]])

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})
    ncol = max(len(v) for v in S.values())
    n = len(S)
    H = 5.8 * n + 3.9                                   # inches: a map row per sample, the summary row, the header
    fig = plt.figure(figsize=(14, H), facecolor=SURF)
    gs = fig.add_gridspec(n + 1, 1, height_ratios=[1] * n + [0.78], left=0.09, right=0.87, top=1 - 1.95 / H,
                          bottom=0.78 / H, hspace=0.5)
    # maps: one row per sample, shared mm extent so sizes compare directly
    zmax = max(sp["z"].max() for s in S.values() for sp in s.values()) + 0.3
    tmax = max(sp["th"].max() for s in S.values() for sp in s.values()) + 0.3
    for r, (s, tests) in enumerate(S.items()):
        sub = gs[r].subgridspec(1, ncol, wspace=0.1)
        for c, (t, sp) in enumerate(tests.items()):
            ax = fig.add_subplot(sub[c])
            ax.set_facecolor(SURF)
            im = ax.pcolormesh(sp["z"], sp["th"], sp["Jn"], cmap=PJ.CMAP, vmin=a.vmin, vmax=a.vmax,
                               shading="nearest", rasterized=True)
            Lz, Lt = sp["spans"]
            ax.add_patch(Rectangle((edge["z"], edge["theta"]), Lz - 2 * edge["z"], Lt - 2 * edge["theta"], fill=False,
                                   ec=MUTED, lw=1, ls="--"))
            z0, z1, t0, t1 = sp["outline"]
            ax.add_patch(Rectangle((z0, t0), z1 - z0, t1 - t0, fill=False, ec=INK, lw=1.8))
            lo = next(x["low_spot_J_over_mean"] for x in sp["rows"] if x["protocol"] == a.ref)
            sd = next(x["interior_sd_J_over_mean"] for x in sp["rows"] if x["protocol"] == a.ref)
            ax.text((z0 + z1) / 2, t1 + 0.12, f"{lo:.2f}", ha="center", va="bottom", fontsize=10.5, weight="bold",
                    bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.85))
            ax.set_xlim(-0.3, zmax); ax.set_ylim(-0.3, tmax); ax.set_aspect("equal")
            ax.set_xticks(np.arange(0, zmax, 2)); ax.set_yticks(np.arange(0, tmax, 2))
            ax.tick_params(length=2, labelsize=8.5)
            ptxt = f"   programmed {sp['prog'][0]:g}% / {sp['prog'][1]:g}%" if sp["prog"] else ""
            ax.text(0, 1.075, f"{s} {t}", transform=ax.transAxes, fontsize=12, va="bottom", weight="bold")
            ax.text(0, 1.015, f"λ_z {sp['lz']:.2f}  λ_θ {sp['lt']:.2f}   SD {sd:.3f}{ptxt}", transform=ax.transAxes,
                    fontsize=9, color=MUTED, va="bottom")
            if c == 0:
                ax.set_ylabel("θ, circumferential (mm)", fontsize=9.5)
            else:
                ax.set_yticklabels([])
            ax.set_xlabel("z, axial (mm)", fontsize=9)
    maps = fig.axes
    y_lo, y_hi = min(m.get_position().y0 for m in maps), max(m.get_position().y1 for m in maps)
    cax = fig.add_axes([0.89, y_lo + 0.1 * (y_hi - y_lo), 0.013, 0.8 * (y_hi - y_lo)])
    cb = fig.colorbar(im, cax=cax, extend="both")
    cb.outline.set_visible(False)
    cb.set_label("areal stretch J / specimen mean", color=INK, fontsize=10)
    cb.set_ticks(np.round(np.linspace(a.vmin, a.vmax, 5), 2))
    cb.ax.tick_params(labelsize=9, length=2)
    cax.text(0.5, -0.04, "less\ndeformation", transform=cax.transAxes, fontsize=8.5, color=MUTED, va="top", ha="center")
    cax.text(0.5, 1.04, "more\ndeformation", transform=cax.transAxes, fontsize=8.5, color=MUTED, va="bottom", ha="center")

    # summary: each test's range over Def11-19 (line) and its reference protocol (dot)
    metrics = [("interior_sd_J_over_mean", "spread: SD of J/mean", "{:.3f}"),
               ("interior_pct_below_0p9", "interior more than 10% below\nthe mean (%)", "{:.1f}"),
               ("low_spot_J_over_mean", "lowest 3 × 2 mm region\n(J/mean)", "{:.2f}")]
    sub = gs[n].subgridspec(1, 3, wspace=0.3)
    labels = [(s, t) for s, ts in S.items() for t in ts]
    xs = np.arange(len(labels)) + np.repeat(np.arange(len(S)), [len(v) for v in S.values()]) * 0.6
    for m, (col, title, fmt) in enumerate(metrics):
        ax = fig.add_subplot(sub[m])
        ax.set_facecolor(SURF)
        ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
        for x, (s, t) in zip(xs, labels):
            d = df[(df["sample"] == s) & (df["test"] == t)]
            color = SAMPLE_COLOR[a.samples.index(s)]
            ax.plot([x, x], [d[col].min(), d[col].max()], color=color, lw=2, solid_capstyle="round")
            v = float(d.loc[d.protocol == a.ref, col].iloc[0])
            ax.plot(x, v, "o", ms=8, color=color, mec=SURF, mew=2, zorder=3)
        ax.set_xticks(xs); ax.set_xticklabels([t.replace("test", "") for _, t in labels], fontsize=8.5)
        for i, s in enumerate(S):
            xi = xs[[k for k, (ss, _) in enumerate(labels) if ss == s]]
            ax.text(xi.mean(), -0.16, s, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=10,
                    color=INK, weight="bold")
        ax.set_xlim(xs[0] - 0.6, xs[-1] + 0.6)
        ax.tick_params(length=2, labelsize=8.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title, fontsize=10.5, color=INK, loc="left")
    y1 = max(x.get_position().y1 for x in fig.axes[-3:])
    fig.text(0.09, y1 + 0.62 / H, f"Per-test summary (test number under each sample): dot = {a.ref}, line = range "
             "over Def11–Def19", fontsize=11, color=INK)

    fig.text(0.09, 1 - 0.43 / H, f"Local areal stretch at peak: {' vs '.join(a.samples)}", fontsize=15, weight="bold")
    fig.text(0.09, 1 - 0.74 / H, f"J = λ_1·λ_2 divided by the interior mean at {a.ref}'s peak, {a.window_mm:g} mm Gaussian strain "
                         f"window for every test.\nOutlined: the lowest {a.box[0]:g} × {a.box[1]:g} mm interior region "
                         f"(its mean J/mean above it). Dashed: interior, {edge['z']:g} mm from the left/right and {edge['theta']:g} mm "
                         "from the top/bottom tine lines.\nGray: untracked points. Same mm scale in every map. Programmed: axial / "
                         f"circumferential rake stretch at {a.ref}.",
             fontsize=10, color=MUTED, va="top")
    tag = "_".join(a.samples) + f"_w{a.window_mm:g}mm" + (f"_edge{edge['z']:g}x{edge['theta']:g}mm" if a.edge_mm else "")
    out = os.path.join(testcfg.ROOT, "comparisons", f"peak_J_{tag}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out + ".png", dpi=130, facecolor=SURF)
    df.to_csv(out + ".csv", index=False, float_format="%.4f")
    with pd.option_context("display.width", 160):
        print(df.groupby(["sample", "test"], sort=False).agg(
            sd_min=("interior_sd_J_over_mean", "min"), sd_max=("interior_sd_J_over_mean", "max"),
            below_max=("interior_pct_below_0p9", "max"), low_min=("low_spot_J_over_mean", "min"),
            low_max=("low_spot_J_over_mean", "max")).round(3))
    print("wrote", out + ".png and .csv")


if __name__ == "__main__":
    main()
