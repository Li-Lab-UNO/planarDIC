"""
Local areal stretch at each protocol's peak, as small multiples: J = lambda_1 lambda_2 divided by
the specimen's interior mean, for Def11-Def19 of one test, with the specimen's low spot outlined.

Computed from the full-domain tracks (track_protocols.py), no exports needed: at the frame of
maximum mean axial stretch, F from a Gaussian strain window of window_mm (analysis_settings.json,
or --window-mm), J = det F, divided by the mean over the interior (grid points more than edge_mm
from the tine-tip lines: edge_mm z from the left/right lines, theta from the top/bottom lines;
outlined dashed). Points valid in fewer than half the frames are untracked: shown gray
and left out of the mean. --box-px 75 gives the earlier 75 px box window instead.

Low spot: the lowest-mean box of --box Z TH mm (default low_spot_box_mm) lying fully in the interior at
--ref (Def11), outlined in every panel with its mean J/mean. Rows group the protocols by
loading family; ratios are axial : circumferential rake stretch (same design in every test).
The color scale is fixed (--vmin/--vmax, default 0.8-1.2) so specimens can be compared.

    .venv/bin/python strain_tracking/peak_j_figure.py --test N1690/test1
    .venv/bin/python strain_tracking/peak_j_figure.py --test N1739/test1
    .venv/bin/python strain_tracking/peak_j_figure.py --test N1690/test1 --window-mm 1.4   (another physical window)
    .venv/bin/python strain_tracking/peak_j_figure.py --test N1739/test3 --edge-mm 1.6 2.0   (another interior)
Output: <test>/DIC_results/peak_J_over_mean_Def11-19[_w<L>mm | _box<N>px][_edge<Z>x<TH>mm].png (and a .csv of the
per-protocol numbers, with the interior edges used); no suffix = the analysis_settings.json window and interior
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

INK, MUTED, SURF, GRID = "#1f1f1e", "#6b6b67", "#fcfcfb", "#d9d9d6"
# diverging: blue (less deformation) - neutral gray at the mean - red (more deformation)
CMAP = LinearSegmentedColormap.from_list(
    "blue_gray_red", ["#104281", "#2a78d6", "#9ec5f4", "#f0efec", "#f3b3ad", "#e34948", "#8f1d1d"])
CMAP.set_bad(GRID)
ROWS = [("Equibiaxial", ["Def11", "Def15", "Def19"]),
        ("Axial reduced", ["Def12", "Def13", "Def14"]),
        ("Circumferential reduced", ["Def16", "Def17", "Def18"])]
RATIO = {"Def11": "1 : 1", "Def15": "1 : 1", "Def19": "1 : 1 (save)",
         "Def12": "0.75 : 1", "Def13": "0.5 : 1", "Def14": "0.25 : 1",
         "Def16": "1 : 0.75", "Def17": "1 : 0.5", "Def18": "1 : 0.25"}


def peak_J(T, defn, window_mm=None, box_px=None, edge_mm=None, R=None):
    """J/mean at the peak of protocol defn. Defaults (window_mm, edge_mm) from analysis_settings.json;
    box_px uses the square box window of that many px instead. R: a track dict to use instead of
    the test's saved full-domain track."""
    S = testcfg.analysis()
    window_mm = window_mm or S["window_mm"]
    edge_mm = testcfg.edge_mm(edge_mm)
    if R is None:
        R, _ = D.load_track(T.results, T.folder(defn), "tine_region")
    X0, sh = R["X0"], tuple(R["shape"])
    mpp = D.load_scale(T.scale_file)["mm_per_px"]
    best = (-1.0, 0)
    for k in range(len(R["times"])):                      # frame of maximum mean axial stretch
        v = R["valid"][k]
        if v.sum() < 10:
            continue
        A = D.affine_fit(X0[v], R["pos"][k][v])
        lz = np.sqrt((A.T @ A)[0, 0])
        if lz > best[0]:
            best = (lz, k)
    k = best[1]
    if box_px:
        F, ok = D.deformation_gradient_field(X0, R["pos"][k], R["valid"][k], sh,
                                             D.window_points(float(R["step"]), box_px))
    else:
        F, ok = D.deformation_gradient_field_mm(X0, R["pos"][k], R["valid"][k], sh, window_mm, mpp)
    F = D.to_specimen(F)
    with np.errstate(invalid="ignore"):
        J = np.where(ok, np.linalg.det(F), np.nan)
    C = np.einsum("...ki,...kj->...ij", F, F)
    inner = D.interior_mask(X0, sh, T.setting("tine_tips_px"), mpp, edge_mm)
    J[R["valid"].mean(0).reshape(sh) < 0.5] = np.nan       # untracked: not shown, not in the mean
    Jn = J / np.nanmean(J[inner])
    return dict(Jn=Jn, inner=inner, X=X0[:, 0].reshape(sh), Y=X0[:, 1].reshape(sh), t=float(R["times"][k]),
                lz=float(np.nanmean(np.sqrt(C[..., 0, 0])[ok])), lt=float(np.nanmean(np.sqrt(C[..., 1, 1])[ok])))


def low_spot(Jn, z, th, box_mm=(3.0, 2.0), inner=None):
    """Lowest-mean box of box_mm (z, theta) lying fully in the interior mask `inner` (default: all but
    the outer 2 grid rows/cols).
    Returns (mask, (z0, z1, th0, th1) outline in mm, mean J/mean)."""
    step = abs(z[0, 1] - z[0, 0])
    nz, nt = int(round(box_mm[0] / step)) + 1, int(round(box_mm[1] / step)) + 1
    best = (np.inf, None)
    if inner is None:
        inner = np.zeros(Jn.shape, bool)
        inner[2:-2, 2:-2] = True
    for i in range(0, Jn.shape[0] - nt + 1):
        for j in range(0, Jn.shape[1] - nz + 1):
            w = Jn[i:i + nt, j:j + nz]
            if inner[i:i + nt, j:j + nz].all() and np.isfinite(w).mean() >= 0.8 and np.nanmean(w) < best[0]:
                best = (np.nanmean(w), (i, j))
    i0, j0 = best[1]
    mask = np.zeros(Jn.shape, bool)
    mask[i0:i0 + nt, j0:j0 + nz] = True
    outline = (z[0, j0] - step / 2, z[0, j0 + nz - 1] + step / 2, th[i0 + nt - 1, 0] - step / 2, th[i0, 0] + step / 2)
    return mask, outline, float(best[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="Def11", help="protocol whose lowest box defines the low spot")
    ap.add_argument("--box", nargs=2, type=float, default=None, metavar=("Z_MM", "THETA_MM"),
                    help="low-spot box (default: low_spot_box_mm of analysis_settings.json)")
    ap.add_argument("--exclude", nargs=4, type=float, action="append", default=[],
                    metavar=("Z0", "Z1", "TH0", "TH1"),
                    help="region (mm from the tine tips) to leave out of the maps and the low-spot search, "
                         "e.g. a known tracking artifact; repeatable")
    ap.add_argument("--window-mm", type=float, default=None,
                    help="strain window (Gaussian, equivalent box width in mm); default: window_mm of "
                         "analysis_settings.json")
    ap.add_argument("--box-px", type=int, default=None, help="square box window of this many px instead (e.g. 75)")
    ap.add_argument("--edge-mm", nargs=2, type=float, default=None, metavar=("Z_MM", "THETA_MM"),
                    help="interior edge distances from the left/right and top/bottom tine lines; default: "
                         "edge_mm of analysis_settings.json")
    ap.add_argument("--vmin", type=float, default=0.80)
    ap.add_argument("--vmax", type=float, default=1.20)
    a = ap.parse_args()
    T = testcfg.get()
    mpp = D.load_scale()["mm_per_px"]
    tips, cal = T.setting("tine_tips_px"), T.setting("calibration_tines_px")
    S = testcfg.analysis()
    a.box = a.box or S["low_spot_box_mm"]
    em = testcfg.edge_mm(a.edge_mm)
    P = {d: peak_J(T, d, a.window_mm, a.box_px, em) for _, defs in ROWS for d in defs}
    s_ = D.load_scale()["px_per_mm"]
    wtxt = (f"{a.box_px} px box ({a.box_px / s_:.2f} mm)" if a.box_px
            else f"{a.window_mm or S['window_mm']:g} mm Gaussian")
    ref = P[a.ref]
    z = (ref["X"] - tips["z"][0]) * mpp
    th = (tips["theta"][1] - ref["Y"]) * mpp
    step = abs(z[0, 1] - z[0, 0])
    for z0, z1, t0, t1 in a.exclude:                       # known artifacts: shown gray, never the low spot
        m = (z >= z0) & (z <= z1) & (th >= t0) & (th <= t1)
        for p in P.values():
            p["Jn"][m] = np.nan

    box, (bz0, bz1, bt0, bt1), lo = low_spot(ref["Jn"], z, th, a.box, ref["inner"])
    ez, et = em["z"], em["theta"]
    Lz, Lt = (tips["z"][1] - tips["z"][0]) * mpp, (tips["theta"][1] - tips["theta"][0]) * mpp
    deficit = 100 * (1 - lo)

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})
    fig, ax = plt.subplots(3, 3, figsize=(13.5, 14.2), facecolor=SURF)
    fig.subplots_adjust(left=0.1, right=0.86, top=0.855, bottom=0.06, wspace=0.12, hspace=0.36)
    ext = (z.min() - step / 2, z.max() + step / 2, th.min() - step / 2, th.max() + step / 2)
    rows = []
    for r, (rname, defs) in enumerate(ROWS):
        for c, d in enumerate(defs):
            ax_, p = ax[r, c], P[d]
            ax_.set_facecolor(SURF)
            im = ax_.pcolormesh(z, th, p["Jn"], cmap=CMAP, vmin=a.vmin, vmax=a.vmax, shading="nearest",
                                rasterized=True)
            ax_.add_patch(Rectangle((ez, et), Lz - 2 * ez, Lt - 2 * et, fill=False, ec=MUTED, lw=1, ls="--"))
            ax_.add_patch(Rectangle((bz0, bt0), bz1 - bz0, bt1 - bt0, fill=False, ec=INK, lw=1.8))
            val = float(np.nanmean(p["Jn"][box]))
            tx, ha = (bz1 + 0.15, "left") if bz1 + 1.2 < ext[1] else (bz0 - 0.15, "right")
            ax_.text(tx, (bt0 + bt1) / 2, f"{val:.2f}", va="center", ha=ha, fontsize=10.5, color=INK, weight="bold",
                     bbox=dict(boxstyle="round,pad=0.18", fc=SURF, ec="none", alpha=0.85))
            for side, spec in cal.items():                  # rest-state tine positions on the edges
                for cc, _ in spec["tips"]:
                    if spec["horizontal"]:
                        ax_.plot(ext[0] if spec["out"] < 0 else ext[1], (tips["theta"][1] - cc) * mpp,
                                 marker=">" if spec["out"] < 0 else "<", ms=5, color=MUTED, clip_on=False)
                    else:
                        ax_.plot((cc - tips["z"][0]) * mpp, ext[3] if spec["out"] < 0 else ext[2],
                                 marker="v" if spec["out"] < 0 else "^", ms=5, color=MUTED, clip_on=False)
            ax_.set_xlim(ext[0], ext[1]); ax_.set_ylim(ext[2], ext[3]); ax_.set_aspect("equal")
            prog = D.programmed_stretch(T.root, T.folder(d))
            ptxt = f"   {prog[0]:g}% / {prog[1]:g}%" if prog else ""
            ax_.text(0, 1.085, f"{d}  {RATIO[d]}{ptxt}", transform=ax_.transAxes, fontsize=11.5, color=INK,
                     va="bottom")
            ax_.text(0, 1.02, f"λ_z {p['lz']:.2f}   λ_θ {p['lt']:.2f}", transform=ax_.transAxes, fontsize=9,
                     color=MUTED, va="bottom")
            ticks = np.arange(0, ext[1], 2)
            ax_.set_xticks(ticks); ax_.set_yticks(np.arange(0, ext[3], 2))
            ax_.tick_params(length=2, labelsize=8.5)
            if c == 0:
                ax_.set_ylabel("θ, circumferential (mm)", fontsize=9.5)
                ax_.text(-0.32, 0.5, rname, transform=ax_.transAxes, rotation=90, va="center", ha="center",
                         fontsize=12, color=INK, weight="bold")
            else:
                ax_.set_yticklabels([])
            if r == 2:
                ax_.set_xlabel("z, axial (mm)", fontsize=9.5)
            inn = p["inner"] & np.isfinite(p["Jn"])
            rows.append(dict(protocol=d, ratio_axial_circ=RATIO[d],
                             programmed_axial_pct=prog[0] if prog else np.nan,
                             programmed_circ_pct=prog[1] if prog else np.nan, t_s=p["t"], mean_lambda_z=p["lz"],
                             mean_lambda_theta=p["lt"], low_spot_J_over_mean=val,
                             interior_sd_J_over_mean=float(np.nanstd(p["Jn"][inn])),
                             interior_pct_below_0p9=float(np.mean(p["Jn"][inn] < 0.9) * 100),
                             interior_edge_z_mm=ez, interior_edge_theta_mm=et))
    ax_ = ax[0, 0]
    ax_.annotate(f"low spot: lowest {a.box[0]:g} × {a.box[1]:g} mm region\n(J about {deficit:.0f}% below the mean)",
                 xy=((bz0 + bz1) / 2, bt0 if bt0 > ext[2] + 2.5 else bt1),
                 xytext=(0.25, 1.6 if bt0 > ext[2] + 2.5 else ext[3] - 1.2), fontsize=9.5, color=INK, ha="left",
                 arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2),
                 bbox=dict(boxstyle="round,pad=0.3", fc=SURF, ec=GRID))
    cax = fig.add_axes([0.885, 0.30, 0.015, 0.45])
    cb = fig.colorbar(im, cax=cax, extend="both")
    cb.outline.set_visible(False)
    cb.set_label("areal stretch J / specimen mean", color=INK, fontsize=10)
    cb.set_ticks(np.round(np.linspace(a.vmin, a.vmax, 5), 2))
    cb.ax.tick_params(labelsize=9, length=2)
    cax.text(2.6, 0.02, "less\ndeformation", transform=cax.transAxes, fontsize=8.5, color=MUTED, va="bottom")
    cax.text(2.6, 0.98, "more\ndeformation", transform=cax.transAxes, fontsize=8.5, color=MUTED, va="top")
    sample = T.config.get("sample", f"{T.sample} {T.name}")
    fig.text(0.1, 0.965, f"{sample}: local areal stretch at each protocol's peak", fontsize=14,
             weight="bold", color=INK)
    fig.text(0.1, 0.943, f"J = λ_1·λ_2 divided by the interior mean. Outlined: the lowest {a.box[0]:g} × {a.box[1]:g} mm "
             f"region at {a.ref} (number = its mean J/mean in each protocol).\n"
             f"Ratios and programmed rake stretches (%) are axial : circumferential. Strain window: {wtxt}. "
             f"Gray: untracked points" + (" and excluded artifact regions" if a.exclude else "")
             + f".\nDashed: interior, {ez:g} mm from the left/right and {et:g} mm from the top/bottom tine lines. "
             "Arrowheads: rest-state tine positions.",
             fontsize=9.5, color=MUTED, va="top")
    tag = f"_box{a.box_px}px" if a.box_px else (f"_w{a.window_mm:g}mm" if a.window_mm else "")
    tag += f"_edge{ez:g}x{et:g}mm" if a.edge_mm else ""
    out = os.path.join(T.results, "peak_J_over_mean_Def11-19" + tag + ".png")
    os.makedirs(T.results, exist_ok=True)
    fig.savefig(out, dpi=110, facecolor=SURF)
    df = pd.DataFrame(rows)
    df.to_csv(out[:-4] + ".csv", index=False)
    print(f"low spot center z {(bz0 + bz1) / 2:.1f} mm, theta {(bt0 + bt1) / 2:.1f} mm "
          f"(z {bz0:.1f}-{bz1:.1f}, theta {bt0:.1f}-{bt1:.1f})")
    print(df.round(3).to_string(index=False))
    print("wrote", out, "and .csv")


if __name__ == "__main__":
    main()
