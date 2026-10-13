"""
Animate DIC fields over a tracked deformation (MP4): displacement, stretch, strain.

One MP4 per quantity. Top row: three maps drawn on the deformed image, with
color limits fixed for the whole sequence. Bottom row: the region-of-interest
(ROI) mean over time and the rake force against the matching measure, each with a
marker for the current frame. Outputs go next to the .npz as
<tag>_<quantity>_animation.mp4.

    displacement : u_z, u_theta, |u|  (mm if scale_calibration.json exists, else px);
                   lower panels: engineering strain (change in length / starting length)
                   of the tracked region and of the rake separation, each over its own length
    stretch      : lambda_z, lambda_theta, J = det F
    strain       : E_zz, E_thetatheta, E_ztheta (Green–Lagrange)

Axes: image x = axial (z); circumferential theta increases upward (theta = -image y).
Maps are drawn on the photo, so the picture keeps its image orientation.

    ../.venv/bin/python animate_fields.py ../N1739_LFPA/test1_08252026_PLC25_38-20/Images/Def11_Equibiaxial/strain_results/Def11_Equibiaxial_tine_region.npz
    ../.venv/bin/python animate_fields.py <npz> --quantity stretch strain --fps 30
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import cv2  # noqa: E402
import imageio.v2 as imageio  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

QUANTITIES = ("displacement", "stretch", "strain")
SPEC = {
    "displacement": dict(
        labels=[r"$u_z$, axial (px)", r"$u_\theta$, circumferential (px)", r"$|u|$, magnitude (px)"],
        short=[r"$|u_z|$", r"$|u_\theta|$", r"magnitude $|u|$"], cmaps=["diverging", "diverging", "viridis"], base=[None, None, 0.0],
        ylabel="ROI mean of absolute displacement (px)", mean_abs=True),
    "stretch": dict(
        labels=[r"$\lambda_z$, axial stretch", r"$\lambda_\theta$, circumferential stretch", "J = det F, area ratio"],
        short=[r"$\lambda_z$", r"$\lambda_\theta$", "J"], cmaps=["viridis", "viridis", "viridis"], base=[1.0, 1.0, 1.0],
        ylabel="stretch / area ratio", mean_abs=False),
    "strain": dict(
        labels=[r"$E_{zz}$, axial strain", r"$E_{\theta\theta}$, circumferential strain", r"$E_{z\theta}$, shear strain"],
        short=[r"$E_{zz}$", r"$E_{\theta\theta}$", r"$E_{z\theta}$"], cmaps=["viridis", "viridis", "diverging"], base=[0.0, 0.0, None],
        ylabel="Green–Lagrange strain", mean_abs=False),
}
COLS = ["#1f77b4", "#d62728", "#2ca02c"]
ALPHA = .72              # map opacity over the photo

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("npz")
ap.add_argument("--quantity", nargs="+", default=["all"], choices=QUANTITIES + ("all",))
ap.add_argument("--win", type=int, default=None,
                help="square strain window half-width in grid points (default: the Gaussian window_mm "
                     "of analysis_settings.json)")
ap.add_argument("--fps", type=int, default=20, help="MP4 frame rate (data are 5 Hz)")
a = ap.parse_args()
quantities = list(QUANTITIES) if "all" in a.quantity else a.quantity

npz = os.path.abspath(a.npz)
TEST_ROOT = testcfg.get(npz).root
out_dir, tag = os.path.dirname(npz), os.path.splitext(os.path.basename(npz))[0]
R = np.load(npz, allow_pickle=True)
if "branch" in R and str(R["branch"]) == "loading":
    print("track holds the loading branch only (track_protocols.py --loading-only)")
X0, shape, times = R["X0"], tuple(R["shape"]), R["times"]
pos, valid = R["pos"], R["valid"]
T = len(times)
img_dir = os.path.dirname(out_dir)                      # Images/DefNN_*/
img_paths = [os.path.join(img_dir, p) if os.path.exists(os.path.join(img_dir, p))
             else glob.glob(os.path.join(TEST_ROOT, "Images", "**", p), recursive=True)[0]
             for p in R["paths"]]

SCALE = D.load_scale()                                   # from calibrate_scale.py
MPP = SCALE["mm_per_px"] if SCALE else 1.0
UNIT = "mm" if SCALE else "px"
SPEC["displacement"]["labels"] = [rf"$u_z$, axial displacement ({UNIT})",
                                  rf"$u_\theta$, circumferential displacement ({UNIT})",
                                  rf"$|u|$, displacement magnitude ({UNIT})"]
PROTO = D.protocol_title(tag, TEST_ROOT)
TP = times - times[0]                                    # protocol time (s)

# change in ROI axial length (width) / circumferential length (height), from the mean
# of the edge columns / rows, and change in rake separation
P = pos.reshape((T,) + shape + (2,))
dW = (np.nanmean(P[:, :, -1, 0], 1) - np.nanmean(P[:, :, 0, 0], 1)) * MPP
dH = (np.nanmean(P[:, -1, :, 1], 1) - np.nanmean(P[:, 0, :, 1], 1)) * MPP
W0, H0 = dW[0], dH[0]                                   # tracked-region length at the protocol's first frame
dW, dH = dW - W0, dH - H0
Dd = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "*Data.csv"))[0], encoding="latin1")
Dd = Dd.drop_duplicates("Time_S").set_index("Time_S").reindex(np.round(times, 1))
rake_dX = ((Dd.X_Size_um - Dd.X_Size_um.iloc[0]) / 1000).to_numpy(float)   # mm
rake_dY = ((Dd.Y_Size_um - Dd.Y_Size_um.iloc[0]) / 1000).to_numpy(float)
RX0, RY0 = Dd.X_Size_um.iloc[0] / 1000, Dd.Y_Size_um.iloc[0] / 1000         # rake separation at start (mm)
# engineering strain dL / L0, each over its own starting length, so tissue and rakes compare directly
eW, eH, eRX, eRY = dW / W0, dH / H0, rake_dX / RX0, rake_dY / RY0
Dall = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "*Data.csv"))[0], encoding="latin1").drop_duplicates("Time_S")
_lx, _ly = D.rake_stretch(Dall.set_index("Time_S"), times, D.reference_time(R, times))
J_RAKE = _lx * _ly                                       # J of a uniform specimen following the rakes

# ---------------------------------------------------------------- fields
Q = {q: np.full((T,) + shape + (3,), np.nan) for q in QUANTITIES}
bar = D.Progress(T, "fields")
for k in range(T):
    u = (pos[k] - X0).reshape(shape + (2,)) * MPP
    u[..., 1] *= D.THETA_SIGN
    Q["displacement"][k] = np.stack([u[..., 0], u[..., 1], np.hypot(u[..., 0], u[..., 1])], -1)
    F, ok, WLABEL = D.strain_window_F(X0, pos[k], valid[k], shape, float(R["step"]), MPP, a.win)
    F = D.to_specimen(F)
    C = np.einsum("...ki,...kj->...ij", F, F)
    Q["stretch"][k] = np.stack([np.sqrt(C[..., 0, 0]), np.sqrt(C[..., 1, 1]), np.linalg.det(F)], -1)
    E = 0.5 * (C - np.eye(2))
    Q["strain"][k] = np.stack([E[..., 0, 0], E[..., 1, 1], E[..., 0, 1]], -1)
    Q["stretch"][k][~ok] = np.nan
    Q["strain"][k][~ok] = np.nan
    bar.update(k + 1)
bar.close()
k_peak = int(np.nanargmax(np.nanmean(Q["strain"][..., 0].reshape(T, -1), 1)))

# CellScale 4-point and force traces, if compare_cellscale.py has been run
cmp_csv = os.path.join(out_dir, f"{tag}_comparison.csv")
CS = pd.read_csv(cmp_csv) if os.path.exists(cmp_csv) else None
FZ, FT = D.FORCE_KEYS
has_f = CS is not None and FZ in CS

# fixed view: union of deformed grid over all frames
pad = 70
cx0, cy0 = int(max(np.nanmin(pos[..., 0]) - pad, 0)), int(max(np.nanmin(pos[..., 1]) - pad, 0))
cx1, cy1 = int(np.nanmax(pos[..., 0]) + pad), int(np.nanmax(pos[..., 1]) + pad)


def build(q):
    """Create the figure for quantity q; returns render(k) -> RGB uint8 array."""
    sp, V = SPEC[q], Q[q]
    flat = V.reshape(T, -1, 3)
    mean = np.nanmean(np.abs(flat) if sp["mean_abs"] else flat, 1)
    lims = []
    for i in range(3):
        lo, hi = np.nanpercentile(V[k_peak, ..., i], [1, 99])
        if sp["base"][i] is None:                      # diverging, centered on 0
            m = max(abs(lo), abs(hi)); lims.append((-m, m))
        else:
            lims.append((min(sp["base"][i], lo), hi))
    cm = [D.diverging_cmap() if c == "diverging" else plt.get_cmap(c) for c in sp["cmaps"]]


    fig = plt.figure(figsize=(15, 9.2), dpi=80)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.9, 1], hspace=0.12, wspace=0.22,
                          left=0.05, right=0.97, top=0.94, bottom=0.07)
    axm = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_t = fig.add_subplot(gs[1, :2])
    ax_f = fig.add_subplot(gs[1, 2])

    # bottom left: ROI means vs time (+ CellScale 4-point where comparable);
    # displacement: change in ROI width/height vs change in rake separation
    if q == "displacement":
        series = [eW, eH]
        ax_t.plot(TP, eW, color=COLS[0], lw=2, label=f"tracked region, axial ($L_0$ {W0:.2f} {UNIT})")
        ax_t.plot(TP, eH, color=COLS[1], lw=2, label=f"tracked region, circumferential ($L_0$ {H0:.2f} {UNIT})")
        ax_t.plot(TP, eRX, color=COLS[0], lw=1.2, ls="--", label=f"rakes, axial ($L_0$ {RX0:.2f} mm)")
        ax_t.plot(TP, eRY, color=COLS[1], lw=1.2, ls="--", label=f"rakes, circumferential ($L_0$ {RY0:.2f} mm)")
        sp = dict(sp, ylabel=r"engineering strain $\Delta L / L_0$")
    else:
        series = [mean[:, i] for i in range(3)]
        for i in range(3):
            ax_t.plot(TP, mean[:, i], color=COLS[i], lw=2, label=f"DIC ROI mean {sp['short'][i]}")
        if q == "stretch":
            ax_t.plot(TP, J_RAKE, color=COLS[2], lw=1, ls=":", label="J, uniform rake")
    if CS is not None and q in ("strain", "stretch"):
        for i, c in enumerate(D.STRAIN_KEYS[:3 if q == "strain" else 2]):
            e = CS[f"cs4_{c}"]
            ax_t.plot(CS.time_s - times[0], e if q == "strain" else np.sqrt(1 + 2 * e), color=COLS[i],
                      lw=1, ls="--", alpha=.7)
        ax_t.plot([], [], "k--", lw=1, label="CellScale 4-point")
    ax_t.set_xlabel("time (s)"); ax_t.set_ylabel(sp["ylabel"])
    ax_t.grid(alpha=.3)
    ax_t.legend(loc="upper left" if q == "displacement" else "upper right", fontsize=9,
                ncol=1 if q == "displacement" else 2)
    cur_t = ax_t.axvline(TP[0], color="k", lw=1)
    dots_t = [ax_t.plot([], [], "o", color=COLS[i], ms=7)[0] for i in range(len(series))]

    # bottom right: force vs the matching measure (vs time for displacement)
    fx = fy = None
    if has_f:
        if q == "displacement":
            fx, fy = eW, eH
            ax_f.plot(eW, CS[FZ], color=COLS[0], lw=1.2, alpha=.7, label="axial, tracked region")
            ax_f.plot(eH, CS[FT], color=COLS[1], lw=1.2, alpha=.7, label="circ., tracked region")
            ax_f.plot(eRX, CS[FZ], color=COLS[0], lw=1, ls="--", alpha=.7, label="axial, rakes")
            ax_f.plot(eRY, CS[FT], color=COLS[1], lw=1, ls="--", alpha=.7, label="circ., rakes")
            ax_f.set_xlabel(r"engineering strain $\Delta L / L_0$")
        else:
            fx, fy = mean[:, 0], mean[:, 1]
            ax_f.plot(fx, CS[FZ], color=COLS[0], lw=1.2, alpha=.6, label=f"axial force vs {sp['short'][0]}")
            ax_f.plot(fy, CS[FT], color=COLS[1], lw=1.2, alpha=.6, label=f"circ. force vs {sp['short'][1]}")
            ax_f.set_xlabel(f"DIC ROI mean {q}")
        dot_fx, = ax_f.plot([], [], "o", color=COLS[0], ms=8)
        dot_fy, = ax_f.plot([], [], "o", color=COLS[1], ms=8)
        ax_f.set_ylabel("rake force (mN)"); ax_f.grid(alpha=.3); ax_f.legend(fontsize=9, loc="upper left")
    else:
        ax_f.axis("off")

    meshes, ims, cbs = [], [], []
    for i, ax in enumerate(axm):
        ims.append(ax.imshow(np.zeros((cy1 - cy0, cx1 - cx0, 3), np.uint8),
                             extent=(cx0, cx1, cy1, cy0)))
        m = ax.pcolormesh(pos[0, :, 0].reshape(shape), pos[0, :, 1].reshape(shape), V[0, ..., i],
                          cmap=cm[i], vmin=lims[i][0], vmax=lims[i][1], alpha=ALPHA, shading="nearest")
        meshes.append(m)
        cbs.append(fig.colorbar(m, ax=ax, shrink=.85, pad=.02))
        ax.set_xlim(cx0, cx1); ax.set_ylim(cy1, cy0); ax.set_xticks([]); ax.set_yticks([])
    title = fig.suptitle("", fontsize=14)
    map_titles = [ax.set_title("", fontsize=11) for ax in axm]
    fmt = ("{:.3f}" if SCALE else "{:.1f}") if q == "displacement" else ("{:.3f}" if q == "stretch" else "{:+.3f}")

    def render(k):
        img = cv2.cvtColor(cv2.imread(img_paths[k]), cv2.COLOR_BGR2RGB)[cy0:cy1, cx0:cx1]
        for i, ax in enumerate(axm):
            ims[i].set_data(img)
            meshes[i].remove()
            meshes[i] = ax.pcolormesh(pos[k, :, 0].reshape(shape), pos[k, :, 1].reshape(shape),
                                      V[k, ..., i], cmap=cm[i], vmin=lims[i][0], vmax=lims[i][1],
                                      alpha=ALPHA, shading="nearest")
            cbs[i].update_normal(meshes[i])
            stat = "mean |·|" if sp["mean_abs"] and i < 2 else "mean"
            map_titles[i].set_text(f"{sp['labels'][i]}\n{stat} {fmt.format(mean[k, i])}")
        title.set_text(f"{PROTO} · {q} · t = {TP[k]:5.1f} s · {D.settings_label(R, MPP, WLABEL)}")
        cur_t.set_xdata([TP[k], TP[k]])
        for i, sv in enumerate(series):
            dots_t[i].set_data([TP[k]], [sv[k]])
        if has_f:
            dot_fx.set_data([fx[k]], [CS[FZ][k]])
            dot_fy.set_data([fy[k]], [CS[FT][k]])
        fig.canvas.draw()
        return np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()

    return render


for q in quantities:
    render = build(q)
    mp4 = os.path.join(out_dir, f"{tag}_{q}_animation.mp4")
    writer = imageio.get_writer(mp4, fps=a.fps, codec="libx264", quality=8,
                                pixelformat="yuv420p", macro_block_size=16)
    bar = D.Progress(T, q)
    for k in range(T):
        writer.append_data(render(k))
        bar.update(k + 1)
    writer.close()
    bar.close()
    plt.close("all")
    print("wrote", mp4, f"({os.path.getsize(mp4)/1e6:.1f} MB)")
