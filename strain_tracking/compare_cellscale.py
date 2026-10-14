"""
Compare DIC tracking (track.py output) against the CellScale 4-point tracking
in Coords1.csv, and against rake (grip) strain from the test Data.csv.

Strain measures compared (all Green-Lagrange; image x = tester X = axial z,
circumferential theta = -image y, i.e. increasing upward; dic2d.to_specimen):
  cellscale_4pt : bilinear F at element center from the 4 CellScale points
  dic_4pt       : same formula, applied to the same 4 material points as
                  tracked by DIC (DIC displacement interpolated at the
                  CellScale reference positions)
  dic_field_quad: mean of the local DIC strain field inside the 4-point quad
  rake          : X_Size / X_Size0 and Y_Size / Y_Size0 from Data.csv

CSV columns: <method>_E_zz, <method>_E_thetatheta, <method>_E_ztheta for methods
cs4, dic4, dicquad, dicroi, rake; ptdiff_<point>_px; force_z_mN, force_theta_mN.

Example:
    .venv/bin/python strain_tracking/compare_cellscale.py N1739_LFPA/test1_08252026_PLC25_38-20/DIC_results/Def11_Equibiaxial.npz

Outputs are written next to the .npz file.
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
from matplotlib.path import Path  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("npz")
ap.add_argument("--coords", default=None, help="CellScale Coords1.csv (default: the test's)")
ap.add_argument("--data", default=None, help="CellScale *Data.csv (auto-detected)")
ap.add_argument("--win", type=int, default=None,
                help="square strain window half-width in grid points (default: the Gaussian window_mm "
                     "of analysis_settings.json)")
a = ap.parse_args()
TEST_ROOT = testcfg.get(a.npz).root
a.coords = a.coords or os.path.join(TEST_ROOT, "Coords1.csv")
out_dir = os.path.dirname(os.path.abspath(a.npz))
tag = os.path.splitext(os.path.basename(a.npz))[0]

R = np.load(a.npz, allow_pickle=True)
LOADING_ONLY = "branch" in R and str(R["branch"]) == "loading"
if LOADING_ONLY:
    print("track holds the loading branch only (track_protocols.py --loading-only)")
MPP_ = D.load_scale()["mm_per_px"]
X0, shape, times = R["X0"], tuple(R["shape"]), R["times"]
pos, valid, zncc = R["pos"], R["valid"], R["zncc"]
T = len(times)

# ---------------------------------------------------------------- CellScale
C = pd.read_csv(a.coords, skiprows=1)
C = C.set_index("FileName").loc[list(R["paths"])]
P_cs = C[["X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4"]].to_numpy().reshape(T, 4, 2)
ORDER = [0, 2, 3, 1]          # points 1,3,4,2 = TL, TR, BR, BL (counter-clockwise)
P_cs = P_cs[:, ORDER]

data_csv = a.data or glob.glob(os.path.join(TEST_ROOT, "*Data.csv"))[0]
Dd = pd.read_csv(data_csv, encoding="latin1").drop_duplicates("Time_S").set_index("Time_S")
Dd = Dd.reindex(np.round(times, 1))
lam_rake = np.stack([Dd.X_Size_um / Dd.X_Size_um.iloc[0],
                     Dd.Y_Size_um / Dd.Y_Size_um.iloc[0]], 1)

# ----------------------------------------------------- DIC at the 4 points
xs = X0[:, 0].reshape(shape)[0]
ys = X0[:, 1].reshape(shape)[:, 0]
P_dic = np.zeros_like(P_cs)
for k in range(T):
    for c in range(2):
        f = RegularGridInterpolator((ys, xs), pos[k, :, c].reshape(shape))
        P_dic[k, :, c] = f(P_cs[0][:, ::-1])
P_dic[0] = P_cs[0]

def E_of(F):
    E = D.green_lagrange(D.to_specimen(F))
    return np.array([E[0, 0], E[1, 1], E[0, 1]])

E_cs = np.array([E_of(D.bilinear_F_center(P_cs[0], P_cs[k])) for k in range(T)])
E_d4 = np.array([E_of(D.bilinear_F_center(P_dic[0], P_dic[k])) for k in range(T)])
E_rk = np.stack([0.5 * (lam_rake[:, 0] ** 2 - 1), 0.5 * (lam_rake[:, 1] ** 2 - 1),
                 np.zeros(T)], 1)

# ------------------------------------------------------- DIC strain field
inside = Path(P_cs[0]).contains_points(X0).reshape(shape)
E_fq = np.full((T, 3), np.nan)
E_roi = np.full((T, 3), np.nan)
fields = {}
k_peak = int(np.nanargmax(E_cs[:, 0]))
for k in range(T):
    F, ok, WLABEL = D.strain_window_F(X0, pos[k], valid[k], shape, float(R["step"]), MPP_, a.win)
    E = D.green_lagrange(D.to_specimen(F))
    comp = np.stack([E[..., 0, 0], E[..., 1, 1], E[..., 0, 1]], -1)
    E_fq[k] = np.nanmean(comp[inside & ok], 0)
    E_roi[k] = np.nanmean(comp[ok], 0)
    if k == k_peak:
        fields = dict(E=comp, ok=ok)

# -------------------------------------------------------------- summaries
names = list(D.STRAIN_KEYS)                           # E_zz, E_thetatheta, E_ztheta
labels = [r"$E_{zz}$, axial strain", r"$E_{\theta\theta}$, circumferential strain", r"$E_{z\theta}$, shear strain"]
tips = testcfg.get().setting("tine_tips_px")
proto = D.protocol_title(tag, TEST_ROOT)
run = D.settings_label(R, MPP_, WLABEL)
tp = times - times[0]                                  # protocol time (s)
d_pts = np.linalg.norm(P_dic - P_cs, axis=2)          # (T,4) px
rows = []
for lab, E in [("CellScale 4-pt", E_cs), ("DIC 4-pt (same points)", E_d4),
               ("DIC field mean in quad", E_fq), ("DIC field mean, whole ROI", E_roi),
               ("Rake (grip)", E_rk)]:
    rows.append([lab, *E[k_peak]])
summary = pd.DataFrame(rows, columns=["method", *names]).set_index("method")
diff = pd.DataFrame({
    "DIC4pt - CS (RMS over time)": np.sqrt(np.nanmean((E_d4 - E_cs) ** 2, 0)),
    "DIC4pt - CS (max |diff|)": np.nanmax(np.abs(E_d4 - E_cs), 0),
    "Field-quad - CS (RMS)": np.sqrt(np.nanmean((E_fq - E_cs) ** 2, 0)),
}, index=names).T

closure_cs = np.linalg.norm(P_cs[-1] - P_cs[0], axis=1)
closure_dic = np.linalg.norm(P_dic[-1] - P_dic[0], axis=1)
pd.set_option("display.width", 140)
print(f"Peak frame (max CellScale E_zz): t = {times[k_peak]:.1f} s")
print(summary.round(4).to_string())
print("\nAgreement over the full sequence (strain):")
print(diff.round(4).to_string())
print(f"\nPoint-position difference DIC vs CellScale (px): "
      f"median {np.median(d_pts):.2f}, 95th pct {np.percentile(d_pts, 95):.2f}, "
      f"max {d_pts.max():.2f}")
if LOADING_ONLY:
    print("Closure error: not computed (the track ends at the rake peak, not back at rest)")
else:
    print("Closure error, last vs first frame (px), points TL,TR,BR,BL:")
    print("  CellScale", np.round(closure_cs, 2), "  DIC", np.round(closure_dic, 2))
print(f"DIC valid points per frame: median {np.median(valid.mean(1))*100:.1f}%, "
      f"min {valid.mean(1).min()*100:.1f}%; keyframes {len(R['keyframes'])}")

ts = pd.DataFrame({"time_s": times})
for lab, E in [("cs4", E_cs), ("dic4", E_d4), ("dicquad", E_fq),
               ("dicroi", E_roi), ("rake", E_rk)]:
    for i, n in enumerate(D.STRAIN_KEYS):
        ts[f"{lab}_{n}"] = E[:, i]
for j, n in enumerate(["TL", "TR", "BR", "BL"]):
    ts[f"ptdiff_{n}_px"] = d_pts[:, j]
ts[D.FORCE_KEYS[0]], ts[D.FORCE_KEYS[1]] = Dd.X_Force_mN.values, Dd.Y_Force_mN.values
summary.to_csv(os.path.join(out_dir, f"{tag}_peak_summary.csv"))

# ------------------------------------------------------------------ plots
ts.insert(1, "t_protocol_s", tp)
ts.to_csv(os.path.join(out_dir, f"{tag}_comparison.csv"), index=False)
cols = {"CellScale 4-pt": "#d62728", "DIC 4-pt": "#1f77b4",
        "DIC field (quad)": "#2ca02c", "Rake": "#7f7f7f"}
fig, ax = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
for i in range(3):
    a0 = ax[0, i]
    a0.plot(tp, E_rk[:, i], color=cols["Rake"], ls=":", label="Rake")
    a0.plot(tp, E_cs[:, i], color=cols["CellScale 4-pt"], lw=2, label="CellScale 4-pt")
    a0.plot(tp, E_d4[:, i], color=cols["DIC 4-pt"], lw=1.2, ls="--", label="DIC 4-pt")
    a0.plot(tp, E_fq[:, i], color=cols["DIC field (quad)"], lw=1, label="DIC field mean (quad)")
    a0.set_title(labels[i]); a0.set_xlabel("time (s)"); a0.grid(alpha=.3)
    a1 = ax[1, i]
    a1.plot(tp, E_d4[:, i] - E_cs[:, i], color=cols["DIC 4-pt"], label="DIC 4-pt − CellScale")
    a1.plot(tp, E_fq[:, i] - E_cs[:, i], color=cols["DIC field (quad)"], label="DIC field − CellScale")
    a1.axhline(0, color="k", lw=.5); a1.set_xlabel("time (s)"); a1.set_ylabel("Δ strain")
    a1.grid(alpha=.3)
ax[0, 0].legend(fontsize=8); ax[1, 0].legend(fontsize=8)
fig.suptitle(f"{proto} · Green–Lagrange strain, DIC vs CellScale 4-point tracking · {run}", fontsize=13)
fig.savefig(os.path.join(out_dir, f"{tag}_strain_vs_time.png"), dpi=130)

fig, ax = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
# d_pts columns are TL, TR, BR, BL = CellScale points 1, 3, 4, 2; plot in point order 1-4
for j, n in [(0, "Point 1, top left"), (3, "Point 2, bottom left"), (1, "Point 3, top right"),
             (2, "Point 4, bottom right")]:
    ax[0].plot(tp, d_pts[:, j] * MPP_ * 1000, label=n)
ax[0].set_ylabel("|x_DIC − x_CellScale| (µm)"); ax[0].set_xlabel("time (s)")
ax[0].set_title("Same material points: position disagreement"); ax[0].legend(); ax[0].grid(alpha=.3)
ax[1].plot(tp, valid.mean(1) * 100, label="valid points (%)")
ax2 = ax[1].twinx(); ax2.plot(tp, np.nanmedian(zncc, 1), color="C1", label="median ZNCC")
for kf in R["keyframes"]:
    ax[1].axvline(tp[kf], color="k", lw=.4, alpha=.4)
ax[1].set_ylabel("valid (%)"); ax2.set_ylabel("median ZNCC", color="C1")
ax[1].set_xlabel("time (s)"); ax[1].set_title("DIC quality (vertical lines = keyframes)")
fig.suptitle(f"{proto} · tracking quality · {run}", fontsize=12)
fig.savefig(os.path.join(out_dir, f"{tag}_points_quality.png"), dpi=130)

# strain maps at peak, drawn on the deformed image
import cv2  # noqa: E402
img_path = [p for p in glob.glob(os.path.join(TEST_ROOT, "Images", "**", R["paths"][k_peak]),
                                 recursive=True)][0]
img = cv2.cvtColor(cv2.imread(img_path), cv2.COLOR_BGR2RGB)
xk = pos[k_peak]
pad = 80 * MPP_
Xd, Yd = D.specimen_mm(xk[:, 0].reshape(shape), xk[:, 1].reshape(shape), tips, MPP_)   # origin: rest-state tine-tip corner
x0, x1 = np.nanmin(Xd) - pad, np.nanmax(Xd) + pad
y0, y1 = np.nanmin(Yd) - pad, np.nanmax(Yd) + pad
h, w = img.shape[:2]
(ex0, ex1), (ey1, ey0) = D.specimen_mm([0, w], [0, h], tips, MPP_)   # image extent in specimen mm
qz, qt = D.specimen_mm(*np.vstack([P_cs[k_peak], P_cs[k_peak][:1]]).T, tips, MPP_)


def cell_corners(A):
    """Corners (ny+1, nx+1) of the cells centered on a curvilinear grid A (ny, nx): the borders extrapolated
    linearly, then the mean of each 2 x 2 block (pcolormesh needs corners on a deformed grid)."""
    P = np.empty((A.shape[0] + 2, A.shape[1] + 2))
    P[1:-1, 1:-1] = A
    P[0, 1:-1], P[-1, 1:-1] = 2 * A[0] - A[1], 2 * A[-1] - A[-2]
    P[:, 0], P[:, -1] = 2 * P[:, 1] - P[:, 2], 2 * P[:, -2] - P[:, -3]
    return 0.25 * (P[:-1, :-1] + P[1:, :-1] + P[:-1, 1:] + P[1:, 1:])


Xe, Ye = cell_corners(Xd), cell_corners(Yd)

fig, ax = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
for i in range(3):
    Ei = np.where(fields["ok"], fields["E"][..., i], np.nan)
    lo, hi = np.nanpercentile(Ei, [2, 98])
    if i == 2:
        hi = max(abs(lo), abs(hi)); lo = -hi
    ax[i].imshow(img, extent=(ex0, ex1, ey0, ey1))
    m = ax[i].pcolormesh(Xe, Ye, Ei, cmap=D.diverging_cmap() if i == 2 else "viridis",
                         vmin=lo, vmax=hi, shading="flat", edgecolors="face", linewidth=0.3)
    ax[i].plot(qz, qt, "r-o", ms=4, lw=1)
    # CellScale point numbers; corners TL, TR, BR, BL are points 1, 3, 4, 2, labeled outside the quad
    for c, (num, dz, dt, ha, va) in enumerate([(1, -0.25, 0.25, "right", "bottom"), (3, 0.25, 0.25, "left", "bottom"),
                                               (4, 0.25, -0.25, "left", "top"), (2, -0.25, -0.25, "right", "top")]):
        ax[i].text(qz[c] + dz, qt[c] + dt, str(num), color="r", fontsize=11, weight="bold", ha=ha, va=va,
                   bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    ax[i].set_xlim(x0, x1); ax[i].set_ylim(y0, y1); ax[i].set_aspect("equal")
    ax[i].set_xlabel("z, axial (mm)"); ax[i].set_ylabel(r"$\theta$, circumferential (mm)")
    ax[i].set_title(f"{labels[i]} (CellScale 4-pt: {E_cs[k_peak, i]:.3f})", fontsize=11)
    fig.colorbar(m, ax=ax[i], shrink=.8)
fig.suptitle(f"{proto} · peak at t = {tp[k_peak]:.1f} s, deformed configuration · {run} · "
             "red: CellScale points 1–4", fontsize=13)
fig.savefig(os.path.join(out_dir, f"{tag}_strain_maps_peak.png"), dpi=110)
print("\nwrote", out_dir)
