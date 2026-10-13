"""
Checks that the E_zz bands seen in the DIC strain maps are not produced by the code.

T1  Null test: real frame-0 texture, *uniform* stretch equal to the measured peak
    strain, ramped over the frames to the rake peak at the real 5 Hz increment, JPEG + noise.
    The full pipeline must return a flat strain map.
T2  Known-band test: same, plus imposed vertical E_zz bands and a horizontal E_θθ
    band at known positions (different from the real ones). They must be recovered
    in the right place and orientation.
T3  Strain operator: deformation_gradient_field vs brute-force per-point lstsq.
T4  Independent tracker: cv2.matchTemplate from frame 0 straight to the peak frame
    (after unwarping with the CellScale 4-point affine only), strain by finite
    differences. Shares no code with dic2d.

Runs on Def01 of the selected test: reference = its first frame, peak = the frame at the rake
peak; grid = the test's roi_px.standard.

    ../.venv/bin/python validate_bands.py --test N1739/test1   (needs Def01 tracked first: track.py --def Def01)
"""
import os
import sys
import tempfile

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.special import erf  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

TC = testcfg.get()
TEST_ROOT = TC.root
DEF1 = os.path.join(TC.images, TC.folder("Def01"))
OUT = os.path.join(DEF1, "strain_results")
DEF1_TRACK, DEF1_NPZ = D.load_track(TC.images, os.path.basename(DEF1), "standard")   # sliced from a full-domain track if needed
_t, _p = D.find_frames([DEF1])
_Dd = pd.read_csv(TC.first("*Data.csv"), encoding="latin1").drop_duplicates("Time_S").set_index("Time_S")
_k_peak = int(np.nanargmax(_Dd.X_Size_um.reindex(np.round(_t, 1)).to_numpy(float)))   # rake peak
REF, PEAK = _p[0], _p[_k_peak]
ROI, STEP = tuple(TC.setting("roi_px", "standard")), D.DEFAULT_STEP
WIN = D.window_points(STEP)
F_PEAK = np.array([[1.371, 0.0], [0.0, 1.177]])   # ~ measured peak E_zz 0.44, E_θθ 0.19
C0 = np.array([1020.0, 1060.0])
NFR = _k_peak + 1                                   # first frame to the rake peak
rng = np.random.default_rng(1)


def band_u(X, centers, sig, amp):
    """Displacement whose derivative is amp * sum of Gaussians (exact integral)."""
    return sum(amp * sig * np.sqrt(np.pi / 2) * erf((X - c) / (np.sqrt(2) * sig))
               for c in centers)


def band_du(X, centers, sig, amp):
    return sum(amp * np.exp(-0.5 * ((X - c) / sig) ** 2) for c in centers)


def make_sequence(outdir, bands):
    ref = cv2.imread(REF)
    H, W = ref.shape[:2]
    x1, y1 = np.arange(W, dtype=np.float64), np.arange(H, dtype=np.float64)
    paths = []
    for k in range(NFR):
        s = k / (NFR - 1)
        fx = 1 + s * (F_PEAK[0, 0] - 1)
        fy = 1 + s * (F_PEAK[1, 1] - 1)
        # forward map (F_PEAK is diagonal, so x and y separate):
        #   x = C0x + fx (X - C0x) + s*ub(X),  y = C0y + fy (Y - C0y) + s*vb(Y)
        # inverted per axis by fixed-point iteration
        X, Y = x1.copy(), y1.copy()
        for _ in range(25):
            bx = s * band_u(X, *bands["x"]) if bands else 0
            by = s * band_u(Y, *bands["y"]) if bands else 0
            X = C0[0] + (x1 - C0[0] - bx) / fx
            Y = C0[1] + (y1 - C0[1] - by) / fy
        MX, MY = np.meshgrid(X, Y)
        im = cv2.remap(ref, MX.astype(np.float32), MY.astype(np.float32), cv2.INTER_CUBIC)
        im = np.clip(im + rng.normal(0, 2, im.shape), 0, 255).astype(np.uint8)
        p = os.path.join(outdir, f"s.{k:06d}.000.jpg")
        cv2.imwrite(p, im, [cv2.IMWRITE_JPEG_QUALITY, 92])
        paths.append(p)
    return paths


def run_pipeline(paths):
    X0, shape = D.make_grid(ROI, STEP)
    res = D.track_sequence(paths, X0, D.TrackParams(), verbose=False)
    F, ok = D.deformation_gradient_field(X0, res["pos"][-1], res["valid"][-1], shape, WIN)
    E = D.green_lagrange(F)
    return X0, shape, np.where(ok, E[..., 0, 0], np.nan), np.where(ok, E[..., 1, 1], np.nan), res


def colprof(E):
    return np.nanmean(E, 0)


def real_map():
    R = DEF1_TRACK
    sh = tuple(R["shape"])
    k = int(np.argmin(abs(R["times"] - _t[_k_peak])))
    F, ok = D.deformation_gradient_field(R["X0"], R["pos"][k], R["valid"][k], sh, WIN)
    E = D.green_lagrange(F)
    return np.where(ok, E[..., 0, 0], np.nan), np.where(ok, E[..., 1, 1], np.nan)


X0, shape = D.make_grid(ROI, STEP)
xs, ys = X0[:, 0].reshape(shape)[0], X0[:, 1].reshape(shape)[:, 0]
E_zz_real, E_tt_real = real_map()
pr_real = colprof(E_zz_real)
print(f"REAL data      : E_zz column-profile range {np.nanmin(pr_real):.3f}-{np.nanmax(pr_real):.3f} "
      f"(p2p {np.nanmax(pr_real)-np.nanmin(pr_real):.3f}), map sd {np.nanstd(E_zz_real):.4f}")

# ---------------------------------------------------------------- T1
with tempfile.TemporaryDirectory() as td:
    _, _, E_zz_t1, E_tt_t1, r1 = run_pipeline(make_sequence(td, None))
Et = D.green_lagrange(F_PEAK)
p1 = colprof(E_zz_t1)
print(f"T1 null test   : true E_zz {Et[0,0]:.4f} E_θθ {Et[1,1]:.4f} | recovered mean "
      f"{np.nanmean(E_zz_t1):.4f} / {np.nanmean(E_tt_t1):.4f}; map sd {np.nanstd(E_zz_t1):.4f} / "
      f"{np.nanstd(E_tt_t1):.4f}; column-profile p2p {np.nanmax(p1)-np.nanmin(p1):.4f}; "
      f"keyframes {len(r1['keyframes'])}; valid {r1['valid'][-1].mean()*100:.0f}%")
v = np.isfinite(E_zz_t1) & np.isfinite(E_zz_real)
print(f"               corr(null-test map, real map) r = {np.corrcoef(E_zz_t1[v], E_zz_real[v])[0,1]:.2f}")

# ---------------------------------------------------------------- T2
bands = {"x": ([950.0, 1170.0], 25.0, 0.10), "y": ([1100.0], 25.0, 0.06)}
with tempfile.TemporaryDirectory() as td:
    _, _, E_zz_t2, E_tt_t2, _ = run_pipeline(make_sequence(td, bands))
Xg, Yg = X0[:, 0].reshape(shape), X0[:, 1].reshape(shape)
Ftrue = np.zeros(shape + (2, 2))
Ftrue[..., 0, 0] = F_PEAK[0, 0] + band_du(Xg, *bands["x"])
Ftrue[..., 1, 1] = F_PEAK[1, 1] + band_du(Yg, *bands["y"])
Etrue = D.green_lagrange(Ftrue)
# what an ideal measurement with this window would see: box-average F over the strain window
from scipy import ndimage  # noqa: E402
Fs = Ftrue.copy()
for i in range(2):
    Fs[..., i, i] = ndimage.uniform_filter(Ftrue[..., i, i], 2 * WIN + 1, mode="nearest")
Esm = D.green_lagrange(Fs)
pc2, pt2 = colprof(E_zz_t2), colprof(Etrue[..., 0, 0])
rp2, rt2 = np.nanmean(E_tt_t2, 1), np.nanmean(Etrue[..., 1, 1], 1)
pk = lambda prof, ax: sorted(ax[np.argsort(-np.nan_to_num(prof))[:1]].tolist())
v = np.isfinite(E_zz_t2)
print(f"T2 known bands : E_zz peaks true x = {bands['x'][0]}, recovered column maxima near "
      f"{[xs[np.nanargmax(np.where(abs(xs-c)<60, pc2, -1))] for c in bands['x'][0]]}; "
      f"E_θθ band true y = 1100, recovered row max at {ys[np.nanargmax(rp2)]}")
print(f"               map r vs window-smoothed truth: E_zz {np.corrcoef(E_zz_t2[v], Esm[...,0,0][v])[0,1]:.3f}, "
      f"E_θθ {np.corrcoef(E_tt_t2[v], Esm[...,1,1][v])[0,1]:.3f}; RMS error E_zz "
      f"{np.sqrt(np.nanmean((E_zz_t2-Esm[...,0,0])**2)):.4f}, E_θθ {np.sqrt(np.nanmean((E_tt_t2-Esm[...,1,1])**2)):.4f}")
print(f"               no leakage: E_zz row-profile p2p {np.nanmax(np.nanmean(E_zz_t2,1))-np.nanmin(np.nanmean(E_zz_t2,1)):.4f}, "
      f"E_θθ column-profile p2p {np.nanmax(np.nanmean(E_tt_t2,0))-np.nanmin(np.nanmean(E_tt_t2,0)):.4f}")

# ---------------------------------------------------------------- T3
R = DEF1_TRACK
k = int(np.argmin(abs(R["times"] - _t[_k_peak])))
pos, val = R["pos"][k], R["valid"][k]
F, ok = D.deformation_gradient_field(R["X0"], pos, val, shape, WIN)
U = (pos - R["X0"]).reshape(shape + (2,))
Vm = val.reshape(shape)
errs = []
for (i, j) in zip(rng.integers(0, shape[0], 400), rng.integers(0, shape[1], 400)):
    if not ok[i, j]:
        continue
    sl = (slice(max(i - WIN, 0), i + WIN + 1), slice(max(j - WIN, 0), j + WIN + 1))
    m = Vm[sl]
    A = np.stack([np.ones(m.sum()), Xg[sl][m], Yg[sl][m]], 1)
    sol, *_ = np.linalg.lstsq(A, U[sl][m], rcond=None)      # rows: const, d/dX, d/dY
    Fb = np.eye(2) + sol[1:].T
    errs.append(np.abs(Fb - F[i, j]).max())
print(f"T3 strain op   : max |F_fast - F_bruteforce| over {len(errs)} points = {max(errs):.2e}")

# ---------------------------------------------------------------- T4
C = pd.read_csv(os.path.join(TEST_ROOT, "Coords1.csv"), skiprows=1).set_index("FileName")
P0 = C.loc[os.path.basename(REF)].iloc[1:].to_numpy(float).reshape(4, 2)
Pk = C.loc[os.path.basename(PEAK)].iloc[1:].to_numpy(float).reshape(4, 2)
A_aff, _ = cv2.estimateAffine2D(P0, Pk, ransacReprojThreshold=1e9)   # 4-point LSQ affine
ref = cv2.imread(REF, cv2.IMREAD_GRAYSCALE).astype(np.float32)
peak = cv2.imread(PEAK, cv2.IMREAD_GRAYSCALE).astype(np.float32)
unw = cv2.warpAffine(peak, A_aff, (ref.shape[1], ref.shape[0]),
                     flags=cv2.INTER_CUBIC | cv2.WARP_INVERSE_MAP)  # peak mapped into ref frame
r, srch = D.TrackParams.radius, 14
res = np.full(shape + (2,), np.nan)
cc = np.full(shape, np.nan)
for i in range(shape[0]):
    for j in range(shape[1]):
        x, y = int(Xg[i, j]), int(Yg[i, j])
        tpl = ref[y - r:y + r + 1, x - r:x + r + 1]
        win = unw[y - r - srch:y + r + srch + 1, x - r - srch:x + r + srch + 1]
        M = cv2.matchTemplate(win, tpl, cv2.TM_CCOEFF_NORMED)
        _, mx, _, (bx, by) = cv2.minMaxLoc(M)
        if 0 < bx < M.shape[1] - 1 and 0 < by < M.shape[0] - 1:
            dxs = 0.5 * (M[by, bx-1] - M[by, bx+1]) / (M[by, bx-1] - 2*M[by, bx] + M[by, bx+1])
            dys = 0.5 * (M[by-1, bx] - M[by+1, bx]) / (M[by-1, bx] - 2*M[by, bx] + M[by+1, bx])
            res[i, j] = (bx - srch + dxs, by - srch + dys)
            cc[i, j] = mx
good = cc > 0.6
res[~good] = np.nan
# residual displacement in unwarped (= reference) frame -> F = A_lin (I + grad r)
gx = np.gradient(res[..., 0], STEP, axis=1); gy = np.gradient(res[..., 0], STEP, axis=0)
hx = np.gradient(res[..., 1], STEP, axis=1); hy = np.gradient(res[..., 1], STEP, axis=0)
G = np.stack([np.stack([1 + gx, gy], -1), np.stack([hx, 1 + hy], -1)], -2)
F4 = np.einsum("ij,...jk->...ik", A_aff[:, :2], G)
E4 = D.green_lagrange(F4)
E_zz_t4 = E4[..., 0, 0]
# smooth with the same window size so it is comparable to the DIC map
def nanbox(a, n):
    m = np.isfinite(a)
    num = ndimage.uniform_filter(np.where(m, a, 0), n, mode="constant")
    den = ndimage.uniform_filter(m.astype(float), n, mode="constant")
    return np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan)
E_zz_t4s = nanbox(E_zz_t4, 2 * WIN + 1)
v = np.isfinite(E_zz_t4s) & np.isfinite(E_zz_real)
p4 = colprof(E_zz_t4s)
print(f"T4 independent : matchTemplate valid {good.mean()*100:.0f}% | map r vs DIC = "
      f"{np.corrcoef(E_zz_t4s[v], E_zz_real[v])[0,1]:.2f}; column-profile r = "
      f"{np.corrcoef(p4[(f := np.isfinite(p4) & np.isfinite(pr_real))], pr_real[f])[0,1]:.2f}; "
      f"profile maxima x = {xs[np.nanargmax(np.where(xs<980,p4,-1))]}, {xs[np.nanargmax(np.where(xs>980,p4,-1))]}")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
kw = dict(cmap="viridis", shading="auto")
for a, E, t, lims in ((ax[0, 0], E_zz_real, "REAL: DIC E_zz at peak", (.3, .55)),
                      (ax[0, 1], E_zz_t1, "T1 null test: uniform stretch → DIC E_zz", (.3, .55)),
                      (ax[0, 2], E_zz_t4s, "T4 independent matchTemplate E_zz", (.3, .55)),
                      (ax[1, 0], E_zz_t2, "T2 imposed bands at x=950,1170 → DIC E_zz", None),
                      (ax[1, 1], E_tt_t2, "T2 imposed band at y=1100 → DIC E_θθ", None)):
    m = a.pcolormesh(xs, ys, E, vmin=lims[0] if lims else None, vmax=lims[1] if lims else None, **kw)
    a.invert_yaxis(); a.set_aspect("equal"); a.set_title(t); fig.colorbar(m, ax=a, shrink=.8)
a = ax[1, 2]
a.plot(xs, pr_real, lw=2, label="real, DIC")
a.plot(xs, p4, lw=1.5, label="real, independent matchTemplate")
a.plot(xs, p1, lw=1.5, label="T1 null test (uniform)")
a.axhline(Et[0, 0], color="k", lw=.6, ls=":")
a.set_xlabel("X (ref px)"); a.set_ylabel("column-mean E_zz"); a.legend(); a.grid(alpha=.3)
a.set_title("E_zz averaged over Y")
os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, "bands_validation.png"), dpi=95)
print("wrote", os.path.join(OUT, "bands_validation.png"))
