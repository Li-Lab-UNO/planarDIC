"""
Incremental 2D subset DIC for CellScale biaxial image sequences.

Method
------
* Subset matching: inverse-compositional Gauss-Newton (IC-GN) with a first-order
  (affine) shape function and the ZNSSD criterion (Pan et al., Exp Mech 2013;
  same formulation as Ncorr/DICe). ZNSSD is invariant to per-subset affine
  intensity changes, so vignetting / the central hot spot need no flat-fielding.
* Tracking: Lagrangian. Each grid point is a material point. Every frame is
  matched against a *keyframe* (the current reference), seeded from the previous
  frame's solution. The keyframe is only replaced when the deformation since it
  becomes large or correlation drops, so error accumulates over a handful of
  reference updates rather than every frame.
* Strain: deformed positions x(X) on the reference grid -> masked least-squares
  plane fit of displacement over a square strain window -> F = I + grad(u)
  -> Green-Lagrange E = 0.5 (F^T F - I).

All vectorized with numpy + cv2.remap; ~0.1-0.3 s per 2048x2048 frame.
"""
from __future__ import annotations

import os
import re
import shutil
import sys
import time
import warnings
from dataclasses import dataclass

import cv2
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

# OpenCV starts one thread per logical CPU by default; parallel jobs (run_pipeline.sh) share the CPUs
# through DIC_CV_THREADS. (macOS builds use Grand Central Dispatch, which only honors 0 = no threading.)
if os.environ.get("DIC_CV_THREADS", "").isdigit():
    cv2.setNumThreads(int(os.environ["DIC_CV_THREADS"]))

FNAME_RE = re.compile(r"\.(\d{6})\.(\d{3})\.jpg$")

# Specimen axes. Image x (tester X) is the artery's axial direction z, increasing to the
# right; the circumferential direction theta runs along image y (tester Y) and increases
# UPWARD in the image, so theta = -image y. Tracking and deformation_gradient_field work in
# image coordinates; to_specimen() converts F (and u_theta = THETA_SIGN * u_y) to (z, theta).
# Only the signs of F_ztheta, F_thetaz, E_ztheta, u_theta and angles depend on this;
# stretches and J do not. Export names are ASCII (valid MATLAB identifiers).
THETA_SIGN = -1
AXIS_CONVENTION = ("image x = tester X = axial (z), increasing right; "
                   "theta = circumferential, increasing upward in the image (theta = -image y)")
STRAIN_KEYS = ("E_zz", "E_thetatheta", "E_ztheta")
STRAIN_LABELS = (r"$E_{zz}$", r"$E_{\theta\theta}$", r"$E_{z\theta}$")
STRETCH_KEYS = ("lambda_z", "lambda_theta")
STRETCH_LABELS = (r"$\lambda_z$", r"$\lambda_\theta$")
DISP_KEYS = ("u_z", "u_theta")
DISP_LABELS = (r"$u_z$", r"$u_\theta$")
F_KEYS = (("F_zz", "F_ztheta"), ("F_thetaz", "F_thetatheta"))   # F[i][j] = d x_i / d X_j
FORCE_KEYS = ("force_z_mN", "force_theta_mN")


# ----------------------------------------------------------------------------
# Terminal progress
# ----------------------------------------------------------------------------
def _hms(s):
    s = int(round(s))
    return f"{s // 3600}h{s // 60 % 60:02d}m" if s >= 3600 else f"{s // 60}m{s % 60:02d}s"


class Progress:
    """
    Progress of a loop over `total` items: count, percent, elapsed time and ETA.

    In a terminal the line updates in place. When output is redirected (nohup, `> log.txt`,
    `| tee`) it prints a full line every `every` seconds instead, so log files stay readable.
    Set DIC_PROGRESS=0 to silence it, or DIC_PROGRESS=lines to force full lines; DIC_PROGRESS_EVERY sets
    the interval of full lines in seconds.

        bar = Progress(len(paths), "Def13")
        for k in range(len(paths)):
            ...
            bar.update(k + 1, f"valid {v:.1f}%")
        bar.close()
    """

    def __init__(self, total, label="", every=None, stream=None):
        if every is None:
            every = float(os.environ.get("DIC_PROGRESS_EVERY", 30))
        self.total, self.label, self.every = max(int(total), 1), label, every
        self.stream = stream or sys.stdout
        mode = os.environ.get("DIC_PROGRESS", "").lower()
        self.enabled = mode not in ("0", "off", "false", "no")
        self.inplace = self.enabled and mode != "lines" and self.stream.isatty()
        self.t0 = self.t_last = time.time()
        self.done, self.extra = 0, ""

    def _line(self):
        el = time.time() - self.t0
        eta = el / self.done * (self.total - self.done) if self.done else float("nan")
        n = len(str(self.total))
        s = (f"  {self.label}  {self.done:{n}d}/{self.total}  {100 * self.done / self.total:3.0f}%  "
             f"elapsed {_hms(el)}" + (f"  ETA {_hms(eta)}" if self.done < self.total and self.done else "")
             + (f"  {self.extra}" if self.extra else ""))
        return s

    def update(self, done, extra=""):
        self.done, self.extra = done, extra
        if not self.enabled:
            return
        now = time.time()
        if self.inplace:
            if now - self.t_last >= 0.25 or done >= self.total:
                width = shutil.get_terminal_size((100, 20)).columns - 1
                self.stream.write("\r" + self._line()[:width] + "\033[K")
                self.stream.flush()
                self.t_last = now
        elif now - self.t_last >= self.every:
            print(self._line(), file=self.stream, flush=True)
            self.t_last = now

    def close(self, extra=None):
        if extra is not None:
            self.extra = extra
        if not self.enabled:
            return
        if self.inplace:
            width = shutil.get_terminal_size((100, 20)).columns - 1
            self.stream.write("\r" + self._line()[:width] + "\033[K\n")
            self.stream.flush()
        else:
            print(self._line(), file=self.stream, flush=True)


# ----------------------------------------------------------------------------
# Image I/O
# ----------------------------------------------------------------------------
def frame_time(fname: str) -> float:
    m = FNAME_RE.search(fname)
    return int(m.group(1)) + int(m.group(2)) / 1000.0


def find_frames(dirs, t_min=-np.inf, t_max=np.inf, recursive=False):
    """Return (times, paths) for all CellScale jpgs in `dirs`, sorted by time."""
    out = {}
    for d in dirs:
        walk = os.walk(d) if recursive else [(d, None, os.listdir(d))]
        for root, _, names in walk:
            for f in names:
                if FNAME_RE.search(f):
                    t = frame_time(f)
                    if t_min <= t <= t_max:
                        out[t] = os.path.join(root, f)
    ts = np.array(sorted(out))
    return ts, [out[t] for t in ts]


def load_scale(path=None):
    """
    Image scale from calibrate_scale.py: dict with px_per_mm, mm_per_px, ...
    or None if the calibration has not been run. Default: the selected test's
    scale_calibration.json (testcfg).
    """
    import json
    if path is None:
        import testcfg
        path = testcfg.get().scale_file
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def load_gray(path, blur_sigma=0.8, resample=1.0):
    """Grayscale float32, resampled by `resample` (area averaging) if not 1. A light Gaussian
    prefilter reduces interpolation bias."""
    img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise IOError(path)
    img = img.astype(np.float32)
    if resample != 1.0:
        img = cv2.resize(img, None, fx=resample, fy=resample, interpolation=cv2.INTER_AREA)
    if blur_sigma > 0:
        img = cv2.GaussianBlur(img, (0, 0), blur_sigma)
    return img


def usable_mask(path, level=240, dilate=2, resample=1.0):
    """
    float32 image, 1 where a pixel can be used for matching and 0 on saturated glints:
    raw gray >= level, grown by `dilate` px to cover the blur halo. Specular glints (bubbles,
    wet highlights) do not move with the tissue and change with load, so subsets that
    include them fail to match.
    """
    raw = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    if raw is None:
        raise IOError(path)
    sat = (raw >= level).astype(np.uint8)
    if dilate > 0:
        sat = cv2.dilate(sat, np.ones((2 * dilate + 1, 2 * dilate + 1), np.uint8))
    use = (1 - sat).astype(np.float32)
    if resample != 1.0:          # detected at native resolution; a resampled pixel is usable only if fully clear
        use = (cv2.resize(use, None, fx=resample, fy=resample, interpolation=cv2.INTER_AREA) > 0.999)
        use = use.astype(np.float32)
    return use


def load_frame(path, prm, rects=None):
    """
    (image, usable mask or None) for TrackParams `prm`. `rects` are excluded rectangles
    (x0, y0, x1, y1) in this image's coordinates (prm.exclude_px moved with the tissue).
    """
    img = load_gray(path, prm.blur_sigma, prm.resample)
    use = usable_mask(path, prm.sat_level, prm.sat_dilate, prm.resample) if prm.mask_saturated else None
    if rects:
        use = np.ones(img.shape, np.float32) if use is None else use.copy()
        for x0, y0, x1, y1 in rects:
            use[max(0, int(np.floor(y0))):int(np.ceil(y1)) + 1, max(0, int(np.floor(x0))):int(np.ceil(x1)) + 1] = 0
    return img, use


def displacement_at(X0, x, valid, pts, k=8):
    """Displacement (n,2) at reference points `pts`, inverse-distance weighted from the k nearest valid points."""
    pts = np.atleast_2d(np.asarray(pts, float))
    tree = cKDTree(X0[valid])
    d, j = tree.query(pts, k=min(k, int(valid.sum())))
    w = 1.0 / np.maximum(d, 1e-6)
    u = (x[valid] - X0[valid])[j]
    return (u * w[..., None]).sum(1) / w.sum(1)[:, None]


def shift_rects(rects, X0, x, valid):
    """Move rectangles (x0, y0, x1, y1), given at reference positions X0, with the tissue (positions x)."""
    if not rects:
        return []
    c = np.array([[(r[0] + r[2]) / 2, (r[1] + r[3]) / 2] for r in rects])
    u = displacement_at(X0, x, valid, c)
    return [(r[0] + du[0], r[1] + du[1], r[2] + du[0], r[3] + du[1]) for r, du in zip(rects, u)]


# ----------------------------------------------------------------------------
# Grid
# ----------------------------------------------------------------------------
DEFAULT_STEP = 15          # grid spacing (px)
# Track files per grid: the full domain (roi_px.tine_region, the default) and the standard region
# (roi_px.standard = the full domain without its outer STANDARD_TRIM rows and columns).
GRID_SUFFIX = {"tine_region": "_tine_region", "standard": ""}
STANDARD_TRIM = 2
STRAIN_WINDOW_PX = 75      # side of the square strain-fit window (px)


def window_points(step, win_px=STRAIN_WINDOW_PX):
    """Strain-window half-width in grid points for a grid step (75 px -> 2 at 15 px)."""
    return max(1, int(round((win_px / float(step) - 1) / 2)))
# Loading protocols (same design in every test); ratios are axial : circumferential rake stretch.
PROTOCOLS = {11: "Equibiaxial", 12: "Axial reduced, 0.75 : 1", 13: "Axial reduced, 0.5 : 1",
             14: "Axial reduced, 0.25 : 1", 15: "Equibiaxial, repeat", 16: "Circumferential reduced, 1 : 0.75",
             17: "Circumferential reduced, 1 : 0.5", 18: "Circumferential reduced, 1 : 0.25",
             19: "Equibiaxial, final"}


DIVERGING_COLORS = ["#104281", "#2a78d6", "#9ec5f4", "#ffffff", "#f3b3ad", "#e34948", "#8f1d1d"]


def diverging_cmap():
    """Blue - white - red colormap for signed quantities (white at the center value)."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("blue_white_red", DIVERGING_COLORS)


def specimen_mm(x, y, tips, mm_per_px):
    """Image px -> specimen mm: z from the left tine tips, theta from the bottom tine tips (upward)."""
    return (np.asarray(x) - tips["z"][0]) * mm_per_px, (tips["theta"][1] - np.asarray(y)) * mm_per_px


def reference_time(R, times):
    """Time (s, test clock) of the track's reference frame: the common reference if recorded."""
    m = re.search(r"t=([\d.]+)s", str(R["reference"])) if "reference" in R else None
    return float(m.group(1)) if m else float(times[0])


def rake_stretch(data, times, t_ref):
    """Rake stretches (axial X, circumferential Y) at `times` relative to t_ref, from a *Data.csv
    table indexed by Time_S."""
    idx = data.index.to_numpy(float)
    ref = data.iloc[int(np.argmin(np.abs(idx - t_ref)))]
    d = data.reindex(np.round(times, 1))
    return (d.X_Size_um / ref.X_Size_um).to_numpy(float), (d.Y_Size_um / ref.Y_Size_um).to_numpy(float)


def settings_label(R, mm_per_px, window_label):
    """'1.0 mm subset, 1.5 mm strain window' for figure titles."""
    if "subset_mm" in R:
        subset = f"{float(R['subset_mm']):.1f} mm subset"
    else:
        n = 2 * int(R["radius"]) + 1
        subset = f"{n} px subset ({n * mm_per_px:.2f} mm)"
    return f"{subset}, {window_label.replace(' Gaussian', '')} strain window"


LOADING_SLACK_S = 2.0   # s past the logged Recover start: the logged boundary can lag the rake motion


def loading_end_time(test_root, folder, slack=LOADING_SLACK_S):
    """Time (s) up to which deformation `folder` (DefNN_...) is loading: the logged Recover start of its
    cycle plus `slack`. DefNN is the NN-th Stretch phase of Logs/*_PhaseTiming.csv (sort_images.py); the
    rake peak that ends the loading branch for the fits lies in this window (compare_stretch_inputs.loading_end)."""
    import csv
    import glob
    m = re.match(r"(?i)def(\d+)", os.path.basename(str(folder)))
    logs = sorted(glob.glob(os.path.join(test_root, "Logs", "*PhaseTiming.csv")))
    if not m or not logs:
        raise ValueError(f"no PhaseTiming log or DefNN folder name for {folder}")
    rows = list(csv.DictReader(open(logs[0])))
    stretch = [r for r in rows if r["Phase"] == "Stretch"]
    st = stretch[int(m.group(1)) - 1]
    rec = [r for r in rows if r["Phase"] == "Recover" and r["Set"] == st["Set"] and r["Cycle"] == st["Cycle"]]
    if not rec:
        raise ValueError(f"no Recover phase for {folder} (Set {st['Set']}, cycle {st['Cycle']})")
    return float(rec[0]["StartTime(s)"]) + slack


def programmed_stretch(test_root, folder):
    """Programmed rake stretch (axial %, circumferential %) of a deformation folder, read from the
    protocol table of the test's .tst file (X = axial, Y = circumferential; NumReps expands into
    DefNN folders in order). None if the file or row is missing or not in %."""
    import glob
    tst = sorted(glob.glob(os.path.join(test_root, "*.tst")))
    m = re.match(r"(?i)def(\d+)", os.path.basename(str(folder)))
    if not tst or not m:
        return None
    lines = open(tst[0], encoding="utf-8", errors="replace").read().splitlines()
    head = next((i for i, l in enumerate(lines) if l.startswith("Name\tXMode")), None)
    if head is None:
        return None
    cols, rows = lines[head].split("\t"), []
    for line in lines[head + 1:]:
        f = line.split("\t")
        if len(f) < 18:
            break
        r = dict(zip(cols, f))
        rows += [r] * int(float(r.get("NumReps") or 1))
    n = int(m.group(1))
    if n > len(rows) or rows[n - 1].get("XUnits") != "%" or rows[n - 1].get("YUnits") != "%":
        return None
    return float(rows[n - 1]["XMagnitude"]), float(rows[n - 1]["YMagnitude"])


def protocol_title(folder, test_root):
    """'Equibiaxial (Def11), programmed 38% axial / 20% circ.' for figure titles."""
    p = programmed_stretch(test_root, folder)
    return protocol_label(folder) + (f", programmed {p[0]:g}% axial / {p[1]:g}% circ." if p else "")


def protocol_label(folder):
    """Readable protocol name for a deformation folder, e.g. 'Equibiaxial (Def11)'."""
    m = re.match(r"(?i)def(\d+)", os.path.basename(str(folder)))
    if not m:
        return str(folder)
    n = int(m.group(1))
    name = PROTOCOLS.get(n, f"Preconditioning cycle {n}" if n <= 10 else f"Def{n:02d}")
    return f"{name} (Def{n:02d})"


def to_resampled(xy, f):
    """Native image px -> px of the image resampled by f (cv2.resize pixel-center convention)."""
    return (np.asarray(xy, np.float64) + 0.5) * f - 0.5


def to_native(xy, f):
    """Px of the image resampled by f -> native image px."""
    return (np.asarray(xy, np.float64) + 0.5) / f - 0.5


def resample_factor(px_per_mm, s=None):
    """Factor the test's images are resampled by before tracking: target_px_per_mm of
    analysis_settings.json over the native scale, never above 1 (images are only downsampled)."""
    import testcfg
    t = (s or testcfg.analysis()).get("target_px_per_mm")
    return min(1.0, t / px_per_mm) if t else 1.0


def settings_px(px_per_mm, s=None):
    """analysis_settings.json (physical sizes, mm) in px at `px_per_mm` (the tracking scale):
    subset radius, grid step, plus the mm values themselves. Subset = 2*radius + 1 px."""
    import testcfg
    s = dict(s or testcfg.analysis())
    # rounded to 1e-6 first, so a scale of 49.99999... px/mm (from resampling) converts like 50
    s["radius"] = int(np.floor(round((s["subset_mm"] * px_per_mm - 1) / 2 + 0.5, 6)))
    s["step"] = int(np.floor(round(s["step_mm"] * px_per_mm + 0.5, 6)))
    return s


def tip_lattice_roi(tips, step, margin=5):
    """Full-domain grid extents [x0, y0, x1, y1] (px): a lattice of `step` centered between the
    rest-state tine-tip lines (tine_tips_px), at least `margin` px inside them."""
    def axis(lo, hi):
        a, b = lo + margin, hi - margin
        n = int((b - a) // step)
        x0 = int(round(a + ((b - a) - n * step) / 2))
        return x0, x0 + n * step
    (x0, x1), (y0, y1) = axis(*tips["z"]), axis(*tips["theta"])
    return [x0, y0, x1, y1]


def interior_mask(X0, shape, tips, mm_per_px, edge_mm):
    """Grid points far enough from the rest-state tine-tip lines (tine_tips_px): the interior used for
    J/mean and spread statistics. edge_mm is one distance for all edges, or {"z": ..., "theta": ...}:
    z = distance from the left/right (axial-pulling) tine lines, theta = from the top/bottom ones.
    Defined in mm, so it is the same region on any grid."""
    ez, et = (edge_mm["z"], edge_mm["theta"]) if isinstance(edge_mm, dict) else (edge_mm, edge_mm)
    ez, et = ez / mm_per_px, et / mm_per_px
    x, y = X0[:, 0].reshape(shape), X0[:, 1].reshape(shape)
    return ((x >= tips["z"][0] + ez) & (x <= tips["z"][1] - ez)
            & (y >= tips["theta"][0] + et) & (y <= tips["theta"][1] - et))


def make_grid(roi, step):
    """roi = (x0, y0, x1, y1) in px. Returns (points (N,2), shape (ny, nx))."""
    x0, y0, x1, y1 = roi
    xs = np.arange(x0, x1 + 1e-9, step, dtype=np.float64)
    ys = np.arange(y0, y1 + 1e-9, step, dtype=np.float64)
    X, Y = np.meshgrid(xs, ys)
    return np.stack([X.ravel(), Y.ravel()], 1), X.shape


def standard_from_full(R, trim=STANDARD_TRIM):
    """
    Standard-region track sliced from a full-domain track: rows and columns trim:-trim of the
    grid. Per-point arrays are cut; everything else is kept. A separately tracked standard grid
    agrees with this to about 0.03 px (same reference, settings and lattice).
    """
    R = dict(R)
    ny, nx = (int(v) for v in R["shape"])
    keep = np.zeros((ny, nx), bool)
    keep[trim:ny - trim, trim:nx - trim] = True
    keep = keep.ravel()
    out = {}
    for k, v in R.items():
        v = np.asarray(v)
        if v.ndim >= 1 and v.shape[0] == ny * nx and k in ("X0", "X_start", "bridge_zncc"):
            out[k] = v[keep]
        elif v.ndim >= 2 and v.shape[1] == ny * nx:
            out[k] = v[:, keep]
        else:
            out[k] = v
    out["shape"] = np.array([ny - 2 * trim, nx - 2 * trim])
    step = float(R["step"])
    out["roi"] = np.asarray(R["roi"], float) + trim * step * np.array([1, 1, -1, -1])
    out["derived_from"] = f"full-domain track, rows/cols {trim}:-{trim}"
    return out


def load_track(images, folder, grid="tine_region"):
    """
    (track dict, path) for deformation `folder` on `grid` ("tine_region" or "standard").
    A standard track that was not tracked separately is sliced from the full-domain one.
    """
    base = os.path.join(images, folder, "strain_results", folder)
    path = base + GRID_SUFFIX[grid] + ".npz"
    if os.path.exists(path):
        return dict(np.load(path, allow_pickle=True)), path
    full = base + GRID_SUFFIX["tine_region"] + ".npz"
    if grid == "standard" and os.path.exists(full):
        return standard_from_full(np.load(full, allow_pickle=True)), full
    raise FileNotFoundError(path)


# ----------------------------------------------------------------------------
# IC-GN subset matching
# ----------------------------------------------------------------------------
def _affine(p):
    """p = (u, ux, uy, v, vx, vy) -> (N,3,3) warp matrices acting on local coords."""
    n = p.shape[0]
    A = np.zeros((n, 3, 3))
    A[:, 0, 0] = 1 + p[:, 1]
    A[:, 0, 1] = p[:, 2]
    A[:, 0, 2] = p[:, 0]
    A[:, 1, 0] = p[:, 4]
    A[:, 1, 1] = 1 + p[:, 5]
    A[:, 1, 2] = p[:, 3]
    A[:, 2, 2] = 1
    return A


def _params(A):
    return np.stack([A[:, 0, 2], A[:, 0, 0] - 1, A[:, 0, 1],
                     A[:, 1, 2], A[:, 1, 0], A[:, 1, 1] - 1], 1)


def _sample(img, x, y):
    """Bicubic sampling of img at float arrays x, y (any 2D shape <= 32767); NaN where the 4 x 4
    bicubic stencil leaves the image. (A NaN borderValue in cv2.remap leaks NaN into some interior
    samples at certain map widths, e.g. 63 or 65 px subsets, so the border is masked explicitly.)"""
    xf, yf = x.astype(np.float32), y.astype(np.float32)
    h, w = img.shape[:2]
    out = cv2.remap(img, xf, yf, cv2.INTER_CUBIC, borderMode=cv2.BORDER_REPLICATE).astype(np.float32)
    with np.errstate(invalid="ignore"):
        out[~((xf >= 1) & (yf >= 1) & (xf <= w - 3) & (yf <= h - 3))] = np.nan
    return out


class Reference:
    """
    Precomputed IC-GN quantities for every subset. Each subset has its own
    reference (keyframe image + center); update() re-references a chosen subset of
    points to a new image, so points can sit on different keyframes.
    """

    def __init__(self, centers, radius):
        self.centers = np.asarray(centers, np.float64).copy()   # (N,2)
        self.radius = radius
        self.w = None          # per-pixel weights (N,S) when glints are masked, else None
        d = np.arange(-radius, radius + 1, dtype=np.float64)
        DX, DY = np.meshgrid(d, d)
        self.dx, self.dy = DX.ravel(), DY.ravel()                 # (S,)
        self.dx32, self.dy32 = self.dx.astype(np.float32), self.dy.astype(np.float32)
        self.B3 = np.stack([np.ones_like(self.dx), self.dx, self.dy], 1)   # (S,3): [1, dx, dy]
        self.corners = np.array([[-radius, radius, -radius, radius],
                                 [-radius, -radius, radius, radius]], np.float64)

    def _subset_data(self, img, centers, use=None):
        cx, cy = centers[:, :1], centers[:, 1:]
        X, Y = cx + self.dx, cy + self.dy                         # (n,S)
        gy, gx = np.gradient(img)
        f = _sample(img, X, Y).astype(np.float64)
        fx = _sample(gx.astype(np.float32), X, Y).astype(np.float64)
        fy = _sample(gy.astype(np.float32), X, Y).astype(np.float64)
        ok = np.isfinite(f).all(1)
        f = np.nan_to_num(f)
        fx, fy = np.nan_to_num(fx), np.nan_to_num(fy)
        w = None
        if use is not None:
            # glint mask: weight 0 on saturated pixels; statistics over the usable ones
            w = (cv2.remap(use, X.astype(np.float32), Y.astype(np.float32), cv2.INTER_NEAREST,
                           borderMode=cv2.BORDER_CONSTANT, borderValue=0) > 0.5).astype(np.float64)
            ok &= w.mean(1) >= 0.5
            fx, fy = fx * w, fy * w
            f_raw = f
            f = (f - (w * f).sum(1, keepdims=True) / np.maximum(w.sum(1, keepdims=True), 1)) * w
        else:
            f_raw = None
            f = f - f.mean(1, keepdims=True)
        df = np.sqrt((f ** 2).sum(1)) + 1e-12
        J = np.stack([fx, fx * self.dx, fx * self.dy,
                      fy, fy * self.dx, fy * self.dy], 2)
        H = np.einsum("nsi,nsj->nij", J, J)
        # Weak/no texture -> ill-conditioned H; flag via condition number.
        ok &= np.linalg.cond(H) < 1e8
        H[~ok] = np.eye(6)
        return f, df, J, np.linalg.inv(H), ok, w, f_raw

    def build(self, img, use=None):
        """Reference every subset to `img` at self.centers (`use`: usable-pixel mask, or None)."""
        (self.f, self.df, self.J, self.Hinv, self.ok,
         self.w, self.f_raw) = self._subset_data(img, self.centers, use)
        # contiguous copies of the gradients for icgn's steepest-descent sums
        self.fx, self.fy = np.ascontiguousarray(self.J[:, :, 0]), np.ascontiguousarray(self.J[:, :, 3])
        return self

    def update(self, img, idx, centers, use=None):
        """Re-reference subsets `idx` to `img`, centered at `centers` (len(idx),2)."""
        if len(idx) == 0:
            return self
        f, df, J, Hinv, ok, w, f_raw = self._subset_data(img, np.asarray(centers, np.float64), use)
        self.centers[idx] = centers
        self.f[idx], self.df[idx], self.J[idx], self.Hinv[idx], self.ok[idx] = f, df, J, Hinv, ok
        self.fx[idx], self.fy[idx] = J[:, :, 0], J[:, :, 3]
        if self.w is not None:
            self.w[idx], self.f_raw[idx] = w, f_raw
        return self


def _weighted_pair(ref, idx, g, wc):
    """
    Masked ZNSSD terms for subsets idx: the usable pixels are those usable in both the
    reference and the current image. Returns (f~, |f~|, g~, |g~|), zero on unusable pixels.
    """
    w = ref.w[idx] if wc is None else ref.w[idx] * wc
    n = np.maximum(w.sum(1, keepdims=True), 1)
    fr = ref.f_raw[idx]
    fc = (fr - (w * fr).sum(1, keepdims=True) / n) * w
    gc = (g - (w * g).sum(1, keepdims=True) / n) * w
    return fc, np.sqrt((fc ** 2).sum(1)) + 1e-12, gc, np.sqrt((gc ** 2).sum(1)) + 1e-12


def _usable_at(use, xw, yw):
    if use is None:
        return None
    return (cv2.remap(use, xw.astype(np.float32), yw.astype(np.float32), cv2.INTER_NEAREST,
                      borderMode=cv2.BORDER_CONSTANT, borderValue=0) > 0.5).astype(np.float64)


def _warp_coords(ref, idx, p):
    """Warped pixel coordinates of subsets idx under affine parameters p, in float32 (the precision
    cv2.remap samples at): each subset's center offset is rounded once, then the local terms added."""
    cx = (ref.centers[idx, 0] + p[:, 0]).astype(np.float32)[:, None]
    cy = (ref.centers[idx, 1] + p[:, 3]).astype(np.float32)[:, None]
    xw = cx + (1 + p[:, 1:2]).astype(np.float32) * ref.dx32 + p[:, 2:3].astype(np.float32) * ref.dy32
    yw = cy + p[:, 4:5].astype(np.float32) * ref.dx32 + (1 + p[:, 5:6]).astype(np.float32) * ref.dy32
    return xw, yw


def _leaves_image(ref, idx, p, shape):
    """True for subsets whose bicubic stencil leaves the image; the warp is affine, so its 4 corners decide."""
    h, w = shape[:2]
    cdx, cdy = ref.corners
    x = ref.centers[idx, :1] + p[:, :1] + (1 + p[:, 1:2]) * cdx + p[:, 2:3] * cdy
    y = ref.centers[idx, 1:] + p[:, 3:4] + p[:, 4:5] * cdx + (1 + p[:, 5:6]) * cdy
    with np.errstate(invalid="ignore"):
        return ~((x >= 1) & (y >= 1) & (x <= w - 3) & (y <= h - 3)).all(1)


def icgn(ref: Reference, img, p0, max_iter=100, tol=1e-3, use=None):
    """
    Match all subsets of `ref` into `img` starting from p0 (N,6).
    Returns p (N,6), zncc (N,), converged (N,).
    With a glint-masked reference (ref.w), pixels that are saturated in the reference or
    (via `use`, the current image's usable mask) in the current image are left out.
    """
    p = p0.copy()
    n = p.shape[0]
    active = ref.ok.copy()
    conv = np.zeros(n, bool)
    r = ref.radius
    for _ in range(max_iter):
        idx = np.nonzero(active)[0]
        if idx.size == 0:
            break
        pa = p[idx]
        xw, yw = _warp_coords(ref, idx, pa)
        bad = _leaves_image(ref, idx, pa, img.shape)
        g = cv2.remap(img, xw, yw, cv2.INTER_CUBIC, borderMode=cv2.BORDER_REPLICATE).astype(np.float64)
        g[bad] = 0
        if ref.w is None:
            g -= g.mean(1, keepdims=True)
            dg = np.sqrt(np.einsum("ns,ns->n", g, g)) + 1e-12
            e = ref.f[idx] - (ref.df[idx] / dg)[:, None] * g
        elif use is None:
            # usable pixels are the reference's, so ref.f and ref.df already hold the masked, centered terms
            w = ref.w[idx]
            g -= (np.einsum("ns,ns->n", w, g) / np.maximum(w.sum(1), 1))[:, None]
            g *= w
            dg = np.sqrt(np.einsum("ns,ns->n", g, g)) + 1e-12
            e = ref.f[idx] - (ref.df[idx] / dg)[:, None] * g
        else:
            fc, df_, g, dg = _weighted_pair(ref, idx, g, _usable_at(use, xw, yw))
            e = fc - (df_ / dg)[:, None] * g
        # b = J^T e with J = [fx, fx dx, fx dy, fy, fy dx, fy dy], as two (n,S) x (S,3) products
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            # (numpy 2.0 with Apple Accelerate flags these on finite inputs; a NaN would still propagate to p)
            b = np.concatenate([(ref.fx[idx] * e) @ ref.B3, (ref.fy[idx] * e) @ ref.B3], 1)
        dp = -np.einsum("nij,nj->ni", ref.Hinv[idx], b)
        A = _affine(pa) @ np.linalg.inv(_affine(dp))
        p[idx] = _params(A)
        step = np.sqrt(dp[:, 0] ** 2 + dp[:, 3] ** 2
                       + ((r * dp[:, [1, 2, 4, 5]]) ** 2).sum(1))
        done = (step < tol) | bad
        conv[idx[done & ~bad]] = True
        active[idx[done]] = False
        active[idx[bad]] = False
        p[idx[bad]] = np.nan
    # Final correlation coefficient at the converged warp
    cx, cy = ref.centers[:, :1], ref.centers[:, 1:]
    pp = np.nan_to_num(p)
    xw = cx + pp[:, :1] + (1 + pp[:, 1:2]) * ref.dx + pp[:, 2:3] * ref.dy
    yw = cy + pp[:, 3:4] + pp[:, 4:5] * ref.dx + (1 + pp[:, 5:6]) * ref.dy
    g = _sample(img, xw, yw).astype(np.float64)
    if ref.w is None:
        g = g - np.nanmean(g, 1, keepdims=True)
        dg = np.sqrt(np.nansum(g ** 2, 1)) + 1e-12
        zncc = np.nansum(ref.f * g, 1) / (ref.df * dg)
    else:
        finite = np.isfinite(g)
        wc = _usable_at(use, xw, yw)
        wc = finite.astype(np.float64) if wc is None else wc * finite
        fc, df_, gc, dg = _weighted_pair(ref, np.arange(n), np.nan_to_num(g), wc)
        zncc = (fc * gc).sum(1) / (df_ * dg)
    zncc[~ref.ok | ~np.isfinite(p).all(1)] = np.nan
    return p, zncc, conv


def _fill_from_neighbors(pts, p, good, k=8):
    """Replace p[~good] with inverse-distance-weighted mean of k good neighbors."""
    if good.all() or good.sum() < 3:
        return p
    tree = cKDTree(pts[good])
    d, j = tree.query(pts[~good], k=min(k, good.sum()))
    w = 1.0 / np.maximum(d, 1e-6)
    pg = p[good]
    p = p.copy()
    p[~good] = (pg[j] * w[..., None]).sum(1) / w.sum(1)[:, None]
    return p


# ----------------------------------------------------------------------------
# Sequence tracker
# ----------------------------------------------------------------------------
@dataclass
class TrackParams:
    radius: int = 30              # subset half-width (px); subset = 2r+1 = 61 px
    zncc_min: float = 0.75        # below -> point invalid in that frame
    median_thr: float = 3.0       # normalized median test threshold (frame increments)
    median_eps: float = 0.1       # px; noise floor in the median test
    update_strain: float = 0.04   # max |displacement gradient| before new keyframe
    update_zncc: float = 0.93     # median ZNCC below -> new keyframe
    max_gap: int = 60             # frames between forced keyframe updates
    max_iter: int = 100           # IC-GN iterations
    blur_sigma: float = 0.8
    mask_saturated: bool = False  # leave saturated glint pixels out of the matching
    sat_level: float = 240        # raw gray level counted as saturated
    sat_dilate: int = 2           # px the saturated area is grown by
    sat_current: bool = False     # also mask saturated pixels of the current image (else reference only)
    exclude_px: tuple = ()        # rectangles (x0, y0, x1, y1) in the first frame's px, left out of the matching
                                  # (e.g. a speckle particle that slides over the tissue); they move with the tissue
    resample: float = 1.0         # images are resampled by this factor (area averaging) before matching; all px
                                  # quantities above are then in resampled px


def median_test(du, nb, good, eps, thr):
    """
    Normalized median test (Westerweel & Scarano 2005) on displacement vectors du
    (N,2) using neighbor indices nb (N,K). Only `good` neighbors are used; points
    with fewer than 3 good neighbors pass. Returns a boolean pass mask.
    """
    G = good[nb]
    U = np.where(G[..., None], du[nb], np.nan)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(U, 1)
        res = np.nanmedian(np.abs(U - med[:, None, :]), 1)
        r = np.abs(du - med) / (res + eps)
    return (G.sum(1) < 3) | (np.nanmax(r, 1) <= thr)


def track_sequence(paths, X0, prm: TrackParams, verbose=True, times=None, label="frames"):
    """
    Track material points X0 (N,2) through image `paths`.
    Returns dict with pos (T,N,2), zncc (T,N), valid (T,N), keyframes list.
    verbose shows a Progress line labeled `label`.

    A point is valid in a frame when IC-GN converges, ZNCC >= zncc_min and its
    frame-to-frame displacement passes a normalized median test against its
    neighbors. Failed points get a position interpolated from valid neighbors
    (flagged invalid) and keep their last good reference, so they rejoin as real
    measurements as soon as they match again. At a keyframe update only currently
    valid points are re-referenced to the new image.
    """
    T, N = len(paths), len(X0)
    pos = np.full((T, N, 2), np.nan)
    zn = np.full((T, N), np.nan)
    valid = np.zeros((T, N), bool)
    # why a point was rejected, per frame (bitmask): 1 IC-GN not converged,
    # 2 ZNCC < zncc_min, 4 failed median test, 8 unusable reference subset
    fail = np.zeros((T, N), np.uint8)
    rects0 = [tuple(r) for r in prm.exclude_px]
    img, use = load_frame(paths[0], prm, rects0)
    pos[0], zn[0] = X0, 1.0
    ref = Reference(X0, prm.radius).build(img, use)
    valid[0] = ref.ok
    nb = cKDTree(X0).query(X0, k=9)[1][:, 1:]          # 8 nearest grid neighbors
    k_ref, keyframes = 0, [0]
    p = np.zeros((N, 6))
    p_prev = p.copy()
    bar = Progress(T, label) if verbose else None
    for k in range(1, T):
        img, use = load_frame(paths[k], prm)
        # constant-velocity predictor, scaled by the actual time step (frame gaps)
        r_dt = 1.0
        if times is not None and k >= 2:
            r_dt = (times[k] - times[k - 1]) / max(times[k - 1] - times[k - 2], 1e-9)
        p_init = p + (p - p_prev) * r_dt
        p_new, c, conv = icgn(ref, img, p_init, max_iter=prm.max_iter, use=use if prm.sat_current else None)
        good = conv & (c >= prm.zncc_min) & ref.ok & np.isfinite(p_new).all(1)
        fail[k] = (1 * ~conv) | (2 * ~(c >= prm.zncc_min)) | (8 * ~ref.ok)
        p_new = np.where(np.isfinite(p_new), p_new, p)
        x = ref.centers + p_new[:, [0, 3]]
        mt = median_test(x - pos[k - 1], nb, good, prm.median_eps, prm.median_thr)
        fail[k] |= (4 * (good & ~mt)).astype(np.uint8)
        good &= mt
        bad = ~good
        if bad.any() and good.sum() >= 3:
            # position from neighbors' current displacement; keep own reference
            u = _fill_from_neighbors(X0, x - X0, good)
            x[bad] = X0[bad] + u[bad]
            p_new[bad, 0] = x[bad, 0] - ref.centers[bad, 0]
            p_new[bad, 3] = x[bad, 1] - ref.centers[bad, 1]
            p_new[np.ix_(bad, [1, 2, 4, 5])] = p[np.ix_(bad, [1, 2, 4, 5])]
        pos[k], zn[k], valid[k] = x, c, good
        p_prev, p = p, p_new

        grad = np.abs(p_new[good][:, [1, 2, 4, 5]]).max(1) if good.any() else [0]
        need = (np.median(grad) > prm.update_strain
                or np.nanmedian(c) < prm.update_zncc
                or k - k_ref >= prm.max_gap)
        if need and k < T - 1:
            # New keyframe for currently valid points only, re-centered on their
            # current positions. Failing points keep their last good reference.
            idx = np.nonzero(good)[0]
            vel = pos[k, idx] - pos[k - 1, idx]
            if rects0:            # excluded rectangles follow the tissue into the new keyframe
                _, use = load_frame(paths[k], prm, shift_rects(rects0, X0, pos[k], good))
            ref.update(img, idx, pos[k, idx], use)
            p[idx] = 0
            p_prev[idx] = 0
            p_prev[idx, 0], p_prev[idx, 3] = -vel[:, 0], -vel[:, 1]
            k_ref = k
            keyframes.append(k)
        if bar:
            bar.update(k + 1, f"valid {good.mean()*100:5.1f}%  ZNCC {np.nanmedian(c):.4f}  "
                              f"keyframes {len(keyframes)}")
    if bar:
        bar.close()
    return dict(pos=pos, zncc=zn, valid=valid, fail=fail, keyframes=np.array(keyframes))


# ----------------------------------------------------------------------------
# Strain
# ----------------------------------------------------------------------------
def deformation_gradient_field(X0, x, valid, shape, win, min_frac=0.5, fill_holes=True):
    """
    Masked least-squares plane fit of displacement over a (2*win+1)^2 window of
    grid points. Returns F (ny,nx,2,2) and a mask of points with enough support.

    A point needs >= min_frac of its window valid. With fill_holes, a point that is
    itself invalid still gets F from its valid neighbors, provided they surround
    it (support centroid within one grid step), i.e. interpolation, never
    extrapolation. Use min_frac=0.6, fill_holes=False for the original behavior.
    """
    ny, nx = shape
    Xg = X0[:, 0].reshape(shape)
    Yg = X0[:, 1].reshape(shape)
    U = (x[:, 0] - X0[:, 0]).reshape(shape)
    V = (x[:, 1] - X0[:, 1]).reshape(shape)
    w = valid.reshape(shape).astype(np.float64)
    U, V = np.where(w > 0, U, 0), np.where(w > 0, V, 0)
    size = 2 * win + 1
    S = lambda a: ndimage.uniform_filter(a * w, size, mode="constant") * size ** 2
    # local coords relative to each window center
    n = S(np.ones(shape))
    sx, sy = S(Xg), S(Yg)
    sxx, syy, sxy = S(Xg * Xg), S(Yg * Yg), S(Xg * Yg)
    mx, my = sx / np.maximum(n, 1), sy / np.maximum(n, 1)
    cxx = sxx - n * mx * mx
    cyy = syy - n * my * my
    cxy = sxy - n * mx * my
    det = cxx * cyy - cxy ** 2
    F = np.zeros((ny, nx, 2, 2))
    for i, Q in enumerate((U, V)):
        sq = S(Q)
        sqx = S(Q * Xg) - mx * sq
        sqy = S(Q * Yg) - my * sq
        F[..., i, 0] = (cyy * sqx - cxy * sqy) / np.where(det == 0, np.nan, det)
        F[..., i, 1] = (cxx * sqy - cxy * sqx) / np.where(det == 0, np.nan, det)
    F[..., 0, 0] += 1
    F[..., 1, 1] += 1
    ok = n >= min_frac * size ** 2
    if fill_holes:
        step = max(abs(Xg[0, 1] - Xg[0, 0]), abs(Yg[1, 0] - Yg[0, 0])) if min(shape) > 1 else 1
        ok &= np.hypot(mx - Xg, my - Yg) <= step
    else:
        ok &= valid.reshape(shape)
    F[~ok] = np.nan
    return F, ok


def strain_window_F(X0, x, valid, shape, step, mm_per_px, win=None):
    """F field with the project's strain window: the Gaussian of window_mm (analysis_settings.json),
    or, with `win` (half-width in grid points), the square box window. Returns (F, ok, label)."""
    if win is not None:
        F, ok = deformation_gradient_field(X0, x, valid, shape, win)
        return F, ok, f"{(2 * win + 1) * step:g} px box"
    import testcfg
    w = testcfg.analysis()["window_mm"]
    F, ok = deformation_gradient_field_mm(X0, x, valid, shape, w, mm_per_px)
    return F, ok, f"{w:g} mm Gaussian"


def deformation_gradient_field_mm(X0, x, valid, shape, window_mm, mm_per_px, min_frac=0.5, fill_holes=True,
                                  truncate=2.5):
    """
    Like deformation_gradient_field, but the plane fit is Gaussian-weighted with a width set in mm,
    so specimens imaged at different zoom can be compared at the same physical averaging scale.
    `window_mm` is the equivalent box width L: the Gaussian has the same spread, sigma = L / sqrt(12)
    (per axis), so --window-mm 1.19 matches the 75 px box at 62.8 px/mm. Support rule: at least
    min_frac of the full kernel weight must be valid and inside the grid.
    """
    ny, nx = shape
    Xg = X0[:, 0].reshape(shape)
    Yg = X0[:, 1].reshape(shape)
    step = max(abs(Xg[0, 1] - Xg[0, 0]), abs(Yg[1, 0] - Yg[0, 0])) if min(shape) > 1 else 1
    sigma = window_mm / np.sqrt(12) / (step * mm_per_px)            # in grid points
    U = (x[:, 0] - X0[:, 0]).reshape(shape)
    V = (x[:, 1] - X0[:, 1]).reshape(shape)
    w = valid.reshape(shape).astype(np.float64)
    U, V = np.where(w > 0, U, 0), np.where(w > 0, V, 0)
    S = lambda a: ndimage.gaussian_filter(a * w, sigma, mode="constant", truncate=truncate)   # noqa: E731
    n = S(np.ones(shape))                                              # valid fraction of the kernel weight
    sx, sy = S(Xg), S(Yg)
    sxx, syy, sxy = S(Xg * Xg), S(Yg * Yg), S(Xg * Yg)
    mx, my = sx / np.maximum(n, 1e-12), sy / np.maximum(n, 1e-12)
    cxx = sxx - n * mx * mx
    cyy = syy - n * my * my
    cxy = sxy - n * mx * my
    det = cxx * cyy - cxy ** 2
    F = np.zeros((ny, nx, 2, 2))
    for i, Q in enumerate((U, V)):
        sq = S(Q)
        sqx = S(Q * Xg) - mx * sq
        sqy = S(Q * Yg) - my * sq
        F[..., i, 0] = (cyy * sqx - cxy * sqy) / np.where(det <= 0, np.nan, det)
        F[..., i, 1] = (cxx * sqy - cxy * sqx) / np.where(det <= 0, np.nan, det)
    F[..., 0, 0] += 1
    F[..., 1, 1] += 1
    ok = n >= min_frac
    if fill_holes:
        ok &= np.hypot(mx - Xg, my - Yg) <= step
    else:
        ok &= valid.reshape(shape)
    F[~ok] = np.nan
    return F, ok


def to_specimen(F):
    """Image-frame F (..., 2, 2) -> specimen frame (z, theta): S F S with S = diag(1, THETA_SIGN)."""
    S = np.diag([1.0, float(THETA_SIGN)])
    return S @ np.asarray(F) @ S


def green_lagrange(F):
    C = np.einsum("...ki,...kj->...ij", F, F)
    return 0.5 * (C - np.eye(2))


def affine_fit(X, x):
    """Least-squares homogeneous F mapping X -> x (both (N,2))."""
    Xc, xc = X - X.mean(0), x - x.mean(0)
    sol, *_ = np.linalg.lstsq(Xc, xc, rcond=None)
    return sol.T


def bilinear_F_center(P0, P):
    """
    F at the centroid of a 4-node bilinear element (Sacks 2000 style).
    P0, P: (4,2) nodes ordered counter-clockwise in the reference config.
    """
    xi = np.array([-1, 1, 1, -1.0])
    eta = np.array([-1, -1, 1, 1.0])
    dN = 0.25 * np.stack([xi, eta], 1)          # dN/d(xi,eta) at center
    Jx0 = P0.T @ dN                              # dX/dxi
    Jx = P.T @ dN                                # dx/dxi
    return Jx @ np.linalg.inv(Jx0)
