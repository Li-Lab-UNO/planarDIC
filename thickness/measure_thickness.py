"""
Measure the wall thickness of FPA strips photographed on a measuring mat or ruler.

Each strip lies on its cut face, so its wall thickness is the width of the
tissue band in the photo. For each image the script:

  1. Scale. Detects the background and reads the image scale from it.
     * grid:  blue cutting mat, 0.5 in grid lines. The two line families are found
       by projection (Radon), the spacing by a harmonic comb on the autocorrelation,
       then refined from the individual line positions in a window around the strip
       (so perspective across the frame does not bias it).
     * ruler: gray steel ruler. The tick period is read along both long edges; the mm
       edge is told from the inch edge (1/16 or 1/32 in) by the ratio of the periods.
  2. Strip. Tissue is warm (high Lab b*) against the blue mat and the white board, and brighter than
     the dark steel of the ruler (there brightness alone decides: wet films on the steel look warm).
     The mask is opened with a disk of OPEN_MM radius, which strips off loose fibers,
     Petri-dish rim highlights and glare. Strips cut by the image border are skipped.
  3. Thickness. The centerline is the longest path through the mask skeleton,
     smoothed. At stations every STEP_MM along it (ends trimmed) the chord along the
     local normal is measured, sub-pixel, on the mask. Reports mean ± SD, median, IQR.

    .venv/bin/python thickness/measure_thickness.py "N1690_RFPA/Thickness, 1690/test 1"
    .venv/bin/python thickness/measure_thickness.py IMG.JPG --background ruler --ruler-tick-mm 1
    .venv/bin/python thickness/measure_thickness.py "N1690_RFPA/Thickness, 1690" --recursive

Writes the results of a photo folder <SAMPLE>/Thickness, <ID>/test <N> to <SAMPLE>/Thickness, <ID>/thickness_auto/test<N>/
(any other image folder: <image folder>/thickness_auto/): an overlay PNG per image, a per-station CSV per strip, and
summary.csv (one row per strip).

--write puts the result into the test's test_config.json as thickness_mm, for each photo folder
<SAMPLE>/Thickness, <ID>/test <N>: H and H_sd are the mean and SD of all stations of the selected
strips of the folder's photos whose scale is reliable (as compare_manual.py pools them); one such photo
is enough, and the output lists the photos used and every photo left out, with the reason (excluded,
failed, no strip found, scale uncertain). It replaces a workbook value (make_config.py), which is kept in
thickness_mm.workbook for comparison. With --keep-workbook, a workbook value stays in use (and one kept
in thickness_mm.workbook is put back); the photos then fill only tests without one. When no photo gives
a thickness, the config keeps what it has.

    .venv/bin/python thickness/measure_thickness.py --recursive --write "N1757_LFPA/Thickness, 1757"
    .venv/bin/python thickness/measure_thickness.py --recursive --write --keep-workbook "N1739_LFPA/Thickness, 1739"

`selected` marks the strips each photo is for: every complete strip (the strips in a test's photos are
pieces of its specimen), unless <image folder>/strip_selection.json names one by position
({"DSC00503.jpg": "lower"}; upper, lower, left, right or center; "exclude" leaves a photo out).
A strip much fainter against its background than the photo's strongest strip (STRIP_CONTRAST_FRAC) is
dropped as stained mat or a film, and the drop is printed.
"""
import argparse
import glob
import json
import os
import re
import sys

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402
from scipy.signal import find_peaks, savgol_filter  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import dijkstra  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from skimage.morphology import reconstruction, skeletonize  # noqa: E402
from skimage.segmentation import watershed  # noqa: E402

MM_PER_IN = 25.4
DS = 4              # downsampling for scale detection
WORK_PX_MM = 40.0   # strip images are resampled to this scale for segmentation
OPEN_MM = 0.25      # opening radius: removes fibers/rims thinner than ~2*OPEN_MM
CLOSE_MM = 0.4      # closing radius: joins the medial and adventitial layers across the
                    # darker seam between them (gaps up to ~2*CLOSE_MM)
MIN_AREA_MM2 = 3.0  # smallest tissue blob kept
RULER_TISSUE_PCT = 75   # on the ruler: tissue brightness = this percentile of the pixels brighter than the steel
RULER_EDGE_FRAC = 0.6   # ... and the tissue edge this fraction of the way from the steel's brightness to it
# On the mat the threshold mask (warm and bright) holds the confident tissue; each strip is then grown into the rest
# of its tissue (pale layers, a darker or shaded side, blurred edges), judged against the mat around it:
STRIP_CONTRAST_FRAC = 0.5   # a strip stands out from its local background at least this fraction as much as the photo's
                            # strongest strip (fainter ones are stained mat or films, dropped and reported)
SEED_WARM_MAX = 25   # on the mat, confident tissue is at most this much warmer (b*) than the mat need be
GROW_FRAC = 0.5      # brighter pixels: differ from the local mat by at least this fraction of the strip core's contrast
GROW_MM = 0.75       # growth reaches at most this far from the confident tissue
SEED_MIN_MM2 = 0.3   # confident pieces this large take part (smaller ones are specks)
BRIDGE_MM = 3.0      # pieces of one strip this close are joined when tissue-colored all along the gap ...
BRIDGE_FRAC = 0.5    # ... (every point within 0.3 mm of the path at least this fraction of the core contrast)
CORRIDOR_MM = 1.0    # ... by growing inside a corridor of this half-width
RING_MM = (1.5, 3.0) # the local mat: a ring this far around each piece, away from other tissue
LINE_MM = 2.5        # printed grid lines: warm straight runs at least this long along a grid direction
CHUNK_MM, CUT_CHUNKS = 1.0, 3   # growth may not join parts of a piece more than CUT_CHUNKS x CHUNK_MM apart along it
STEP_MM = 0.1       # station spacing along the centerline
TRIM_FRAC = 0.10    # fraction of centerline length dropped at each end
SMOOTH_MM = 1.0     # centerline smoothing window, at least ...
SMOOTH_THICK = 2.0  # ... this many strip thicknesses (bends sharper than that fan the normals)


# ---------------------------------------------------------------- scale: helpers

def _rotate(a, deg):
    h, w = a.shape
    M = cv2.getRotationMatrix2D((w / 2, h / 2), deg, 1)
    r = cv2.warpAffine(a, M, (w, h), flags=cv2.INTER_LINEAR, borderValue=0)
    m = cv2.warpAffine(np.ones_like(a), M, (w, h), flags=cv2.INTER_NEAREST, borderValue=0)
    return r, m, M


def _profile(a, deg, mask=None, hp=30):
    """Mean of `a` along columns after rotating by deg, high-passed. Lines at angle deg
    in the original become columns; returns profile and the rotation matrix."""
    w = np.ones_like(a) if mask is None else mask.astype(np.float32)
    r, m, M = _rotate(a * w, deg)
    rw, _, _ = _rotate(w, deg)
    s, c = r.sum(0), rw.sum(0)
    p = np.where(c > 20, s / np.maximum(c, 1e-6), 0)
    p = p - cv2.GaussianBlur(p.reshape(1, -1).astype(np.float64), (0, 0), hp).ravel()
    p[c <= 20] = 0
    return p, M


def _norm_autocorr(p):
    """Autocorrelation normalized over the overlap at each lag (Pearson), so lags of
    half the profile length are not penalised; zeros in p (masked) are ignored."""
    n = len(p)
    w = (p != 0).astype(float)
    nf = 1 << int(np.ceil(np.log2(2 * n)))
    F = lambda x: np.fft.rfft(x, nf)
    I = lambda X: np.fft.irfft(X, nf)[:n]
    Fp, Fw, Fp2 = F(p), F(w), F(p * p)
    num = I(Fp * np.conj(Fp))
    cnt = I(Fw * np.conj(Fw))
    e1 = I(Fp2 * np.conj(Fw))          # sum p_i^2 over overlap (leading)
    e2 = I(Fw * np.conj(Fp2))          # sum p_{i+l}^2 over overlap (trailing)
    ac = num / np.sqrt(np.maximum(e1 * e2, 1e-12))
    return ac, cnt


def _comb_period(p, tmin, tmax):
    """Fundamental period of a line profile from the mean autocorrelation at T, 2T, 3T, ..."""
    p = p - p[p != 0].mean() * (p != 0)
    # clip at a robust noise level so a few strong edges (dish rim, a thick line) do not
    # outweigh the regular lines
    sig = 1.4826 * np.median(np.abs(p[p != 0])) + 1e-9
    p = np.clip(p, -6 * sig, 6 * sig)
    ac, cnt = _norm_autocorr(p)
    n = len(p)
    nvalid = cnt[0]
    lags = np.arange(n)
    Ts = np.arange(tmin, tmax, 0.25)
    sc = np.full(len(Ts), -1.0)
    for i, T in enumerate(Ts):
        ks = np.arange(1, 8)
        ks = ks[np.interp(ks * T, lags, cnt) >= max(1.0 * T, 0.25 * nvalid)]  # enough overlap
        if len(ks):
            # Pearson values over a small overlap are noisy: weight by the overlap fraction
            wgt = np.sqrt(np.interp(ks * T, lags, cnt) / nvalid)
            sc[i] = np.mean(np.interp(ks * T, lags, ac) * wgt)
    if sc.max() <= -1:
        return np.nan, 0.0
    best = sc.max()
    # smallest lag that is a strong local maximum: a multiple of the period scores about as
    # well as the period itself, a half-period (peaks averaged with troughs) scores far lower
    pk, _ = find_peaks(np.r_[-1, sc, -1])
    pk = pk - 1
    strong = pk[sc[pk] >= 0.6 * best]
    j = int(strong.min()) if len(strong) else int(np.argmax(sc))
    return Ts[j], sc[j]


def _best_angle(R, mask, center, span, step):
    angs = np.arange(center - span, center + span + 1e-9, step)
    sc = [np.var(_profile(R, a, mask)[0]) for a in angs]
    return angs[int(np.argmax(sc))]


# ---------------------------------------------------------------- scale: grid

def grid_scale(img, focus_xy=None, grid_in=0.5, roi_cells=2.5, angles=None):
    """px per mm from the blue mat grid. focus_xy (full-res px) sets where the local
    spacing is read; None uses the whole frame. angles: line directions from an earlier call.
    When a Petri dish is in the frame, the grid is read outside it first (its rim and the
    worn, refracted grid under it confuse the spacing); the whole frame is the fallback."""
    small = cv2.resize(img, None, fx=1 / DS, fy=1 / DS, interpolation=cv2.INTER_AREA)
    lab = cv2.cvtColor(small, cv2.COLOR_BGR2LAB).astype(np.float32)
    g = lab[..., 2]                                   # grid lines are yellowish: bright in b*
    R = g - cv2.GaussianBlur(g, (0, 0), 6)
    L = lab[..., 0]
    # keep the strip out of the projections. Over the whole frame anything strongly yellow
    # goes (strip, dust smears, protractor lines); outside a dish only warm *and* bright
    # pixels go, since there the clean grid lines are themselves strongly yellow
    warm = g > np.median(g) + 25
    k = np.ones((15, 15), np.uint8)
    valid = cv2.dilate(warm.astype(np.uint8), k) == 0
    valid_bright = cv2.dilate((warm & (L > np.median(L) + 30)).astype(np.uint8), k) == 0
    h, w = g.shape

    def measure(valid, angles, focus_xy):
        if angles is None:
            a1 = _best_angle(R, valid, 90, 90, 1.0)
            a1 = _best_angle(R, valid, a1, 1.5, 0.1)
            a2 = _best_angle(R, valid, a1 + 90, 4, 0.1)
        else:
            a1, a2 = angles

        tmin, tmax = 300 / DS, 4000 / DS
        fams = []
        for a in (a1, a2):
            p, _ = _profile(R, a, valid)
            T, s = _comb_period(p, tmin, tmax)
            fams.append(dict(angle=a % 180, T=T, score=s))
        # both families must agree (square grid); if one is a multiple of the other, the
        # family with the stronger comb decides
        Ta, Tb = fams[0]["T"], fams[1]["T"]
        for k in (2, 3):
            if abs(Ta / Tb - k) < 0.08 * k or abs(Tb / Ta - k) < 0.08 * k:
                strong = max(fams, key=lambda f: f["score"])
                for f_ in fams:
                    f_["T"] = strong["T"]

        # local period around the strip (perspective changes the spacing across the frame):
        # same comb, on a window of +-roi_cells cells, searched within 15% of the global period
        for f in fams:
            f["T_local"] = f["T"]
        if focus_xy is not None:
            cx, cy = focus_xy[0] / DS, focus_xy[1] / DS
            for f in fams:
                rad = roi_cells * f["T"]
                yy, xx = np.mgrid[0:h, 0:w]
                local = valid & (np.abs(xx - cx) < rad) & (np.abs(yy - cy) < rad)
                if local.sum() < (1.5 * f["T"]) ** 2:
                    continue
                p, _ = _profile(R, f["angle"], local)
                nz = np.nonzero(p)[0]
                p = p[nz.min():nz.max() + 1]
                Tl, sl = _comb_period(p, 0.85 * f["T"], 1.15 * f["T"])
                if np.isfinite(Tl) and sl > 0.3 * f["score"]:
                    f["T_local"] = Tl
        px_per_cell = DS * np.mean([f["T_local"] for f in fams])
        px_mm = px_per_cell / (grid_in * MM_PER_IN)
        aniso = float(fams[0]["T_local"] / fams[1]["T_local"])
        return dict(method="grid", px_per_mm=px_mm, px_per_unit=px_per_cell, unit=f"{grid_in} in",
                    angle_deg=fams[0]["angle"], aniso=aniso,
                    # the two line families of a square grid must give the same spacing
                    scale_ok=bool(abs(aniso - 1) < 0.05), detail=fams)

    whole = measure(valid, angles, focus_xy)
    whole.update(dish_px=None, region="whole frame")
    dish = find_dish(img)
    if dish is None:
        return whole
    whole["dish_px"] = list(map(float, dish))
    # line directions from the whole frame, spacing from the mat outside the dish
    yy, xx = np.mgrid[0:h, 0:w]
    x0, y0, r = (v / DS for v in dish)
    outside = valid_bright & (np.hypot(xx - x0, yy - y0) > 1.05 * r)
    if outside.sum() < 0.2 * h * w:
        return whole
    # the strips lie in the dish, so there is no local window to read: one spacing outside it
    out = measure(outside, [d["angle"] for d in whole["detail"]], None)
    if not out["scale_ok"]:
        return whole
    out.update(dish_px=whole["dish_px"], region="outside dish")
    return out


def find_dish(img, min_score=0.18):
    """Petri dish rim as (x, y, r) in full-res px, or None. Hough circles are checked against
    the rim itself: a bright, neutral (not tissue-colored) line along the in-frame arc."""
    d = 8
    s = cv2.resize(img, None, fx=1 / d, fy=1 / d, interpolation=cv2.INTER_AREA)
    lab = cv2.cvtColor(s, cv2.COLOR_BGR2LAB).astype(np.float32)
    L, bs = lab[..., 0], lab[..., 2]
    h, w = L.shape
    g = cv2.GaussianBlur(L.astype(np.uint8), (5, 5), 1.5)
    cs = cv2.HoughCircles(g, cv2.HOUGH_GRADIENT, dp=2, minDist=10, param1=80, param2=40,
                          minRadius=int(0.25 * h), maxRadius=int(0.75 * w))
    if cs is None:
        return None
    rim = (cv2.morphologyEx(L, cv2.MORPH_TOPHAT, np.ones((9, 9), np.uint8)) > 12) & (bs < np.median(bs) + 12)
    rim = cv2.dilate(rim.astype(np.uint8), np.ones((3, 3), np.uint8))
    th = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    best, best_sc = None, min_score
    for x, y, r in cs[0][:10]:
        px, py = x + r * np.cos(th), y + r * np.sin(th)
        ins = (px >= 0) & (px < w - 1) & (py >= 0) & (py < h - 1)
        if ins.sum() < 0.2 * len(th):
            continue
        sc = rim[py[ins].astype(int), px[ins].astype(int)].mean()
        if sc > best_sc:
            best, best_sc = (x * d, y * d, r * d), sc
    return best


# ---------------------------------------------------------------- scale: ruler

def ruler_scale(img, tick_mm=None):
    """px per mm from a steel ruler: tick period along each long edge."""
    small = cv2.resize(img, None, fx=1 / DS, fy=1 / DS, interpolation=cv2.INTER_AREA)
    gray = cv2.cvtColor(small, cv2.COLOR_BGR2GRAY).astype(np.float32)
    # ruler = large dark region against a light background
    t, _ = cv2.threshold(gray.astype(np.uint8), 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    dark = (gray < t).astype(np.uint8)
    dark = cv2.morphologyEx(dark, cv2.MORPH_CLOSE, np.ones((25, 25), np.uint8))
    n, lab, st, _ = cv2.connectedComponentsWithStats(dark)
    k = 1 + int(np.argmax(st[1:, cv2.CC_STAT_AREA]))
    ys, xs = np.nonzero(lab == k)
    (cx, cy), (rw, rh), ang = cv2.minAreaRect(np.c_[xs, ys].astype(np.float32))
    # long axis direction
    if rw < rh:
        ang += 90
        rw, rh = rh, rw
    width = rh                                    # ruler width (small px)
    th = np.radians(ang)
    u = np.array([np.cos(th), np.sin(th)])        # along ruler
    v = np.array([-np.sin(th), np.cos(th)])       # across ruler
    # sample bands just inside each long edge, along the ruler
    hp = gray - cv2.GaussianBlur(gray, (0, 0), 4)
    L = int(rw * 0.95)
    s = np.arange(-L / 2, L / 2)
    edges = []
    for side in (-1, 1):
        offs = np.linspace(0.30, 0.44, 12) * width * side   # 6-20% of width in from the edge
        prof = []
        for o in offs:
            X = cx + s * u[0] + o * v[0]
            Y = cy + s * u[1] + o * v[1]
            prof.append(ndi.map_coordinates(hp, [Y, X], order=1, mode="nearest"))
        p = -np.mean(prof, 0)                     # ticks are dark
        p = p - cv2.GaussianBlur(p.reshape(1, -1).astype(np.float64), (0, 0), 20).ravel()
        T, sc = _comb_period(p, width / 60, width / 8)
        edges.append(dict(side=side, T=T * DS, score=sc))
    if tick_mm is not None:
        e = max(edges, key=lambda e: e["score"])
        mm = e
        mm["tick_mm"] = tick_mm
    else:
        # classify the edges by their period ratio: mm vs 1/16 in (1.5875 mm) or 1/32 in (0.794 mm)
        a, b = edges
        r = b["T"] / a["T"]
        cands = {"b=1/16in": (1.5875, a, 1.0), "b=1/32in": (0.79375, a, 1.0),
                 "a=1/16in": (1 / 1.5875, b, 1.0), "a=1/32in": (1 / 0.79375, b, 1.0)}
        name, (ratio, mm, tmm) = min(cands.items(), key=lambda kv: abs(np.log(r / kv[1][0])))
        if abs(np.log(r / ratio)) > 0.05:
            raise RuntimeError(f"ruler edges have period ratio {r:.3f}; cannot tell the mm edge. "
                               "Pass --ruler-tick-mm")
        mm["tick_mm"] = tmm
    px_mm = mm["T"] / mm["tick_mm"]
    return dict(method="ruler", px_per_mm=px_mm, px_per_unit=px_mm * 10, unit="1 cm",
                angle_deg=ang % 180, aniso=1.0, scale_ok=True, detail=edges,
                ruler_poly=(cv2.boxPoints(((cx, cy), (rw, rh), ang)) * DS).tolist())


def detect_background(img):
    """'grid' (blue mat) or 'ruler' (gray/neutral)."""
    small = cv2.resize(img, None, fx=1 / 16, fy=1 / 16, interpolation=cv2.INTER_AREA)
    lab = cv2.cvtColor(small, cv2.COLOR_BGR2LAB).astype(np.float32)
    blue = np.mean(lab[..., 2] < 115)
    return "grid" if blue > 0.3 else "ruler"


# ---------------------------------------------------------------- segmentation

def tissue_mask(img, px_mm, use_L=True, ruler_poly=None, grid_angles=None):
    """Binary tissue mask at WORK_PX_MM px/mm, and the resampling factor. use_L adds the
    brightness test (blue mat); on the ruler the white board is as bright as tissue, so
    tissue is found by b* off the ruler and by brightness alone on it (ruler_poly, full-res px): the edge
    on the steel is RULER_EDGE_FRAC of the way from the steel's brightness to the tissue's. On the mat
    (use_L) the result is refined by grow_on_mat (grid_angles: the grid's line directions from grid_scale)."""
    f = WORK_PX_MM / px_mm
    im = cv2.resize(img, None, fx=f, fy=f, interpolation=cv2.INTER_AREA)
    lab = cv2.cvtColor(im, cv2.COLOR_BGR2LAB).astype(np.float32)
    bst = cv2.GaussianBlur(lab[..., 2], (0, 0), 1.0)
    L = cv2.GaussianBlur(lab[..., 0], (0, 0), 1.0)
    # tissue is both warmer (b*) and brighter (L) than the background; grid lines and dust
    # are warm but not brighter, dish rims and glare are bright but not warm
    bg_b, bg_L = np.median(bst), np.median(L)
    sb = 1.4826 * np.median(np.abs(bst - bg_b)) + 1
    sL = 1.4826 * np.median(np.abs(L - bg_L)) + 1
    # warmer than the mat by 4 of its spreads, at most SEED_WARM_MAX on the mat: a stained, patchy mat has a wide
    # spread, and the brightness test keeps the (warm, not bright) printed lines out
    warm_min = max(15, 4 * sb) if not use_L else min(max(15, 4 * sb), SEED_WARM_MAX)
    seed = (bst > bg_b + warm_min) & ((not use_L) | (L > bg_L + max(20, 3 * sL)))
    if seed.sum() < 0.5 * WORK_PX_MM ** 2:
        raise RuntimeError("no tissue found")
    # tissue level from the seeds, leaving out near-saturated glare, which would set the
    # brightness threshold above the (darker, pinker) medial layer
    seed &= L < 235
    tb, tL = np.percentile(bst[seed], 75), np.median(L[seed])
    thr = (bg_b + 0.5 * (tb - bg_b), bg_L + 0.4 * (tL - bg_L))
    m = ((bst > thr[0]) & ((not use_L) | (L > thr[1]))).astype(np.uint8)
    if ruler_poly is not None:
        on = np.zeros(m.shape, np.uint8)
        cv2.fillPoly(on, [np.round(np.array(ruler_poly) * f).astype(np.int32)], 1)
        on = cv2.erode(on, np.ones((5, 5), np.uint8)).astype(bool)
        t, _ = cv2.threshold(L[on].astype(np.uint8).reshape(-1, 1), 0, 255,
                             cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        # on the dark steel, brightness alone decides: fluid films there can be warm (b*) but are dark, and
        # the strip's rounded side shows as a gray band between steel and tissue; both stay below the edge
        lo, hi = np.percentile(L[on], 50), np.percentile(L[on & (L > t)], RULER_TISSUE_PCT)
        m[on] = (L[on] > lo + RULER_EDGE_FRAC * (hi - lo)).astype(np.uint8)
    r = max(1, int(round(OPEN_MM * WORK_PX_MM)))
    disk = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
    m = cv2.morphologyEx(m, cv2.MORPH_OPEN, disk)
    rc = max(1, int(round(CLOSE_MM * WORK_PX_MM)))
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE,
                         cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * rc + 1, 2 * rc + 1)))
    # fill small holes (layer gaps), not the inside of C-shaped strips
    holes = ndi.binary_fill_holes(m) & ~m.astype(bool)
    hl, nh = ndi.label(holes)
    if nh:
        sizes = ndi.sum(holes, hl, range(1, nh + 1))
        small_h = np.isin(hl, 1 + np.nonzero(sizes < 1.0 * WORK_PX_MM ** 2)[0])
        m[small_h] = 1
    if use_L:
        m = grow_on_mat(m, cv2.cvtColor(im, cv2.COLOR_BGR2LAB).astype(np.float32), grid_angles or [])
    return m, f, bst, thr


def _disk(r):
    r = max(1, int(round(r)))
    return cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))


def _tissue_like(lab, bg, C, core):
    """Pixels that look like a strip's tissue against the local mat (bg): brighter than the mat and differing from it
    by at least GROW_FRAC of the core's contrast C (the edge half-way through a blurred step), or darker than the mat
    but warm (b* at least half-way from the mat's to the core's): darker layers and the strip's shaded side. Printed
    grid lines are not darker than the mat (or are found by _grid_lines), and shadows on the blue mat are not warm."""
    dE = np.linalg.norm(lab - bg, axis=-1)
    bright = (lab[..., 0] > bg[0]) & (dE >= GROW_FRAC * C)
    dark_warm = (lab[..., 0] <= bg[0]) & (lab[..., 2] >= bg[2] + 0.5 * max(core[2] - bg[2], 10))
    return bright | dark_warm


def _grid_lines(lab, angles, seeds):
    """The mat's printed grid lines near the strips: warm pixels in straight runs of at least LINE_MM along a grid
    direction, found where no tissue is (beyond the seeds' growth zone) and continued along the line into that zone."""
    W = WORK_PX_MM
    warm = (lab[..., 2] > np.median(lab[..., 2]) + 8).astype(np.uint8)
    warm[cv2.dilate(seeds, _disk((GROW_MM + 0.25) * W)) > 0] = 0
    out = np.zeros(warm.shape, np.uint8)
    for a in angles:
        def kernel(mm):
            n = int(mm * W) | 1
            k = np.zeros((n, n), np.uint8)
            c = n // 2
            dx, dy = np.cos(np.radians(a)) * c, np.sin(np.radians(a)) * c
            cv2.line(k, (int(round(c - dx)), int(round(c - dy))), (int(round(c + dx)), int(round(c + dy))), 1, 1)
            return k
        runs = cv2.morphologyEx(warm, cv2.MORPH_OPEN, kernel(LINE_MM))
        out |= cv2.dilate(runs, kernel(2 * (GROW_MM + 0.5)))
    return cv2.dilate(out, _disk(2)) > 0


def _keep_topology(out, lc, pieces):
    """Cut the mask where growth from parts of one piece that lie far apart along it met (e.g. the two ends of a
    C-shaped strip): seed pixels are labeled by CHUNK_MM chunks of their piece's centerline, every mask pixel takes
    the label of the nearest seed (geodesically), and pixels between labels of one piece more than CUT_CHUNKS apart
    are removed (a 0.2 mm cut)."""
    W = WORK_PX_MM
    markers = np.zeros(out.shape, np.int32)
    piece, pos = [0], [0]
    for k in pieces:
        comp = lc == k
        if comp.sum() < 50:
            continue
        cl = centerline(comp)
        s_arc = np.r_[0, np.cumsum(np.hypot(*np.diff(cl, axis=0).T))]
        ys, xs = np.nonzero(comp)
        _, j = cKDTree(cl).query(np.c_[ys, xs])
        chunk = (s_arc[j] // (CHUNK_MM * W)).astype(int)
        for c in np.unique(chunk):
            markers[ys[chunk == c], xs[chunk == c]] = len(piece)
            piece.append(k)
            pos.append(c)
    if len(piece) == 1:
        return out
    piece, pos = np.array(piece), np.array(pos)
    lab = watershed(np.zeros(out.shape, np.float32), markers, mask=out.astype(bool))
    cut = np.zeros(out.shape, bool)
    for dy, dx in ((0, 1), (1, 0)):
        a, b = lab[:lab.shape[0] - dy, :lab.shape[1] - dx], lab[dy:, dx:]
        far = (a > 0) & (b > 0) & (piece[a] == piece[b]) & (np.abs(pos[a] - pos[b]) > CUT_CHUNKS)
        cut[:lab.shape[0] - dy, :lab.shape[1] - dx] |= far
        cut[dy:, dx:] |= far
    cut = cv2.dilate(cut.astype(np.uint8), _disk(0.1 * W)) > 0
    out = out.copy()
    out[cut & (markers == 0)] = 0
    return out


def grow_on_mat(m0, lab, grid_angles):
    """Grow the mat's threshold mask m0 (the confident tissue) into the rest of each strip's tissue.
    Pieces of at least MIN_AREA_MM2, and smaller pieces bridged to them, are strip pieces; isolated small pieces
    (printed numerals, ruler ticks, specks) stay as they are. Each strip piece is judged against the mat in a ring
    around it (RING_MM) and its own core color: it grows up to GROW_MM into connected _tissue_like pixels, staying out
    of printed grid lines. Pieces within BRIDGE_MM whose gap is tissue-colored all along are joined inside a
    corridor (a strip whose faint middle the threshold mask lost). Where a grid line crosses a strip, the line is
    filled between tissue on both sides; growth never joins far-apart parts of one piece (_keep_topology)."""
    W = WORK_PX_MM
    lab = np.stack([cv2.GaussianBlur(c, (0, 0), 1.0) for c in cv2.split(lab)], -1)
    lc, n = ndi.label(m0)
    if n == 0:
        return m0
    sizes = ndi.sum(m0, lc, range(1, n + 1))
    seeds = np.isin(lc, 1 + np.nonzero(sizes >= SEED_MIN_MM2 * W * W)[0]).astype(np.uint8)
    near_any = cv2.dilate(seeds, _disk(RING_MM[0] * W)) > 0
    lines = _grid_lines(lab, grid_angles, seeds) if grid_angles else np.zeros(m0.shape, bool)
    lc, n = ndi.label(seeds)
    info, box, pts = {}, {}, {}
    for k in range(1, n + 1):
        ys, xs = np.nonzero(lc == k)
        pad = int((RING_MM[1] + 1) * W)
        y0, y1, x0, x1 = max(ys.min() - pad, 0), ys.max() + pad, max(xs.min() - pad, 0), xs.max() + pad
        cc = (lc[y0:y1, x0:x1] == k).astype(np.uint8)
        Lw = lab[y0:y1, x0:x1]
        ring = (cv2.dilate(cc, _disk(RING_MM[1] * W)) > 0) & ~near_any[y0:y1, x0:x1]
        if ring.sum() < 100:
            continue
        bg = np.median(Lw[ring], 0)
        er = cv2.erode(cc, _disk(3)) > 0
        core = np.median(Lw[er if er.any() else cc > 0], 0)
        info[k], box[k], pts[k] = (bg, float(np.linalg.norm(core - bg)), core), (y0, y1, x0, x1), np.c_[ys, xs]
    big = {k for k in info if len(pts[k]) >= MIN_AREA_MM2 * W * W}
    bridges, strip_pieces = [], set(big)
    keys = sorted(info)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if a not in big and b not in big:
                continue
            A, B = pts[a][::max(1, len(pts[a]) // 400)], pts[b][::max(1, len(pts[b]) // 400)]
            dd = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1)
            ia, ib = np.unravel_index(np.argmin(dd), dd.shape)
            if dd[ia, ib] > BRIDGE_MM * W:
                continue
            pa, pb = A[ia].astype(float), B[ib].astype(float)
            bg, C = (info[a][0] + info[b][0]) / 2, min(info[a][1], info[b][1])
            core = (info[a][2] + info[b][2]) / 2
            r = int(0.3 * W)
            ok = all(np.linalg.norm(lab[max(y - r, 0):y + r + 1, max(x - r, 0):x + r + 1] - bg, axis=-1).max()
                     >= BRIDGE_FRAC * C
                     for y, x in (np.round(pa + t * (pb - pa)).astype(int) for t in np.linspace(0, 1, int(dd[ia, ib]) + 1)))
            if ok:
                bridges.append((pa, pb, bg, C, core))
                strip_pieces |= {a, b}
    out = m0.copy()
    for k in strip_pieces:
        y0, y1, x0, x1 = box[k]
        cc = (lc[y0:y1, x0:x1] == k).astype(np.uint8)
        bg, C, core = info[k]
        allowed = (_tissue_like(lab[y0:y1, x0:x1], bg, C, core) & ~lines[y0:y1, x0:x1]
                   & (cv2.dilate(cc, _disk(GROW_MM * W)) > 0)) | (cc > 0)
        out[y0:y1, x0:x1] |= reconstruction(cc, allowed.astype(np.uint8), method="dilation").astype(np.uint8)
    for pa, pb, bg, C, core in bridges:
        corridor = np.zeros(out.shape, np.uint8)
        cv2.line(corridor, (int(pa[1]), int(pa[0])), (int(pb[1]), int(pb[0])), 1, int(2 * CORRIDOR_MM * W))
        allowed = (_tissue_like(lab, bg, C, core) & ~lines & (corridor > 0)) | (out > 0)
        out = reconstruction(out, allowed.astype(np.uint8), method="dilation").astype(np.uint8)
    out |= ((cv2.morphologyEx(out, cv2.MORPH_CLOSE, _disk(0.6 * W)) > 0) & lines).astype(np.uint8)
    out = cv2.morphologyEx(out, cv2.MORPH_OPEN, _disk(OPEN_MM * W))
    out = cv2.morphologyEx(out, cv2.MORPH_CLOSE, _disk(CLOSE_MM * W))
    holes = ndi.binary_fill_holes(out) & ~out.astype(bool)
    hl, nh = ndi.label(holes)
    if nh:
        hs = ndi.sum(holes, hl, range(1, nh + 1))
        out[np.isin(hl, 1 + np.nonzero(hs < 1.0 * W ** 2)[0])] = 1
    return _keep_topology(out, lc, strip_pieces)


def strip_components(m):
    n, lab, st, cen = cv2.connectedComponentsWithStats(m)
    h, w = m.shape
    out = []
    for i in range(1, n):
        x, y, bw, bh, area = st[i]
        if area < MIN_AREA_MM2 * WORK_PX_MM ** 2:
            continue
        touches = x <= 1 or y <= 1 or x + bw >= w - 1 or y + bh >= h - 1
        out.append(dict(label=i, area_mm2=area / WORK_PX_MM ** 2, touches_border=touches,
                        centroid=cen[i], bbox=(x, y, bw, bh)))
    return out, lab


# ---------------------------------------------------------------- thickness

def centerline(comp):
    """Longest geodesic path through the skeleton of a binary component (row, col)."""
    sk = skeletonize(comp.astype(bool))
    ys, xs = np.nonzero(sk)
    idx = -np.ones(sk.shape, int)
    idx[ys, xs] = np.arange(len(ys))
    rows, cols, wts = [], [], []
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        y2, x2 = ys + dy, xs + dx
        ok = (y2 >= 0) & (y2 < sk.shape[0]) & (x2 >= 0) & (x2 < sk.shape[1])
        ok[ok] = sk[y2[ok], x2[ok]]
        rows += list(idx[ys[ok], xs[ok]])
        cols += list(idx[y2[ok], x2[ok]])
        wts += [np.hypot(dy, dx)] * int(ok.sum())
    G = coo_matrix((wts, (rows, cols)), shape=(len(ys),) * 2).tocsr()
    G = G + G.T
    d0 = dijkstra(G, indices=0)
    d0[~np.isfinite(d0)] = -1
    a = int(np.argmax(d0))
    da, pred = dijkstra(G, indices=a, return_predecessors=True)
    da[~np.isfinite(da)] = -1
    b = int(np.argmax(da))
    path = [b]
    while path[-1] != a:
        path.append(pred[path[-1]])
    path = np.array(path[::-1])
    return np.c_[ys[path], xs[path]].astype(float)


def _first_crossing(prof, step):
    """Distance from index 0 to the first sub-pixel 0.5 crossing of a falling profile."""
    below = np.nonzero(prof < 0.5)[0]
    if len(below) == 0 or below[0] == 0:
        return np.nan
    i = below[0]
    p0, p1 = prof[i - 1], prof[i]
    return (i - 1 + (p0 - 0.5) / (p0 - p1)) * step


def measure_strip(comp):
    """Chord thickness along the normal at stations along the centerline (work px)."""
    P = centerline(comp)
    # typical thickness (inscribed-circle diameter on the medial axis) sets the smoothing
    dt = ndi.distance_transform_edt(comp)
    D = 2 * np.median(dt[P[:, 0].astype(int), P[:, 1].astype(int)])
    seg = np.hypot(*np.diff(P, axis=0).T)
    s = np.r_[0, np.cumsum(seg)]
    Ltot = s[-1]
    ds = STEP_MM * WORK_PX_MM
    ss = np.arange(0, Ltot, ds / 4)
    y = np.interp(ss, s, P[:, 0])
    x = np.interp(ss, s, P[:, 1])
    win = int(max(SMOOTH_MM * WORK_PX_MM, SMOOTH_THICK * D) / (ds / 4)) | 1
    if len(ss) > win + 2:
        y = savgol_filter(y, win, 2)
        x = savgol_filter(x, win, 2)
        ty = savgol_filter(y, win, 2, deriv=1)
        tx = savgol_filter(x, win, 2, deriv=1)
    else:
        ty, tx = np.gradient(y), np.gradient(x)
    nrm = np.hypot(tx, ty)
    ty, tx = ty / nrm, tx / nrm
    ny, nx = tx, -ty
    soft = cv2.GaussianBlur(comp.astype(np.float32), (0, 0), 0.7)
    rstep = 0.25
    rr = np.arange(0, 0.6 * max(comp.shape), rstep)
    keep = (ss >= TRIM_FRAC * Ltot) & (ss <= (1 - TRIM_FRAC) * Ltot)
    sel = np.nonzero(keep)[0][::4]
    rows = []
    for i in sel:
        sides = []
        for sgn in (1, -1):
            Y = y[i] + sgn * rr * ny[i]
            X = x[i] + sgn * rr * nx[i]
            pr = ndi.map_coordinates(soft, [Y, X], order=1, mode="constant", cval=0)
            sides.append(_first_crossing(pr, rstep))
        rows.append(dict(s_px=ss[i], y=y[i], x=x[i], ny=ny[i], nx=nx[i],
                         d_plus=sides[0], d_minus=sides[1], chord_px=sides[0] + sides[1]))
    st = pd.DataFrame(rows)
    # the centerline is the medial axis, so the two half-chords should match; a large
    # mismatch means the normal ran out through a ragged edge or along the strip
    asym = np.abs(st.d_plus - st.d_minus) / st.chord_px
    st["ok"] = np.isfinite(st.chord_px) & (asym < 0.5)
    st["crosses"] = _crossing_chords(st)
    st["ok"] &= ~st.crosses
    return st, dict(y=y, x=x), Ltot


def _crossing_chords(st):
    """True for stations whose chord crosses another station's chord. Chords cross where the
    strip bends more sharply than its thickness, or where the normal leaves through a lobe;
    such a chord is not a thickness across the wall."""
    ok = np.isfinite(st.chord_px.values)
    A = np.c_[st.x + st.d_plus * st.nx, st.y + st.d_plus * st.ny]
    B = np.c_[st.x - st.d_minus * st.nx, st.y - st.d_minus * st.ny]
    n = len(st)
    cross = np.zeros(n, bool)

    def orient(p, q, r):
        return np.sign((q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1])
                       - (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0]))
    idx = np.nonzero(ok)[0]
    for i in idx:
        j = idx[idx != i]
        a, b = A[i], B[i]
        c, d = A[j], B[j]
        hit = (orient(a, b, c) * orient(a, b, d) < 0) & (orient(c, d, a) * orient(c, d, b) < 0)
        cross[i] = hit.any()
    return cross


# ---------------------------------------------------------------- driver

def process(path, background="auto", grid_in=0.5, ruler_tick_mm=None, px_per_mm=None, outdir=None,
            pick=None):
    """Measure every complete strip in one photo. pick chooses the strips the photo is
    for: all (default; the strips are pieces of the test's specimen), or one by position: center
    (nearest the frame center), upper, lower, left or right."""
    img = cv2.imread(path)
    if img is None:
        raise RuntimeError(f"cannot read {path}")
    bg = detect_background(img) if background == "auto" else background
    if px_per_mm is not None:
        scale = dict(method="manual", px_per_mm=px_per_mm, px_per_unit=10 * px_per_mm, unit="1 cm",
                     aniso=1.0, scale_ok=True)
    elif bg == "grid":
        scale = grid_scale(img, None, grid_in)
    else:
        scale = ruler_scale(img, ruler_tick_mm)

    def segment():
        m, f, _, _ = tissue_mask(img, scale["px_per_mm"], use_L=(bg == "grid"),
                                 ruler_poly=scale.get("ruler_poly"),
                                 grid_angles=[d["angle"] for d in scale.get("detail", []) if "angle" in d]
                                 if bg == "grid" else None)
        comps, lab = strip_components(m)
        return m, f, comps, lab

    m, f, comps, lab = segment()
    if scale["method"] == "grid" and any(not c["touches_border"] for c in comps):
        # re-read the grid spacing locally, around the strips
        cx, cy = np.mean([c["centroid"] for c in comps if not c["touches_border"]], 0) / f
        scale = grid_scale(img, (cx, cy), grid_in, angles=[d["angle"] for d in scale["detail"]])
        m, f, comps, lab = segment()

    dish = scale.get("dish_px")
    if dish is not None:
        # blobs lying along the dish rim are rim highlights, not tissue
        yy, xx = np.nonzero(lab)
        rr = np.hypot(xx / f - dish[0], yy / f - dish[1])
        on_rim = np.abs(rr - dish[2]) < 0.06 * dish[2]
        frac = ndi.mean(on_rim, lab[yy, xx], [c["label"] for c in comps])
        for c, fr in zip(comps, frac):
            c["on_rim"] = fr > 0.5
    strips = []
    for c in comps:
        if c["touches_border"] or c.get("on_rim", False):
            continue
        st, cl, Lpx = measure_strip((lab == c["label"]).astype(np.uint8))
        st["thickness_mm"] = st.chord_px / WORK_PX_MM
        st["s_mm"] = st.s_px / WORK_PX_MM
        t = st.loc[st.ok, "thickness_mm"]
        # a strip is long and thin; round blobs (glare, debris) are not measured
        c["is_strip"] = len(t) >= 10 and Lpx / WORK_PX_MM >= 3 * t.median()
        if c["is_strip"]:
            c.update(stations=st, cl=cl, Lpx=Lpx)
            strips.append(c)
    if len(strips) > 1:
        # a strip stands out from the background around it; much fainter "strips" are stains or films
        im_lab = cv2.cvtColor(cv2.resize(img, None, fx=f, fy=f, interpolation=cv2.INTER_AREA),
                              cv2.COLOR_BGR2LAB).astype(np.float32)
        near = cv2.dilate(m, _disk(RING_MM[0] * WORK_PX_MM)) > 0
        for c in strips:
            comp = (lab == c["label"]).astype(np.uint8)
            ring = (cv2.dilate(comp, _disk(RING_MM[1] * WORK_PX_MM)) > 0) & ~near
            er = cv2.erode(comp, _disk(3)) > 0
            core = np.median(im_lab[er if er.any() else comp > 0], 0)
            c["contrast"] = float(np.linalg.norm(core - np.median(im_lab[ring], 0))) if ring.any() else np.inf
        best = max(c["contrast"] for c in strips)
        for c in strips:
            if c["contrast"] < STRIP_CONTRAST_FRAC * best:
                print(f"{os.path.basename(path)}: dropped a strip-shaped region at ({c['centroid'][0] / WORK_PX_MM:.1f}, "
                      f"{c['centroid'][1] / WORK_PX_MM:.1f}) mm: contrast {c['contrast']:.0f} against {best:.0f} for the "
                      "strongest strip (stained mat or a film)")
        strips = [c for c in strips if c["contrast"] >= STRIP_CONTRAST_FRAC * best]
    strips.sort(key=lambda c: -c["area_mm2"])
    h, w = m.shape
    rule = {"centre": "center"}.get(pick, pick) or "all"      # older strip_selection.json files say "centre"
    picker = {"all": None, "center": lambda c: np.hypot(*(c["centroid"] - [w / 2, h / 2])),
              "upper": lambda c: c["centroid"][1], "lower": lambda c: -c["centroid"][1],
              "left": lambda c: c["centroid"][0], "right": lambda c: -c["centroid"][0]}
    if rule not in picker:
        raise RuntimeError(f"unknown strip selection {rule!r} (use {', '.join(picker)})")
    chosen = ({c["label"] for c in strips} if rule == "all"
              else {min(strips, key=picker[rule])["label"]} if strips else set())

    stem = os.path.splitext(os.path.basename(path))[0]
    results = []
    for k, c in enumerate(strips):
        t = c["stations"].loc[c["stations"].ok, "thickness_mm"]
        q1, q3 = np.percentile(t, [25, 75])
        res = dict(image=os.path.basename(path), strip=k + 1, selected=c["label"] in chosen,
                   selected_by=rule,
                   background=bg, scale_method=scale["method"], px_per_mm=scale["px_per_mm"],
                   scale_unit=scale["unit"], px_per_unit=scale["px_per_unit"],
                   grid_aniso=scale["aniso"], scale_ok=scale["scale_ok"],
                   length_mm=c["Lpx"] / WORK_PX_MM, area_mm2=c["area_mm2"],
                   n_stations=len(t), mean_mm=t.mean(), sd_mm=t.std(ddof=1),
                   median_mm=t.median(), q1_mm=q1, q3_mm=q3, min_mm=t.min(), max_mm=t.max())
        results.append(res)
        c["res"] = res
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        for c in strips:
            c["stations"].assign(x_full=c["stations"].x / f, y_full=c["stations"].y / f).to_csv(
                os.path.join(outdir, f"{stem}_strip{c['res']['strip']}_stations.csv"), index=False)
        overlay(img, f, m, comps, strips, scale, os.path.join(outdir, f"{stem}_overlay.png"))
    return results, scale


def overlay(img, f, m, comps, strips, scale, out):
    im = cv2.resize(img, None, fx=f, fy=f, interpolation=cv2.INTER_AREA)[..., ::-1]
    h, w = m.shape
    fig, axes = plt.subplots(1, 1 + len(strips), figsize=(6 * (1 + len(strips)), 5),
                             gridspec_kw=dict(width_ratios=[1.6] + [1] * len(strips)))
    axes = np.atleast_1d(axes)
    ax = axes[0]
    ax.imshow(im)
    ax.contour(m, [0.5], colors="yellow", linewidths=0.6)
    for c in comps:
        x, y, bw, bh = c["bbox"]
        if c["touches_border"] or not c.get("is_strip", False):
            ax.add_patch(plt.Rectangle((x, y), bw, bh, fill=False, ec="gray", ls="--", lw=0.8))
    for c in strips:
        sel = c["res"]["selected"]
        ax.text(*c["centroid"], f" {c['res']['strip']}" + (f" ({c['res']['selected_by']})" if sel else ""),
                color="yellow" if sel else "white", fontsize=12, weight="bold")
    # scale bar: one grid cell / 1 cm
    L = scale["px_per_unit"] * f
    ax.plot([w * 0.05, w * 0.05 + L], [h * 0.95] * 2, color="w", lw=3)
    ax.text(w * 0.05, h * 0.93, f"{scale['unit']} = {scale['px_per_unit']:.0f} px "
            f"({scale['px_per_mm']:.1f} px/mm)", color="w", fontsize=9)
    ax.set_title(os.path.basename(out).replace("_overlay.png", "") + f"  [{scale['method']}]"
                 + ("" if scale["scale_ok"] else "  SCALE UNCERTAIN"), color="k" if scale["scale_ok"] else "r")
    ax.axis("off")
    for ax, c in zip(axes[1:], strips):
        st, r = c["stations"], c["res"]
        x0, y0, bw, bh = c["bbox"]
        pad = int(1.5 * WORK_PX_MM)
        ys, xs = slice(max(0, y0 - pad), y0 + bh + pad), slice(max(0, x0 - pad), x0 + bw + pad)
        ax.imshow(im[ys, xs], extent=(xs.start, xs.start + im[ys, xs].shape[1],
                                      ys.start + im[ys, xs].shape[0], ys.start))
        ax.contour(np.arange(xs.start, xs.start + m[ys, xs].shape[1]),
                   np.arange(ys.start, ys.start + m[ys, xs].shape[0]),
                   m[ys, xs], [0.5], colors="yellow", linewidths=0.8)
        ax.plot(c["cl"]["x"], c["cl"]["y"], "c-", lw=0.8)
        for _, sr in st.iloc[::3].iterrows():
            col = "lime" if sr.ok else "red"
            if np.isfinite(sr.chord_px):
                ax.plot([sr.x + sr.d_plus * sr.nx, sr.x - sr.d_minus * sr.nx],
                        [sr.y + sr.d_plus * sr.ny, sr.y - sr.d_minus * sr.ny], color=col, lw=0.6)
        ax.set_title(f"strip {r['strip']}: {r['mean_mm']:.2f} ± {r['sd_mm']:.2f} mm "
                     f"(median {r['median_mm']:.2f}, n={r['n_stations']})", fontsize=10)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def load_selection(folder):
    """{photo name: rule} from <folder>/strip_selection.json, e.g. {"DSC00503.jpg": "lower"}; rules: all, center,
    upper, lower, left, right, or exclude (the photo is left out). Photos not listed use all their strips."""
    f = os.path.join(folder, "strip_selection.json")
    if not os.path.exists(f):
        return {}
    with open(f) as fh:
        return {k: v for k, v in json.load(fh).items() if not k.startswith("_")}


AUTO_SOURCE = "measure_thickness.py"    # start of the source text of an automatic thickness_mm


def output_location(folder):
    """Output folder for a photo folder: <SAMPLE>/Thickness, <ID>/thickness_auto/test<N>/ for
    <SAMPLE>/Thickness, <ID>/test <N> (or test<N>), else <folder>/thickness_auto/."""
    d = os.path.abspath(folder)
    m = re.fullmatch(r"(?i)test\s*(\d+)", os.path.basename(d))
    if m and os.path.basename(os.path.dirname(d)).lower().startswith("thickness"):
        return os.path.join(os.path.dirname(d), "thickness_auto", f"test{int(m[1])}")
    return os.path.join(d, "thickness_auto")


def test_of_photo_folder(folder):
    """The test folder for a photo folder <SAMPLE>/Thickness, <ID>/test <N> (or test<N>), else None."""
    d = os.path.abspath(folder)
    m = re.fullmatch(r"(?i)test\s*(\d+)", os.path.basename(d))
    if not m or not os.path.basename(os.path.dirname(d)).lower().startswith("thickness"):
        return None
    hits = sorted(glob.glob(os.path.join(os.path.dirname(os.path.dirname(d)), f"test{int(m[1])}_*")))
    return hits[0] if len(hits) == 1 and os.path.isdir(hits[0]) else None


def pooled_thickness(rows, notes):
    """thickness_mm from result rows of one photo folder: mean and SD of all valid stations of the selected strips
    of the photos whose scale is reliable (scale_ok). notes ({photo: why it was left out}) is extended. Returns
    (thickness_mm or None, photos used)."""
    vals, photos = [], []
    for r in rows:
        if not r["selected"]:
            continue
        if not r["scale_ok"]:
            notes[r["image"]] = "scale uncertain"
            continue
        st = pd.read_csv(r["stations_csv"])
        vals.append(st.loc[st.ok, "thickness_mm"].to_numpy())
        if r["image"] not in photos:
            photos.append(r["image"])
    v = np.concatenate(vals) if vals else np.array([])
    if len(v) < 2:
        return None, photos
    return ({"H": round(float(v.mean()), 5), "H_sd": round(float(v.std(ddof=1)), 4),
             "source": f"{AUTO_SOURCE}: mean and SD of {len(v)} stations of the selected strips in " + ", ".join(photos)},
            photos)


def write_thickness(folder, th, keep_workbook=False):
    """Put th into the test_config.json of the photo folder's test. A workbook value it replaces is kept in
    thickness_mm.workbook; with keep_workbook the workbook value stays (or is put back) instead."""
    root = test_of_photo_folder(folder)
    if root is None:
        print(f"{folder}: not <SAMPLE>/Thickness, <ID>/test <N> with a matching test folder; nothing written")
        return
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "strain_tracking"))
    import testcfg
    T = testcfg.Test(root)
    name = f"{os.path.basename(os.path.dirname(root))}/{os.path.basename(root)}"
    if not T.config:
        print(f"{name}: no test_config.json yet (run make_config.py first); nothing written")
        return
    old = T.config.get("thickness_mm") or {}
    auto = str(old.get("source", "")).startswith(AUTO_SOURCE)
    book = old.get("workbook") if auto else ({k: old[k] for k in ("H", "H_sd", "source") if k in old} if old else None)
    if keep_workbook and book:
        T.config["thickness_mm"] = book
        T.write_config(T.config)
        print(f"{name}: workbook thickness_mm H = {book.get('H')} mm in use ({book.get('source', 'source not given')}); "
              f"photos {th['H']} mm")
        return
    if book:
        th = dict(th, workbook=book)
    T.config["thickness_mm"] = th
    T.write_config(T.config)
    print(f"{name}: thickness_mm H = {th['H']} ± {th['H_sd']} mm written"
          + (f" (workbook {book.get('H')} mm, {(th['H'] / book['H'] - 1) * 100:+.1f}%)" if book and book.get("H") else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="image files or folders")
    ap.add_argument("--recursive", action="store_true", help="search folders recursively")
    ap.add_argument("--background", choices=["auto", "grid", "ruler"], default="auto")
    ap.add_argument("--grid-in", type=float, default=0.5, help="grid spacing (in)")
    ap.add_argument("--ruler-tick-mm", type=float, help="tick spacing on the ruler edge to use "
                    "(skips the automatic mm/inch edge classification)")
    ap.add_argument("--px-per-mm", type=float, help="use this scale instead of detecting it")
    ap.add_argument("--out", help="output folder (default: <SAMPLE>/Thickness, <ID>/thickness_auto/test<N>)")
    ap.add_argument("--write", action="store_true",
                    help="put the result into each test's test_config.json as thickness_mm, in place of a workbook "
                         "value (kept in thickness_mm.workbook)")
    ap.add_argument("--keep-workbook", action="store_true",
                    help="with --write: keep (or put back) a workbook value; fill only tests without one")
    a = ap.parse_args()

    files = []
    for p in a.inputs:
        if os.path.isdir(p):
            pat = os.path.join(p, "**", "*") if a.recursive else os.path.join(p, "*")
            files += [f for f in glob.glob(pat, recursive=a.recursive)
                      if f.lower().endswith((".jpg", ".jpeg", ".png", ".tif", ".tiff"))
                      and "thickness_auto" not in f]
        else:
            files.append(p)
    rows, notes = {}, {}          # notes: {photo folder: {photo: why it gave no thickness}}
    for fpath in sorted(files):
        outdir = a.out or output_location(os.path.dirname(fpath))
        folder, photo = os.path.dirname(fpath), os.path.basename(fpath)
        notes.setdefault(folder, {})
        pick = load_selection(folder).get(photo)
        for old in glob.glob(os.path.join(glob.escape(outdir), os.path.splitext(photo)[0] + "_*")):
            os.remove(old)            # an earlier run's overlay and stations of this photo
        if pick == "exclude":
            print(f"{photo:16s} excluded in strip_selection.json")
            notes[folder][photo] = "excluded in strip_selection.json"
            continue
        try:
            res, sc = process(fpath, a.background, a.grid_in, a.ruler_tick_mm, a.px_per_mm, outdir, pick)
        except Exception as e:  # keep going through a folder
            print(f"{fpath}: FAILED ({e})")
            notes[folder][photo] = f"failed: {e}"
            continue
        for r in res:
            r["folder"] = os.path.dirname(fpath)
            r["stations_csv"] = os.path.join(outdir, f"{os.path.splitext(r['image'])[0]}_strip{r['strip']}_stations.csv")
            print(f"{os.path.basename(fpath):16s} strip {r['strip']}{'*' if r['selected'] else ' '} "
                  f"{sc['method']:6s} {sc['px_per_mm']:6.1f} px/mm{'' if sc['scale_ok'] else ' (?)'}  "
                  f"{r['mean_mm']:.3f} ± {r['sd_mm']:.3f} mm "
                  f"(median {r['median_mm']:.3f}, n={r['n_stations']}, L={r['length_mm']:.1f} mm)")
        if not res:
            print(f"{os.path.basename(fpath):16s} no complete strip found")
            notes[folder][photo] = "no complete strip found"
        rows.setdefault(outdir, []).extend(res)
    for outdir, rr in rows.items():
        pd.DataFrame(rr).drop(columns="stations_csv").to_csv(os.path.join(outdir, "summary.csv"), index=False)
    if a.write:
        by_folder = {}
        for rr in rows.values():
            for r in rr:
                by_folder.setdefault(r["folder"], []).append(r)
        for folder in sorted(notes):
            th, used = pooled_thickness(by_folder.get(folder, []), notes[folder])
            left = "; ".join(f"{p} ({why})" for p, why in sorted(notes[folder].items())) or "none"
            print(f"{folder}: thickness from {len(used)} photo(s): {', '.join(used) or 'none'}; left out: {left}")
            if th is None:
                print(f"{folder}: NO THICKNESS: no photo with a measured strip at a reliable scale; "
                      "the config keeps what it has")
            else:
                write_thickness(folder, th, a.keep_workbook)


if __name__ == "__main__":
    main()
