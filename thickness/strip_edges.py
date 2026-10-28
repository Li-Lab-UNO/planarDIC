"""
Wall thickness of a strip lying on its cut face in a photo: the width of the tissue band, measured across a
centerline the user clicks (the seed). The Python version of the MATLAB version's thicknessCore, step for
step; measure_interactive.py is the window around it.

Along the smoothed seed, every STEP_MM, the image is read along the local perpendicular (up to SEARCH_MM
either side). Each pixel gets a tissue score from its color (CIE Lab) against the tissue's and the local
background's (crossing): 1 is tissue, 0 background, on the blue mat and on the steel ruler alike. The edge
on each side is where the score first drops below one half, to a fraction of a pixel. The centerline is then
moved to the midpoints between the edges and the strip measured again, with its normals smoothed over at
least 2 wall thicknesses, so a rough seed does not widen the chords. As in measure_thickness.py, positions
are not used where an edge is missing, the chord is lopsided (one half 3 times the other) or crosses a
neighbor's, or within TRIM_MM of the strip's ends; positions far from the strip's median (3 robust SDs: a
loose fiber) are flagged as outliers. reason says why a position is not used.

Coordinates are image px, 0-based (pixel centers at integers), x to the right and y down.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d
from scipy.ndimage import map_coordinates

STEP_MM = 0.1          # spacing of the positions along the strip
SEARCH_MM = 3.0        # read this far either side of the centerline
TRIM_MM = 0.5          # positions this close to the strip's ends are not used


def measure_strip(rgb, seed, pxmm, step_mm=STEP_MM, search_mm=SEARCH_MM, trim_mm=TRIM_MM):
    """rgb: the photo (uint8, H x W x 3, RGB); seed: n x 2 clicked points (x, y); pxmm: image scale.
    Returns a dict of arrays, one value per position: x, y (center), nx, ny (unit normal), s (arc length),
    rp, rm (distance to the edge on the + and - side, px), thickness_mm, ok, outlier, reason; and method."""
    seed = np.asarray(seed, float)
    R = search_mm * pxmm
    lab, x0, y0 = lab_crop(rgb, seed, R + 1.5 * pxmm)
    C = curve(seed, step_mm * pxmm, pxmm)
    for p in range(2):
        S = edges(lab, x0, y0, C, R, pxmm)
        if p == 0 and S["ok"].sum() >= 3:                     # recenter between the edges
            m = (S["rp"] - S["rm"]) / 2
            ok = S["ok"]
            if (~ok).any():
                m[~ok] = interp1d(np.nonzero(ok)[0], m[ok], fill_value="extrapolate")(np.nonzero(~ok)[0])
            P = np.c_[S["x"] + m * S["nx"], S["y"] + m * S["ny"]]
            # as in measure_thickness.py, the normals are smoothed over at least 2 wall thicknesses
            D = np.median(S["rp"][ok] + S["rm"][ok])
            C = curve(P, step_mm * pxmm, pxmm, max(pxmm, 2 * D))
    # as in measure_thickness.py: a chord much longer on one side than the other ran out through a ragged
    # edge, and chords that cross a neighbor's (tight bends) measure the bend, not the wall
    reason = np.array([""] * len(S["ok"]), dtype=object)
    reason[~S["ok"]] = "no edge"
    with np.errstate(invalid="ignore"):
        lop = S["ok"] & (np.abs(S["rp"] - S["rm"]) / (S["rp"] + S["rm"]) >= 0.5)
    reason[lop] = "lopsided"
    S["ok"] = S["ok"] & ~lop
    crs = crossing_chords(S) & S["ok"]
    reason[crs] = "crosses a neighbor"
    S["ok"] = S["ok"] & ~crs
    S["thickness_mm"] = np.where(S["ok"], (S["rp"] + S["rm"]) / pxmm, np.nan)
    trim = (S["s"] < trim_mm * pxmm) | (S["s"] > S["s"][-1] - trim_mm * pxmm)
    reason[trim & S["ok"]] = "strip end"
    S["ok"] = S["ok"] & ~trim
    S["outlier"] = outliers(S["thickness_mm"], S["ok"])
    S["reason"] = reason
    S["method"] = "auto"
    return S


def measure_manual(edge_a, edge_b, pxmm, step_mm=STEP_MM):
    """Thickness between two clicked edges (n x 2 each): positions along their middle line, chord along its
    normal to each edge."""
    A, B = resample(np.asarray(edge_a, float), 200), resample(np.asarray(edge_b, float), 200)
    if (np.linalg.norm(A[0] - B[0]) + np.linalg.norm(A[-1] - B[-1])
            > np.linalg.norm(A[0] - B[-1]) + np.linalg.norm(A[-1] - B[0])):
        B = B[::-1]                                           # clicked in opposite directions
    C = curve((A + B) / 2, step_mm * pxmm, pxmm)
    k = len(C["x"])
    S = dict(C, rp=np.full(k, np.nan), rm=np.full(k, np.nan))
    for i in range(k):
        p, d = np.array([C["x"][i], C["y"][i]]), np.array([C["nx"][i], C["ny"][i]])
        t = np.r_[hit_polyline(p, d, edge_a), hit_polyline(p, d, edge_b)]
        if (t > 0).any() and (t < 0).any():
            S["rp"][i], S["rm"][i] = t[t > 0].min(), -t[t < 0].max()
    S["ok"] = np.isfinite(S["rp"]) & np.isfinite(S["rm"])
    S["reason"] = np.where(S["ok"], "", "no edge").astype(object)
    S["thickness_mm"] = np.where(S["ok"], (S["rp"] + S["rm"]) / pxmm, np.nan)
    S["outlier"] = outliers(S["thickness_mm"], S["ok"])
    S["method"] = "manual"
    return S


def scale_from_lines(L1, L2, size_mm):
    """Scale from a line across one square (length L1 px) and a long line (L2 px) between grid lines:
    (px/mm, N whole squares in the long line, ratio L2 / L1, which should be near N)."""
    ratio = L2 / L1
    N = max(1, int(round(ratio)))
    return L2 / (N * size_mm), N, ratio


def outliers(t, ok):
    """Positions more than 3 robust SDs (1.4826 x MAD) from the median."""
    out = np.zeros(len(t), bool)
    v = t[ok]
    if len(v) < 5:
        return out
    med = np.median(v)
    mad = 1.4826 * np.median(np.abs(v - med))
    with np.errstate(invalid="ignore"):
        return ok & (np.abs(t - med) > 3 * max(mad, 1e-6))


# ------------------------------------------------------------------ internals
def movmean(x, w):
    """Moving mean over w points, shrinking at the ends (MATLAB's smoothdata 'movmean')."""
    kb = w // 2
    kf = w - 1 - kb
    c = np.r_[0.0, np.cumsum(x)]
    i = np.arange(len(x))
    lo, hi = np.maximum(i - kb, 0), np.minimum(i + kf + 1, len(x))
    return (c[hi] - c[lo]) / (hi - lo)


def curve(P, step, pxmm, smooth_px=None):
    """Smooth curve through points P (n x 2), sampled every `step` px along its length and smoothed over
    smooth_px (default 1 mm), with unit normals and arc length s."""
    smooth_px = pxmm if smooth_px is None else smooth_px
    P = np.asarray(P, float)[:, :2]
    keep = np.r_[True, np.any(np.diff(P, axis=0) != 0, axis=1)]
    P = P[keep]
    s = np.r_[0.0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    si = np.arange(0, s[-1] + 1e-9, step)
    if len(P) >= 3:
        x, y = PchipInterpolator(s, P[:, 0])(si), PchipInterpolator(s, P[:, 1])(si)
    else:
        x, y = np.interp(si, s, P[:, 0]), np.interp(si, s, P[:, 1])
    w = max(3, int(round(smooth_px / step)))
    x, y = movmean(x, w), movmean(y, w)
    tx, ty = np.gradient(x), np.gradient(y)
    L = np.hypot(tx, ty)
    return dict(x=x, y=y, nx=-ty / L, ny=tx / L, s=si)


def edges(lab, x0, y0, C, R, pxmm):
    """Edge on each side of every position of curve C (see the module help)."""
    k = len(C["x"])
    r = np.arange(-R, R + 1e-9, 0.5)                          # along the normal (px)
    off = np.array([-0.05, 0, 0.05]) * pxmm                   # averaged across 0.1 mm along the strip
    prof = np.zeros((k, len(r), 3))
    for j in range(3):
        acc = 0
        for o in off:
            X = C["x"][:, None] + o * (-C["ny"][:, None]) + r[None, :] * C["nx"][:, None] - x0
            Y = C["y"][:, None] + o * C["nx"][:, None] + r[None, :] * C["ny"][:, None] - y0
            acc = acc + map_coordinates(lab[:, :, j], [Y.ravel(), X.ravel()], order=1, mode="constant",
                                        cval=np.nan).reshape(X.shape)
        prof[:, :, j] = acc / len(off)
    mid = np.abs(r) <= 0.15 * pxmm
    tis = np.nanmedian(prof[:, mid, :].reshape(-1, 3), axis=0)          # tissue color
    far = np.abs(r) >= R - 0.5 * pxmm
    bg_all = np.nanmedian(prof[:, far, :].reshape(-1, 3), axis=0)       # the strip's background
    C_all = np.linalg.norm(tis - bg_all)
    S = dict(C, rp=np.full(k, np.nan), rm=np.full(k, np.nan))
    for i in range(k):
        p = prof[i]
        for side in (1, -1):
            sel = side * r >= 0
            o = np.argsort(side * r[sel], kind="stable")
            rr, pp = (side * r[sel])[o], p[sel][o]
            # background: at the end of this side, unless that looks like tissue (the other arm of a tightly
            # curled strip), then the strip's; then refined next to the edge found, likewise
            with np.errstate(all="ignore"):
                bg = np.nanmedian(p[far & (side * r > 0)], axis=0)
            if not (np.linalg.norm(bg - tis) >= 0.5 * C_all):
                bg = bg_all
            e = crossing(rr, pp, tis, bg)
            if np.isfinite(e):
                nb = (rr >= e + 0.3 * pxmm) & (rr <= e + 1.0 * pxmm)
                with np.errstate(all="ignore"):
                    b2 = np.nanmedian(pp[nb], axis=0)
                if nb.sum() >= 3 and np.linalg.norm(b2 - tis) >= 0.5 * C_all:
                    e = crossing(rr, pp, tis, b2)
            if side > 0:
                S["rp"][i] = e
            else:
                S["rm"][i] = e
    S["ok"] = np.isfinite(S["rp"]) & np.isfinite(S["rm"])
    return S


def crossing(rr, pp, tis, bg):
    """First distance where the tissue score drops below 1/2 (NaN if tissue and background hardly differ, the
    start is not tissue, or there is no drop). As in measure_thickness.py, a pixel brighter than the
    background scores its color distance from the background over the tissue's; a darker one scores its
    warmth (b*) over the tissue's (darker layers, the strip's shaded side), so dark shadows on the mat do not
    count."""
    C = np.linalg.norm(tis - bg)
    if not np.all(np.isfinite(bg)) or C < 6:                  # Lab difference too small to tell
        return np.nan
    dE = np.sqrt(((pp - bg) ** 2).sum(1))
    warm = (pp[:, 2] - bg[2]) / max(tis[2] - bg[2], 10)
    f = dE / C
    dark = pp[:, 0] <= bg[0]
    f[dark] = warm[dark]
    if not (f[0] > 0.5):
        return np.nan
    j = np.nonzero(f < 0.5)[0]
    if not len(j):
        return np.nan
    j = j[0]
    return rr[j - 1] + (0.5 - f[j - 1]) / (f[j] - f[j - 1]) * (rr[j] - rr[j - 1])


def lab_crop(rgb, P, margin):
    """CIE Lab (D65) of the part of rgb around the points P, with `margin` px around them, and its corner."""
    h, w = rgb.shape[:2]
    x0, x1 = max(0, int(np.floor(P[:, 0].min() - margin))), min(w - 1, int(np.ceil(P[:, 0].max() + margin)))
    y0, y1 = max(0, int(np.floor(P[:, 1].min() - margin))), min(h - 1, int(np.ceil(P[:, 1].max() + margin)))
    return srgb_to_lab(rgb[y0:y1 + 1, x0:x1 + 1]), x0, y0


def srgb_to_lab(rgb):
    """sRGB (uint8) to CIE Lab, D65 white."""
    v = rgb.astype(np.float32) / 255
    lin = np.where(v > 0.04045, ((v + 0.055) / 1.055) ** 2.4, v / 12.92)
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]], np.float32)
    XYZ = lin @ M.T / np.array([0.95047, 1, 1.08883], np.float32)
    f = np.where(XYZ <= 0.008856, 7.787 * XYZ + 16 / 116, np.cbrt(XYZ))
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]), 200 * (f[..., 1] - f[..., 2])], -1)


def crossing_chords(S):
    """Positions whose chord (edge to edge) crosses the chord of another position nearby."""
    k = len(S["x"])
    c = np.zeros(k, bool)
    A = np.c_[S["x"] - S["rm"] * S["nx"], S["y"] - S["rm"] * S["ny"]]
    B = np.c_[S["x"] + S["rp"] * S["nx"], S["y"] + S["rp"] * S["ny"]]
    ok = np.nonzero(S["ok"])[0]
    for i in ok:
        reach = (S["rp"][i] + S["rm"][i]) * 2                 # neighbors within two thicknesses
        for j in ok[(np.abs(S["s"][ok] - S["s"][i]) <= reach) & (ok != i)]:
            if segments_cross(A[i], B[i], A[j], B[j]):
                c[i] = True
                break
    return c


def segments_cross(p1, p2, q1, q2):
    def o(a, b, c):
        return np.sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
    return o(p1, p2, q1) * o(p1, p2, q2) < 0 and o(q1, q2, p1) * o(q1, q2, p2) < 0


def resample(P, n):
    """n points evenly spaced along polyline P."""
    s = np.r_[0.0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    s, u = np.unique(s, return_index=True)
    P = P[u]
    t = np.linspace(0, s[-1], n)
    return np.c_[np.interp(t, s, P[:, 0]), np.interp(t, s, P[:, 1])]


def hit_polyline(p, d, Pl):
    """Signed distances along direction d from p to every crossing of polyline Pl."""
    Pl = np.asarray(Pl, float)
    t = []
    for a, b in zip(Pl[:-1], Pl[1:]):
        e = b - a
        M = np.array([[d[0], -e[0]], [d[1], -e[1]]])
        if abs(np.linalg.det(M)) < 1e-9:
            continue
        q = np.linalg.solve(M, a - p)
        if 0 <= q[1] <= 1:
            t.append(q[0])
    return np.array(t)
