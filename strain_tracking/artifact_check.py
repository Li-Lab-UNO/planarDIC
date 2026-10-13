"""
Find local tracking errors in a test's tracks and exclude their cause from the matching.

Some image features move differently from the tissue: a speckle particle that slides over it, or a
saturated glare patch whose edges shift with load. Subsets that overlap one follow it, which shows up
as a red/blue pair in J. Here they are found and confirmed without a person looking:

1. Candidates. At the peak of the reference protocol (Def11), each grid point's displacement is compared
   with the median of its 9 x 9 neighborhood. Clusters of at least MIN_POINTS connected points above
   max(RES_MIN_UM, RES_FACTOR x the median residual), with a peak of at least PEAK_UM, at least EDGE_MM from
   every tine-tip line (tine contact zones have steep real gradients), are candidates.
2. Confirmation. Small textured patches of the rest frame are matched into the peak frame by plain
   normalized cross-correlation, independently of the tracker. A candidate is confirmed when the tracker
   disagrees with these matches by at least CONFIRM_PX (median) inside it and at least CONFIRM_RATIO times
   more than in a ring 1-2 mm around it.
3. Fix. Each confirmed cluster's bounding box (reference-frame px, grown by half a grid step) is added to
   dic.exclude_px in test_config.json, with a note. Those pixels get zero weight in the matching and the
   neighboring subsets track on the rest (see README, Excluded regions). The test then has to be tracked
   again (run_pipeline.sh does this).

    ../.venv/bin/python artifact_check.py --test N1739/test3
    ../.venv/bin/python artifact_check.py --test N1739/test3 --report-only
    ../.venv/bin/python artifact_check.py --test N1739/test3 --track other_Def11_tine_region.npz --report-only
Exit status: 0 nothing to change, 3 test_config.json updated (track the test again).
Output: <test>/Images/<reference>/strain_results/artifact_check.png
"""
import argparse
import os
import sys

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402
from scipy.ndimage import generic_filter, label  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

RES_MIN_UM = 60.0      # residual threshold floor (um) ...
RES_FACTOR = 8.0       # ... or this many times the median residual, whichever is larger
MIN_POINTS = 6         # connected grid points in a candidate cluster
PEAK_UM = 100.0        # a candidate's largest residual (um)
EDGE_MM = 1.0          # candidates must be at least this far from every tine-tip line
CONFIRM_PX = 3.0       # tracker vs patch matching, median disagreement inside the cluster (px) ...
CONFIRM_RATIO = 3.0    # ... and this many times the disagreement in the surrounding ring
INK, MUTED = "#1f1f1e", "#6b6b67"


def peak_frame(R):
    u = R["pos"] - R["X0"][None]
    return int(np.nanargmax([np.nanmean(np.linalg.norm(u[k][R["valid"][k]], axis=1)) for k in range(len(u))]))


def residual_um(R, k, mpp):
    """|displacement - 9 x 9 neighborhood median| per grid point (um), on the grid's shape."""
    sh = tuple(R["shape"])
    u = np.where(R["valid"][k][:, None], R["pos"][k] - R["X0"], np.nan)
    res = np.zeros(sh)
    for c in range(2):
        U = u[:, c].reshape(sh)
        res += (U - generic_filter(U, np.nanmedian, size=9, mode="nearest")) ** 2
    return np.sqrt(res) * mpp * 1000


def patch_disagreement(g0, gk, R, k, pts):
    """For textured 25 px patches centered at pts (reference px): |patch match - tracked displacement| (px)."""
    v = R["valid"][k]
    u = R["pos"][k] - R["X0"]
    out = []
    for x, y in pts:
        x, y = int(round(x)), int(round(y))
        tpl = g0[y - 12:y + 13, x - 12:x + 13]
        if tpl.shape != (25, 25) or tpl.std() < 10:
            continue
        ug = griddata(R["X0"][v], u[v], (x, y), method="linear")
        if not np.all(np.isfinite(ug)):
            continue
        gx, gy, S = int(x + ug[0]), int(y + ug[1]), 45
        win = gk[gy - 12 - S:gy + 13 + S, gx - 12 - S:gx + 13 + S]
        if win.shape != (25 + 2 * S, 25 + 2 * S):
            continue
        r = cv2.matchTemplate(win, tpl, cv2.TM_CCOEFF_NORMED)
        _, mx, _, loc = cv2.minMaxLoc(r)
        if mx < 0.8:
            continue
        out.append(np.hypot(gx - S + loc[0] - (x + ug[0]), gy - S + loc[1] - (y + ug[1])))
    return np.array(out)


def check(T, R, paths, mpp):
    """Candidates and confirmations. Returns (list of dicts, residual map, z, theta)."""
    sh = tuple(R["shape"])
    tips = T.setting("tine_tips_px")
    k = peak_frame(R)
    res = residual_um(R, k, mpp)
    z, th = D.specimen_mm(R["X0"][:, 0].reshape(sh), R["X0"][:, 1].reshape(sh), tips, mpp)
    Lz, Lt = (tips["z"][1] - tips["z"][0]) * mpp, (tips["theta"][1] - tips["theta"][0]) * mpp
    edge = np.minimum.reduce([z, Lz - z, th, Lt - th])
    bg = float(np.nanmedian(res))
    thr = max(RES_MIN_UM, RES_FACTOR * bg)
    lab, n = label(np.nan_to_num(res) > thr)
    g0 = cv2.GaussianBlur(cv2.imread(paths[0], 0).astype(np.float32), (0, 0), 1)
    gk = cv2.GaussianBlur(cv2.imread(paths[k], 0).astype(np.float32), (0, 0), 1)
    X0 = R["X0"].reshape(sh + (2,))
    step = float(R["step"])
    found = []
    for i in range(1, n + 1):
        m = lab == i
        cz, ct = float(z[m].mean()), float(th[m].mean())
        c = dict(z=cz, theta=ct, points=int(m.sum()), peak_um=float(np.nanmax(res[m])),
                 edge_mm=float(edge[m].min()), confirmed=False, why="")
        if c["points"] < MIN_POINTS or c["peak_um"] < PEAK_UM:
            continue                                    # small or weak: tracking noise
        if c["edge_mm"] < EDGE_MM:
            c["why"] = f"within {EDGE_MM:g} mm of a tine line"
            found.append(c)
            continue
        xy = X0[m]
        x0, y0 = xy.min(0) - step / 2
        x1, y1 = xy.max(0) + step / 2
        inside = [(x, y) for x in np.arange(x0, x1, 15) for y in np.arange(y0, y1, 15)]
        cx, cy = xy.mean(0)
        ring = [(cx + r * np.cos(a), cy + r * np.sin(a)) for r in np.arange(1.0, 2.01, 0.25) / mpp
                for a in np.linspace(0, 2 * np.pi, 24, endpoint=False)]
        d_in, d_ring = patch_disagreement(g0, gk, R, k, inside), patch_disagreement(g0, gk, R, k, ring)
        c.update(match_in_px=float(np.median(d_in)) if d_in.size else np.nan,
                 match_ring_px=float(np.median(d_ring)) if d_ring.size else np.nan,
                 n_in=int(d_in.size), n_ring=int(d_ring.size),
                 rect=[int(np.floor(x0)), int(np.floor(y0)), int(np.ceil(x1)), int(np.ceil(y1))])
        if d_in.size < 3 or d_ring.size < 5:
            c["why"] = "too few textured patches to confirm"
        elif c["match_in_px"] >= CONFIRM_PX and c["match_in_px"] >= CONFIRM_RATIO * max(c["match_ring_px"], 0.5):
            c["confirmed"] = True
            c["why"] = "tracker disagrees with independent patch matching"
        else:
            c["why"] = "patch matching agrees with the tracker"
        found.append(c)
    return found, res, z, th, dict(k=k, bg=bg, thr=thr)


def figure(path, res, z, th, found, info, title):
    fig, ax = plt.subplots(figsize=(7.2, 7.4), constrained_layout=True)
    im = ax.pcolormesh(z, th, res, cmap="magma", vmin=0, vmax=max(2 * info["thr"], 100), shading="nearest")
    fig.colorbar(im, ax=ax, shrink=0.8, label="|displacement − neighborhood median| at peak (µm)")
    for c in found:
        col = "#2ecc71" if c["confirmed"] else "#9ec5f4"
        ax.plot(c["z"], c["theta"], "o", ms=18, mfc="none", mec=col, mew=2)
        ax.text(c["z"] + 0.35, c["theta"], f"{c['points']} pts, {c['peak_um']:.0f} µm\n{c['why']}", color=col,
                fontsize=8, va="center")
    ax.set_aspect("equal"); ax.set_xlabel("z, axial (mm)"); ax.set_ylabel("θ, circumferential (mm)")
    ax.set_title(title, fontsize=10, loc="left")
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--track", default=None, help="track .npz to check (default: the reference protocol's track)")
    ap.add_argument("--report-only", action="store_true", help="report, don't change test_config.json")
    a = ap.parse_args()
    T = testcfg.get()
    ref = T.config.get("reference_deformation", "Def11")
    folder = T.folder(ref)
    R = dict(np.load(a.track, allow_pickle=True)) if a.track else D.load_track(T.images, folder, "tine_region")[0]
    _, paths = D.find_frames([os.path.join(T.images, folder)])
    mpp = D.load_scale(T.scale_file)["mm_per_px"]
    found, res, z, th, info = check(T, R, paths, mpp)
    print(f"{ref} peak frame {info['k']}: median residual {info['bg']:.1f} um, threshold {info['thr']:.0f} um")
    for c in found:
        extra = (f"; patch matching vs tracker {c['match_in_px']:.1f} px inside (n={c['n_in']}), "
                 f"{c['match_ring_px']:.1f} px around (n={c['n_ring']})") if "match_in_px" in c else ""
        print(f"  ({c['z']:.1f}, {c['theta']:.1f}) mm: {c['points']} points, peak {c['peak_um']:.0f} um: "
              f"{'CONFIRMED' if c['confirmed'] else 'not confirmed'}, {c['why']}{extra}")
    conf = [c for c in found if c["confirmed"]]
    if not found:
        print("  no candidate clusters")
    out = os.path.join(T.images, folder, "strain_results", "artifact_check.png")
    figure(out, res, z, th, found, info, f"{T.config.get('sample', T.name)}: tracking artifact check ({ref} peak)")
    print("wrote", out)
    if not conf or a.report_only:
        return 0
    cfg = dict(T.config)
    dic = dict(cfg.get("dic", {}))
    rects = [list(r) for r in dic.get("exclude_px", [])]
    added = []
    for c in conf:
        r = c["rect"]
        overlap = any(min(r[2], q[2]) > max(r[0], q[0]) and min(r[3], q[3]) > max(r[1], q[1]) for q in rects)
        if not overlap:
            rects.append(r)
            added.append(f"{r} at ({c['z']:.1f}, {c['theta']:.1f}) mm")
    if not added:
        print("  confirmed clusters already lie in excluded rectangles: nothing to change")
        return 0
    dic["exclude_px"] = rects
    note = f"artifact_check.py excluded {'; '.join(added)} (tracker disagreed with independent patch matching)"
    dic["note"] = (dic["note"] + "; " + note) if dic.get("note") else note
    cfg["dic"] = dic
    T.write_config(cfg)
    print(f"  added to dic.exclude_px: {'; '.join(added)}\n  wrote {T.config_file}: track this test again")
    return 3


if __name__ == "__main__":
    sys.exit(main())
