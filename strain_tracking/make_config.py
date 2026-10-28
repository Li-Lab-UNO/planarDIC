"""
Draft a test's test_config.json from its Def11 rest frame: the 20 tine tips, the tine-tip lines, the ROIs,
the sample label and a glint check.

Rods: each rake's five rods are found as long straight dark lines in the background band between the
image edge and the tissue. In every column (left/right rakes) or row (top/bottom rakes) of that band
the profile across the rods is searched for narrow dark dips; dips that line up across the band are
fitted with a straight line per rod, so speckle and bubbles cannot pass for a rod.
Tips: each rod is followed inward from the outer end of the band along its fitted line (re-centered on
its dark core, at most 6 px off the line). The rod's center (the darkest point within 6 px of the
centerline, which skips the bright highlight on the ball tip) stays near the rod's gray level up to the
tip and rises to tissue brightness past it. The largest sustained step up locates the tip's
neighborhood; the tip is the end of the continuous rod-dark stretch there. A tip more than 20 px from
its rake's median is searched again near the median with longer gaps (long ball-tip highlight) and a
looser darkness limit (rod crossing the dark tissue margin), keeping the result nearest the median; if
it is still off, the calibration's tip finder walks in from around the median tip (a bent or shortened
tine), keeping the clear rod-to-tissue edge (contrast > 30) nearest the median.

Self-check (nothing is written if it fails, unless --allow-warnings): five rods found on every rake;
each rake's tips within 35 px of their median (the tines are the same length); a clear brightness step
at every tip. It catches gross errors; tips 10-25 px off (most often on a top rake, past the rod end)
can pass, so check the overlay.

Written: sample, roi_px, tine_tips_px (median tip per rake), calibration_tines_px, reference_deformation;
dic.mask_saturated (on when the tissue in the rest frame has saturated glints, off when it has none or is
overexposed, >= 10% saturated) unless the config sets it. widths_mm is left to measure_tines.py --write
and thickness_mm to the thickness step (thickness/measure_interactive.py); when a config exists, --update
replaces only the keys written here and keeps the rest, and a rake that fails the self-check starts from
its existing tines. Finally every tip is refined with calibrate_scale.find_tip (the tip finder the scale
calibration applies frame by frame) and tine_tips_px is the median of the refined tips per rake.

With --review (runDIC.sh: review = yes), the tips are shown and the question is whether every red circle
is on the tip of a tine; the rakes answered wrong are clicked (5 tips each) and refined the same way. A
rake that fails the self-check, with no tines in the config to start from, is clicked too.

    .venv/bin/python strain_tracking/make_config.py --test N1757/test1 --label "labeled healthy"
    .venv/bin/python strain_tracking/make_config.py --test N1739/test1 --out /tmp/check.json       (compare, no overwrite)
Output: <test>/test_config.json (or --out) and <test>/DIC_results/config_tine_tips.png
"""
import argparse
import os
import re
import sys

import cv2
import numpy as np
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import review  # noqa: E402
import testcfg  # noqa: E402

N_RODS = 5
BANDS = [(0.10, 0.27), (0.04, 0.18), (0.08, 0.20), (0.16, 0.30)]   # rod search bands, fraction of the image from its edge
SIGMAS = [4.0, 8.0]              # px: smoothing across the rods (wider for wide rods at high zoom)
MAX_SHIFT = 6                    # px the followed centerline may leave the fitted rod line
SAT_LEVEL = 240                  # gray level the tracker treats as saturated (dic2d)
SAT_MAX = 10                     # % of the tine region: above this, glint masking stays off (overexposed)
RAKE_TOL = 20                    # px a tip may differ from its rake's median before it is searched again
SETTLE_PX = 3                    # refine(): the last two find_tip passes must agree within this (px)
MAX_GAP = 12                     # px of non-rod bridged along a rod (ball-tip highlight)
RETRY = [(0.3, 20), (0.5, 12), (0.5, 20)]  # (dark, max_gap) tried for a tip off its rake's median
SPREAD_TOL = 35                  # px: self-check limit on a rake's tips around their median (tines bend)
STEP_PX = 30                     # px averaged on each side of a candidate tip
# rake -> (rods run along x?, direction from tip out to the image edge)
RAKES = {"left": (True, -1), "right": (True, 1), "top": (False, -1), "bottom": (False, 1)}


def sample_line(g, along, across, horizontal):
    """Gray values at (along, across) for a rod running along x (horizontal) or along y."""
    X, Y = (along, across) if horizontal else (across, along)
    return cv2.remap(g, np.float32(X), np.float32(Y), cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)


def find_rods(g, horizontal, out, band, sigma=4.0):
    """Straight-line fits (c = a + b * s) of the rake's rods inside the band of s (coordinate along the rods):
    the N_RODS lines with the most total dip contrast."""
    H, W = g.shape
    n_across = H if horizontal else W
    lo, hi = int(0.2 * n_across), int(0.8 * n_across)
    hits = []                                          # (s, c) of narrow dark dips
    for s in range(band[0], band[1], 6):
        prof = g[lo:hi, s - 2:s + 3].mean(1) if horizontal else g[s - 2:s + 3, lo:hi].mean(0)
        prof = cv2.GaussianBlur(prof[None, :].astype(np.float32), (0, 0), sigma)[0]   # one dip per rod (shiny
        base = cv2.blur(prof[None, :], (int(15 * sigma) | 1, 1))[0]                  # rods have a bright stripe)
        dip = base - prof
        pk, pr = find_peaks(dip, prominence=20, width=(4, 40), distance=25)
        hits += [(s, lo + p, q) for p, q in zip(pk, pr["prominences"])]
    if not hits:
        return []
    h = np.array(hits, float)
    # group dips into rods: seed on nearby c values, then refit along the line (rods tilt slightly)
    lines, used = [], np.zeros(len(h), bool)
    for c0 in np.sort(h[:, 1]):
        m = ~used & (np.abs(h[:, 1] - c0) < 10)
        if m.sum() < 0.5 * (band[1] - band[0]) / 6:
            continue
        for _ in range(3):                             # refit and regather along the fitted line
            b, a = np.polyfit(h[m, 0], h[m, 1], 1)
            m = ~used & (np.abs(h[:, 1] - (a + b * h[:, 0])) < 4)
        if m.sum() >= 0.6 * (band[1] - band[0]) / 6:
            used |= m
            b, a = np.polyfit(h[m, 0], h[m, 1], 1)
            lines.append((a, b, float(h[m, 2].sum())))       # score: total dip contrast along the line
    lines.sort(key=lambda L: -L[2])
    return sorted(lines[:N_RODS], key=lambda L: L[0] + L[1] * np.mean(band))


def rodness(g, s, c, horizontal):
    """Center darkness relative to both sides, along the rod at positions s with centerline c(s)."""
    off = np.arange(-14, 15, dtype=float)
    P = sample_line(g, np.repeat(s[:, None], len(off), 1), c[:, None] + off[None, :], horizontal)
    center = P[:, np.abs(off) <= 3].min(1)
    sides = 0.5 * (P[:, (off >= -14) & (off <= -9)].mean(1) + P[:, (off >= 9) & (off <= 14)].mean(1))
    return sides - center


def find_tip(g, line, horizontal, out, s_start, s_stop, dark=0.3, shape=False, max_gap=None, last_step=False,
             levels=None):
    """Follow the rod inward from s_start (outer end of the band) toward s_stop. The centerline follows
    the rod's dark core but stays within MAX_SHIFT px of the fitted line. The rod's center brightness
    (darkest point within 6 px of the centerline: the rod's dark edges stay dark under its bright center
    stripe) stays near the rod's gray level up to the tip and rises to tissue brightness past it. The tip
    is the largest sustained step up: mean brightness over the next STEP_PX px minus that over the previous
    STEP_PX px, where the previous stretch is still rod-dark; this finds the tip's neighborhood without
    being misled by bright spots along the rod. Within 40 px of it, the tip is placed at the end of the
    continuous rod-dark stretch (gaps up to MAX_GAP px bridged): no brighter than `dark` of the way from
    the rod's gray level to the background. 0.3 keeps the dark ring of tissue around the ball end off the
    rod. A tip off its rake's median is searched again near it with a longer gap (a long highlight on the
    ball tip) or dark=0.5 (a rod crossing the dark tissue margin), keeping the result nearest the median.
    With `shape`, the rod must also be darker than its sides.
    Returns (tip, centerline at the tip, step, rod gray, background gray)."""
    a, b = line[:2]
    off = np.arange(-15, 16, dtype=float)
    core, side = np.abs(off) <= 6, np.abs(off) >= 10
    s_all = np.arange(s_start, s_stop, -out, dtype=float)
    shift, cs, cen, sid = 0.0, [], [], []
    for s in s_all:
        c = a + b * s + shift
        P = sample_line(g, np.full((1, len(off)), s), (c + off)[None, :], horizontal)[0]
        P = np.convolve(P, np.ones(3) / 3, mode="same")
        i = int(np.argmin(np.where(np.abs(off) <= 3, P, np.inf)))
        shift = float(np.clip(shift + 0.3 * off[i], -MAX_SHIFT, MAX_SHIFT))
        cs.append(c); cen.append(P[core].min()); sid.append(P[side].mean())
    cen, sid = np.array(cen), np.array(sid)
    # rod and background gray from the start of the walk, or `levels` (rod_gray, bg) from an earlier walk
    rod_gray, bg = levels or (float(np.median(cen[:60])), float(np.median(sid[:60])))
    # contrast scale: rod vs background at the start of the walk, or half of rod vs the brightest tissue
    # along it when the start lies in the image's dark corner vignette
    contrast = max(bg - rod_gray, 0.5 * (float(np.nanpercentile(cen, 90)) - rod_gray), 30.0)
    k = np.arange(STEP_PX + 5, len(cen) - STEP_PX - 5)
    before = np.array([np.median(cen[q - STEP_PX - 3:q - 3]) for q in k])
    after = np.array([np.median(cen[q + 3:q + STEP_PX + 3]) for q in k])
    score = np.where(before < rod_gray + 0.35 * contrast, after - before, -np.inf)
    step = float(np.max(score))
    # last_step: the farthest-in strong step (a bright highlight on the ball tip can outscore the rod's end)
    q = int(k[np.nonzero(score >= 0.6 * step)[0].max()] if last_step else k[np.argmax(score)])
    # place the tip: end of the continuous rod-dark stretch through the step region (gaps up to MAX_GAP px
    # bridged); the step alone lands past the tip, in the dark ring of tissue around the ball end
    on = cen < rod_gray + dark * max(contrast, step)
    if shape:                                           # narrow dark core with brighter sides, not a broad dark fold
        on &= (sid - cen) > 0.25 * contrast
    last, gap = max(q - 40, 0), 0
    for i in range(max(q - 40, 0), min(q + 40, len(cen))):
        if on[i]:
            last, gap = i, 0
        else:
            gap += 1
            if gap > (max_gap or MAX_GAP):
                break
    return float(s_all[last]), float(cs[last]), step, rod_gray, bg


def tip_from(g, c, t, horizontal, out, passes=6):
    """Tip found by calibrate_scale.find_tip started at tip coordinate t on the rod at c: each pass moves it
    up to about 15 px toward the rod end. Returns (tip, rod-to-tissue contrast)."""
    from calibrate_scale import find_tip as cal_find_tip
    cc, tt, con = float(c), float(t), 0.0
    for _ in range(passes):
        cc, tn, con = cal_find_tip(g, cc, tt, horizontal, out)
        if not np.isfinite(tn):
            return np.nan, con
        tt = tn
    return tt, con


def detect(g):
    """Rods, tips and self-check problems per rake."""
    H, W = g.shape
    rakes, problems = {}, []
    for side, (hor, out) in RAKES.items():
        n = W if hor else H
        best = None             # the first band and smoothing with five rods (later bands avoid a large specimen;
        for sig in SIGMAS:      # wider smoothing merges the two dark edges of a wide shiny rod), else the most
            for f0, f1 in BANDS:
                bd = (int(f0 * n), int(f1 * n)) if out < 0 else (int((1 - f1) * n), int((1 - f0) * n))
                ls = find_rods(g, hor, out, bd, sig)
                if best is None or len(ls) > len(best[1]):
                    best = (bd, ls)
                if len(ls) == N_RODS:
                    break
            if len(best[1]) == N_RODS:
                break
        band, lines = best
        if len(lines) != N_RODS:
            problems.append(f"{side} rake: found {len(lines)} rods, expected {N_RODS}")
            continue
        start = band[0] if out < 0 else band[1]         # outer end of the band: tips may lie inside it
        found = [find_tip(g, L, hor, out, start, n // 2) for L in lines]
        med = float(np.median([f[0] for f in found]))
        tips, redo = [], []
        for L, f in zip(lines, found):
            if abs(f[0] - med) > RAKE_TOL:              # tines are the same length: search again near the rake's tips
                alts = [find_tip(g, L, hor, out, med + out * 80, med - out * 80, dark=d, max_gap=m, last_step=ls,
                                  levels=f[3:5]) for d, m in RETRY for ls in (False, True)]              # longer ball-tip highlight, or rod crossing the dark margin
                f = min(alts, key=lambda x: abs(x[0] - med))
                redo.append(f"{f[1]:.0f}")
            t, c, step, rod_gray, bg = f
            weak = step < max(25.0, 0.4 * (bg - rod_gray))
            if abs(t - med) > RAKE_TOL or weak:         # last resort: walk the calibration's tip finder in from
                cands = [tip_from(g, c, med + dk, hor, out) for dk in (-20, -10, 0, 10, 20)]   # the rake's median
                cands = [x for x in cands if np.isfinite(x[0]) and x[1] > 30 and abs(x[0] - med) <= SPREAD_TOL]
                tn, con = min(cands, key=lambda x: abs(x[0] - med)) if cands else (np.nan, 0.0)  # tip (a bent or
                if np.isfinite(tn):                                                              # shortened tine)
                    t, weak = tn, False
                    redo.append(f"{c:.0f} (from the median)")
            if weak:
                problems.append(f"{side} rake, rod at {c:.0f} px: weak tip step ({step:.0f} gray levels)")
            tips.append([round(c), round(t)])
        if redo:
            print(f"  {side} rake: searched again near the rake's median tip for rods at {', '.join(redo)} px")
        spread = max(abs(t - np.median([u for _, u in tips])) for _, t in tips)
        if spread > SPREAD_TOL:
            problems.append(f"{side} rake: tips spread {spread:.0f} px from their median (tines are the same length)")
        rakes[side] = dict(horizontal=hor, out=out, tips=tips)
    return rakes, problems


def sample_label(T, label):
    m = re.match(r"(test\d+)_(\d\d)(\d\d)(\d{4})", T.name)
    date = f" ({m[2]}/{m[3]}/{m[4]})" if m else ""
    test = m[1] if m else T.name
    return f"{T.sample.replace('_', ' ')}, {test}{date}" + (f", {label}" if label else "")


def overlay(rgb, rakes, tips_px, path):
    for side, r in rakes.items():
        for c, t in r["tips"]:
            x, y = (t, c) if r["horizontal"] else (c, t)
            cv2.circle(rgb, (int(x), int(y)), 8, (0, 0, 255), 2)
    (x0, x1), (y0, y1) = tips_px["z"], tips_px["theta"]
    cv2.rectangle(rgb, (x0, y0), (x1, y1), (0, 255, 255), 2)
    pad = 150
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cv2.imwrite(path, rgb[max(0, y0 - pad):y1 + pad, max(0, x0 - pad):x1 + pad])


def show_tips(gray8, rakes, tips_px):
    """The tips over the rest frame, with the question and its answers in the same window: "Continue" or
    "Click <RAKE> again" (closing the window answers "Continue")."""
    fig = review.figure("runDIC: tine tips")
    ax = fig.add_axes([0.02, 0.15, 0.96, 0.78])
    ax.imshow(gray8, cmap="gray", vmin=0, vmax=255)
    for r in rakes.values():
        xy = np.array([(t, c) if r["horizontal"] else (c, t) for c, t in r["tips"]], float)
        ax.plot(xy[:, 0], xy[:, 1], "o", ms=10, mfc="none", mec="r", mew=1.5)
    (x0, x1), (y0, y1) = tips_px["z"], tips_px["theta"]
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], "y-", lw=1.5)
    pad = 150
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y1 + pad, y0 - pad)
    ax.axis("off")
    ax.set_title("Tine tips (red circles) and tine-tip lines (yellow): every circle should sit on a tine tip",
                 fontsize=11)
    fig.canvas.mpl_connect("scroll_event", lambda e: review.zoom(ax, e))
    return review.ask(fig, "Is every red circle on the tip of a tine? (The yellow rectangle joins the tips of "
                      "each rake.) If a rake's circles are off, click that rake again.",
                      ["Continue"] + [f"Click {s.upper()} again" for s in RAKES])


def click_rake(gray8, side, near=None):
    """The 5 tips of one rake, clicked on the rest frame (None if the window is closed). The view starts on
    the rake: around its tips `near` ([centerline, tip] pairs) when given, else the half of the frame it
    comes from; the scroll wheel zooms."""
    hor, out = RAKES[side]
    fig = review.figure(f"Click the 5 tine tips of the {side.upper()} rake")
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.9])
    ax.imshow(gray8, cmap="gray", vmin=0, vmax=255)
    ax.axis("off")
    H, W = gray8.shape
    if near:
        xy = np.array([(t, c) if hor else (c, t) for c, t in near], float)
        (x0, y0), (x1, y1) = xy.min(0), xy.max(0)
        pad = max(60.0, 0.5 * max(x1 - x0, y1 - y0))
        x0, x1, y0, y1 = x0 - pad, x1 + pad, y0 - pad, y1 + pad
    else:                                                  # the half of the frame the rake comes from
        x0, x1, y0, y1 = 0, W, 0, H
        if hor:
            x0, x1 = (0, W / 2) if out < 0 else (W / 2, W)
        else:
            y0, y1 = (0, H / 2) if out < 0 else (H / 2, H)
    ax.set_xlim(max(x0, -0.5), min(x1, W - 0.5))
    ax.set_ylim(min(y1, H - 0.5), max(y0, -0.5))
    pts = review.click_points(fig, ax, N_RODS, f"Click the tip of each of the 5 tines of the {side.upper()} rake "
                              "(the end of the rod, where it meets the tissue)")
    if pts is None:
        return None
    ct = sorted((y, x) if hor else (x, y) for x, y in pts)       # (centerline, tip), image px
    return dict(horizontal=hor, out=out, tips=[[round(c), round(t)] for c, t in ct])


def refine(g, rakes):
    """Refine every tip with calibrate_scale.find_tip, the tip finder the scale calibration uses frame by
    frame (four passes, as there). Only the tip coordinate is taken: find_tip's centerline is the rod's
    darkest line, which on these shiny rods is an edge, about half a rod width off the axis (as in
    measure_tines.py). A refinement is kept when its last two passes agree within SETTLE_PX (on some rods
    find_tip alternates between the two edges, and the tip with them), its rod-to-tissue contrast exceeds 30
    gray levels and it moves the tip by at most 25 px (centerline at most 10 px); otherwise the detected tip
    is kept. Clicked tips are not refined (they stay where they were clicked). Returns the refined rakes and
    the number of tips kept unrefined."""
    from calibrate_scale import find_tip as cal_find_tip
    out, kept = {}, 0
    for side, r in rakes.items():
        tips = []
        for c0, t0 in r["tips"]:
            c, t, con, ts = float(c0), float(t0), 0.0, []
            for _ in range(4):
                c, tn, con = cal_find_tip(g, c, t, r["horizontal"], r["out"])
                if np.isfinite(tn):
                    t = tn
                    ts.append(tn)
            settled = len(ts) >= 2 and abs(ts[-1] - ts[-2]) <= SETTLE_PX
            if settled and con > 30 and abs(t - t0) <= 25 and abs(c - c0) <= 10:
                tips.append([c0, round(t)])            # tip only: find_tip's centerline sits on a rod edge
            else:
                tips.append([c0, t0]); kept += 1
        out[side] = dict(r, tips=tips)
    return out, kept


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", default="", help='appended to the sample name, e.g. "labeled healthy"')
    ap.add_argument("--reference", default="Def11", help="reference deformation (rest frame read here)")
    ap.add_argument("--out", default=None, help="write here instead of the test's test_config.json")
    ap.add_argument("--allow-warnings", action="store_true", help="write even if the self-check reports problems")
    ap.add_argument("--update", action="store_true",
                    help="update an existing test_config.json: replace the keys written here, keep the rest; a rake "
                         "that fails the self-check starts from its existing tines instead")
    ap.add_argument("--review", action="store_true", help="show the tips and ask whether they are right")
    a = ap.parse_args()
    if a.review:
        review.enable()
    T = testcfg.get()
    out_path = a.out or T.config_file
    if a.out is None and os.path.exists(T.config_file) and not a.update:
        sys.exit(f"{T.config_file} exists: use --update to update it, or --out to write elsewhere")

    folder = T.folder(a.reference)
    _, paths = D.find_frames([os.path.join(T.images, folder)])
    rgb = cv2.imread(paths[0])
    gray8 = cv2.cvtColor(rgb, cv2.COLOR_BGR2GRAY)
    g = cv2.GaussianBlur(gray8.astype(np.float32), (0, 0), 1.0)
    rakes, problems = detect(g)
    failed = sorted({p.split()[0] for p in problems})
    clicked = []                                        # rakes whose tips were clicked: kept as clicked
    if problems:
        print("SELF-CHECK FAILED:\n  " + "\n  ".join(problems))
        existing = T.config.get("calibration_tines_px", {}) if a.update else {}
        unresolved = []
        for sd in failed:
            if sd in existing:                          # start the rake from its existing tines
                rakes[sd] = {k: existing[sd][k] for k in ("horizontal", "out", "tips")}
                print(f"  {sd} rake: using its existing tines as starting points")
            elif review.enabled() and (r := click_rake(gray8, sd, rakes.get(sd, {}).get("tips"))):
                rakes[sd] = r
                clicked.append(sd)
                print(f"  {sd} rake: using the tips clicked")
            else:
                unresolved.append(sd)
        if any(sd not in rakes for sd in unresolved) or (unresolved and not a.allow_warnings):
            sys.exit("nothing written; check the rest frame, run with --review to click the tips, or pass "
                     "--allow-warnings to write anyway")
    g08 = cv2.GaussianBlur(gray8.astype(np.float32), (0, 0), 0.8)
    detected = {sd: [list(t) for t in r["tips"]] for sd, r in rakes.items()}
    refined, kept = refine(g08, {sd: r for sd, r in rakes.items() if sd not in clicked})
    rakes = {sd: refined.get(sd, rakes[sd]) for sd in rakes}
    kept += N_RODS * len(clicked)
    moved = [abs(t - d[1]) for sd in rakes for (c, t), d in zip(rakes[sd]["tips"], detected[sd])]
    print(f"refined with calibrate_scale.find_tip: {20 - kept} of 20 tips (median move {np.median(moved):.0f} px, "
          f"max {max(moved):.0f} px); {kept} kept as detected or clicked")
    med = lambda side: int(round(float(np.median([t for _, t in rakes[side]["tips"]]))))  # noqa: E731
    tips_px = {"z": [med("left"), med("right")], "theta": [med("top"), med("bottom")]}
    while review.enabled() and (answer := show_tips(gray8, rakes, tips_px)) != "Continue":
        sd = answer.split()[1].lower()                     # "Click LEFT again"
        r = click_rake(gray8, sd, rakes[sd]["tips"])
        if r:
            rakes[sd] = r
            print(f"  {sd} rake: tips clicked (kept as clicked)")
        tips_px = {"z": [med("left"), med("right")], "theta": [med("top"), med("bottom")]}
    (z0, z1), (t0, t1) = tips_px["z"], tips_px["theta"]
    sat = 100 * float(np.mean(gray8[t0:t1, z0:z1] >= SAT_LEVEL))

    cfg = dict(T.config) if (a.update and a.out is None) else {}
    cfg.update({
        "sample": sample_label(T, a.label) if (a.label or "sample" not in cfg) else cfg["sample"],
        "roi_px": {"standard": [z0 + 40, t0 + 40, z1 - 40, t1 - 40], "tine_region": [z0 + 10, t0 + 10, z1 - 10, t1 - 10]},
        "tine_tips_px": tips_px,
        "calibration_tines_px": rakes,
        "reference_deformation": cfg.get("reference_deformation", a.reference),
    })
    dic = dict(cfg.get("dic", {}))
    if "mask_saturated" not in dic:       # glints: mask them, unless they cover so much that masking removes the speckle
        dic["mask_saturated"] = bool(0 < sat < SAT_MAX)
        why = ("no saturated pixels" if sat == 0 else "glints masked" if sat < SAT_MAX
               else "overexposed: masking would remove too much of the speckle")
        dic.setdefault("note", f"rest frame: {sat:.3f}% of the tine region at gray >= {SAT_LEVEL} ({why})")
    cfg["dic"] = dic
    Tw = testcfg.Test(T.root)
    Tw.config_file = out_path
    Tw.write_config(cfg)
    png = os.path.join(T.results, "config_tine_tips.png")
    overlay(rgb, rakes, tips_px, png)
    print(f"tine_tips_px {tips_px}; saturated {sat:.2f}% of the tine region")
    for side, r in rakes.items():
        print(f"  {side:6s} " + "  ".join(f"[{c}, {t}]" for c, t in r["tips"]))
    print("wrote", out_path, "and", png)


if __name__ == "__main__":
    main()
