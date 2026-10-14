"""
Measure each rake's tine span in the rest frame of the test's reference deformation, and set the
tine-span widths as the mean of opposing rakes.

For every tine (starting from calibration_tines_px in test_config.json):
  * tip: refined with calibrate_scale.find_tip;
  * centerline: at 5 px steps from 15 to 150 px behind the tip, the rod center is the midpoint of
    its two edges (best dark-band-between-bright-flanks match, band width 10-26 px, so the rod's
    specular highlight does not bias it); a robust line through these is evaluated 15 px behind
    the tip, where the tines hold the tissue.
A rake's span = (outer-to-outer centerline distance) x 5/4, i.e. five tine spacings.
The tines are bent and applied unevenly, so opposing rakes differ; both carry the same force,
so the width is their mean:
    w_theta (divides F_z)     = mean(left, right)
    w_z     (divides F_theta) = mean(top, bottom)

    .venv/bin/python strain_tracking/measure_tines.py --test N1739/test2            print
    .venv/bin/python strain_tracking/measure_tines.py --test N1739/test2 --write    also update test_config.json
--write sets widths_mm.tine_span and widths_mm.tine_span_by_rake (needs scale_calibration.json).
"""
import argparse
import json
import os
import sys

import cv2
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402
from calibrate_scale import find_tip  # noqa: E402

N_TINES_SPAN = 5 / 4          # span = outer distance x n / (n - 1), n = 5 tines


def rod_center(p, across):
    """Center of a dark rod in an across-profile: midpoint of the best dark band between bright flanks."""
    cs = np.cumsum(np.r_[0.0, p])
    mean = lambda a, b: (cs[b] - cs[a]) / (b - a)  # noqa: E731
    best, pos = -np.inf, None
    for w in range(10, 27):
        for i in range(6, len(p) - w - 6):
            score = 0.5 * (mean(i - 6, i) + mean(i + w, i + w + 6)) - mean(i, i + w)
            if score > best:
                best, pos = score, i + (w - 1) / 2
    return across[0] + pos


def centerline(g, c, t, horizontal, out, at=15):
    """Rod centerline position `at` px behind the tip (robust line fit over 15-150 px)."""
    ds = np.arange(15, 151, 5)
    across = np.arange(-35, 36)
    cs = []
    for d in ds:
        a = t + out * d + out * np.arange(-2, 3)
        if horizontal:
            X, Y = np.meshgrid(a, c + across)
        else:
            Y, X = np.meshgrid(a, c + across)
        prof = cv2.remap(g, X.astype(np.float32), Y.astype(np.float32), cv2.INTER_LINEAR).mean(1)
        cs.append(c + rod_center(prof, across))
    cs = np.array(cs)
    A = np.c_[np.ones_like(ds), ds].astype(float)
    w = np.ones(len(ds))
    for _ in range(3):
        coef = np.linalg.lstsq(A * w[:, None], cs * w, rcond=None)[0]
        r = cs - A @ coef
        s = 1.4826 * np.median(np.abs(r[w > 0])) + 0.3
        w = (np.abs(r) < 3 * s).astype(float)
    return float(coef[0] + coef[1] * at), float(np.median(np.abs(r)))


def measure(T):
    folder = T.folder(T.setting("reference_deformation"))
    _, paths = D.find_frames([os.path.join(T.images, folder)])
    g = cv2.GaussianBlur(cv2.imread(paths[0], cv2.IMREAD_GRAYSCALE).astype(np.float32), (0, 0), 1.0)
    res = {}
    for side, spec in T.setting("calibration_tines_px").items():
        hz, out = spec["horizontal"], spec["out"]
        tines = []
        for c0, t in spec["tips"]:
            c, t = float(c0), float(t)
            for _ in range(4):                        # tip position (find_tip's centerline is not used)
                c2, t2, _ = find_tip(g, c, t, hz, out)
                if not np.isfinite(t2):
                    break
                c, t = c2, t2
            cl, resid = centerline(g, float(c0), t, hz, out)
            tines.append((cl, t, resid))
        res[side] = tines
    return folder, os.path.basename(paths[0]), res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true", help="update widths_mm in test_config.json")
    a = ap.parse_args()
    T = testcfg.get()
    scale = D.load_scale()
    if scale is None:
        sys.exit(f"no {T.scale_file}; run calibrate_scale.py first")
    s = scale["px_per_mm"]
    folder, frame, res = measure(T)
    print(f"{T.name}: rest frame {frame} ({folder}), {s:.2f} px/mm")
    span = {}
    for side, tines in res.items():
        cl = np.array([x[0] for x in tines])
        span[side] = round(float((cl[-1] - cl[0]) * N_TINES_SPAN / s), 2)
        print(f"  {side:6s} centerlines {np.round(cl, 1).tolist()} px, spacing {np.round(np.diff(cl) / s, 2).tolist()} mm, "
              f"span {span[side]:.2f} mm (fit residual median {np.round([x[2] for x in tines], 1).tolist()} px)")
    w = [round((span["left"] + span["right"]) / 2, 2), round((span["top"] + span["bottom"]) / 2, 2)]
    print(f"  w_theta = mean(left, right) = {w[0]:.2f} mm; w_z = mean(top, bottom) = {w[1]:.2f} mm")
    if a.write:
        cfg = json.load(open(T.config_file))
        widths = cfg.get("widths_mm", {})
        widths["tine_span"] = w
        widths["tine_span_by_rake"] = span
        cfg["widths_mm"] = widths
        T.write_config(cfg)
        print("wrote", T.config_file)


if __name__ == "__main__":
    main()
