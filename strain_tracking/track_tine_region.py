"""
Full-domain tracking (the default grid): the same as `track_protocols.py` with --grid tine_region.

The grid is the test's roi_px.tine_region (test_config.json), on the 15 px lattice and extended
to just inside the rest-state tine tips (tine_tips_px). For N1739 test1: ROI [760, 755, 1285, 1365],
ending 5-15 px inside the tips; the standard 37 x 32 grid is rows 2:-2, cols 2:-2 of this one.
Other scripts import ROI, TINE_TIPS and SUFFIX from here.

    ../.venv/bin/python track_tine_region.py --test N1739/test1
Outputs: <test>/DIC_results/DefNN_*.npz
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402
import track_protocols as TP  # noqa: E402

T = testcfg.get()
IMAGES = T.images
ROI = T.setting("roi_px", "tine_region")
TINE_TIPS = {k: tuple(v) for k, v in T.setting("tine_tips_px").items()}   # rest-state tip lines (px)
SUFFIX = D.GRID_SUFFIX["tine_region"]
folder = T.folder


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--defs", nargs="+", default=TP.DEFAULT_DEFS)
    ap.add_argument("--mask-saturated", action=argparse.BooleanOptionalAction,
                    default=T.config.get("dic", {}).get("mask_saturated", False),
                    help="leave saturated glint pixels out of the matching (default: the test's dic.mask_saturated)")
    a = ap.parse_args()
    TP.run(T, a.defs, "tine_region", mask_saturated=a.mask_saturated)


if __name__ == "__main__":
    main()
