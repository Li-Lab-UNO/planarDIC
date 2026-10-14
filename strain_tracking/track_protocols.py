"""
Track several deformations against ONE common reference configuration.

The reference is the first frame of the test's reference_deformation (Def11). For each
deformation, its first frame is matched directly to the reference frame (IC-GN, rest
state to rest state, so the "bridge" displacement is small). The deformation is then
tracked from those positions. Results are saved in the track.py format, with X0 = the common
reference grid, so export_fields.py / compare_cellscale.py / animate_fields.py give F relative
to the common reference.

Images are resampled (area averaging) to target_px_per_mm of analysis_settings.json before
matching (downsampling only), and subset size and grid step (subset_mm, step_mm) are converted to
px at that scale, so every test is tracked with the same px and mm settings. Tracks are saved in
native image px; radius is in the resampled px, and the npz records resample and track_px_per_mm.
--radius / --step (resampled px) override them.

Grid (--grid):
  tine_region (default)  the full domain: a lattice of the grid step centered between the rest-state
                         tine tips (tine_tips_px) -> DefNN_*.npz
  standard               the full domain without its outer two rows and columns -> DefNN_*.npz.
                         Usually not needed: dic2d.load_track(..., "standard") slices it from the
                         full-domain track.

    ../.venv/bin/python track_protocols.py --test N1739/test1
    ../.venv/bin/python track_protocols.py --test N1690/test1 --defs Def11
    ../.venv/bin/python track_protocols.py --test N1739/test1 --loading-only

--loading-only stops each deformation at the end of loading, the frame where the rakes peak
(dic2d.loading_end_time; the fits end the loading branch there too), about half the frames. Tracking runs forward, so the loading frames are the
same as in a full-cycle track; the unloading frames are left out (animations and CellScale comparisons
then cover loading only). The npz records branch = "loading" or "full".
Outputs: <test>/DIC_results/DefNN_*[_standard].npz
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

DEFAULT_DEFS = [f"Def{i}" for i in range(11, 20)]


def run(T, defs=DEFAULT_DEFS, grid="tine_region", roi=None, step=None, radius=None, mask_saturated=None,
        out_dir=None, loading_only=False):
    """Track deformations `defs` of test T on `grid` against the common reference.
    step / radius (px) default to analysis_settings.json converted with the test's scale.
    out_dir: write the .npz files there instead of the test's DIC_results/ (for trials)."""
    scale = D.load_scale(T.scale_file)
    if scale is None:
        sys.exit(f"{T.scale_file} missing: run calibrate_scale.py first")
    f = D.resample_factor(scale["px_per_mm"])          # images resampled to target_px_per_mm before matching
    ppm = scale["px_per_mm"] * f                       # tracking scale (px/mm of the resampled images)
    S = D.settings_px(ppm)
    step = S["step"] if step is None else step         # resampled px
    radius = S["radius"] if radius is None else radius
    if roi is None:                                    # in resampled px
        tips = {k: list(D.to_resampled(v, f)) for k, v in T.setting("tine_tips_px").items()}
        roi = D.tip_lattice_roi(tips, step)
        if grid == "standard":
            roi = [roi[0] + 2 * step, roi[1] + 2 * step, roi[2] - 2 * step, roi[3] - 2 * step]
    else:
        roi = list(D.to_resampled(roi, f))
    roi = list(roi)
    if mask_saturated is None:
        mask_saturated = T.config.get("dic", {}).get("mask_saturated", False)
    folders = [T.folder(d) for d in defs]
    X0ref, shape = D.make_grid(roi, step)
    excl = [tuple(r) for r in T.config.get("dic", {}).get("exclude_px", [])]   # reference frame, native px
    excl_r = [tuple(D.to_resampled(r, f)) for r in excl]
    prm = D.TrackParams(radius=radius, mask_saturated=mask_saturated, resample=f)
    ref_fd = T.folder(T.setting("reference_deformation"))
    t_ref, p_ref = D.find_frames([os.path.join(T.images, ref_fd)])
    ref = D.Reference(X0ref, prm.radius).build(*D.load_frame(p_ref[0], prm, excl_r))
    print(f"reference: {ref_fd} t = {t_ref[0]:.1f} s; {grid} grid {shape[1]} x {shape[0]}, ROI {roi}, "
          f"subset {2 * radius + 1} px" + ("; saturated glints masked" if prm.mask_saturated else "")
          + (f"; excluded {excl}" if excl else "")
          + (f"; resampled x{f:.3f} to {ppm:.1f} px/mm" if f != 1 else "")
          + f"; {(2 * radius + 1) / ppm:.2f} mm subset, {step / ppm:.3f} mm step",
          flush=True)
    t_all = time.time()
    for i, fd in enumerate(folders, 1):
        times, paths = D.find_frames([os.path.join(T.images, fd)])
        n_all, branch = len(paths), "full"
        if loading_only:
            keep = times <= D.loading_end_time(T.root, fd) + 0.05      # through the rake-peak frame
            times, paths, branch = times[keep], [q for q, k in zip(paths, keep) if k], "loading"
        # bridge: reference frame -> this deformation's first frame
        img0, use0 = D.load_frame(paths[0], prm)
        p, c, conv = D.icgn(ref, img0, np.zeros((len(X0ref), 6)), use=use0 if prm.sat_current else None)
        good = conv & (c >= prm.zncc_min)
        p = D._fill_from_neighbors(X0ref, np.nan_to_num(p), good)
        X_start = X0ref + p[:, [0, 3]]
        u = X_start - X0ref
        prm.exclude_px = tuple(D.shift_rects(excl_r, X0ref, X_start, good))   # into this protocol's first frame
        print(f"[{i}/{len(folders)}] {fd}: {len(paths)} frames"
              + (f" (loading branch of {n_all}, to t = {times[-1]:.1f} s)" if loading_only else "") + f"; bridge valid {good.mean()*100:.1f}%, median ZNCC {np.nanmedian(c):.4f}, "
              f"rest-state offset median {np.median(np.hypot(*u.T)):.2f} px", flush=True)
        t = time.time()
        res = D.track_sequence(paths, X_start, prm, times=times, label=fd.split("_")[0])
        res["valid"] &= good[None, :]
        out = os.path.join(out_dir or T.results, fd + D.GRID_SUFFIX[grid] + ".npz")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        # positions are saved in native image px (so every later script, the configs and the photos agree);
        # radius stays in the resampled px it was matched with
        res["pos"] = D.to_native(res["pos"], f)
        np.savez_compressed(out, X0=D.to_native(X0ref, f), shape=np.array(shape), times=times,
                            paths=np.array([os.path.basename(q) for q in paths]),
                            radius=radius, step=float(step / f), roi=D.to_native(roi, f), grid=grid,
                            mask_saturated=prm.mask_saturated, sat_level=prm.sat_level, sat_dilate=prm.sat_dilate,
                            exclude_px=np.array(excl, float).reshape(-1, 4),
                            reference=f"{ref_fd} t={t_ref[0]:.1f}s", X_start=D.to_native(X_start, f), bridge_zncc=c,
                            px_per_mm=scale["px_per_mm"], resample=f, track_px_per_mm=ppm,
                            subset_mm=(2 * radius + 1) / ppm, step_mm=step / ppm, branch=branch, **res)
        print(f"   tracked in {time.time()-t:.0f} s; valid per frame median {np.median(res['valid'].mean(1))*100:.1f}%, "
              f"min {res['valid'].mean(1).min()*100:.1f}%; saved {out}", flush=True)
    print(f"all {len(folders)} deformations tracked in {D._hms(time.time() - t_all)}", flush=True)


def main():
    T = testcfg.get()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--defs", nargs="+", default=DEFAULT_DEFS)
    ap.add_argument("--grid", choices=tuple(D.GRID_SUFFIX), default="tine_region")
    ap.add_argument("--roi", nargs=4, type=float, default=None, help="override the grid's ROI (px)")
    ap.add_argument("--step", type=int, default=None, help="grid step (px); default: step_mm of analysis_settings.json")
    ap.add_argument("--radius", type=int, default=None,
                    help="subset half-width (px); default: subset_mm of analysis_settings.json")
    ap.add_argument("--mask-saturated", action=argparse.BooleanOptionalAction,
                    default=T.config.get("dic", {}).get("mask_saturated", False),
                    help="leave saturated glint pixels out of the matching (default: the test's dic.mask_saturated)")
    ap.add_argument("--loading-only", action="store_true",
                    help="track each deformation to the end of loading only (the rake peak)")
    a = ap.parse_args()
    run(T, a.defs, a.grid, a.roi, a.step, a.radius, a.mask_saturated, loading_only=a.loading_only)


if __name__ == "__main__":
    main()
