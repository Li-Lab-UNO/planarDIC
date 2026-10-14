"""
Export full-field kinematics for a tracked deformation, per frame, on the
reference grid: displacement, deformation gradient, stretches and Green-Lagrange
strain, together with the rake forces.

Writes next to the tracking .npz:
    <tag>_fields.npz   (Python:  np.load(...))
    <tag>_fields.mat   (MATLAB:  load(...))
    <tag>_fields_peak.png

    .venv/bin/python strain_tracking/export_fields.py N1739_LFPA/test1_08252026_PLC25_38-20/DIC_results/Def11_Equibiaxial.npz

Axes: image x = tester X = artery axial direction z (increasing right); circumferential
direction theta increases UPWARD in the image (theta = -image y; dic2d.THETA_SIGN). Arrays are
[frame, row, col] on the reference grid; lengths in px unless suffixed _mm.

    X, Y                      reference image coordinates of grid points (px)   [row, col]
    x, y                      current image coordinates (px)
    u_z, u_theta              displacement, axial / circumferential (px; u_theta > 0 = upward)
    X_mm, Y_mm, x_mm, y_mm, u_z_mm, u_theta_mm
                              the same in mm, if scale_calibration.json exists
                              (calibrate_scale.py); px_per_mm, mm_per_px stored alongside
    valid_track               point measured in that frame (False = interpolated from neighbors)
    F_zz, F_ztheta,           deformation gradient F[i][j] = d x_i / d X_j, i, j in (z, theta)
    F_thetaz, F_thetatheta    (least-squares fit over the strain window)
    valid_field               F (and everything derived from it) available
    J                         det F = area ratio
    lambda_z, lambda_theta    stretch of material lines initially axial / circumferential
                              (sqrt C_zz, sqrt C_thetatheta)
    lambda_1, lambda_2        principal stretches (lambda_1 >= lambda_2)
    phi_1                     angle of the lambda_1 direction in the reference configuration,
                              degrees from the axial (z) direction towards +theta
    E_zz, E_thetatheta, E_ztheta
                              Green-Lagrange strain, E = (F^T F - I) / 2
    force_z_mN, force_theta_mN  rake forces (tester X / Y)
    times, frame_files, axis_convention, subset_px, step_px, strain_window
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.io import savemat  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

UZ, UT = D.DISP_KEYS
DIVERGING = D.diverging_cmap()             # white at the center value (blue below, red above)
EZZ, ETT, EZT = D.STRAIN_KEYS
LZ, LT = D.STRETCH_KEYS


def kinematics(F):
    """Stretches and strains from F (..., 2, 2); index 0 = z, 1 = theta."""
    C = np.einsum("...ki,...kj->...ij", F, F)
    czz, ctt, czt = C[..., 0, 0], C[..., 1, 1], C[..., 0, 1]
    half_tr = 0.5 * (czz + ctt)
    disc = np.sqrt(np.maximum(half_tr ** 2 - (czz * ctt - czt ** 2), 0))
    E = 0.5 * (C - np.eye(2))
    return {
        "J": np.linalg.det(F),
        LZ: np.sqrt(czz), LT: np.sqrt(ctt),
        "lambda_1": np.sqrt(half_tr + disc), "lambda_2": np.sqrt(np.maximum(half_tr - disc, 0)),
        "phi_1": np.degrees(0.5 * np.arctan2(2 * czt, czz - ctt)),
        EZZ: E[..., 0, 0], ETT: E[..., 1, 1], EZT: E[..., 0, 1],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--win", type=int, default=None,
                    help="square strain window half-width in grid points (default: the Gaussian window_mm "
                         "of analysis_settings.json)")
    ap.add_argument("--no-mat", action="store_true")
    a = ap.parse_args()

    npz = os.path.abspath(a.npz)
    testcfg.get(npz)
    out_dir, tag = os.path.dirname(npz), os.path.splitext(os.path.basename(npz))[0]
    R = np.load(npz, allow_pickle=True)
    if "branch" in R and str(R["branch"]) == "loading":
        print("track holds the loading branch only (track_protocols.py --loading-only)")
    X0, shape, times = R["X0"], tuple(R["shape"]), R["times"]
    pos, valid = R["pos"], R["valid"]
    step = float(R["step"])
    mpp = D.load_scale()["mm_per_px"]
    T = len(times)

    fkeys = [k for row in D.F_KEYS for k in row]
    out = {k: np.full((T,) + shape, np.nan, np.float32) for k in
           ["x", "y", UZ, UT, *fkeys, "J", LZ, LT, "lambda_1", "lambda_2", "phi_1", EZZ, ETT, EZT]}
    out["valid_track"] = valid.reshape((T,) + shape)
    out["valid_field"] = np.zeros((T,) + shape, bool)
    bar = D.Progress(T, "fields")
    for k in range(T):
        F, ok, wlabel = D.strain_window_F(X0, pos[k], valid[k], shape, step, mpp, a.win)
        F = D.to_specimen(F)
        out["x"][k] = pos[k, :, 0].reshape(shape)
        out["y"][k] = pos[k, :, 1].reshape(shape)
        out[UZ][k] = (pos[k, :, 0] - X0[:, 0]).reshape(shape)
        out[UT][k] = D.THETA_SIGN * (pos[k, :, 1] - X0[:, 1]).reshape(shape)
        for i in range(2):
            for j in range(2):
                out[D.F_KEYS[i][j]][k] = F[..., i, j]
        for name, val in kinematics(F).items():
            out[name][k] = val
        out["valid_field"][k] = ok
        bar.update(k + 1)
    bar.close()

    # rake forces at the frame times
    data_csv = testcfg.get().first("*Data.csv")
    Dall = pd.read_csv(data_csv, encoding="latin1").drop_duplicates("Time_S").set_index("Time_S")
    Dd = Dall.reindex(np.round(times, 1))

    out.update({
        "X": X0[:, 0].reshape(shape).astype(np.float32), "Y": X0[:, 1].reshape(shape).astype(np.float32),
        "times": times,
        D.FORCE_KEYS[0]: Dd.X_Force_mN.to_numpy(float), D.FORCE_KEYS[1]: Dd.Y_Force_mN.to_numpy(float),
        "frame_files": np.array([str(p) for p in R["paths"]], dtype=object),
        "axis_convention": D.AXIS_CONVENTION,
        "subset_px": int(round((2 * int(R["radius"]) + 1) / float(R["resample"] if "resample" in R else 1))), "step_px": step, "strain_window": wlabel,
    })
    scale = D.load_scale()
    if scale is not None:
        mpp = scale["mm_per_px"]
        for k_ in ("X", "Y", "x", "y", UZ, UT):
            out[f"{k_}_mm"] = (out[k_] * mpp).astype(np.float32)
        out.update(px_per_mm=scale["px_per_mm"], mm_per_px=mpp, scale_method=scale["method"])
        print(f"mm fields added ({scale['px_per_mm']:.2f} px/mm)")
    else:
        print("no scale_calibration.json: fields in px only (run calibrate_scale.py for mm)")
    np.savez_compressed(os.path.join(out_dir, f"{tag}_fields.npz"), **out)
    print("wrote", os.path.join(out_dir, f"{tag}_fields.npz"))
    if not a.no_mat:
        savemat(os.path.join(out_dir, f"{tag}_fields.mat"),
                {k: (np.asarray(v, dtype=object) if k == "frame_files" else v) for k, v in out.items()},
                do_compression=True)
        print("wrote", os.path.join(out_dir, f"{tag}_fields.mat"))

    # peak-frame overview on the reference configuration
    k = int(np.nanargmax(np.nanmean(out[EZZ].reshape(T, -1), 1)))
    mm = f"{UZ}_mm" in out
    unit = "mm" if mm else "px"
    # J the rakes would give on a uniform specimen, relative to the same reference frame
    lx, ly = D.rake_stretch(Dall, times, D.reference_time(R, times))
    J_rake = float(lx[k] * ly[k])
    panels = [(f"{UZ}_mm" if mm else UZ, rf"$u_z$, axial displacement ({unit})", "viridis"),
              (f"{UT}_mm" if mm else UT, rf"$u_\theta$, circumferential displacement ({unit})", "viridis"),
              ("J", f"J = det F, area ratio (white: uniform rake J = {J_rake:.3f})", DIVERGING),
              (LZ, r"$\lambda_z$, axial stretch", "viridis"), (LT, r"$\lambda_\theta$, circumferential stretch", "viridis"),
              ("phi_1", r"$\varphi_1$, major principal direction from axial (°)", DIVERGING),
              (EZZ, r"$E_{zz}$, axial Green–Lagrange strain", "viridis"),
              (ETT, r"$E_{\theta\theta}$, circumferential Green–Lagrange strain", "viridis"),
              (EZT, r"$E_{z\theta}$, shear Green–Lagrange strain", DIVERGING)]
    fig, ax = plt.subplots(3, 3, figsize=(14, 13), constrained_layout=True)
    if mm:   # z from the left tine tips, theta from the bottom tine tips, both in mm and increasing up/right
        tips = testcfg.get().setting("tine_tips_px")
        Xg, Yg = D.specimen_mm(out["X"], out["Y"], tips, mpp)
        xlab, ylab = "z, axial (mm)", r"$\theta$, circumferential (mm)"
    else:
        Xg, Yg = out["X"], D.THETA_SIGN * out["Y"]
        xlab, ylab = "axial z (image x, px)", r"circumferential $\theta$ (= −image y, px)"
    for a_, (key, label, cmap) in zip(ax.ravel(), panels):
        val = out[key][k].astype(float)
        if not key.startswith("u_"):
            val = np.where(out["valid_field"][k], val, np.nan)
        lo, hi = np.nanpercentile(val, [1, 99])
        center = {"J": J_rake, "phi_1": 0.0, EZT: 0.0}.get(key)
        if center is not None:                       # diverging: white at the center value
            half = max(abs(lo - center), abs(hi - center))
            lo, hi = center - half, center + half
        m = a_.pcolormesh(Xg, Yg, val, cmap=cmap, vmin=lo, vmax=hi, shading="auto")
        if mm:                                   # specimen domain from (0, 0) at the bottom-left tine tips
            a_.set_xlim(0, (tips["z"][1] - tips["z"][0]) * mpp)
            a_.set_ylim(0, (tips["theta"][1] - tips["theta"][0]) * mpp)
        a_.set_aspect("equal"); a_.set_title(label, fontsize=11)
        a_.set_xlabel(xlab); a_.set_ylabel(ylab)
        fig.colorbar(m, ax=a_, shrink=.8)
    fig.suptitle(f"{D.protocol_title(tag, testcfg.get().root)} · peak at t = {times[k] - times[0]:.1f} s · "
                 f"{D.settings_label(R, mpp, out['strain_window'])}", fontsize=14)
    fig.savefig(os.path.join(out_dir, f"{tag}_fields_peak.png"), dpi=90)
    print("wrote", os.path.join(out_dir, f"{tag}_fields_peak.png"))


if __name__ == "__main__":
    main()
