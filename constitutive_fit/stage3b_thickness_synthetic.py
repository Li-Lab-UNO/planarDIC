"""
Stage 3b: does thickness variation in a 3D slab look, to a membrane fit, like stiffness that changes with loading ratio?

Synthetic data from a 3D slab (solid_fe.py) of uniform material, the test's bulk fit (DIC affine F), whose thickness
varies over the grid: H(z, theta) = H exp(pattern), top face (the imaged face) flat, bottom face carrying the
variation (--faces bottom) or both faces with the mid-plane flat (--faces both: no load off the mid-plane). The
slab is loaded at the peak of every protocol Def11-Def19 by the measured boundary positions, applied
through the full thickness. Its top-surface in-plane positions, with Gaussian noise of --noise-px on every point, are
saved for stage2_alpha_field.py --data, which fits them with the membrane stiffness factor alpha(x) exactly as the
measured data. In a membrane, thickness and alpha are the same thing, so a membrane fit of membrane data recovers
log(H/H) with one map for every loading ratio (stage 2 synthetic tests); any loading-ratio dependence in the fit of
these data comes from the slab being 3D. The membrane-equivalent truth, log(H(x)/H) at the membrane Gauss points, is
saved with the positions.

Patterns (log of the thickness ratio):
    bands     the stage 2 'bands' pattern (horizontal bands, stiff at theta ~ 5 mm, and an oblique streak) times --amp;
              --amp 1 gives an interior SD of log H of about 0.11 (test1's thickness strip: CV 10%), --amp 3 about
              0.34, the size of the measured alpha maps
    random    a Gaussian random field with correlation length --corr-mm (default 2.5 mm, the strips) and interior SD
              0.10 x --amp, the same in z and theta

    .venv/bin/python constitutive_fit/stage3b_thickness_synthetic.py --test N1739/test1 --pattern bands --amp 1
    .venv/bin/python constitutive_fit/stage2_alpha_field.py --test N1739/test1 \\
        --data <test>/constitutive_fit/local/stage3b_bands_amp1_data.npz --groups --bases 2x6
Outputs: <test>/constitutive_fit/local/stage3b_<pattern>_amp<amp>[_sym]_data.npz (positions, thickness, truth, solver info)
and _thickness.png (the thickness map and the top-surface E_thetatheta of three protocols)
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import ndimage  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import constitutive as M  # noqa: E402
import dic2d as D  # noqa: E402
from solid_fe import Slab  # noqa: E402
from stage1_uniform_fe import INPUT, strain  # noqa: E402
from stage1b_best_uniform import Problem  # noqa: E402
from stage2_alpha_field import interior_gp, mesh_coords, truth  # noqa: E402


def log_thickness(pattern, amp, corr_mm, pb, seed=3):
    """log(H/H) at the nodes (ny, nx) and at the membrane Gauss points (ne_r, ne_c, 4)."""
    Xg = pb.mesh.Xg
    g, s, t = mesh_coords(pb)
    if pattern == "bands":
        z0, z1, t0, t1 = Xg[..., 0].min(), Xg[..., 0].max(), Xg[..., 1].min(), Xg[..., 1].max()
        sn = 2 * (Xg[..., 0] - z0) / (z1 - z0) - 1
        tn = 2 * (Xg[..., 1] - t0) / (t1 - t0) - 1
        return amp * truth("bands", Xg, sn, tn), amp * truth("bands", g, s, t)
    if pattern == "random":
        rng = np.random.default_rng(seed)
        h = abs(Xg[0, 1, 0] - Xg[0, 0, 0])
        sig = corr_mm / 2 / h                                   # Gaussian kernel: correlation exp(-d^2 / (4 sig^2))
        pad = int(np.ceil(3 * sig))
        wn = rng.normal(size=(Xg.shape[0] + 2 * pad, Xg.shape[1] + 2 * pad))
        fld = ndimage.gaussian_filter(wn, sig)[pad:-pad, pad:-pad]
        gin = interior_gp(pb, g)
        interp = RegularGridInterpolator((Xg[:, 0, 1], Xg[0, :, 0]), fld)
        fg = interp(np.stack([g[..., 1], g[..., 0]], -1))
        scale = 0.10 * amp / fg[gin].std()
        return (fld - fg[gin].mean()) * scale, (fg - fg[gin].mean()) * scale
    raise ValueError(pattern)


_W = {}


def _init(test, params, layers, log_nodes, log_gp, faces):
    np.seterr(over="ignore", invalid="ignore")
    testcfg.select(test)
    T = testcfg.get()
    pb = Problem(T, pd.read_csv(os.path.join(T.fit_dir, "stretch_inputs_data.csv")))
    H = T.config["thickness_mm"]["H"]
    _W.update(pb=pb, params=params, slab=Slab(pb.mesh.Xg, H * np.exp(log_nodes), layers=layers, faces=faces),
              alpha_gp=np.exp(log_gp))


def _solve(folder):
    pb, slab, P = _W["pb"], _W["slab"], _W["params"]
    t0 = time.time()
    xm, _ = pb.solve(folder, P, warm=False, alpha=_W["alpha_gp"])     # membrane with alpha = H/H: the start
    x, info = slab.solve(pb.prot[folder]["xb"], P, x0=slab.from_membrane(xm))
    info["seconds_total"] = time.time() - t0
    return folder, slab.top(x)[..., :2], info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", choices=("bands", "random"), default="bands")
    ap.add_argument("--amp", type=float, default=1.0)
    ap.add_argument("--corr-mm", type=float, default=2.5)
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--faces", choices=("bottom", "both"), default="bottom",
                    help="where the thickness change sits: bottom face (top flat) or both faces (mid-plane flat)")
    ap.add_argument("--noise-px", type=float, default=0.06)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    T = testcfg.get()
    test = os.path.relpath(T.root, testcfg.ROOT)
    np.seterr(over="ignore", invalid="ignore")
    os.makedirs(T.local_dir, exist_ok=True)
    stem = os.path.join(T.local_dir, f"stage3b_{a.pattern}_amp{a.amp:g}" + ("_sym" if a.faces == "both" else ""))
    df = pd.read_csv(os.path.join(T.fit_dir, "stretch_inputs_data.csv"))
    params = pd.read_csv(os.path.join(T.fit_dir, "params_by_stretch_input.csv")).set_index("input") \
        .loc[INPUT, list(M.PARAM_NAMES)].to_numpy(float)
    pb = Problem(T, df)
    folders = list(pb.prot)
    log_nodes, log_gp = log_thickness(a.pattern, a.amp, a.corr_mm, pb)
    g, _, _ = mesh_coords(pb)
    gin = interior_gp(pb, g)
    H = T.config["thickness_mm"]["H"]
    print(f"{T.name}: 3D slab, H = {H:.2f} mm x exp({a.pattern}, amp {a.amp:g}); interior SD of log H "
          f"{log_gp[gin].std():.3f}, thickness {H * np.exp(log_nodes.min()):.2f}-{H * np.exp(log_nodes.max()):.2f} mm; "
          f"{a.layers} layers; thickness change on {'both faces' if a.faces == 'both' else 'the bottom face'}; "
          f"{len(folders)} protocol peaks on {a.workers} workers", flush=True)
    t0 = time.time()
    rng = np.random.default_rng(11)
    pos, infos = {}, {}
    with ProcessPoolExecutor(min(a.workers, len(folders)), initializer=_init,
                             initargs=(test, params, a.layers, log_nodes, log_gp, a.faces)) as pool:
        for folder, top, info in pool.map(_solve, folders):
            pos[folder] = pb.to_px(top) + rng.normal(0, a.noise_px, (len(pb.X0), 2))
            infos[folder] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in info.items()}
            print(f"  {folder}: {'ok' if info['success'] else 'FAILED'} in {info['seconds_total']:.0f} s "
                  f"({info['newton']} Newton iterations) | {time.time() - t0:.0f} s elapsed", flush=True)
    np.savez_compressed(stem + "_data.npz", truth_log_alpha=log_gp, log_thickness_nodes=log_nodes, H_mm=H,
                        noise_px=a.noise_px, layers=a.layers, info=json.dumps(infos),
                        **{f"pos_{f}": v for f, v in pos.items()})

    # figure: thickness map and top-surface E_thetatheta (noisy synthetic) of three protocols
    show = [f for f in folders if f[:5] in ("Def11", "Def14", "Def18")]
    fig, axes = plt.subplots(1, 1 + len(show), figsize=(4.4 * (1 + len(show)), 4.4), constrained_layout=True)
    Xg = pb.mesh.Xg
    pm = axes[0].pcolormesh(D.cell_corners(Xg[..., 0]), D.cell_corners(Xg[..., 1]), H * np.exp(log_nodes),
                            cmap="viridis", shading="flat", edgecolors="face")
    fig.colorbar(pm, ax=axes[0], shrink=0.8, label="reference thickness (mm)")
    axes[0].set_aspect("equal"); axes[0].set_title(f"thickness ({a.pattern}, amp {a.amp:g})", fontsize=9)
    Es = {f: strain(pb.X0, pos[f], np.ones(len(pb.X0), bool), pb.shape, pb.window_mm, pb.mpp)[0][..., 1, 1]
          for f in show}
    for ax, f in zip(axes[1:], show):
        z, th = D.specimen_mm(pos[f][:, 0], pos[f][:, 1], pb.tips, pb.mpp)
        Z, Tt = z.reshape(pb.shape), th.reshape(pb.shape)
        v = Es[f][pb.interior]
        pm = ax.pcolormesh(D.cell_corners(Z), D.cell_corners(Tt), Es[f], cmap="viridis",
                           vmin=np.nanpercentile(v, 1), vmax=np.nanpercentile(v, 99), shading="flat", edgecolors="face")
        ax.contour(Z, Tt, pb.interior.astype(float), levels=[0.5], colors="w", linewidths=0.7, linestyles="--")
        fig.colorbar(pm, ax=ax, shrink=0.8)
        ax.set_aspect("equal"); ax.set_title(f"{f}: E_thetatheta, 3D top surface", fontsize=9)
    for ax in axes:
        ax.set_xlabel("z (mm)", fontsize=8); ax.set_ylabel("θ (mm)", fontsize=8); ax.tick_params(labelsize=7)
    fig.suptitle(f"{T.name}: stage 3b synthetic data, uniform material, thickness varying on "
                 f"{'both faces (mid-plane flat)' if a.faces == 'both' else 'the bottom face (top flat)'}", fontsize=10)
    fig.savefig(stem + "_thickness.png", dpi=110)
    print(f"wrote {stem}_data.npz, _thickness.png ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
