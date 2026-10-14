"""
Synthetic test of candidate 1: the stiffening rate (exponents) varies by patch while the single-factor Stage 2 fit
(ladder (a)) assumes uniform exponents. Question: how much virtual-field misfit does that leave, and does the
misfit grow with load the way the real-data misfit does (VF 0.048 -> 0.099 from the lowest to the highest load
quarter)?

Same forward model and pipeline as stage2_synthetic.py (FE with measured smoothed boundary displacements, outer
nodes from the smoothed tracks, 45 px window, 8 x 8 virtual-field modes, all blocks together, neighbor penalty
with lambda by leave-one-protocol-out CV). The exponents enter the FE through per-Gauss-point parameter arrays;
the fit uses the bulk parameters. Pattern p = Stage 1 interior log alpha scaled to [-1, 1] (ring patches p = 0).
True cases (C1 constants uniform, alpha = 1):
  fib20   C2_smc, C2_col x (1 + 0.2 p)
  fib40   C2_smc, C2_col x (1 + 0.4 p)
  all20   C2_el, C2_smc, C2_col x (1 + 0.2 p)
The single-factor truth from stage2_synthetic.py (alpha = Stage 1 pattern, fe_stage1.npz) is scored the same way
as the reference.

    .venv/bin/python constitutive_fit/stage2_synthetic_c.py --test N1739/test1 [--noise 0 0.03]
Outputs: <test>/constitutive_fit/local/stage2_synthetic/ (fe_exp_*.npz, stage2_synthetic_c.*)
"""
import argparse
import json
import os
import sys
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import baseline as B  # noqa: E402
import constitutive as M  # noqa: E402
from membrane_fe import Mesh  # noqa: E402
from stage2_setup import LOCAL, assemble, geometry, measured_positions, window_F  # noqa: E402
from stage2_synthetic import difference_matrix, neighbor_pairs  # noqa: E402
from stage2_synthetic_b import cv_fit, scaled_blocks  # noqa: E402

LOADS = ((0, .25), (.25, .5), (.5, .75), (.75, 1.01))
REAL_FIT = {"0-25%": 0.048, "25-50%": 0.054, "50-75%": 0.074, "75-100%": 0.099, "all": 0.084}


def gp_params(P, E, cols, gi, gj):
    """Parameter vector with the exponents in `cols` scaled by the patch field E at every Gauss point."""
    Pg = np.empty(len(P), object)
    Pg[:] = list(P)
    for c in cols:
        Pg[c] = P[c] * E[gi, gj]
    return Pg


def rel_vf(Gv2, x, m=None):
    r, s = Gv2 @ x, np.abs(Gv2) @ np.abs(x)
    if m is not None:
        r, s = r[m], s[m]
    return float(np.sqrt(np.mean(r ** 2)) / np.sqrt(np.mean(s ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--every", type=int, default=8)
    ap.add_argument("--noise", type=float, nargs="+", default=[0.0, 0.03])
    ap.add_argument("--win", type=int, default=1)
    ap.add_argument("--modes", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(LOCAL, "stage2_synthetic"))
    a = ap.parse_args()
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    d = np.load(os.path.join(LOCAL, "local_data.npz"))
    base, df = B.load(str(d["baseline"]))
    df = df.reset_index(drop=True)
    st = base["settings"]
    H, w_t, w_z = st["H_mm"], st["w_theta_mm"], st["w_z_mm"]
    P = np.array([base["inputs"]["dic_affine"]["params"][n] for n in M.PARAM_NAMES])
    mpp = float(d["mm_per_px"])
    g = geometry(d, 1, a.modes)
    inner, npat, nt, nz, pz, pt = g["inner"], g["npat"], g["nt"], g["nz"], g["pz"], g["pt"]
    Xg = np.stack([d["z_mm"], d["theta_mm"]], -1)
    m_ = 2
    sl = (slice(m_, -m_), slice(m_, -m_))
    mesh = Mesh(Xg[sl])
    gp = mesh.gauss_points
    gi = np.clip(np.searchsorted(pt, gp[..., 1]) - 1, 0, nt - 1)
    gj = np.clip(np.searchsorted(pz, gp[..., 0]) - 1, 0, nz - 1)
    frames = np.arange(0, len(df), a.every)
    sets = df.set.to_numpy()[frames]
    Lt_in, Lz_in = (mesh.ny - 1) * abs(mesh.hy), (mesh.nx - 1) * abs(mesh.hx)
    xm = measured_positions(df, frames, mpp)
    xsm = Xg + np.stack([gaussian_filter(xm[..., i] - Xg[..., i], (0, 1.0, 1.0), mode="nearest") for i in range(2)], -1)

    s1 = pd.read_csv(os.path.join(LOCAL, "stage1_interior1", "stage1_patches.csv"))
    A1 = np.ones((nt, nz))
    A1[s1.i_theta, s1.j_z] = s1.alpha
    p = np.zeros((nt, nz))
    p[inner] = np.log(A1[inner]) / np.abs(np.log(A1[inner])).max()
    cases = {"fib20": ((4, 6), 0.2), "fib40": ((4, 6), 0.4), "all20": ((2, 4, 6), 0.2)}

    R1 = difference_matrix(neighbor_pairs(inner), npat)
    lams = np.r_[0.0, np.logspace(-6, 4, 41)]
    ones_gp = np.ones(gp.shape[:3])
    sims = {}
    for name, (cols, amp) in cases.items():
        E = 1 + amp * p
        cache = os.path.join(a.out, f"fe_exp_{name}.npz")
        if os.path.exists(cache) and np.array_equal(np.load(cache)["frames"], frames):
            C = np.load(cache)
            X, Fz, Ft = C["x"], C["Fz"], C["Ft"]
            print(f"[{name}] FE loaded", flush=True)
        else:
            Pg = gp_params(P, E, cols, gi, gj)
            X, Fz, Ft, conv = [], [], [], []
            t0 = time.time()
            for n in range(len(frames)):
                xb = xsm[n][sl]
                x, info = mesh.solve(xb, ones_gp, Pg, H, xb)
                r = mesh.reactions(x, ones_gp, Pg, H)
                X.append(x); conv.append(info["max_free_force"])
                Fz.append(0.5 * (r[:, -1, 0].sum() - r[:, 0, 0].sum()) * w_t / Lt_in)
                Ft.append(0.5 * (r[-1, :, 1].sum() - r[0, :, 1].sum()) * w_z / Lz_in)
                if n % 60 == 0:
                    print(f"  [{name}] frame {n}/{len(frames)}, {time.time() - t0:.0f} s", flush=True)
            X, Fz, Ft, conv = map(np.array, (X, Fz, Ft, conv))
            np.savez_compressed(cache, frames=frames, x=X, Fz=Fz, Ft=Ft, max_free_force=conv, exp_factor=E,
                                exp_cols=np.array(cols), inset=m_)
            print(f"[{name}] FE done in {time.time() - t0:.0f} s; max free force {conv.max():.1e} mN", flush=True)
        sims[name] = (X, Fz, Ft, E, cols)
    C = np.load(os.path.join(a.out, "fe_stage1.npz"))
    sims["alpha_only"] = (C["x"], C["Fz"], C["Ft"], None, ())

    summary, figdata = {}, {}
    for name, (X, Fz, Ft, E, cols) in sims.items():
        load = np.zeros(len(frames))
        for s in np.unique(sets):
            m = sets == s
            load[m] = np.maximum(Fz[m] / Fz[m].max(), Ft[m] / Ft[m].max())
        for sig in a.noise:
            rng = np.random.default_rng(1)
            xs = xsm.copy()
            xs[:, m_:-m_, m_:-m_] = X
            xs = xs + rng.normal(0, sig * mpp, xs.shape)
            Fg = window_F(Xg, xs, np.ones(xs.shape[:3], bool), a.win)
            Gz, Gt, Gv, bz, bt = assemble(g, Fg, Fz, Ft, P, H, w_t, w_z)
            _, _, _, sysa = scaled_blocks(((Gz, bz), (Gt, bt), (Gv, np.zeros(Gv.shape[:2]))), sets, npat)
            xa, lam = cv_fit(*sysa, R1, lams)
            Gv2 = Gv.reshape(-1, npat)
            nvf = Gv.shape[1]
            e = dict(lam=lam, alpha_min=float(xa.min()), alpha_max=float(xa.max()), alpha_mean=float(xa.mean()),
                     vf={"all": rel_vf(Gv2, xa)})
            for lo, hi in LOADS:
                e["vf"][f"{int(lo * 100)}-{min(int(hi * 100), 100)}%"] = rel_vf(Gv2, xa, np.repeat((load >= lo) & (load < hi), nvf))
            if name == "alpha_only":
                e["alpha_rms"] = float(np.sqrt(np.mean((xa - A1[inner]) ** 2)))
            else:
                e["r_alpha_exp"] = float(np.corrcoef(xa, E[inner])[0, 1])
            key = f"{name}, noise {sig:g} px"
            summary[key] = e
            figdata[key] = (E, xa)
            v = e["vf"]
            print(f"[{key}] alpha {e['alpha_min']:.2f}-{e['alpha_max']:.2f} (mean {e['alpha_mean']:.2f}), lambda "
                  f"{lam:.2g} | vf all {v['all']:.4f}; by load " + ", ".join(f"{k} {v[k]:.4f}" for k in list(v)[1:]),
                  flush=True)

    with open(os.path.join(a.out, "stage2_synthetic_c.json"), "w") as f:
        json.dump(dict(description="Synthetic test of candidate 1: exponents vary by patch, fitted with ladder (a)",
                       frames=int(len(frames)), window_px=(2 * a.win + 1) * 15, modes=a.modes, real_data_fit=REAL_FIT,
                       results=summary), f, indent=2)

    # figure: exponent pattern, recovered alpha per case (lowest noise), and VF residual vs load
    names = [n for n in sims if n != "alpha_only"]
    sig0 = a.noise[0]
    fig, ax = plt.subplots(2, len(names) + 1, figsize=(5.2 * (len(names) + 1), 8.4), constrained_layout=True)
    for c_, n_ in enumerate(names):
        E, xa = figdata[f"{n_}, noise {sig0:g} px"]
        for r_, (v, t_) in enumerate(((E[inner], f"{n_}: true exponent factor"),
                                      (xa, f"{n_}: fitted α (ladder a)"))):
            G_ = np.full((nt, nz), np.nan)
            G_[inner] = v
            span = max(np.max(np.abs(v - 1)), 0.05)
            a_ = ax[r_, c_]
            a_.set_facecolor("0.85")
            mm = a_.pcolormesh(pz, pt, G_, cmap="RdBu_r", norm=TwoSlopeNorm(1.0, 1 - span, 1 + span))
            fig.colorbar(mm, ax=a_, shrink=.8)
            a_.set_xlim(0, pz[-1]); a_.set_ylim(0, pt[-1]); a_.set_aspect("equal")
            a_.set_xlabel("z (mm)"); a_.set_ylabel("θ (mm)")
            a_.set_title(t_, fontsize=9)
    lab = [k for k in REAL_FIT if k != "all"]
    xx = np.arange(len(lab))
    for r_, sig in enumerate(a.noise[:2]):
        a_ = ax[r_, -1]
        a_.plot(xx, [REAL_FIT[k] for k in lab], "k-o", lw=2, label="real data, Stage 2 fit")
        for n_ in sims:
            v = summary[f"{n_}, noise {sig:g} px"]["vf"]
            a_.plot(xx, [v[k] for k in lab], "-o", label=n_)
        a_.set_xticks(xx, lab); a_.set_xlabel("load level (fraction of protocol peak)")
        a_.set_ylabel("relative VF residual"); a_.set_ylim(0, None)
        a_.set_title(f"VF residual by load, noise {sig:g} px", fontsize=9)
        a_.legend(fontsize=8)
    fig.suptitle("Candidate 1: exponents varying by patch, fitted with a single stiffness factor per patch")
    fig.savefig(os.path.join(a.out, "stage2_synthetic_c.png"), dpi=75)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
