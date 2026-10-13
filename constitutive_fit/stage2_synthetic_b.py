"""
Synthetic test of Stage 2 at ladder level (b): two factors per patch, beta_k for the matrix and elastin
(Cgr, C1_el) and phi_k for the fiber families (C1_smc, C1_col), exponents and gamma fixed. The membrane
tension is linear in the four constants, so each coefficient column of the Stage 2 equations splits into a
matrix part and a fiber part (2 x 64 unknowns) and the fit stays a linear least-squares problem.

Same forward model and pipeline as stage2_synthetic.py (FE with measured smoothed boundary displacements,
outer nodes from the smoothed tracks, 45 px window, 8 x 8 virtual-field modes, all blocks together, neighbor
penalty on both factor maps with lambda by leave-one-protocol-out CV). True patterns:
  uniform     beta = phi = 1                                  (FE reused from stage2_synthetic)
  same        beta = phi = Stage 1 interior alpha              (FE reused; single-factor truth)
  fiber_only  beta = 1, phi = Stage 1 interior alpha
  split       beta = Stage 1 alpha, phi = the same map mirrored left to right
Ring patches have beta = phi = 1. For each case the single-factor fit (ladder (a)) is also run, to compare
residuals.

    ../.venv/bin/python stage2_synthetic_b.py --test N1739/test1 [--noise 0 0.03]
Outputs: Images/constitutive_fit/local/stage2_synthetic/ (fe_fiber_only.npz, fe_split.npz, stage2_synthetic_b.*)
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


def split_params(P):
    Pm, Pf = P.copy(), P.copy()
    Pm[[3, 5]] = 0.0                       # matrix + elastin: Cgr, C1_el
    Pf[[0, 1]] = 0.0                       # fibers: C1_smc, C1_col
    return Pm, Pf


def scaled_blocks(blocks, sets, ncol):
    """Rows of each block scaled to unit norm and each block to equal total weight; per-protocol normal equations."""
    rows, rhs, rsets = [], [], []
    for G, b in blocks:
        G2, b2 = G.reshape(-1, ncol), b.ravel()
        nr = np.linalg.norm(G2, axis=1)
        kk = nr > 0
        sc = 1 / (nr[kk] * np.sqrt(kk.sum()))
        rows.append(G2[kk] * sc[:, None]); rhs.append(b2[kk] * sc)
        rsets.append(np.repeat(sets, G.shape[1])[kk])
    A, b, rs = np.vstack(rows), np.concatenate(rhs), np.concatenate(rsets)
    uset = np.unique(sets)
    return A, b, rs, ({s: A[rs == s].T @ A[rs == s] for s in uset}, {s: A[rs == s].T @ b[rs == s] for s in uset},
                      {s: b[rs == s] @ b[rs == s] for s in uset})


def cv_fit(AtA, Atb, btb, R, lams):
    uset = list(AtA)
    AA, Ab = sum(AtA.values()), sum(Atb.values())
    cv = []
    for lam in lams:
        e = 0.0
        for s in uset:
            x = np.linalg.solve(AA - AtA[s] + lam * R, Ab - Atb[s])
            e += x @ AtA[s] @ x - 2 * x @ Atb[s] + btb[s]
        cv.append(e)
    lam = float(lams[int(np.argmin(cv))])
    return np.linalg.solve(AA + lam * R, Ab), lam


def vf_rel(Gv, x):
    return float(np.sqrt(np.mean((Gv @ x) ** 2)) / np.sqrt(np.mean((np.abs(Gv) @ np.abs(x)) ** 2)))


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
    Pm, Pf = split_params(P)
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
    A1m = np.ones((nt, nz))
    A1m[1:-1, 1:-1] = A1[1:-1, 1:-1][:, ::-1]                 # interior mirrored left to right
    truths = {"uniform": (np.ones((nt, nz)), np.ones((nt, nz)), "fe_uniform.npz"),
              "same": (A1, A1, "fe_stage1.npz"),
              "fiber_only": (np.ones((nt, nz)), A1, "fe_fiber_only.npz"),
              "split": (A1, A1m, "fe_split.npz")}

    R1 = difference_matrix(neighbor_pairs(inner), npat)
    R2 = np.kron(np.eye(2), R1)
    lams = np.r_[0.0, np.logspace(-6, 4, 41)]
    summary, figdata = {}, {}
    for pat, (Bt, Ft_, fn) in truths.items():
        cache = os.path.join(a.out, fn)
        if os.path.exists(cache) and np.array_equal(np.load(cache)["frames"], frames):
            C = np.load(cache)
            X, Fz, Ft = C["x"], C["Fz"], C["Ft"]
            print(f"[{pat}] FE loaded from {fn}")
        else:
            agp = np.stack([Bt[gi, gj], Ft_[gi, gj]], -1)
            X, Fz, Ft, conv = [], [], [], []
            t0 = time.time()
            for n in range(len(frames)):
                xb = xsm[n][sl]
                x, info = mesh.solve(xb, agp, [Pm, Pf], H, xb)
                r = mesh.reactions(x, agp, [Pm, Pf], H)
                X.append(x); conv.append(info["max_free_force"])
                Fz.append(0.5 * (r[:, -1, 0].sum() - r[:, 0, 0].sum()) * w_t / Lt_in)
                Ft.append(0.5 * (r[-1, :, 1].sum() - r[0, :, 1].sum()) * w_z / Lz_in)
                if n % 60 == 0:
                    print(f"  [{pat}] frame {n}/{len(frames)}, {time.time() - t0:.0f} s", flush=True)
            X, Fz, Ft, conv = map(np.array, (X, Fz, Ft, conv))
            np.savez_compressed(cache, frames=frames, x=X, Fz=Fz, Ft=Ft, max_free_force=conv,
                                beta_true=Bt, phi_true=Ft_, inset=m_)
            print(f"[{pat}] FE done in {time.time() - t0:.0f} s; max free force {conv.max():.1e} mN")
        bt_true, ft_true = Bt[inner], Ft_[inner]
        for sig in a.noise:
            rng = np.random.default_rng(1)
            xs = xsm.copy()
            xs[:, m_:-m_, m_:-m_] = X
            xs = xs + rng.normal(0, sig * mpp, xs.shape)
            Fg = window_F(Xg, xs, np.ones(xs.shape[:3], bool), a.win)
            Gzm, Gtm, Gvm, bz, bt = assemble(g, Fg, Fz, Ft, Pm, H, w_t, w_z)
            Gzf, Gtf, Gvf, _, _ = assemble(g, Fg, Fz, Ft, Pf, H, w_t, w_z)
            Gz, Gt, Gv = (np.concatenate(p_, -1) for p_ in ((Gzm, Gzf), (Gtm, Gtf), (Gvm, Gvf)))
            # ladder (b)
            Ab, _, _, sysb = scaled_blocks(((Gz, bz), (Gt, bt), (Gv, np.zeros(Gv.shape[:2]))), sets, 2 * npat)
            xb_, lam_b = cv_fit(*sysb, R2, lams)
            sv = np.linalg.svd(Ab, compute_uv=False)
            beta, phi = xb_[:npat], xb_[npat:]
            # ladder (a) on the same data
            Ga = [(G[..., :npat] + G[..., npat:], b) for G, b in ((Gz, bz), (Gt, bt), (Gv, np.zeros(Gv.shape[:2])))]
            _, _, _, sysa = scaled_blocks(Ga, sets, npat)
            xa, lam_a = cv_fit(*sysa, R1, lams)
            Gv2 = Gv.reshape(-1, 2 * npat)
            res_b = vf_rel(Gv2, xb_)
            res_a = vf_rel(Gv2, np.r_[xa, xa])
            res_t = vf_rel(Gv2, np.r_[bt_true, ft_true])
            e = dict(beta_rms=float(np.sqrt(np.mean((beta - bt_true) ** 2))),
                     phi_rms=float(np.sqrt(np.mean((phi - ft_true) ** 2))),
                     beta_r=float(np.corrcoef(beta, bt_true)[0, 1]) if bt_true.std() > 0 else None,
                     phi_r=float(np.corrcoef(phi, ft_true)[0, 1]) if ft_true.std() > 0 else None,
                     beta_mean=float(beta.mean()), phi_mean=float(phi.mean()),
                     lam_b=lam_b, lam_a=lam_a, cond_b=float(sv[0] / sv[-1]),
                     vf_residual=dict(true=res_t, ladder_b=res_b, ladder_a=res_a))
            key = f"{pat}, noise {sig:g} px"
            summary[key] = e
            figdata[key] = (bt_true, ft_true, beta, phi, xa)
            print(f"[{key}] (b): beta RMS {e['beta_rms']:.3f} r {e['beta_r'] if e['beta_r'] is None else round(e['beta_r'], 2)}"
                  f", phi RMS {e['phi_rms']:.3f} r {e['phi_r'] if e['phi_r'] is None else round(e['phi_r'], 2)}, "
                  f"means {e['beta_mean']:.2f}/{e['phi_mean']:.2f}, lambda {lam_b:.2g}, cond {e['cond_b']:.0f} | "
                  f"vf residual: true {res_t:.4f}, (b) {res_b:.4f}, (a) {res_a:.4f}", flush=True)

    with open(os.path.join(a.out, "stage2_synthetic_b.json"), "w") as f:
        json.dump(dict(description="Synthetic test of ladder (b): beta (Cgr, C1_el) and phi (C1_smc, C1_col) per patch",
                       frames=int(len(frames)), window_px=(2 * a.win + 1) * 15, modes=a.modes, results=summary), f,
                  indent=2)
    keys = list(figdata)
    fig, ax = plt.subplots(len(keys), 5, figsize=(22, 3.6 * len(keys)), constrained_layout=True, squeeze=False)
    for r_, key in enumerate(keys):
        vals = figdata[key]
        span = max(max(np.max(np.abs(v - 1)) for v in vals[:4]), 0.1)
        norm = TwoSlopeNorm(1.0, max(1 - span, -0.5), 1 + span)
        for c_, (v, t_) in enumerate(zip(vals, ("true β (matrix/elastin)", "true φ (fibers)", "recovered β",
                                                "recovered φ", "single-factor α (ladder a)"))):
            G_ = np.full((nt, nz), np.nan)
            G_[inner] = v
            a_ = ax[r_, c_]
            a_.set_facecolor("0.85")
            mm = a_.pcolormesh(pz, pt, G_, cmap="RdBu_r", norm=norm)
            fig.colorbar(mm, ax=a_, shrink=.8, extend="both")
            a_.set_xlim(0, pz[-1]); a_.set_ylim(0, pt[-1]); a_.set_aspect("equal")
            s_ = summary[key]
            extra = {2: f" (RMS {s_['beta_rms']:.2f})", 3: f" (RMS {s_['phi_rms']:.2f})"}.get(c_, "")
            a_.set_title(f"{key}: {t_}{extra}", fontsize=8)
    fig.suptitle("Ladder (b) synthetic test: matrix/elastin factor β and fiber factor φ per patch")
    fig.savefig(os.path.join(a.out, "stage2_synthetic_b.png"), dpi=75)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
