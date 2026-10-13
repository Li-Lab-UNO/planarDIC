"""
Synthetic test of the Stage 2 virtual-field equations (stage2_setup.py) with a known patch stiffness map.

Forward model (membrane_fe.py): bilinear elements on the DIC grid, inset by --inset nodes from the
tracked edge (the raw element-scale stretch next to the tine tips, up to lambda_z ~1.6, is far outside
the range the bulk model was fitted on). For each selected loading frame, the measured positions of the
inset boundary nodes, smoothed with a Gaussian of one node, are imposed, and the energy is minimized for
a known alpha map (alpha at each Gauss point = the alpha of its patch). Material: the full-domain bulk
baseline (dic_affine) with Cgr, C1_el, C1_smc, C1_col scaled by alpha.

Synthetic "load cell" forces: the mean of the reactions on the two opposite edges, scaled to the tine-span
width as the real forces are assumed to be (F_z = mean edge reaction x w_theta / inset height, likewise F_theta),
so the cut equations apply the same uniform-flux fractions as for the measured data.

Recovery: node positions (FE solution inside the inset, smoothed tracks on the outer nodes, + Gaussian noise of
--noise px) -> F with the Stage 2 window (--win 1 = 45 px) -> the same interior equations (stage2_setup.assemble). Two methods:

  shape+scale (recommended)  the pattern of alpha from the virtual fields alone (homogeneous equations, so
              alpha is found up to a common factor; normalized to mean 1):
                  min ||A_vf alpha||^2 + (lam / n_pairs) sum_{k~l} (alpha_k - alpha_l)^2   s.t. mean(alpha) = 1
              then one scale factor from all cut equations (least squares). lam by leave-one-protocol-out
              cross-validation of the virtual-field residual.
  combined    all blocks at once (rows unit-normalized, equal block weights) with the same penalty; shown for
              comparison. The cut equations carry 6-15% systematic error from the uniform-flux assumption even
              with perfect data, and the combined fit absorbs it by distorting alpha.

Patterns: 'uniform' (alpha = 1 everywhere: the null test; its recovered RMS error is the detection floor for
patch contrasts) and 'stage1' (interior alpha from Stage 1, outer ring alpha = 1).

    ../.venv/bin/python stage2_synthetic.py --test N1739/test1 [--every 8] [--modes 8] [--win 1] [--noise 0 0.03 0.1]
Outputs: Images/constitutive_fit/local/stage2_synthetic/ (FE solutions are cached as fe_<pattern>.npz)
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
import testcfg  # noqa: E402
import baseline as B  # noqa: E402
import constitutive as M  # noqa: E402
from membrane_fe import Mesh  # noqa: E402
from stage2_setup import assemble, geometry, measured_positions, window_F  # noqa: E402

IMAGES = testcfg.get().images
LOCAL = os.path.join(IMAGES, "constitutive_fit", "local")


def neighbor_pairs(mask):
    nt, nz = mask.shape
    idx = -np.ones(mask.shape, int)
    idx[mask] = np.arange(mask.sum())
    pairs = [(idx[i, j], idx[i + di, j + dj]) for i in range(nt) for j in range(nz) if mask[i, j]
             for di, dj in ((0, 1), (1, 0)) if i + di < nt and j + dj < nz and mask[i + di, j + dj]]
    return np.array(pairs)


def difference_matrix(pairs, n):
    Dm = np.zeros((len(pairs), n))
    Dm[np.arange(len(pairs)), pairs[:, 0]] = 1
    Dm[np.arange(len(pairs)), pairs[:, 1]] = -1
    return Dm.T @ Dm / len(pairs)


def unit_rows(G):
    G = G.reshape(-1, G.shape[-1])
    n = np.linalg.norm(G, axis=1)
    keep = n > 0
    return G[keep] / n[keep, None], keep


def vf_shape(K, lam, DtD):
    """argmin a'(K + lam DtD)a subject to mean(a) = 1."""
    n = K.shape[0]
    Mx = np.zeros((n + 1, n + 1))
    Mx[:n, :n] = 2 * (K + lam * DtD)
    Mx[:n, n] = Mx[n, :n] = 1
    r = np.zeros(n + 1)
    r[n] = n
    return np.linalg.solve(Mx, r)[:n]


def cut_scale(shape, Gz, bz, Gt, bt):
    G = np.vstack([Gz.reshape(-1, len(shape)), Gt.reshape(-1, len(shape))])
    b = np.r_[bz.ravel(), bt.ravel()]
    m = G @ shape
    s = (m @ b) / (m @ m)
    mz, mt = Gz.reshape(-1, len(shape)) @ shape, Gt.reshape(-1, len(shape)) @ shape
    return s, (mz @ bz.ravel()) / (mz @ mz), (mt @ bt.ravel()) / (mt @ mt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--every", type=int, default=8)
    ap.add_argument("--inset", type=int, default=2, help="FE boundary this many nodes inside the tracked edge")
    ap.add_argument("--noise", type=float, nargs="+", default=[0.0, 0.03, 0.1], help="position noise (px)")
    ap.add_argument("--patterns", nargs="+", default=["uniform", "stage1"])
    ap.add_argument("--exclude-ring", type=int, default=1)
    ap.add_argument("--modes", type=int, default=8)
    ap.add_argument("--win", type=int, default=1, help="strain window half-width (1 = 45 px, 2 = 75 px)")
    ap.add_argument("--out", default=os.path.join(LOCAL, "stage2_synthetic"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    d = np.load(os.path.join(LOCAL, "local_data.npz"))
    base, df = B.load(str(d["baseline"]))
    df = df.reset_index(drop=True)
    st = base["settings"]
    H, w_t, w_z = st["H_mm"], st["w_theta_mm"], st["w_z_mm"]
    P = np.array([base["inputs"]["dic_affine"]["params"][n] for n in M.PARAM_NAMES])
    mpp = float(d["mm_per_px"])
    g = geometry(d, a.exclude_ring, a.modes)
    inner, npat, nt, nz, pz, pt = g["inner"], g["npat"], g["nt"], g["nz"], g["pz"], g["pt"]
    Xg = np.stack([d["z_mm"], d["theta_mm"]], -1)
    m_ = a.inset
    sl = (slice(m_, -m_), slice(m_, -m_))
    mesh = Mesh(Xg[sl])
    gp = mesh.gauss_points
    gi = np.clip(np.searchsorted(pt, gp[..., 1]) - 1, 0, nt - 1)
    gj = np.clip(np.searchsorted(pz, gp[..., 0]) - 1, 0, nz - 1)
    frames = np.arange(0, len(df), a.every)
    sets = df.set.to_numpy()[frames]
    uset = np.unique(sets)
    Lt_in = (mesh.ny - 1) * abs(mesh.hy)
    Lz_in = (mesh.nx - 1) * abs(mesh.hx)
    print(f"{len(frames)} frames (every {a.every}); FE mesh {mesh.nx} x {mesh.ny} nodes (inset {m_}), "
          f"{Lz_in:.2f} x {Lt_in:.2f} mm; {npat} interior patches; {2 * a.modes ** 2} virtual fields per frame")

    s1 = pd.read_csv(os.path.join(LOCAL, f"stage1_interior{a.exclude_ring}", "stage1_patches.csv"))
    truths = {"uniform": np.ones((nt, nz))}
    t1 = np.ones((nt, nz))
    t1[s1.i_theta, s1.j_z] = s1.alpha
    truths["stage1"] = t1
    DtD = difference_matrix(neighbor_pairs(inner), npat)
    lams = np.r_[0.0, np.logspace(-6, 2, 33)]
    summary, results = {}, {}
    xm = measured_positions(df, frames, mpp)
    xsm = Xg + np.stack([gaussian_filter(xm[..., i] - Xg[..., i], (0, 1.0, 1.0), mode="nearest")
                         for i in range(2)], -1)

    for pat in a.patterns:
        A_true = truths[pat]
        cache = os.path.join(a.out, f"fe_{pat}.npz")
        if os.path.exists(cache) and np.array_equal(np.load(cache)["frames"], frames) \
                and int(np.load(cache)["inset"]) == m_:
            C = np.load(cache)
            X, Fz, Ft = C["x"], C["Fz"], C["Ft"]
            print(f"[{pat}] FE solutions loaded from cache (max free force {C['max_free_force'].max():.1e} mN)")
        else:
            alpha_gp = A_true[gi, gj]
            X, Fz, Ft, conv, side = [], [], [], [], []
            t0 = time.time()
            for n in range(len(frames)):
                xb = xsm[n][sl]
                x, info = mesh.solve(xb, alpha_gp, P, H, xb)
                r = mesh.reactions(x, alpha_gp, P, H)
                X.append(x); conv.append(info["max_free_force"])
                Fz.append(0.5 * (r[:, -1, 0].sum() - r[:, 0, 0].sum()) * w_t / Lt_in)
                Ft.append(0.5 * (r[-1, :, 1].sum() - r[0, :, 1].sum()) * w_z / Lz_in)
                side.append((r[:, -1, 0].sum(), -r[:, 0, 0].sum(), r[-1, :, 1].sum(), -r[0, :, 1].sum()))
                if n % 40 == 0:
                    print(f"  [{pat}] frame {n}/{len(frames)}  {info['nit']} iters, max free force "
                          f"{info['max_free_force']:.1e} mN, {time.time() - t0:.0f} s", flush=True)
            X, Fz, Ft, conv, side = map(np.array, (X, Fz, Ft, conv, side))
            np.savez_compressed(cache, frames=frames, x=X, Fz=Fz, Ft=Ft, max_free_force=conv, side_reactions=side,
                                alpha_true=A_true, inset=m_)
            print(f"[{pat}] FE done in {time.time() - t0:.0f} s; max free force {conv.max():.1e} mN")

        truth = A_true[inner]
        for sig in a.noise:
            # "measured" positions on the whole grid: the FE solution inside, the smoothed tracks on the outer
            # nodes (DIC sees them too), so every interior cell has a full window as in the real data
            rng = np.random.default_rng(1)
            xs = xsm.copy()
            xs[:, m_:-m_, m_:-m_] = X
            xs = xs + rng.normal(0, sig * mpp, xs.shape)
            Fgrid = window_F(Xg, xs, np.ones(xs.shape[:3], bool), a.win)
            Gz, Gt, Gv, bz, bt = assemble(g, Fgrid, Fz, Ft, P, H, w_t, w_z)
            rmse = lambda al: float(np.sqrt(np.mean((al - truth) ** 2)))                     # noqa: E731
            corr = lambda al: float(np.corrcoef(al, truth)[0, 1]) if truth.std() > 0 else float("nan")  # noqa: E731

            # --- shape + scale
            Av, keep = unit_rows(Gv)
            rs = np.repeat(sets, Gv.shape[1])[keep]
            Kp = {s: Av[rs == s].T @ Av[rs == s] for s in uset}
            K = sum(Kp.values())
            cv = np.array([sum(vf_shape(K - Kp[s], lam, DtD) @ Kp[s] @ vf_shape(K - Kp[s], lam, DtD) for s in uset)
                           for lam in lams])
            m_cv = int(np.argmin(cv))
            ss = {}
            for name, lam in (("lambda=0", 0.0), ("lambda_cv", float(lams[m_cv]))):
                shp = vf_shape(K, lam, DtD)
                s_all, s_z, s_t = cut_scale(shp, Gz, bz, Gt, bt)
                al = s_all * shp
                ss[name] = dict(lam=lam, rms=rmse(al), max=float(np.abs(al - truth).max()), corr=corr(al),
                                mean=float(al.mean()), true_mean=float(truth.mean()),
                                scale_z_cuts=float(s_z), scale_theta_cuts=float(s_t),
                                shape_rms=float(np.sqrt(np.mean((shp - truth / truth.mean()) ** 2))), alpha=al)
            ecurve = np.array([rmse(cut_scale(vf_shape(K, lam, DtD), Gz, bz, Gt, bt)[0] * vf_shape(K, lam, DtD))
                               for lam in lams])

            # --- combined (for comparison)
            rows, rhs, rsets = [], [], []
            for G, b in ((Gz, bz), (Gt, bt), (Gv, np.zeros(Gv.shape[:2]))):
                G2, b2 = G.reshape(-1, npat), b.ravel()
                nr = np.linalg.norm(G2, axis=1)
                kk = nr > 0
                sc = 1 / (nr[kk] * np.sqrt(kk.sum()))
                rows.append(G2[kk] * sc[:, None]); rhs.append(b2[kk] * sc); rsets.append(np.repeat(sets, G.shape[1])[kk])
            Ac, bc, rc = np.vstack(rows), np.concatenate(rhs), np.concatenate(rsets)
            AtA = {s: Ac[rc == s].T @ Ac[rc == s] for s in uset}
            Atb = {s: Ac[rc == s].T @ bc[rc == s] for s in uset}
            btb = {s: bc[rc == s] @ bc[rc == s] for s in uset}
            AA, Ab = sum(AtA.values()), sum(Atb.values())
            cvc = []
            for lam in lams:
                e = 0.0
                for s in uset:
                    al = np.linalg.solve(AA - AtA[s] + lam * DtD, Ab - Atb[s])
                    e += al @ AtA[s] @ al - 2 * al @ Atb[s] + btb[s]
                cvc.append(e)
            lam_c = float(lams[int(np.argmin(cvc))])
            al_c = np.linalg.solve(AA + lam_c * DtD, Ab)
            comb = dict(lam=lam_c, rms=rmse(al_c), corr=corr(al_c), alpha=al_c)

            key = f"{pat}, noise {sig:g} px"
            results[key] = dict(ss=ss, comb=comb, cv=cv, ecurve=ecurve, truth=truth)
            summary[key] = dict(shape_scale={k: {kk: vv for kk, vv in v.items() if kk != "alpha"} for k, v in ss.items()},
                                combined={k: v for k, v in comb.items() if k != "alpha"},
                                best_lambda_on_grid=dict(lam=float(lams[int(np.argmin(ecurve))]), rms=float(ecurve.min())))
            c = ss["lambda_cv"]
            print(f"[{key}] shape+scale: CV lambda {c['lam']:.2g}, RMS error {c['rms']:.3f}, r {c['corr']:.2f}, "
                  f"mean {c['mean']:.3f} (true {c['true_mean']:.3f}; z cuts {c['scale_z_cuts']:.3f}, theta cuts "
                  f"{c['scale_theta_cuts']:.3f}) | lambda=0: RMS {ss['lambda=0']['rms']:.3f} | combined: RMS "
                  f"{comb['rms']:.3f}, r {comb['corr']:.2f}", flush=True)

    with open(os.path.join(a.out, "stage2_synthetic.json"), "w") as f:
        json.dump(dict(description="Synthetic recovery of interior patch alpha with the Stage 2 equations; FE forward "
                       "model with measured (smoothed) boundary displacements",
                       frames=int(len(frames)), every=a.every, inset_nodes=m_, modes=a.modes,
                       window_px=(2 * a.win + 1) * 15,
                       virtual_fields_per_frame=2 * a.modes ** 2, results=summary), f, indent=2)

    # ---------------------------------------------------------------- figure
    keys = list(results)
    fig, ax = plt.subplots(len(keys), 4, figsize=(19, 3.9 * len(keys)), constrained_layout=True, squeeze=False)
    for r_, key in enumerate(keys):
        R_ = results[key]
        maps = [(R_["truth"], "true α"),
                (R_["ss"]["lambda_cv"]["alpha"], f"shape+scale, CV λ = {R_['ss']['lambda_cv']['lam']:.2g} "
                                                f"(RMS {R_['ss']['lambda_cv']['rms']:.3f})"),
                (R_["comb"]["alpha"], f"combined, CV λ = {R_['comb']['lam']:.2g} (RMS {R_['comb']['rms']:.3f})")]
        allv = np.concatenate([v for v, _ in maps[:2]])
        span = max(np.max(np.abs(allv - 1)), 0.05)
        norm = TwoSlopeNorm(1.0, 1 - span, 1 + span)
        for c_, (v, title) in enumerate(maps):
            grid = np.full((nt, nz), np.nan)
            grid[inner] = v
            a_ = ax[r_, c_]
            a_.set_facecolor("0.85")
            mm = a_.pcolormesh(pz, pt, grid, cmap="RdBu_r", norm=norm)
            fig.colorbar(mm, ax=a_, shrink=.8, extend="both")
            a_.set_xlim(0, pz[-1]); a_.set_ylim(0, pt[-1]); a_.set_aspect("equal")
            a_.set_title(f"{key}: {title}", fontsize=9)
        a_ = ax[r_, 3]
        a_.semilogx(lams[1:], R_["ecurve"][1:], "o-", ms=3, color="C1", label="shape+scale: RMS error vs truth")
        a_.axhline(R_["ecurve"][0], color="C1", ls=":", lw=1, label=f"λ = 0: {R_['ecurve'][0]:.3f}")
        a_.axvline(R_["ss"]["lambda_cv"]["lam"] or 1e-6, color="r", ls="--", lw=1, label="CV choice")
        a2 = a_.twinx()
        a2.semilogx(lams[1:], R_["cv"][1:] / R_["cv"].min(), "-", color="C0", lw=1)
        a2.set_ylabel("CV residual / minimum", color="C0")
        a_.set_xlabel("λ"); a_.set_ylabel("RMS error of α"); a_.grid(alpha=.3); a_.legend(fontsize=7)
        a_.set_title("error vs truth (orange), CV residual (blue)", fontsize=9)
    fig.suptitle(f"Stage 2 synthetic test ({len(frames)} frames, {2 * a.modes ** 2} virtual fields per frame, "
                 f"{(2 * a.win + 1) * 15} px window): "
                 "FE data → DIC F pipeline → interior equations")
    fig.savefig(os.path.join(a.out, "stage2_synthetic.png"), dpi=80)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
