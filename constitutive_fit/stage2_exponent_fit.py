"""
Candidate 1 on the data: per patch, a stiffness factor alpha_k (all C1 constants and Cgr, as in ladder (a)) and an
exponent factor eps_k scaling the stiffening exponents (C2_smc, C2_col by default; --exps all adds C2_el).

Column k of the Stage 2 equations depends only on the material of patch k, so the columns are assembled once on a
grid of eps levels (one assembly per level with the exponent scaled everywhere) and interpolated per patch with
Catmull-Rom splines. Unknowns alpha_k and s_k = log eps_k; residual = row-scaled combined system (cuts + virtual
fields, rows scaled to unit norm at eps = 1, blocks equal weight) plus neighbor penalties sqrt(lam_a) D alpha and
sqrt(lam_e) D s. lam_a is the ladder (a) value from stage2_fit.json; lam_e by leave-one-protocol-out CV.
Solved with scipy least_squares (trust region, eps bounded to the level grid), started from ladder (a).

    ../.venv/bin/python stage2_exponent_fit.py --test N1739/test1 [--exps fiber|all] [--synthetic fib40]
--synthetic NAME fits the FE data from stage2_synthetic_c.py (fe_exp_NAME.npz) instead, as a recovery check.
Outputs: Images/constitutive_fit/local/stage2_exponent/ (stage2_exponent_<tag>.json/.npz/.png)
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
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from scipy.ndimage import gaussian_filter  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import baseline as B  # noqa: E402
import constitutive as M  # noqa: E402
from stage2_setup import LOCAL, assemble, geometry, measured_positions, window_F  # noqa: E402
from stage2_synthetic import difference_matrix, neighbor_pairs  # noqa: E402

EXPS = {"fiber": (4, 6), "all": (2, 4, 6)}
LEVELS = np.round(np.arange(0.4, 2.01, 0.1), 2)
LOADS = ((0, .25), (.25, .5), (.5, .75), (.75, 1.01))


def catmull_rom(u):
    """Weights and their derivatives (per unit of u) for the 4 levels around fractional index u."""
    i = np.clip(np.floor(u).astype(int), 1, len(LEVELS) - 3)
    t = u - i
    w = np.stack([(-t ** 3 + 2 * t ** 2 - t), (3 * t ** 3 - 5 * t ** 2 + 2), (-3 * t ** 3 + 4 * t ** 2 + t),
                  (t ** 3 - t ** 2)], 0) / 2
    dw = np.stack([(-3 * t ** 2 + 4 * t - 1), (9 * t ** 2 - 10 * t), (-9 * t ** 2 + 8 * t + 1),
                   (3 * t ** 2 - 2 * t)], 0) / 2
    return i, w, dw


class System:
    """Row-scaled columns on the eps grid: A (levels, rows, npat), b (rows); block, set and frame of each row."""

    def __init__(self, g, F, Fz, Ft, P, H, w_t, w_z, cols, sets):
        npat = g["npat"]
        A, sc, b, blk, rs, rf = [], None, None, None, None, None
        for lev in np.r_[1.0, LEVELS]:                 # eps = 1 first: it sets the row scaling
            Pl = P.copy()
            Pl[list(cols)] *= lev
            Gz, Gt, Gv, bz, bt = assemble(g, F, Fz, Ft, Pl, H, w_t, w_z)
            parts = ((Gz, bz), (Gt, bt), (Gv, np.zeros(Gv.shape[:2])))
            G2 = np.vstack([G_.reshape(-1, npat) for G_, _ in parts])
            if sc is None:
                nr = [np.linalg.norm(G_.reshape(-1, npat), axis=1) for G_, _ in parts]
                keep = np.concatenate([n_ > 0 for n_ in nr])
                sc = np.concatenate([np.where(n_ > 0, 1 / np.where(n_ > 0, n_, 1) / np.sqrt((n_ > 0).sum()), 0)
                                     for n_ in nr])[keep]
                b = np.concatenate([b_.ravel() for _, b_ in parts])[keep] * sc
                blk = np.concatenate([np.full(G_.shape[0] * G_.shape[1], j) for j, (G_, _) in enumerate(parts)])[keep]
                fr = np.concatenate([np.repeat(np.arange(G_.shape[0]), G_.shape[1]) for G_, _ in parts])[keep]
                rs, rf = sets[fr], fr
                continue
            A.append((G2[keep] * sc[:, None]).astype(np.float32))
        self.A, self.sc, self.b, self.blk, self.rs, self.rf = np.stack(A), sc, b, blk, rs, rf
        self.npat = npat

    def columns(self, s):
        """G (rows, npat) at eps = exp(s) and dG/ds."""
        eps = np.exp(s)
        u = (eps - LEVELS[0]) / (LEVELS[1] - LEVELS[0])
        i, w, dw = catmull_rom(u)
        k = np.arange(self.npat)
        G = sum(w[j] * self.A[i - 1 + j, :, k].T for j in range(4))
        dG = sum(dw[j] * self.A[i - 1 + j, :, k].T for j in range(4)) * (eps / (LEVELS[1] - LEVELS[0]))
        return G.astype(float), dG.astype(float)


def fit(S, m, Da, De, lam_a, lam_e, x0, lo, hi):
    npat = S.npat
    b = S.b[m]

    def fun(x):
        al, s = x[:npat], x[npat:]
        G, _ = S.columns(s)
        return np.r_[G[m] @ al - b, np.sqrt(lam_a) * Da @ al, np.sqrt(lam_e) * De @ s]

    def jac(x):
        al, s = x[:npat], x[npat:]
        G, dG = S.columns(s)
        J = np.zeros((m.sum() + 2 * Da.shape[0], 2 * npat))
        J[:m.sum(), :npat] = G[m]
        J[:m.sum(), npat:] = dG[m] * al[None]
        J[m.sum():m.sum() + Da.shape[0], :npat] = np.sqrt(lam_a) * Da
        J[m.sum() + Da.shape[0]:, npat:] = np.sqrt(lam_e) * De
        return J
    r = least_squares(fun, x0, jac=jac, bounds=(lo, hi), method="trf", x_scale="jac", max_nfev=200,
                      ftol=1e-10, xtol=1e-10)
    return r.x, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exps", choices=list(EXPS), default="fiber")
    ap.add_argument("--synthetic", default=None, help="fe_exp_<NAME>.npz from stage2_synthetic_c.py")
    ap.add_argument("--every", type=int, default=4)
    ap.add_argument("--win", type=int, default=1)
    ap.add_argument("--modes", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(LOCAL, "stage2_exponent"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    t0 = time.time()
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
    fitinfo = json.load(open(os.path.join(LOCAL, "stage2_fit", "stage2_fit.json")))
    lam_a = float(np.load(os.path.join(LOCAL, "stage2_fit", "stage2_fit.npz"))["lam"])
    truth = None
    if a.synthetic:
        C = np.load(os.path.join(LOCAL, "stage2_synthetic", f"fe_exp_{a.synthetic}.npz"))
        fidx = C["frames"]
        xm = measured_positions(df, fidx, mpp)
        xs = Xg + np.stack([gaussian_filter(xm[..., i] - Xg[..., i], (0, 1.0, 1.0), mode="nearest") for i in range(2)], -1)
        m_ = int(C["inset"])
        xs[:, m_:-m_, m_:-m_] = C["x"]
        F = window_F(Xg, xs, np.ones(xs.shape[:3], bool), a.win)
        Fz, Ft = C["Fz"], C["Ft"]
        truth = C["exp_factor"][inner]
        if tuple(C["exp_cols"]) != EXPS[a.exps]:
            print(f"note: synthetic exponents {tuple(C['exp_cols'])} differ from the fitted set {EXPS[a.exps]}")
        tag = f"synthetic_{a.synthetic}_{a.exps}"
    else:
        fidx = np.arange(0, len(df), a.every)
        xm, vm = measured_positions(df, fidx, mpp, with_valid=True)
        F = window_F(Xg, xm, vm, a.win)
        Fz, Ft = d["force_z_mN"][fidx], d["force_theta_mN"][fidx]
        tag = f"real_{a.exps}"
    sets = df.set.to_numpy()[fidx]
    S = System(g, F, Fz, Ft, P, H, w_t, w_z, EXPS[a.exps], sets)
    print(f"[{tag}] {len(fidx)} frames, {len(S.b)} rows, {len(LEVELS)} exponent levels; assembled in "
          f"{time.time() - t0:.0f} s; lam_a {lam_a:g}", flush=True)

    Da = np.zeros((0, npat))
    pairs = neighbor_pairs(inner)
    Da = np.zeros((len(pairs), npat))
    for r_, (i_, j_) in enumerate(pairs):
        Da[r_, i_], Da[r_, j_] = 1, -1
    De = Da.copy()
    allm = np.ones(len(S.b), bool)

    # start: ladder (a) at eps = 1 with the same lam_a
    G1, _ = S.columns(np.zeros(npat))
    al0 = np.linalg.solve(G1.T @ G1 + lam_a * Da.T @ Da, G1.T @ S.b)
    x0 = np.r_[al0, np.zeros(npat)]
    lo = np.r_[np.full(npat, -np.inf), np.full(npat, np.log(LEVELS[1] + 1e-6))]
    hi = np.r_[np.full(npat, np.inf), np.full(npat, np.log(LEVELS[-2] - 1e-6))]

    # leave-one-protocol-out CV for lam_e
    lam_es = np.logspace(-4, 1, 6)
    uset = np.unique(sets)
    cv = []
    for le in lam_es:
        e = 0.0
        for s_ in uset:
            m = S.rs != s_
            x, _ = fit(S, m, Da, De, lam_a, le, x0, lo, hi)
            G, _ = S.columns(x[npat:])
            e += float(np.sum((G[~m] @ x[:npat] - S.b[~m]) ** 2))
        cv.append(e)
        print(f"  lam_e {le:.0e}: CV {e:.5g} ({time.time() - t0:.0f} s)", flush=True)
    e0 = sum(float(np.sum((G1[S.rs == s_] @ np.linalg.solve(G1[S.rs != s_].T @ G1[S.rs != s_] + lam_a * Da.T @ Da,
                                                              G1[S.rs != s_].T @ S.b[S.rs != s_]) - S.b[S.rs == s_]) ** 2))
             for s_ in uset)
    print(f"  ladder (a), eps = 1: CV {e0:.5g}")
    lam_e = float(lam_es[int(np.argmin(cv))])
    x, res = fit(S, allm, Da, De, lam_a, lam_e, x0, lo, hi)
    al, eps = x[:npat], np.exp(x[npat:])
    print(f"  final fit: {res.nfev} evaluations, status {res.status}")

    # residuals on the unscaled equations
    load = np.zeros(len(fidx))
    for s_ in uset:
        mm = sets == s_
        load[mm] = np.maximum(Fz[mm] / Fz[mm].max(), Ft[mm] / Ft[mm].max())
    rload = load[S.rf]

    def rel(al_, s_, m=None):
        G, _ = S.columns(s_)
        Gu = G / S.sc[:, None]
        bu = S.b / S.sc
        out = {}
        for j, nm in enumerate(("cut_z", "cut_theta")):
            k = S.blk == j
            r = Gu[k] @ al_ - bu[k]
            out[nm] = float(np.sqrt(np.mean(r ** 2)) / np.sqrt(np.mean(bu[k] ** 2)))
        k = S.blk == 2
        out["vf"] = float(np.sqrt(np.mean((Gu[k] @ al_) ** 2)) / np.sqrt(np.mean((np.abs(Gu[k]) @ np.abs(al_)) ** 2)))
        for lo_, hi_ in LOADS:
            kk = k & (rload >= lo_) & (rload < hi_)
            out[f"vf {int(lo_ * 100)}-{min(int(hi_ * 100), 100)}%"] = float(
                np.sqrt(np.mean((Gu[kk] @ al_) ** 2)) / np.sqrt(np.mean((np.abs(Gu[kk]) @ np.abs(al_)) ** 2)))
        return out
    ra, rb = rel(al0, np.zeros(npat)), rel(al, x[npat:])
    fmt = lambda r: ", ".join(f"{k} {v:.4f}" for k, v in r.items())  # noqa: E731
    print(f"  ladder (a):   {fmt(ra)}")
    print(f"  alpha + eps:  {fmt(rb)}")
    print(f"  lam_e {lam_e:g}; alpha {al.min():.2f}-{al.max():.2f} (mean {al.mean():.2f}); eps {eps.min():.2f}-"
          f"{eps.max():.2f} (mean {eps.mean():.2f}); at bounds {np.sum(np.isclose(x[npat:], lo[npat:]) | np.isclose(x[npat:], hi[npat:]))}; "
          f"r(alpha, alpha_a) {np.corrcoef(al, al0)[0, 1]:.2f}, r(alpha, eps) {np.corrcoef(al, eps)[0, 1]:.2f}")
    summ = dict(tag=tag, exponents=[M.PARAM_NAMES[c] for c in EXPS[a.exps]], frames=int(len(fidx)),
                lam_a=lam_a, lam_e=lam_e, lam_e_grid=lam_es.tolist(), cv=cv, cv_ladder_a=e0,
                residuals_ladder_a=ra, residuals_alpha_eps=rb,
                alpha=dict(min=float(al.min()), max=float(al.max()), mean=float(al.mean())),
                eps=dict(min=float(eps.min()), max=float(eps.max()), mean=float(eps.mean())),
                r_alpha_alpha_a=float(np.corrcoef(al, al0)[0, 1]), r_alpha_eps=float(np.corrcoef(al, eps)[0, 1]))
    if truth is not None:
        summ["eps_rms"] = float(np.sqrt(np.mean((eps - truth) ** 2)))
        summ["eps_r"] = float(np.corrcoef(eps, truth)[0, 1])
        summ["alpha_rms_vs_1"] = float(np.sqrt(np.mean((al - 1) ** 2)))
        summ["alpha_a_rms_vs_1"] = float(np.sqrt(np.mean((al0 - 1) ** 2)))
        print(f"  truth: eps RMS {summ['eps_rms']:.3f}, r {summ['eps_r']:.2f}; alpha RMS vs 1: joint "
              f"{summ['alpha_rms_vs_1']:.3f}, ladder (a) {summ['alpha_a_rms_vs_1']:.3f}")
    with open(os.path.join(a.out, f"stage2_exponent_{tag}.json"), "w") as f:
        json.dump(summ, f, indent=2)
    np.savez_compressed(os.path.join(a.out, f"stage2_exponent_{tag}.npz"), alpha=al, eps=eps, alpha_a=al0,
                        inner=inner, lam_a=lam_a, lam_e=lam_e)

    maps = [(al0, "ladder (a) α"), (al, "joint fit α"), (eps, "joint fit ε (exponent factor)")]
    if truth is not None:
        maps.append((truth, "true ε"))
    fig, ax = plt.subplots(1, len(maps) + 1, figsize=(5.2 * (len(maps) + 1), 4.6), constrained_layout=True)
    amax = max(np.max(np.abs(al0 - 1)), np.max(np.abs(al - 1)))
    for c_, (v, t_) in enumerate(maps):
        G_ = np.full((nt, nz), np.nan)
        G_[inner] = v
        span = amax if c_ < 2 else max(np.max(np.abs(eps - 1)), np.max(np.abs(truth - 1)) if truth is not None else 0)
        a_ = ax[c_]
        a_.set_facecolor("0.85")
        mm = a_.pcolormesh(pz, pt, G_, cmap="RdBu_r", norm=TwoSlopeNorm(1.0, 1 - span, 1 + span))
        fig.colorbar(mm, ax=a_, shrink=.8, extend="both")
        a_.set_xlim(0, pz[-1]); a_.set_ylim(0, pt[-1]); a_.set_aspect("equal")
        a_.set_xlabel("z (mm)"); a_.set_ylabel("θ (mm)"); a_.set_title(t_, fontsize=10)
    a_ = ax[-1]
    lab = [k for k in ra if k.startswith("vf ")]
    xx = np.arange(len(lab))
    a_.plot(xx, [ra[k] for k in lab], "-o", label=f"ladder (a), all {ra['vf']:.3f}")
    a_.plot(xx, [rb[k] for k in lab], "-o", label=f"α + ε, all {rb['vf']:.3f}")
    a_.set_xticks(xx, [k[3:] for k in lab]); a_.set_ylim(0, None)
    a_.set_xlabel("load level"); a_.set_ylabel("relative VF residual"); a_.legend(fontsize=8)
    a_.set_title("VF residual by load", fontsize=10)
    fig.suptitle(f"Per-patch stiffness α and exponent factor ε ({', '.join(summ['exponents'])}): {tag}")
    fig.savefig(os.path.join(a.out, f"stage2_exponent_{tag}.png"), dpi=75)
    print(f"wrote {a.out} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
