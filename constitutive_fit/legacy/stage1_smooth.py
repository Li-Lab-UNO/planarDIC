"""
Stage 1 with a smoothness penalty on log alpha (Tikhonov regularization between neighboring patches).

Same model as stage1_patches.py: patch k's stress is alpha_k * t(F_k; P), with P the full-domain bulk
baseline (dic_affine), F_k the patch-mean F and the GLOBAL stresses as data (uniform-stress
assumption). With beta_k = log alpha_k, the fit minimizes

    sum_k c_k (alpha_k - alpha_hat_k)^2  +  lam * mean_{k~l} (beta_k - beta_l)^2,
    c_k = S2_k / sum_j S2_j

where alpha_hat_k is the unpenalized estimate and S2_k = sum over frames of (t_zz,k^2 + t_tt,k^2) at
alpha = 1. The first term equals the stress misfit up to a constant (the data enter only through
per-patch sums), normalized so lam is dimensionless; k~l runs over patches sharing an edge.

lam is chosen by leave-one-protocol-out cross-validation: fit on eight of Def11-Def19, predict the
held-out protocol's global stresses from every patch, and take the lam with the lowest pooled RMS.
Maps are also shown for the lam values that cut the roughness (RMS of neighbor differences in log
alpha) to 50% and 10% of the unpenalized value, to show what stronger smoothing does.
Uncertainty at the chosen lam: protocol-level bootstrap.

--exclude-ring N leaves the outer N rings of patches (next to the tine tips) out of the fit and the
penalty entirely, so their extreme values cannot pull their neighbors.

    .venv/bin/python constitutive_fit/stage1_smooth.py --test N1739/test1                    -> local/stage1_smooth/
    .venv/bin/python constitutive_fit/stage1_smooth.py --test N1739/test1 --exclude-ring 1   -> local/stage1_smooth_interior1/
Outputs: stage1_smooth.csv / .npz / .json and stage1_smooth.png
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402
from scipy.sparse import coo_matrix, diags, vstack  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.dirname(HERE))                                              # constitutive_fit/
sys.path.insert(2, os.path.join(os.path.dirname(os.path.dirname(HERE)), "strain_tracking"))
from stage1_patches import LOCAL, load_inputs  # noqa: E402


def neighbor_pairs(mask):
    """Index pairs (into the flattened mask) of included patches that share an edge."""
    nt, nz = mask.shape
    idx = -np.ones(mask.shape, int)
    idx[mask] = np.arange(mask.sum())
    pairs = []
    for i in range(nt):
        for j in range(nz):
            if not mask[i, j]:
                continue
            for di, dj in ((0, 1), (1, 0)):
                if i + di < nt and j + dj < nz and mask[i + di, j + dj]:
                    pairs.append((idx[i, j], idx[i + di, j + dj]))
    return np.array(pairs)


def fit(s1, s2, pairs, lam, beta0=None):
    """Penalized estimate of beta = log alpha from per-patch sums s1 = sum(b*y), s2 = sum(b^2)."""
    ah = s1 / s2
    c = s2 / s2.sum()
    n, npair = len(s1), len(pairs)
    D = coo_matrix((np.r_[np.ones(npair), -np.ones(npair)],
                    (np.r_[np.arange(npair), np.arange(npair)], np.r_[pairs[:, 0], pairs[:, 1]])),
                   shape=(npair, n)).tocsr()
    w = np.sqrt(lam / npair)

    def res(b):
        return np.r_[np.sqrt(c) * (np.exp(b) - ah), w * (D @ b)]

    def jac(b):
        return vstack([diags(np.sqrt(c) * np.exp(b)), w * D]).tocsr()
    b0 = np.log(np.clip(ah, 1e-3, None)) if beta0 is None else beta0
    r = least_squares(res, b0, jac=jac, method="trf", x_scale="jac", max_nfev=2000)
    return r.x


def roughness(beta, pairs):
    return np.sqrt(np.mean((beta[pairs[:, 0]] - beta[pairs[:, 1]]) ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exclude-ring", type=int, default=0)
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--lams", default="-5,2,36", help="log10 range and count of the lambda grid")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    k_ring = a.exclude_ring
    a.out = a.out or os.path.join(LOCAL, "stage1_smooth" + (f"_interior{k_ring}" if k_ring else ""))
    os.makedirs(a.out, exist_ok=True)

    I = load_inputs()
    d, sets, tzz, ttt, bz, bt = I["d"], I["sets"], I["tzz"], I["ttt"], I["bz"], I["bt"]
    nt, nz = bz.shape[1:]
    ii, jj = np.meshgrid(np.arange(nt), np.arange(nz), indexing="ij")
    mask = (ii >= k_ring) & (ii < nt - k_ring) & (jj >= k_ring) & (jj < nz - k_ring)
    pairs = neighbor_pairs(mask)
    uset = np.unique(sets)
    yz, yt = tzz[:, None, None], ttt[:, None, None]
    # per-protocol sufficient statistics for the included patches: (protocols, n_patches)
    S1 = np.stack([(bz * yz + bt * yt)[sets == s][:, mask].sum(0) for s in uset])
    S2 = np.stack([(bz ** 2 + bt ** 2)[sets == s][:, mask].sum(0) for s in uset])
    S0 = np.stack([(tzz ** 2 + ttt ** 2)[sets == s].sum() for s in uset])
    nf = np.array([(sets == s).sum() for s in uset])
    npatch = mask.sum()
    lo_, hi_, n_ = (float(v) for v in a.lams.split(","))
    lams = np.r_[0.0, np.logspace(lo_, hi_, int(n_))]
    print(f"{npatch} patches, {len(pairs)} neighbor pairs, {len(uset)} protocols; lambda grid {len(lams)} values")

    # ---------------------------------------------------------------- full-data path and CV
    betas, prev = [], None
    for lam in lams:
        prev = fit(S1.sum(0), S2.sum(0), pairs, lam, prev)
        betas.append(prev)
    betas = np.array(betas)
    rough = np.array([roughness(b, pairs) for b in betas])
    ah = S1.sum(0) / S2.sum(0)
    c = S2.sum(0) / S2.sum()
    misfit = np.array([np.sqrt(np.sum(c * (np.exp(b) - ah) ** 2)) for b in betas])
    cv = np.zeros(len(lams))
    for p in range(len(uset)):
        keep = np.arange(len(uset)) != p
        prev = None
        for m, lam in enumerate(lams):
            prev = fit(S1[keep].sum(0), S2[keep].sum(0), pairs, lam, prev)
            al = np.exp(prev)
            cv[m] += np.sum(al ** 2 * S2[p] - 2 * al * S1[p] + S0[p])   # held-out squared error, all patches
    cv_rms = np.sqrt(cv / (npatch * 2 * nf.sum()))
    m_cv = int(np.argmin(cv_rms))
    m50 = int(np.argmin(np.abs(rough - 0.5 * rough[0])))
    m10 = int(np.argmin(np.abs(rough - 0.1 * rough[0])))
    print(f"CV RMS (held-out protocol, all patches): unpenalized {cv_rms[0]:.3f} kPa; minimum {cv_rms[m_cv]:.3f} kPa "
          f"at lambda = {lams[m_cv]:.3g}")
    for name, m in (("unpenalized", 0), ("CV choice", m_cv), ("roughness 50%", m50), ("roughness 10%", m10)):
        al = np.exp(betas[m])
        print(f"  {name:14s} lambda {lams[m]:9.3g}: alpha {al.min():.2f}-{al.max():.2f}, roughness {rough[m]:.3f}, "
              f"CV RMS {cv_rms[m]:.3f} kPa")

    # ---------------------------------------------------------------- bootstrap at the CV choice
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(a.boot):
        w = np.bincount(rng.integers(0, len(uset), len(uset)), minlength=len(uset)).astype(float)
        boots.append(np.exp(fit(w @ S1, w @ S2, pairs, lams[m_cv], betas[m_cv])))
    lo, hi = np.percentile(np.array(boots), [2.5, 97.5], axis=0)

    # ---------------------------------------------------------------- save
    def grid(v):
        g = np.full((nt, nz), np.nan)
        g[mask] = v
        return g
    A0, Acv, A50, A10 = (grid(np.exp(betas[m])) for m in (0, m_cv, m50, m10))
    ez, et = d["patch_edges_z_mm"], d["patch_edges_theta_mm"]
    rows = [dict(i_theta=i, j_z=j, z0_mm=ez[j], z1_mm=ez[j + 1], theta0_mm=et[i], theta1_mm=et[i + 1],
                 alpha_unpenalized=A0[i, j], alpha_cv=Acv[i, j], alpha_cv_lo95=grid(lo)[i, j],
                 alpha_cv_hi95=grid(hi)[i, j], alpha_rough50=A50[i, j], alpha_rough10=A10[i, j])
            for i in range(nt) for j in range(nz) if mask[i, j]]
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "stage1_smooth.csv"), index=False)
    np.savez_compressed(os.path.join(a.out, "stage1_smooth.npz"), lambdas=lams, alpha_path=np.exp(betas),
                        mask=mask, roughness=rough, misfit=misfit, cv_rms_kPa=cv_rms, lambda_cv=lams[m_cv],
                        alpha_cv=Acv, lo95=grid(lo), hi95=grid(hi),
                        layout="alpha_path columns follow the patches where mask is True, row-major")
    with open(os.path.join(a.out, "stage1_smooth.json"), "w") as f:
        json.dump(dict(description="Stage 1 with a Tikhonov penalty on neighbor differences of log alpha; "
                       "lambda by leave-one-protocol-out cross-validation",
                       exclude_ring=k_ring, n_patches=int(npatch), n_pairs=int(len(pairs)),
                       lambda_cv=float(lams[m_cv]), cv_rms_kPa=dict(unpenalized=float(cv_rms[0]),
                                                                      chosen=float(cv_rms[m_cv])),
                       lambda_rough50=float(lams[m50]), lambda_rough10=float(lams[m10]),
                       roughness=dict(unpenalized=float(rough[0]), cv=float(rough[m_cv]))), f, indent=2)

    # ---------------------------------------------------------------- figure
    zc, tc = 0.5 * (ez[1:] + ez[:-1]), 0.5 * (et[1:] + et[:-1])
    allv = np.concatenate([A0[mask], Acv[mask], A50[mask], A10[mask]])
    mlog = np.nanpercentile(np.abs(np.log(allv)), 98)
    norm = LogNorm(np.exp(-mlog), np.exp(mlog))
    ticks = [v for v in (0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5) if norm.vmin <= v <= norm.vmax]
    fig = plt.figure(figsize=(21, 11), constrained_layout=True)
    gs = fig.add_gridspec(2, 4)
    a_ = fig.add_subplot(gs[0, 0])
    a_.semilogx(lams[1:], cv_rms[1:], "o-", ms=3)
    a_.axhline(cv_rms[0], color="0.5", ls="--", label=f"unpenalized {cv_rms[0]:.3f} kPa")
    a_.plot(lams[m_cv], cv_rms[m_cv], "r*", ms=14, label=f"minimum, λ = {lams[m_cv]:.2g}")
    a_.set_xlabel("λ"); a_.set_ylabel("held-out protocol RMS (kPa)"); a_.grid(alpha=.3); a_.legend(fontsize=8)
    a_.set_title("Leave-one-protocol-out cross-validation", fontsize=10)
    a_ = fig.add_subplot(gs[1, 0])
    a_.plot(rough, misfit, "o-", ms=3)
    for name, m, col in (("CV", m_cv, "r"), ("50%", m50, "C2"), ("10%", m10, "C1")):
        a_.plot(rough[m], misfit[m], "*", ms=13, color=col, label=f"{name}: λ = {lams[m]:.2g}")
    a_.set_xlabel("roughness: RMS neighbor difference of log α"); a_.set_ylabel("weighted RMS of α − α̂")
    a_.grid(alpha=.3); a_.legend(fontsize=8); a_.set_title("Trade-off (L-curve)", fontsize=10)
    panels = [(A0, "unpenalized (λ = 0)"), (Acv, f"cross-validated λ = {lams[m_cv]:.2g}"),
              (A50, f"λ = {lams[m50]:.2g}: roughness 50%"), (A10, f"λ = {lams[m10]:.2g}: roughness 10%")]
    for k_, (V, title) in enumerate(panels):
        a_ = fig.add_subplot(gs[k_ // 2, 1 + k_ % 2] if k_ < 4 else gs[0, 3])
        a_.set_facecolor("0.85")
        m = a_.pcolormesh(ez, et, V, cmap="RdBu_r", norm=norm)
        for i in range(nt):
            for j in range(nz):
                if np.isfinite(V[i, j]):
                    r, g, b_, _ = m.cmap(m.norm(V[i, j]))
                    a_.text(zc[j], tc[i], f"{V[i, j]:.2f}", ha="center", va="center", fontsize=6,
                            color="w" if 0.299 * r + 0.587 * g + 0.114 * b_ < 0.45 else "k")
        cb = fig.colorbar(m, ax=a_, shrink=.85)
        cb.set_ticks(ticks); cb.set_ticklabels([f"{v:g}" for v in ticks]); cb.minorticks_off()
        a_.set_xlim(0, ez[-1]); a_.set_ylim(0, et[-1]); a_.set_aspect("equal")
        a_.set_title(f"α, {title}", fontsize=10)
        a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    a_ = fig.add_subplot(gs[:, 3])
    a_.set_facecolor("0.85")
    wid = 100 * (grid(hi) - grid(lo)) / Acv
    m = a_.pcolormesh(ez, et, wid, cmap="magma_r", vmin=0, vmax=np.nanmax(wid))
    for i in range(nt):
        for j in range(nz):
            if np.isfinite(wid[i, j]):
                a_.text(zc[j], tc[i], f"{wid[i, j]:.0f}", ha="center", va="center", fontsize=6)
    fig.colorbar(m, ax=a_, shrink=.6)
    a_.set_xlim(0, ez[-1]); a_.set_ylim(0, et[-1]); a_.set_aspect("equal")
    a_.set_title(f"cross-validated α: protocol-bootstrap\n95% width (% of α)", fontsize=10)
    a_.set_xlabel("z from left tine tips (mm)")
    fig.suptitle(f"Stage 1 with a smoothness penalty on log α ({nz} × {nt} patches"
                 + (f", outer {k_ring} ring(s) left out of the fit" if k_ring else "") + "); uniform-stress assumption")
    fig.savefig(os.path.join(a.out, "stage1_smooth.png"), dpi=85)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
