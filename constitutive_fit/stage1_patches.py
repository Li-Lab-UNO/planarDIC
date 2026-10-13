"""
Stage 1 of the local estimation (local_estimation_handoff.md): a quick map of relative patch
stiffness from the patch-mean F and the GLOBAL stress (uniform-stress assumption).

For each patch k, the patch-mean deformation gradient F_k (area-weighted over the patch's grid
cells, shear included) goes through the full-domain bulk baseline (dic_affine parameters P).
One factor alpha_k scales the stiffness parameters (Cgr, C1_el, C1_smc, C1_col) with the
exponents and gamma fixed (ladder level (a)). The model is linear in those four, so

    t_model,k = alpha_k * t(F_k; P)

and alpha_k has a closed-form least-squares solution against the global Cauchy stresses
t_zz = F_z lambda_z / (w_theta H), t_thetatheta = F_theta lambda_theta / (w_z H) of the same
frames (all 1492 loading frames of Def11-Def19, tared forces, from the full-domain baseline).

alpha_k is a membrane-stiffness factor (H_k / H x modulus scale) relative to the bulk baseline:
> 1 stiffer, < 1 softer. Uniform stress is a reasonable approximation in theta (the theta
variation acts mostly in series) but biased in z, where a stiff region carries more than the
mean stress; so alpha_z (fit to t_zz only) is expected to underestimate the contrast. alpha_z and
alpha_theta (fit to t_thetatheta only) are reported as diagnostics. Not a final result.

Uncertainty: protocol-level bootstrap (resample whole protocols with replacement).

Each alpha_k is fitted independently, so --exclude-ring N leaves every interior alpha unchanged.
It drops the outer N rings of patches (next to the tine tips) from the color scales, the summary
statistics, the diagnostics and the example fits, so the interior pattern can be read on its own.

    ../.venv/bin/python stage1_patches.py --test N1739/test1                    -> local/stage1/
    ../.venv/bin/python stage1_patches.py --test N1739/test1 --exclude-ring 1   -> local/stage1_interior1/
Outputs: Images/constitutive_fit/local/stage1[_interiorN]/
    stage1_patches.csv     per-patch alpha, alpha_z, alpha_theta, 95% intervals, RMS, coverage
    stage1_results.npz     the same as arrays [theta row, z column], plus the patch-mean F
    stage1_maps.png        alpha, alpha_theta, alpha_z maps, interval width, RMS, alpha_z / alpha_theta
    stage1_diagnostics.png alpha_z vs alpha_theta, alpha vs the strain concentrations
    stage1_fits.png        data and fitted model for the softest, median and stiffest patches
"""
import argparse
import json
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import baseline as B  # noqa: E402
import constitutive as M  # noqa: E402

IMAGES = testcfg.get().images
LOCAL = os.path.join(IMAGES, "constitutive_fit", "local")


def load_inputs(inp="dic_affine"):
    """Patch-mean F, the baseline model stress at it, and the global stresses for every loading frame.
    Returns a dict with d (local_data.npz), base, df, P, sets, tzz, ttt (frames,), Fp (frames, nt, nz, 2, 2),
    bz, bt (frames, nt, nz): t_zz and t_thetatheta of the bulk baseline model at the patch-mean F."""
    d = np.load(os.path.join(LOCAL, "local_data.npz"))
    base, df = B.load(str(d["baseline"]))
    if base["checksums"]["data"]["baseline_data.csv"] != str(d["baseline_data_sha"]):
        sys.exit("local_data.npz was built from a different baseline data set; rerun build_local_data.py")
    st = base["settings"]
    H, w_t, w_z = st["H_mm"], st["w_theta_mm"], st["w_z_mm"]
    P = np.array([base["inputs"][inp]["params"][n] for n in M.PARAM_NAMES])
    tzz = d["force_z_mN"] * df[f"{inp}_z"].to_numpy() / (w_t * H)       # kPa
    ttt = d["force_theta_mN"] * df[f"{inp}_t"].to_numpy() / (w_z * H)
    F = d["F"].astype(np.float64)                                   # (frames, ny, nx, 2, 2)
    W = d["patch_weights_mm2"].astype(np.float64)                   # (nt, nz, ny, nx)
    fin = np.isfinite(F[..., 0, 0])
    Fp = np.einsum("abrc,frcij->fabij", W, np.where(fin[..., None, None], F, 0))
    Fp /= np.einsum("abrc,frc->fab", W, fin.astype(float))[..., None, None]
    tb = M.cauchy(Fp, P)                                            # baseline stress at the patch F
    return dict(d=d, base=base, df=df, P=P, sets=d["set"], tzz=tzz, ttt=ttt, Fp=Fp, W=W,
                bz=tb[..., 0, 0], bt=tb[..., 1, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=200, help="protocol-level bootstrap resamples")
    ap.add_argument("--input", default="dic_affine", help="baseline stretch input supplying P and the global λ")
    ap.add_argument("--exclude-ring", type=int, default=0,
                    help="leave out the outer N rings of patches from scales, statistics and diagnostics")
    ap.add_argument("--out", default=None, help="default local/stage1, or local/stage1_interiorN")
    a = ap.parse_args()
    k_ring = a.exclude_ring
    a.out = a.out or os.path.join(LOCAL, "stage1" + (f"_interior{k_ring}" if k_ring else ""))
    os.makedirs(a.out, exist_ok=True)
    warnings.filterwarnings("ignore", "Mean of empty slice")
    np.seterr(invalid="ignore", divide="ignore")

    I = load_inputs(a.input)
    d, base, df, P, sets = I["d"], I["base"], I["df"], I["P"], I["sets"]
    tzz, ttt, Fp, W, bz, bt = I["tzz"], I["ttt"], I["Fp"], I["W"], I["bz"], I["bt"]
    nt, nz = W.shape[:2]
    print(f"patch-mean F: {Fp.shape[0]} frames x {nt} x {nz} patches")
    yz, yt = tzz[:, None, None], ttt[:, None, None]

    # per-protocol sums, so every estimate (and each bootstrap resample) is a ratio of sums
    uset = np.unique(sets)
    S = {k: np.stack([v[sets == s].sum(0) for s in uset]) for k, v in
         dict(zz=bz * yz, z2=bz ** 2, tt=bt * yt, t2=bt ** 2).items()}

    def estimates(w):                   # w: weight per protocol (bootstrap counts)
        g = {k: np.tensordot(w, v, 1) for k, v in S.items()}
        return ((g["zz"] + g["tt"]) / (g["z2"] + g["t2"]), g["zz"] / g["z2"], g["tt"] / g["t2"])
    alpha, alpha_z, alpha_t = estimates(np.ones(len(uset)))
    rng = np.random.default_rng(0)
    boots = np.array([estimates(np.bincount(rng.integers(0, len(uset), len(uset)), minlength=len(uset)))
                      for _ in range(a.boot)])                      # (boot, 3, nt, nz)
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)

    rz, rt = alpha * bz - yz, alpha * bt - yt
    rms = np.sqrt(np.mean(np.r_[rz, rt] ** 2, 0))
    rms1 = np.sqrt(np.mean(np.r_[bz - yz, bt - yt] ** 2, 0))        # alpha = 1 (bulk parameters)
    y = np.r_[tzz, ttt]
    r2 = 1 - np.sum(np.r_[rz, rt] ** 2, 0) / np.sum((y - y.mean()) ** 2)
    area = W.sum((2, 3))
    ii, jj = np.meshgrid(np.arange(nt), np.arange(nz), indexing="ij")
    inner = (ii >= k_ring) & (ii < nt - k_ring) & (jj >= k_ring) & (jj < nz - k_ring)
    def M_(V):                          # interior only (outer rings -> NaN) for display and statistics
        return np.where(inner, V, np.nan)
    a_mean = np.sum((alpha * area)[inner]) / area[inner].sum()
    scope = f"interior {inner.sum()} patches (outer {k_ring} ring(s) left out)" if k_ring else f"all {inner.size} patches"
    print(f"{scope}:")
    print(f"  alpha: {np.nanmin(M_(alpha)):.2f}-{np.nanmax(M_(alpha)):.2f} (area-weighted mean {a_mean:.3f}); "
          f"alpha_theta {np.nanmin(M_(alpha_t)):.2f}-{np.nanmax(M_(alpha_t)):.2f}; alpha_z {np.nanmin(M_(alpha_z)):.2f}-"
          f"{np.nanmax(M_(alpha_z)):.2f}")
    print(f"  patch RMS {np.nanmin(M_(rms)):.2f}-{np.nanmax(M_(rms)):.2f} kPa (median {np.nanmedian(M_(rms)):.2f}); "
          f"with alpha = 1: median {np.nanmedian(M_(rms1)):.2f} kPa; bulk baseline {base['inputs'][a.input]['rms_kPa']:.3f}")
    sig = ((lo[0] > 1) | (hi[0] < 1)) & inner
    print(f"  patches whose alpha 95% interval excludes 1: {sig.sum()} of {inner.sum()}")
    cc = lambda x, y_: np.corrcoef(x[inner], y_[inner])[0, 1]
    print(f"  r(alpha_z, alpha_theta) {cc(alpha_z, alpha_t):.2f}; r(alpha, E_thetatheta conc) "
          f"{cc(alpha, d['patch_conc_theta']):.2f}; r(alpha, E_zz conc) {cc(alpha, d['patch_conc_z']):.2f}")

    # ---------------------------------------------------------------- save
    ez, et = d["patch_edges_z_mm"], d["patch_edges_theta_mm"]
    rows = []
    for i in range(nt):
        for j in range(nz):
            rows.append(dict(i_theta=i, j_z=j, z0_mm=ez[j], z1_mm=ez[j + 1], theta0_mm=et[i], theta1_mm=et[i + 1],
                             coverage=d["patch_coverage"][i, j], interior=bool(inner[i, j]),
                             alpha=alpha[i, j], alpha_lo95=lo[0, i, j], alpha_hi95=hi[0, i, j],
                             alpha_theta=alpha_t[i, j], alpha_theta_lo95=lo[2, i, j], alpha_theta_hi95=hi[2, i, j],
                             alpha_z=alpha_z[i, j], alpha_z_lo95=lo[1, i, j], alpha_z_hi95=hi[1, i, j],
                             rms_kPa=rms[i, j], rms_alpha1_kPa=rms1[i, j], R2=r2[i, j],
                             conc_theta=d["patch_conc_theta"][i, j], conc_z=d["patch_conc_z"][i, j]))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(a.out, "stage1_patches.csv"), index=False)
    np.savez_compressed(os.path.join(a.out, "stage1_results.npz"), alpha=alpha, alpha_z=alpha_z, alpha_theta=alpha_t,
                        lo95=lo, hi95=hi, boot_order="[alpha, alpha_z, alpha_theta]", rms_kPa=rms,
                        rms_alpha1_kPa=rms1, R2=r2, F_patch=Fp.astype(np.float32),
                        patch_edges_z_mm=ez, patch_edges_theta_mm=et, params_baseline=P,
                        layout="[theta row from the bottom tine tips, z column from the left tine tips]")
    with open(os.path.join(a.out, "stage1_summary.json"), "w") as f:
        json.dump(dict(description="Stage 1: one stiffness factor per patch (ladder (a)) from the patch-mean F "
                       "and the global stress (uniform-stress assumption)",
                       baseline=str(d["baseline"]), input=a.input, params=dict(zip(M.PARAM_NAMES, P.tolist())),
                       n_frames=int(len(sets)), patches=[int(nz), int(nt)], bootstrap=a.boot,
                       exclude_ring=k_ring, statistics_over=scope,
                       alpha_range=[float(np.nanmin(M_(alpha))), float(np.nanmax(M_(alpha)))],
                       alpha_area_weighted_mean=float(a_mean), n_alpha_interval_excludes_1=int(sig.sum()),
                       rms_median_kPa=float(np.nanmedian(M_(rms))),
                       rms_alpha1_median_kPa=float(np.nanmedian(M_(rms1))),
                       r_alpha_z_alpha_theta=float(cc(alpha_z, alpha_t)),
                       r_alpha_conc_theta=float(cc(alpha, d["patch_conc_theta"])),
                       r_alpha_conc_z=float(cc(alpha, d["patch_conc_z"]))),
                  f, indent=2)

    # ---------------------------------------------------------------- figures
    Wz, Wt = ez[-1], et[-1]
    zc, tc = 0.5 * (ez[1:] + ez[:-1]), 0.5 * (et[1:] + et[:-1])

    def pmap(ax_, V, title, cmap, norm, fmt="{:.2f}"):
        m = ax_.pcolormesh(ez, et, V, cmap=cmap, norm=norm)
        ax_.set_facecolor("0.85")                                   # left-out patches show gray
        for i in range(nt):
            for j in range(nz):
                if not np.isfinite(V[i, j]):
                    continue
                r, g, b_, _ = m.cmap(m.norm(V[i, j]))
                ax_.text(zc[j], tc[i], fmt.format(V[i, j]), ha="center", va="center", fontsize=6,
                         color="w" if 0.299 * r + 0.587 * g + 0.114 * b_ < 0.45 else "k")
        ax_.set_xlim(0, Wz); ax_.set_ylim(0, Wt); ax_.set_aspect("equal")
        ax_.set_title(title, fontsize=10); ax_.set_xlabel("z from left tine tips (mm)")
        cb = fig.colorbar(m, ax=ax_, shrink=.85)
        if isinstance(norm, LogNorm):
            t_ = [v for v in (0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10) if norm.vmin <= v <= norm.vmax]
            cb.set_ticks(t_); cb.set_ticklabels([f"{v:g}" for v in t_]); cb.minorticks_off()

    def log_sym(V):                     # log color scale centered on 1, limits from the 2nd-98th percentile
        m = np.nanpercentile(np.abs(np.log(np.asarray(V).ravel())), 98)
        return LogNorm(np.exp(-m), np.exp(m))
    anorm = log_sym([M_(alpha), M_(alpha_t), M_(alpha_z)])          # common limits for the three alpha maps
    fig, ax = plt.subplots(2, 3, figsize=(19, 13), constrained_layout=True)
    pmap(ax[0, 0], M_(alpha), "α: stiffness factor vs bulk baseline\n(fit to t_zz and t_θθ)", "RdBu_r", anorm)
    pmap(ax[0, 1], M_(alpha_t), "α_θ: fit to t_θθ only\n(uniform stress reasonable in θ)", "RdBu_r", anorm)
    pmap(ax[0, 2], M_(alpha_z), "α_z: fit to t_zz only\n(uniform stress biased in z)", "RdBu_r", anorm)
    wid = M_(100 * (hi[0] - lo[0]) / alpha)
    pmap(ax[1, 0], wid, "α: protocol-bootstrap 95% width (% of α)", "magma_r",
         plt.Normalize(0, np.nanmax(wid)), "{:.0f}")
    pmap(ax[1, 1], M_(rms), f"RMS of the α fit (kPa)\n(bulk baseline {base['inputs'][a.input]['rms_kPa']:.2f} kPa)",
         "magma_r", plt.Normalize(0, np.nanmax(M_(rms))))
    rr = M_(alpha_z / alpha_t)
    pmap(ax[1, 2], rr, "α_z / α_θ\n(1 = both directions agree)", "PuOr_r", log_sym(rr))
    for a_ in ax[:, 0]:
        a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    fig.suptitle(f"Stage 1: one stiffness factor per patch ({nz} × {nt}) from the patch-mean F and the global stress "
                 f"(uniform-stress assumption); full-domain baseline, {len(sets)} frames, Def11–Def19"
                 + (f"\nouter {k_ring} ring(s) of patches next to the tine tips left out (gray)" if k_ring else ""))
    fig.savefig(os.path.join(a.out, "stage1_maps.png"), dpi=85)

    fig, ax = plt.subplots(1, 3, figsize=(18, 5.6), constrained_layout=True)
    rowc = plt.get_cmap("viridis")(ii[inner] / max(nt - 1, 1))
    for a_, (x, yv, xl, yl) in zip(ax, ((alpha_t, alpha_z, "α_θ", "α_z"),
                                        (d["patch_conc_theta"], alpha, "E_θθ concentration", "α"),
                                        (d["patch_conc_z"], alpha, "E_zz concentration", "α"))):
        a_.scatter(x[inner], yv[inner], c=rowc, s=22, edgecolor="k", lw=.3)
        a_.set_xlabel(xl); a_.set_ylabel(yl); a_.grid(alpha=.3)
        a_.set_title(f"{yl} vs {xl} (r = {cc(x, yv):.2f})", fontsize=10)
    lim = [min(np.nanmin(M_(alpha_t)), np.nanmin(M_(alpha_z))), max(np.nanmax(M_(alpha_t)), np.nanmax(M_(alpha_z)))]
    ax[0].plot(lim, lim, "k", lw=.8)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, nt - 1))
    fig.colorbar(sm, ax=ax[2], label="patch θ row (from the bottom)")
    fig.suptitle(f"Stage 1 diagnostics (one point per patch; {scope})")
    fig.savefig(os.path.join(a.out, "stage1_diagnostics.png"), dpi=85)

    order = np.argsort(np.where(inner, alpha, np.nan).ravel())[:inner.sum()]   # interior patches, by alpha
    picks = [order[0], order[len(order) // 2], order[-1]]
    fig, ax = plt.subplots(len(picks), 2, figsize=(13, 4 * len(picks)), constrained_layout=True)
    cols = {s: plt.get_cmap("tab10")(n) for n, s in enumerate(uset)}
    for r_, p in enumerate(picks):
        i, j = divmod(int(p), nz)
        C = np.einsum("fki,fkj->fij", Fp[:, i, j], Fp[:, i, j])
        lam = (np.sqrt(C[:, 0, 0]), np.sqrt(C[:, 1, 1]))
        for c_, (y_, b_, lab) in enumerate(((tzz, bz, "zz"), (ttt, bt, "θθ"))):
            for s in uset:
                m_ = sets == s
                ax[r_, c_].plot(lam[c_][m_], y_[m_], ".", ms=2, color=cols[s],
                                label=str(df.folder[m_].iloc[0]).split("_")[0] if r_ == 0 and c_ == 0 else None)
                ax[r_, c_].plot(lam[c_][m_], alpha[i, j] * b_[m_, i, j], "-", lw=1, color=cols[s])
            ax[r_, c_].set_xlabel(f"patch λ_{'z' if c_ == 0 else 'θ'}"); ax[r_, c_].set_ylabel(f"t_{lab} (kPa)")
            ax[r_, c_].grid(alpha=.3)
            ax[r_, c_].set_title(f"patch ({i},{j}), α = {alpha[i, j]:.2f} [{lo[0, i, j]:.2f}, {hi[0, i, j]:.2f}], "
                                 f"RMS {rms[i, j]:.2f} kPa", fontsize=10)
    ax[0, 0].legend(fontsize=7, markerscale=4, ncol=3)
    fig.suptitle("Softest, median and stiffest patch" + (" in the interior" if k_ring else "")
                 + ": global stress (dots) vs α × baseline model at the patch F (lines)")
    fig.savefig(os.path.join(a.out, "stage1_fits.png"), dpi=85)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
