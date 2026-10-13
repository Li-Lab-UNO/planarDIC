"""
Stage 2 fit on the measured data: interior patch stiffness factors from the virtual-field equations.

Model and equations as in stage2_setup.py (ladder (a): alpha_k scales Cgr, C1_el, C1_smc, C1_col of the
full-domain bulk baseline, exponents and gamma fixed; interior patches; cut_z, cut_theta and virtual-field
blocks). Settings follow the synthetic test (stage2_synthetic.py): F from the tracks with a 45 px window
(--win 1), 8 x 8 sine modes per direction.

Methods
  combined     (primary) all blocks at once, rows scaled to unit norm and each block to equal total weight,
               plus (lam / n_pairs) sum_{k~l} (alpha_k - alpha_l)^2; lam by leave-one-protocol-out CV
  shape+scale  the alpha pattern from the virtual fields alone (normalized to mean 1, same penalty), one scale
               factor from all cut equations
Uncertainty: protocol-level bootstrap (resample protocols; per-protocol sums make each resample a 64 x 64 solve).
Synthetic test limits: spurious patch contrast about +-0.07 RMS for a uniform material; overall level about +-5%
(the z and theta cuts disagree by about +-7% because of the uniform-flux assumption).

Checks
  * residual of each block at alpha = 1, at the Stage 1 alpha and at the fitted alpha
  * bulk score: the area-weighted mean alpha applied to the baseline, scored with baseline.evaluate
  * forward check: FE model of the tracked domain (membrane_fe.py; boundary 2 nodes inside the tracked edge,
    smoothed measured boundary displacements) with the fitted alpha (ring patches alpha = 1) vs with alpha = 1
    everywhere; interior node positions compared with the smoothed measured positions
  * window sensitivity: the combined fit repeated with the 75 px window

    ../.venv/bin/python stage2_fit.py --test N1739/test1 [--every 4] [--fe-every 32]
Outputs: Images/constitutive_fit/local/stage2_fit/
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
from stage2_synthetic import cut_scale, difference_matrix, neighbor_pairs, unit_rows, vf_shape  # noqa: E402

FLOOR, LEVEL = 0.07, 0.05


def combined_system(Gz, bz, Gt, bt, Gv, sets, npat):
    """Per-protocol normal-equation pieces of the row-scaled combined system."""
    rows, rhs, rsets = [], [], []
    for G, b in ((Gz, bz), (Gt, bt), (Gv, np.zeros(Gv.shape[:2]))):
        G2, b2 = G.reshape(-1, npat), b.ravel()
        nr = np.linalg.norm(G2, axis=1)
        kk = nr > 0
        sc = 1 / (nr[kk] * np.sqrt(kk.sum()))
        rows.append(G2[kk] * sc[:, None]); rhs.append(b2[kk] * sc); rsets.append(np.repeat(sets, G.shape[1])[kk])
    A, b, rs = np.vstack(rows), np.concatenate(rhs), np.concatenate(rsets)
    uset = np.unique(sets)
    return ({s: A[rs == s].T @ A[rs == s] for s in uset}, {s: A[rs == s].T @ b[rs == s] for s in uset},
            {s: b[rs == s] @ b[rs == s] for s in uset})


def fit_combined(AtA, Atb, btb, DtD, lams):
    uset = list(AtA)
    AA, Ab = sum(AtA.values()), sum(Atb.values())
    cv = []
    for lam in lams:
        e = 0.0
        for s in uset:
            al = np.linalg.solve(AA - AtA[s] + lam * DtD, Ab - Atb[s])
            e += al @ AtA[s] @ al - 2 * al @ Atb[s] + btb[s]
        cv.append(e)
    cv = np.array(cv)
    lam = float(lams[int(np.argmin(cv))])
    return np.linalg.solve(AA + lam * DtD, Ab), lam, cv


def block_residuals(al, Gz, bz, Gt, bt, Gv):
    out = {}
    for name, G, b in (("cut_z", Gz, bz), ("cut_theta", Gt, bt)):
        r = G @ al - b
        out[name] = float(np.sqrt(np.mean(r ** 2)) / np.sqrt(np.mean(b ** 2)))
    out["vf"] = float(np.sqrt(np.mean((Gv @ al) ** 2)) / np.sqrt(np.mean((np.abs(Gv) @ np.abs(al)) ** 2)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--every", type=int, default=4)
    ap.add_argument("--win", type=int, default=1)
    ap.add_argument("--modes", type=int, default=8)
    ap.add_argument("--exclude-ring", type=int, default=1)
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--fe-every", type=int, default=32, help="frames for the forward FE check (every n-th)")
    ap.add_argument("--out", default=os.path.join(LOCAL, "stage2_fit"))
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
    DtD = difference_matrix(neighbor_pairs(inner), npat)
    lams = np.r_[0.0, np.logspace(-6, 4, 41)]
    Xg = np.stack([d["z_mm"], d["theta_mm"]], -1)
    fidx = np.arange(0, len(df), a.every)
    sets = df.set.to_numpy()[fidx]
    uset = np.unique(sets)
    Fz, Ft = d["force_z_mN"][fidx], d["force_theta_mN"][fidx]
    xm, vm = measured_positions(df, fidx, mpp, with_valid=True)
    print(f"{len(fidx)} frames, {npat} patches, {2 * a.modes ** 2} virtual fields per frame")

    fits = {}
    for win in (a.win, 2 if a.win != 2 else 1):
        F = window_F(Xg, xm, vm, win)
        Gz, Gt, Gv, bz, bt = assemble(g, F, Fz, Ft, P, H, w_t, w_z)
        AtA, Atb, btb = combined_system(Gz, bz, Gt, bt, Gv, sets, npat)
        al, lam, cv = fit_combined(AtA, Atb, btb, DtD, lams)
        fits[win] = dict(alpha=al, lam=lam, cv=cv, G=(Gz, bz, Gt, bt, Gv), sys=(AtA, Atb, btb))
        print(f"  window {(2 * win + 1) * 15} px: CV lambda {lam:.3g}, alpha {al.min():.2f}-{al.max():.2f}")

    W0 = fits[a.win]
    Gz, bz, Gt, bt, Gv = W0["G"]
    AtA, Atb, btb = W0["sys"]
    al, lam = W0["alpha"], W0["lam"]
    area = d["patch_weights_mm2"].sum((2, 3))[inner]
    a_mean = float(np.average(al, weights=area))

    # shape + scale
    Av, keep = unit_rows(Gv)
    rs = np.repeat(sets, Gv.shape[1])[keep]
    Kp = {s: Av[rs == s].T @ Av[rs == s] for s in uset}
    K = sum(Kp.values())
    cvs = np.array([sum(vf_shape(K - Kp[s], l_, DtD) @ Kp[s] @ vf_shape(K - Kp[s], l_, DtD) for s in uset) for l_ in lams])
    lam_s = float(lams[int(np.argmin(cvs))])
    shp = vf_shape(K, lam_s, DtD)
    s_all, s_z, s_t = cut_scale(shp, Gz, bz, Gt, bt)
    al_ss = s_all * shp
    print(f"combined: lambda {lam:.3g}, alpha {al.min():.2f}-{al.max():.2f}, mean {a_mean:.3f}; shape+scale: lambda "
          f"{lam_s:.3g}, scale {s_all:.3f} (z cuts {s_z:.3f}, theta cuts {s_t:.3f}), r(combined, shape+scale) "
          f"{np.corrcoef(al, al_ss)[0, 1]:.3f}")

    # bootstrap (combined, at the chosen lambda)
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(a.boot):
        w = np.bincount(rng.integers(0, len(uset), len(uset)), minlength=len(uset))
        AA = sum(w[i] * AtA[s] for i, s in enumerate(uset))
        Ab = sum(w[i] * Atb[s] for i, s in enumerate(uset))
        boots.append(np.linalg.solve(AA + lam * DtD, Ab))
    boots = np.array(boots)
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    rel = boots / boots.mean(1, keepdims=True)                     # pattern relative to each resample's mean
    rlo, rhi = np.percentile(rel, [2.5, 97.5], axis=0)
    # a patch differs from the mean when its relative interval excludes 1 by more than the synthetic floor
    differs = (rlo > 1 + FLOOR) | (rhi < 1 - FLOOR)
    print(f"bootstrap: median 95% width {np.median(hi - lo):.3f}; patches whose relative alpha differs from the mean by "
          f"more than the +-{FLOOR} floor (95%): {differs.sum()} of {npat}")

    # comparisons
    s1 = pd.read_csv(os.path.join(LOCAL, f"stage1_interior{a.exclude_ring}", "stage1_patches.csv"))
    a1 = np.full((nt, nz), np.nan)
    a1[s1.i_theta, s1.j_z] = s1.alpha
    at1 = np.full((nt, nz), np.nan)
    at1[s1.i_theta, s1.j_z] = s1.alpha_theta
    az1 = np.full((nt, nz), np.nan)
    az1[s1.i_theta, s1.j_z] = s1.alpha_z
    corr = {k: float(np.corrcoef(al, v[inner])[0, 1]) for k, v in (("stage1_alpha", a1), ("stage1_alpha_theta", at1),
                                                                    ("stage1_alpha_z", az1))}
    corr["window_75px_vs_45px"] = float(np.corrcoef(fits[a.win]["alpha"], fits[2 if a.win != 2 else 1]["alpha"])[0, 1])
    res = {name: block_residuals(v, Gz, bz, Gt, bt, Gv) for name, v in
           (("alpha = 1", np.ones(npat)), ("Stage 1 alpha", a1[inner]), ("Stage 2 fit", al))}
    for k, v in res.items():
        print(f"  residual, {k:14s}: " + ", ".join(f"{b_} {x:.3f}" for b_, x in v.items()))
    print("  correlations: " + ", ".join(f"{k} {v:.2f}" for k, v in corr.items()))
    Pm = P.copy()
    Pm[[0, 1, 3, 5]] *= a_mean
    bulk = B.evaluate(list(Pm), "dic_affine", verbose=False, tracks="tine_region")
    print(f"  bulk score with the area-weighted mean alpha {a_mean:.3f}: RMS {bulk['rms_kPa']:.3f} kPa "
          f"(baseline {base['inputs']['dic_affine']['rms_kPa']:.3f})")

    # ---------------------------------------------------------------- forward FE check
    m_ = 2
    sl = (slice(m_, -m_), slice(m_, -m_))
    mesh = Mesh(Xg[sl])
    gp = mesh.gauss_points
    gi = np.clip(np.searchsorted(pt, gp[..., 1]) - 1, 0, nt - 1)
    gj = np.clip(np.searchsorted(pz, gp[..., 0]) - 1, 0, nz - 1)
    A_fit = np.ones((nt, nz))
    A_fit[inner] = al
    fe_frames = np.arange(0, len(df), a.fe_every)
    xf = measured_positions(df, fe_frames, mpp)
    xs = Xg + np.stack([gaussian_filter(xf[..., i] - Xg[..., i], (0, 1.0, 1.0), mode="nearest") for i in range(2)], -1)
    dom = (Xg[..., 0] > g["domain"][0]) & (Xg[..., 0] < g["domain"][1]) & \
          (Xg[..., 1] > g["domain"][2]) & (Xg[..., 1] < g["domain"][3])
    err, t0 = {}, time.time()
    for name, Amap in (("alpha = 1", np.ones((nt, nz))), ("Stage 2 fit", A_fit)):
        e = []
        for n in range(len(fe_frames)):
            x, _ = mesh.solve(xs[n][sl], Amap[gi, gj], P, H, xs[n][sl])
            full = xs[n].copy()
            full[sl] = x
            e.append(np.linalg.norm(full - xs[n], axis=-1))
        err[name] = np.array(e)                                     # (frames, ny, nx) mm
    em = {k: float(np.median(v[:, dom]) / mpp) for k, v in err.items()}
    e95 = {k: float(np.percentile(v[:, dom], 95) / mpp) for k, v in err.items()}
    print(f"  forward FE ({len(fe_frames)} frames, {time.time() - t0:.0f} s): interior |x_FE - x_measured| median "
          + ", ".join(f"{k} {em[k]:.2f} px (95th pct {e95[k]:.2f})" for k in em))

    # ---------------------------------------------------------------- save
    rows = []
    for k_, (i, j) in enumerate(np.argwhere(inner)):
        rows.append(dict(i_theta=i, j_z=j, z0_mm=pz[j], z1_mm=pz[j + 1], theta0_mm=pt[i], theta1_mm=pt[i + 1],
                         alpha=al[k_], alpha_lo95=lo[k_], alpha_hi95=hi[k_], relative=al[k_] / a_mean,
                         relative_lo95=rlo[k_], relative_hi95=rhi[k_], differs_from_mean=bool(differs[k_]),
                         alpha_shape_scale=al_ss[k_], alpha_window75=fits[2 if a.win != 2 else 1]["alpha"][k_],
                         alpha_stage1=a1[i, j], alpha_theta_stage1=at1[i, j], alpha_z_stage1=az1[i, j]))
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "stage2_fit.csv"), index=False)
    np.savez_compressed(os.path.join(a.out, "stage2_fit.npz"), alpha=al, lo95=lo, hi95=hi, boots=boots, inner=inner,
                        alpha_shape_scale=al_ss, lambdas=lams, cv=W0["cv"], lam=lam, fe_frames=fe_frames,
                        fe_error_alpha1=err["alpha = 1"].astype(np.float32),
                        fe_error_fit=err["Stage 2 fit"].astype(np.float32))
    with open(os.path.join(a.out, "stage2_fit.json"), "w") as f:
        json.dump(dict(description="Stage 2 fit of interior patch stiffness factors (ladder (a)) from the virtual-field "
                       "equations on the measured data",
                       frames=int(len(fidx)), window_px=(2 * a.win + 1) * 15, modes=a.modes, n_patches=npat,
                       lambda_cv=lam, alpha_range=[float(al.min()), float(al.max())], alpha_area_mean=a_mean,
                       shape_scale=dict(lambda_cv=lam_s, scale=float(s_all), scale_z_cuts=float(s_z),
                                        scale_theta_cuts=float(s_t)),
                       bootstrap=dict(resamples=a.boot, median_width=float(np.median(hi - lo)),
                                      patches_differing_from_mean=int(differs.sum())),
                       synthetic_limits=dict(patch_contrast_floor=FLOOR, level=LEVEL),
                       residuals=res, correlations=corr, bulk_score=dict(alpha_mean=a_mean, rms_kPa=bulk["rms_kPa"]),
                       forward_fe=dict(frames=int(len(fe_frames)), median_px=em, p95_px=e95)), f, indent=2)

    # ---------------------------------------------------------------- figure
    def grid(v):
        G_ = np.full((nt, nz), np.nan)
        G_[inner] = v
        return G_
    zc, tc = 0.5 * (pz[1:] + pz[:-1]), 0.5 * (pt[1:] + pt[:-1])
    fig, ax = plt.subplots(2, 4, figsize=(22, 11), constrained_layout=True)
    span = max(np.max(np.abs(np.r_[al, al_ss, a1[inner]] - 1)), 0.1)
    norm = TwoSlopeNorm(1.0, 1 - span, 1 + span)
    panels = [(ax[0, 0], grid(al), f"Stage 2 α (combined, λ = {lam:.2g})", norm, "RdBu_r", "{:.2f}"),
              (ax[0, 1], grid(al_ss), f"Stage 2 α (virtual-field shape, cut scale {s_all:.3f})", norm, "RdBu_r", "{:.2f}"),
              (ax[0, 2], grid(fits[2 if a.win != 2 else 1]["alpha"]), "Stage 2 α, 75 px window (sensitivity)", norm,
               "RdBu_r", "{:.2f}"),
              (ax[0, 3], a1 * np.where(inner, 1, np.nan), "Stage 1 α (uniform stress), for comparison", norm, "RdBu_r",
               "{:.2f}"),
              (ax[1, 0], grid(hi - lo), "protocol-bootstrap 95% width of α", plt.Normalize(0, np.max(hi - lo)),
               "magma_r", "{:.2f}")]
    rel_ = grid(al / a_mean)
    for a_, V, title, nm, cmap, fmt in panels:
        a_.set_facecolor("0.85")
        mm = a_.pcolormesh(pz, pt, V, cmap=cmap, norm=nm)
        for i in range(nt):
            for j in range(nz):
                if np.isfinite(V[i, j]):
                    r_, g_, b_, _ = mm.cmap(mm.norm(V[i, j]))
                    a_.text(zc[j], tc[i], fmt.format(V[i, j]), ha="center", va="center", fontsize=6,
                            color="w" if 0.299 * r_ + 0.587 * g_ + 0.114 * b_ < 0.45 else "k")
        fig.colorbar(mm, ax=a_, shrink=.8)
        a_.set_xlim(0, pz[-1]); a_.set_ylim(0, pt[-1]); a_.set_aspect("equal"); a_.set_title(title, fontsize=10)
        a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    # significance overlay on the main map
    for k_, (i, j) in enumerate(np.argwhere(inner)):
        if differs[k_]:
            ax[0, 0].add_patch(plt.Rectangle((pz[j], pt[i]), pz[j + 1] - pz[j], pt[i + 1] - pt[i], fill=False,
                                             ec="k", lw=1.6))
    ax[0, 0].set_title(ax[0, 0].get_title() + f"\nbold outline: differs from the mean by > ±{FLOOR} (95%)", fontsize=10)
    a_ = ax[1, 1]
    a_.semilogx(lams[1:], W0["cv"][1:] / W0["cv"].min(), "o-", ms=3)
    a_.axvline(lam or lams[1], color="r", ls="--", lw=1)
    a_.set_xlabel("λ"); a_.set_ylabel("CV residual / minimum"); a_.grid(alpha=.3)
    a_.set_title("leave-one-protocol-out cross-validation (combined)", fontsize=10)
    a_ = ax[1, 2]
    a_.scatter(a1[inner], al, s=20, label=f"vs Stage 1 α (r {corr['stage1_alpha']:.2f})")
    a_.scatter(at1[inner], al, s=20, marker="x", label=f"vs Stage 1 α_θ (r {corr['stage1_alpha_theta']:.2f})")
    a_.plot([0.3, 1.9], [0.3, 1.9], "k", lw=.7)
    a_.set_xlabel("Stage 1"); a_.set_ylabel("Stage 2 α"); a_.legend(fontsize=8); a_.grid(alpha=.3)
    a_.set_title("Stage 2 vs Stage 1 (per patch)", fontsize=10)
    a_ = ax[1, 3]
    a_.axis("off")
    txt = [f"{len(fidx)} frames, {(2 * a.win + 1) * 15} px window, {2 * a.modes ** 2} virtual fields/frame", "",
           "relative residual (cut_z / cut_θ / vf):"]
    txt += [f"  {k:14s} {v['cut_z']:.3f} / {v['cut_theta']:.3f} / {v['vf']:.3f}" for k, v in res.items()]
    txt += ["", f"area-weighted mean α {a_mean:.3f}; bulk RMS with it {bulk['rms_kPa']:.3f} kPa",
            f"  (baseline {base['inputs']['dic_affine']['rms_kPa']:.3f} kPa)",
            f"cut scales: z {s_z:.3f}, θ {s_t:.3f}", "",
            "forward FE, interior |x_FE − x_measured|:"]
    txt += [f"  {k:12s} median {em[k]:.2f} px, 95th pct {e95[k]:.2f} px" for k in em]
    txt += ["", f"75 px vs 45 px window: r {corr['window_75px_vs_45px']:.2f}",
            f"synthetic limits: contrast ±{FLOOR}, level ±{int(LEVEL * 100)}%"]
    a_.text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=9)
    fig.suptitle("Stage 2: interior patch stiffness factors from the virtual-field equations (measured data, Def11–Def19)")
    fig.savefig(os.path.join(a.out, "stage2_fit.png"), dpi=85)

    # forward-check maps
    fig, ax = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
    e1, e2 = np.median(err["alpha = 1"], 0) / mpp, np.median(err["Stage 2 fit"], 0) / mpp
    vmax = np.percentile(np.r_[e1[2:-2, 2:-2].ravel(), e2[2:-2, 2:-2].ravel()], 99)
    ez, et = d["cell_edges_z_mm"], d["cell_edges_theta_mm"]
    for a_, E, title in ((ax[0], e1, "α = 1 everywhere"), (ax[1], e2, "Stage 2 α (ring α = 1)")):
        mm = a_.pcolormesh(ez, et, np.where(dom, E, np.nan), cmap="magma_r", vmin=0, vmax=vmax)
        fig.colorbar(mm, ax=a_, shrink=.8, label="median |x_FE − x_measured| (px)")
        a_.set_aspect("equal"); a_.set_title(f"forward FE, {title}", fontsize=10)
        a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    Dlt = np.where(dom, e2 - e1, np.nan)
    v_ = np.nanmax(np.abs(Dlt))
    mm = ax[2].pcolormesh(ez, et, Dlt, cmap="PuOr_r", vmin=-v_, vmax=v_)
    fig.colorbar(mm, ax=ax[2], shrink=.8, label="change (px); negative = closer to the measurement")
    ax[2].set_aspect("equal"); ax[2].set_title("Stage 2 α minus α = 1", fontsize=10)
    fig.suptitle(f"Forward check: FE with measured boundary displacements, {len(fe_frames)} frames; interior domain only")
    fig.savefig(os.path.join(a.out, "stage2_forward_check.png"), dpi=85)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
