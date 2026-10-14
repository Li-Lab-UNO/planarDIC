"""
Where the Stage 2 fit leaves misfit: residual breakdown of the virtual-field and cut equations at the fitted
alpha (stage2_fit.py), and whether one alpha map can fit every loading ratio.

Breakdown (relative residual = RMS(G alpha - b) / RMS(b) for the cuts, RMS(G alpha) / RMS(|G| |alpha|) for the
virtual fields), at alpha = 1, at the Stage 2 fit, and, for scale, at the true alpha of the synthetic Stage 1
pattern (stage2_synthetic.py, no noise):
  * by protocol (Def11-Def19)
  * by load level (force as a fraction of the protocol's peak)
  * by virtual-field direction (u*_z, u*_theta) and mode order (m + n)

Direction test: the combined fit (same settings and lambda as stage2_fit.py) repeated on three protocol groups,
  theta-dominant  Def12-Def14 (sets 3-5)
  equibiaxial     Def11, Def15, Def19 (sets 2, 6, 10)
  z-dominant      Def16-Def18 (sets 7-9)
then each group's map is scored on the other groups' equations. If one alpha per patch were adequate, the maps
would agree and cross-prediction would cost little.

    .venv/bin/python constitutive_fit/stage2_residuals.py --test N1739/test1
Outputs: <test>/constitutive_fit/local/stage2_fit/stage2_residuals.png / .json
"""
import json
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import baseline as B  # noqa: E402
import constitutive as M  # noqa: E402
from stage2_fit import combined_system  # noqa: E402
from stage2_setup import LOCAL, assemble, geometry, measured_positions, window_F  # noqa: E402
from stage2_synthetic import difference_matrix, neighbor_pairs  # noqa: E402

GROUPS = {"θ-dominant (Def12–14)": (3, 4, 5), "equibiaxial (Def11/15/19)": (2, 6, 10),
          "z-dominant (Def16–18)": (7, 8, 9)}


def rel_cut(G, b, al, m):
    r = (G @ al - b)[m]
    return float(np.sqrt(np.mean(r ** 2)) / np.sqrt(np.mean(b[m] ** 2)))


def rel_vf(G, al, m):
    return float(np.sqrt(np.mean((G @ al)[m] ** 2)) / np.sqrt(np.mean((np.abs(G) @ np.abs(al))[m] ** 2)))


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    out = os.path.join(LOCAL, "stage2_fit")
    S = np.load(os.path.join(out, "stage2_fit.npz"))
    fitinfo = json.load(open(os.path.join(out, "stage2_fit.json")))
    al_fit, lam = S["alpha"], float(S["lam"])
    d = np.load(os.path.join(LOCAL, "local_data.npz"))
    base, df = B.load(str(d["baseline"]))
    df = df.reset_index(drop=True)
    st = base["settings"]
    H, w_t, w_z = st["H_mm"], st["w_theta_mm"], st["w_z_mm"]
    P = np.array([base["inputs"]["dic_affine"]["params"][n] for n in M.PARAM_NAMES])
    mpp = float(d["mm_per_px"])
    modes = fitinfo["modes"]
    win = (fitinfo["window_px"] // 15 - 1) // 2
    g = geometry(d, 1, modes)
    npat, inner = g["npat"], g["inner"]
    every = max(1, round(len(df) / fitinfo["frames"]))
    fidx = np.arange(0, len(df), every)
    sets = df.set.to_numpy()[fidx]
    Fz, Ft = d["force_z_mN"][fidx], d["force_theta_mN"][fidx]
    xm, vm = measured_positions(df, fidx, mpp, with_valid=True)
    Xg = np.stack([d["z_mm"], d["theta_mm"]], -1)
    F = window_F(Xg, xm, vm, win)
    Gz, Gt, Gv, bz, bt = assemble(g, F, Fz, Ft, P, H, w_t, w_z)
    nf = len(fidx)
    # load level per frame: larger of the two forces as a fraction of that protocol's peak
    load = np.zeros(nf)
    for s in np.unique(sets):
        m = sets == s
        load[m] = np.maximum(Fz[m] / max(Fz[m].max(), 1e-9), Ft[m] / max(Ft[m].max(), 1e-9))
    vf_dir = np.array([0 if "u*_z" in n_ else 1 for n_ in g["vf_names"]])
    vf_ord = np.array([int(n_.split("m=")[1].split()[0]) + int(n_.split("n=")[1]) for n_ in g["vf_names"]])

    # synthetic reference at the true alpha (Stage 1 pattern), same window and modes
    syn = None
    fe = os.path.join(LOCAL, "stage2_synthetic", "fe_stage1.npz")
    if os.path.exists(fe):
        from scipy.ndimage import gaussian_filter
        C = np.load(fe)
        sfr = C["frames"]
        sxm = measured_positions(df, sfr, mpp)
        sxs = Xg + np.stack([gaussian_filter(sxm[..., i] - Xg[..., i], (0, 1.0, 1.0), mode="nearest")
                             for i in range(2)], -1)
        m_ = int(C["inset"])
        sxs[:, m_:-m_, m_:-m_] = C["x"]
        sF = window_F(Xg, sxs, np.ones(sxs.shape[:3], bool), win)
        sGz, sGt, sGv, sbz, sbt = assemble(g, sF, C["Fz"], C["Ft"], P, H, w_t, w_z)
        syn = dict(G=(sGz, sbz, sGt, sbt, sGv), al=C["alpha_true"][inner], sets=df.set.to_numpy()[sfr])

    def breakdown(Gz_, bz_, Gt_, bt_, Gv_, al, sets_, load_=None):
        r = {}
        for s in np.unique(sets_):
            m = sets_ == s
            r[f"set {s}"] = dict(cut_z=rel_cut(Gz_.reshape(-1, npat), bz_.ravel(), al, np.repeat(m, Gz_.shape[1])),
                                 cut_theta=rel_cut(Gt_.reshape(-1, npat), bt_.ravel(), al, np.repeat(m, Gt_.shape[1])),
                                 vf=rel_vf(Gv_.reshape(-1, npat), al, np.repeat(m, Gv_.shape[1])))
        if load_ is not None:
            for lo_, hi_ in ((0, .25), (.25, .5), (.5, .75), (.75, 1.01)):
                m = (load_ >= lo_) & (load_ < hi_)
                r[f"load {int(lo_ * 100)}-{min(int(hi_ * 100), 100)}%"] = dict(
                    cut_z=rel_cut(Gz_.reshape(-1, npat), bz_.ravel(), al, np.repeat(m, Gz_.shape[1])),
                    cut_theta=rel_cut(Gt_.reshape(-1, npat), bt_.ravel(), al, np.repeat(m, Gt_.shape[1])),
                    vf=rel_vf(Gv_.reshape(-1, npat), al, np.repeat(m, Gv_.shape[1])))
        for k_, nm in ((0, "u*_z"), (1, "u*_theta")):
            m = np.tile(vf_dir == k_, Gv_.shape[0])
            r[f"vf {nm}"] = dict(vf=rel_vf(Gv_.reshape(-1, npat), al, m))
        for o in sorted(set(vf_ord)):
            m = np.tile(vf_ord == o, Gv_.shape[0])
            r[f"vf order {o}"] = dict(vf=rel_vf(Gv_.reshape(-1, npat), al, m))
        return r
    R = {"alpha = 1": breakdown(Gz, bz, Gt, bt, Gv, np.ones(npat), sets, load),
         "Stage 2 fit": breakdown(Gz, bz, Gt, bt, Gv, al_fit, sets, load)}
    if syn:
        R["synthetic, true alpha"] = breakdown(*syn["G"], syn["al"], syn["sets"])

    # ---------------------------------------------------------------- direction test
    DtD = difference_matrix(neighbor_pairs(inner), npat)
    AtA, Atb, btb = combined_system(Gz, bz, Gt, bt, Gv, sets, npat)
    gfit = {}
    for name, ss in GROUPS.items():
        AA = sum(AtA[s] for s in ss)
        Ab = sum(Atb[s] for s in ss)
        gfit[name] = np.linalg.solve(AA + lam * DtD, Ab)

    def score(al, ss):                      # combined scaled residual on the protocols ss
        return float(sum(al @ AtA[s] @ al - 2 * al @ Atb[s] + btb[s] for s in ss))
    cross = {}
    for fn, al in gfit.items():
        cross[fn] = {tn: score(al, ss) / score(gfit[tn], ss) for tn, ss in GROUPS.items()}
    cross["all-protocol fit"] = {tn: score(al_fit, ss) / score(gfit[tn], ss) for tn, ss in GROUPS.items()}
    names = list(GROUPS)
    cmat = np.array([[np.corrcoef(gfit[a_], gfit[b_])[0, 1] for b_ in names] for a_ in names])

    # ---------------------------------------------------------------- print
    for k, r in R.items():
        print(f"[{k}]")
        for g_, v in r.items():
            print(f"   {g_:14s} " + "  ".join(f"{b_} {x:.3f}" for b_, x in v.items()))
    print("group-fit maps: correlation\n" + "\n".join(f"   {n_:28s} " + " ".join(f"{x:5.2f}" for x in row)
                                                      for n_, row in zip(names, cmat)))
    print("cross-prediction: residual on the target group / residual of the target group's own fit")
    for fn, v in cross.items():
        print(f"   map from {fn:28s}: " + ", ".join(f"{tn.split(' (')[0]} {x:.2f}" for tn, x in v.items()))
    for n_, al in gfit.items():
        print(f"   {n_}: alpha {al.min():.2f}-{al.max():.2f}, mean {al.mean():.2f}")
    with open(os.path.join(out, "stage2_residuals.json"), "w") as f:
        json.dump(dict(breakdown=R, group_map_correlation=dict(groups=names, matrix=cmat.tolist()),
                       cross_prediction=cross,
                       group_alpha={k: v.tolist() for k, v in gfit.items()}), f, indent=2)

    # ---------------------------------------------------------------- figure
    fig = plt.figure(figsize=(18, 11), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    cols = {"alpha = 1": "0.6", "Stage 2 fit": "C0", "synthetic, true alpha": "C2"}
    for c_, (title, keys) in enumerate((("by protocol", [k for k in R["Stage 2 fit"] if k.startswith("set")]),
                                        ("by load level (fraction of protocol peak)",
                                         [k for k in R["Stage 2 fit"] if k.startswith("load")]))):
        for r_, blk in enumerate(("vf", "cut_theta")):
            a_ = fig.add_subplot(gs[r_, c_])
            x = np.arange(len(keys))
            for o_, (nm, r) in enumerate(R.items()):
                vals = [r[k][blk] if k in r else np.nan for k in keys]
                a_.bar(x + (o_ - 1) * 0.27, vals, 0.27, color=cols[nm], label=nm)
            a_.set_xticks(x); a_.set_xticklabels([k.replace("set ", "set\n") for k in keys], fontsize=8)
            a_.set_ylabel(f"relative residual ({blk})"); a_.grid(alpha=.3, axis="y"); a_.legend(fontsize=7)
            a_.set_title(f"{blk} residual {title}", fontsize=10)
    a_ = fig.add_subplot(gs[0, 2])
    keys = [k for k in R["Stage 2 fit"] if k.startswith("vf order")]
    x = np.arange(len(keys))
    for o_, (nm, r) in enumerate(R.items()):
        a_.bar(x + (o_ - 1) * 0.27, [r[k]["vf"] for k in keys], 0.27, color=cols[nm], label=nm)
    a_.set_xticks(x); a_.set_xticklabels([k.split()[-1] for k in keys], fontsize=7)
    a_.set_xlabel("virtual-field mode order m + n"); a_.set_ylabel("relative residual (vf)")
    a_.grid(alpha=.3, axis="y"); a_.legend(fontsize=7); a_.set_title("vf residual by mode order", fontsize=10)
    a_ = fig.add_subplot(gs[1, 2]); a_.axis("off")
    txt = ["Group fits (same λ), map correlation:"] + \
          [f"  {n_.split(' (')[0]:12s} " + " ".join(f"{x:5.2f}" for x in row) for n_, row in zip(names, cmat)] + \
          ["", "Cross-prediction (residual on target /", " target group's own fit):"]
    for fn, v in cross.items():
        txt.append(f"  from {fn.split(' (')[0]:16s}: " + ", ".join(f"{x:.2f}" for x in v.values()))
    txt += ["", "  targets: " + ", ".join(n_.split(" (")[0] for n_ in names)]
    a_.text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=9)
    allv = np.concatenate(list(gfit.values()))
    span = max(np.max(np.abs(allv - 1)), 0.1)
    norm = TwoSlopeNorm(1.0, max(1 - span, 1e-3), 1 + span)
    pz, pt = g["pz"], g["pt"]
    fig.suptitle("Stage 2 residual breakdown at the fitted α (blue), α = 1 (gray) and the synthetic test's true α (green); "
                 "group maps in stage2_group_maps.png")
    fig.savefig(os.path.join(out, "stage2_residuals.png"), dpi=80)

    fig, ax = plt.subplots(1, 4, figsize=(22, 5.6), constrained_layout=True)
    maps = list(gfit.items()) + [("all protocols (stage2_fit)", al_fit)]
    for a_, (n_, al) in zip(ax, maps):
        G_ = np.full(inner.shape, np.nan)
        G_[inner] = al
        a_.set_facecolor("0.85")
        mm = a_.pcolormesh(pz, pt, G_, cmap="RdBu_r", norm=norm)
        fig.colorbar(mm, ax=a_, shrink=.8, extend="both")
        a_.set_xlim(0, pz[-1]); a_.set_ylim(0, pt[-1]); a_.set_aspect("equal")
        a_.set_title(f"α fitted on {n_}", fontsize=10)
        a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    fig.suptitle("Stage 2 fitted separately on each protocol group (same settings and λ)")
    fig.savefig(os.path.join(out, "stage2_group_maps.png"), dpi=85)
    print("wrote", out)


if __name__ == "__main__":
    main()
