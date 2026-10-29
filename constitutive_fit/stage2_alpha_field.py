"""
Stage 2 of the local-property stages: a smooth relative membrane-stiffness field alpha(x) from the interior strains.

Model: the fixed eight-parameter model with the test's bulk parameters (DIC affine F) everywhere, times a stiffness
factor alpha(x) that scales Cgr and every C1 together (exponents and gamma unchanged). alpha is the local membrane
stiffness relative to the bulk: H(x)/H times a modulus factor; thickness and modulus cannot be told apart.
    log alpha(z, theta) = sum_k c_k P_m(s) P_n(t)
with P Legendre polynomials and s, t the mesh coordinates mapped to [-1, 1]. A basis is given as a total degree D
(1 <= m + n <= D, e.g. "4") or as separate maximum degrees in z and theta, "DZxDT" (m <= DZ, n <= DT, (m, n) != (0, 0),
e.g. "2x6"), which resolves horizontal bands with fewer terms than a total degree that reaches the same order in theta. There is no constant term: with
displacements prescribed on the whole boundary the strains do not depend on the stiffness level, so the level is
fixed at the bulk value (mean log alpha ~ 0) and only the pattern is fitted.

FE: the tracked domain (membrane_fe.py) with the measured boundary positions, as in stages 1 and 1b. Objective: at the
peak of the training protocols, every interior point (edge_mm) and strain component (E_zz, E_thetatheta, E_ztheta),
measured minus FE, each block divided by the measured interior SD of that protocol and component. scipy least_squares;
the finite-difference Jacobian is solved in parallel (--workers), each column starting from the current solution.
Training: Def11, Def14, Def18. Held out: the other six peaks, scored with the explained share of stage 1.

--synthetic NAME replaces the measured data by FE data for a known alpha (same boundary and model), with Gaussian
position noise of --noise-px on every point (boundary included) before the strain pipeline. The default 0.06 px is the
high-frequency displacement of the protocols' rest frames (0.04-0.08 px), an upper bound on the DIC noise.
    uniform   alpha = 1: spurious contrast floor
    smooth    in the degree-2 basis, about +-25%
    bands     horizontal bands (stiff at theta ~ 5 mm, soft at ~2.5 and ~8 mm) and an oblique streak, like the
              measured residual pattern; not in any polynomial basis, so the fit can at best recover its projection
    bands_theta   positive control for the direction test: the same pattern (x THETA_GAIN) scales only the SMC term
              (the theta fibers, the larger part of the theta stiffness), so the stiffness differs by direction; the
              fit still uses one alpha on the whole response
    bands_shape   model-error control: the bands scale the whole response, but the data are generated with another
              bulk fit of the same test (SHAPE_INPUT, as good a force fit with a different response shape) and
              fitted with the DIC affine one, so the fit's material law is wrong as it will be on measured data
Recovery is scored on the interior Gauss points after removing the mean of log alpha (the level is not identified):
correlation with the truth, RMS error, and the RMS error of the truth's least-squares projection onto the same basis.

    .venv/bin/python constitutive_fit/stage2_alpha_field.py --test N1739/test1 --synthetic bands --bases 4 1x6 2x6
    .venv/bin/python constitutive_fit/stage2_alpha_field.py --test N1739/test1 --measured --bases 1 2 3 4 2x6

--groups (one basis): direction test. alpha is fitted separately to each loading-ratio group,
    theta_dominant  Def12-Def14 (axial reduced)
    equibiaxial     Def11, Def15, Def19
    z_dominant      Def16-Def18 (circumferential reduced)
and each group's map is scored on every group's protocols. A stiffness factor that scales the response in every
direction gives the same map from every group, and the maps predict each other. Stiffness that differs by direction
keeps the pattern (high correlation between the group maps) but changes its amplitude: a group that loads the stiffer
direction less sees weaker contrast. So the amplitude ratio (least-squares slope of one map's log alpha on the other's,
interior means removed) and the cross-prediction loss are the measures, not the correlation. On synthetic data (alpha uniform in direction by construction)
the disagreement between group maps is the floor that measured disagreement has to exceed.
    .venv/bin/python constitutive_fit/stage2_alpha_field.py --test N1739/test1 --synthetic bands --groups --bases 2x6

--groups sequence: the same test on the three equibiaxial protocols, each alone, in run order: Def11 (first of the
stretch-ratio protocols), Def15 (middle), Def19 (last). The loading is the same, so maps that differ point to the
tissue changing during the test (damage, conditioning, incomplete recovery), which would also make the loading-ratio
groups differ, since the theta-dominant protocols (Def12-14) ran before the z-dominant ones (Def16-18).
    .venv/bin/python constitutive_fit/stage2_alpha_field.py --test N1739/test1 --measured --groups sequence --bases 2x6

--data FILE fits positions generated by another model (stage3b_thickness_synthetic.py: a 3D slab with varying
thickness) with the same membrane alpha, so the direction test asks whether thickness variation in 3D looks to a
membrane fit like stiffness that differs with loading ratio.
Outputs: <test>/constitutive_fit/local/stage2_<synthetic NAME | measured>.json / .npz / .png, with --groups
stage2_<synthetic NAME | measured>_groups_<basis>.json / .npz / .png
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
from numpy.polynomial import legendre  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import constitutive as M  # noqa: E402
import dic2d as D  # noqa: E402
from stage1_uniform_fe import COMPONENTS, INPUT, explained, strain  # noqa: E402
from stage1b_best_uniform import Problem  # noqa: E402

TRAIN = ("Def11", "Def14", "Def18")
THETA_GAIN = 1.5          # bands_theta: amplitude on the SMC term
SHAPE_INPUT = "rake"      # bands_shape: the bulk fit that generates the data
GROUPS = {"theta_dominant": ("Def12", "Def13", "Def14"), "equibiaxial": ("Def11", "Def15", "Def19"),
          "z_dominant": ("Def16", "Def17", "Def18")}
GROUP_SETS = {"loading": GROUPS, "sequence": {"Def11_first": ("Def11",), "Def15_middle": ("Def15",),
                                              "Def19_last": ("Def19",)}}
SHORT = {"theta_dominant": "θ-dom.", "equibiaxial": "equibiax.", "z_dominant": "z-dom.",
         "Def11_first": "Def11", "Def15_middle": "Def15", "Def19_last": "Def19"}
DIFF_STEP = 0.01          # finite-difference step on the coefficients of log alpha


# ------------------------------------------------------------------ basis and truths
def mesh_coords(pb):
    """Gauss point coordinates (ne_r, ne_c, 4, 2) in specimen mm and the [-1, 1] maps of the mesh extent."""
    g = pb.mesh.gauss_points
    z0, z1 = pb.mesh.Xg[..., 0].min(), pb.mesh.Xg[..., 0].max()
    t0, t1 = pb.mesh.Xg[..., 1].min(), pb.mesh.Xg[..., 1].max()
    s = 2 * (g[..., 0] - z0) / (z1 - z0) - 1
    t = 2 * (g[..., 1] - t0) / (t1 - t0) - 1
    return g, s, t


def modes(spec):
    """(m, n) degree pairs of a basis spec: "D" (total degree) or "DZxDT" (maximum degree in z and in theta)."""
    if "x" in spec:
        dz, dt = map(int, spec.split("x"))
        return [(m, n) for m in range(dz + 1) for n in range(dt + 1) if m + n > 0]
    d = int(spec)
    return [(m, k - m) for k in range(1, d + 1) for m in range(k + 1)]


def basis(s, t, spec):
    """(..., K) Legendre products P_m(s) P_n(t) for the basis spec, and their (m, n)."""
    mn = modes(spec)
    B = np.stack([legendre.legval(s, np.eye(m + 1)[m]) * legendre.legval(t, np.eye(n + 1)[n]) for m, n in mn], -1)
    return B, mn


def truth(name, g, s, t):
    z, th = g[..., 0], g[..., 1]
    if name == "uniform":
        return np.zeros_like(z)
    if name == "smooth":
        P1s, P1t, P2t = s, t, 0.5 * (3 * t ** 2 - 1)
        return 0.15 * P1s + 0.12 * P2t - 0.10 * P1s * P1t
    if name in ("bands", "bands_theta", "bands_shape"):
        bands = (0.20 * np.exp(-((th - 5.0) / 1.0) ** 2) - 0.15 * np.exp(-((th - 8.2) / 0.9) ** 2)
                 - 0.12 * np.exp(-((th - 2.5) / 0.8) ** 2))
        streak = 0.12 * np.exp(-(((z - th) / np.sqrt(2) + 0.5) / 0.6) ** 2)      # lower left to upper right
        return (THETA_GAIN if name == "bands_theta" else 1.0) * (bands + streak)
    raise ValueError(name)


def interior_gp(pb, g):
    """Gauss points at least edge_mm from the rest-state tine-tip lines (as dic2d.interior_mask)."""
    ez, et = pb.edge["z"], pb.edge["theta"]
    Zt = (pb.tips["z"][1] - pb.tips["z"][0]) * pb.mpp
    Tt = (pb.tips["theta"][1] - pb.tips["theta"][0]) * pb.mpp
    return (g[..., 0] >= ez) & (g[..., 0] <= Zt - ez) & (g[..., 1] >= et) & (g[..., 1] <= Tt - et)


# ------------------------------------------------------------------ parallel FE workers
_W = {}


def _init(test, xb_override, params):
    np.seterr(over="ignore", invalid="ignore")
    testcfg.select(test)
    T = testcfg.get()
    df = pd.read_csv(os.path.join(T.fit_dir, "stretch_inputs_data.csv"))
    pb = Problem(T, df)
    for f, xb in (xb_override or {}).items():
        pb.prot[f]["xb"] = xb
    _W.update(pb=pb, params=params, basis={})


def _solve(task):
    """task = (folder, degree, c, start) -> (folder, x_FE, interior strains (3, n_mask) in COMPONENTS order)."""
    folder, degree, c, start = task
    pb = _W["pb"]
    if degree not in _W["basis"]:
        _, s, t = mesh_coords(pb)
        _W["basis"][degree] = basis(s, t, degree)[0]
    alpha = np.exp(_W["basis"][degree] @ c) if len(c) else None
    xf, Ef = pb.solve(folder, _W["params"], warm=False, alpha=alpha, start=start)
    m = pb.prot[folder]["mask"]
    return folder, xf, np.stack([Ef[..., i, j][m] for _, i, j in COMPONENTS])


# ------------------------------------------------------------------ fit
class Fit:
    def __init__(self, pool, pb, folders, degree, data):
        self.pool, self.pb, self.folders, self.degree, self.data = pool, pb, folders, degree, data
        self.cache = None                                       # (c, {folder: x_FE}, residual)

    def resid(self, f, E):
        return np.concatenate([(self.data[f]["E"][k] - E[k]) / self.data[f]["sd"][k] for k in range(3)])

    def base(self, c):
        if self.cache is not None and np.array_equal(self.cache[0], c):
            return self.cache
        start = self.cache[1] if self.cache is not None else {f: None for f in self.folders}
        res = list(self.pool.map(_solve, [(f, self.degree, c, start[f]) for f in self.folders]))
        X = {f: x for f, x, _ in res}
        r = np.concatenate([self.resid(f, E) for f, _, E in res])
        self.cache = (c.copy(), X, r)
        return self.cache

    def fun(self, c):
        return self.base(c)[2]

    def jac(self, c):
        _, X, r0 = self.base(c)
        tasks = []
        for k in range(len(c)):
            ck = c.copy()
            ck[k] += DIFF_STEP
            tasks += [(f, self.degree, ck, X[f]) for f in self.folders]
        out = list(self.pool.map(_solve, tasks))
        J = np.empty((len(r0), len(c)))
        for k in range(len(c)):
            rk = np.concatenate([self.resid(f, E) for f, _, E in out[k * len(self.folders):(k + 1) * len(self.folders)]])
            J[:, k] = (rk - r0) / DIFF_STEP
        return J


def plot_groups(pb, gin, true_log, maps, names, cross, base, agree, spec, stem, out, T, synthetic):
    """Group maps (truth first if synthetic) on one log-symmetric color scale, then cross-prediction per component."""
    # figure: the group maps (truth first if synthetic) on one color scale, then cross-prediction per component
    from matplotlib.colors import TwoSlopeNorm
    em = lambda Lg: np.exp(Lg.mean(2))  # noqa: E731
    panels = ([("truth", true_log)] if true_log is not None else []) + [(f"fitted to {n}", maps[n]) for n in names]
    # symmetric in log alpha, set by the interior (outside it the field is unconstrained and saturates)
    span = max(np.exp(max(np.percentile(np.abs(m[gin] - m[gin].mean()), 99) for _, m in panels)), 1.05)
    norm = TwoSlopeNorm(1.0, 1 / span, span)
    Ze, Te = pb.mesh.Xg[..., 0], pb.mesh.Xg[..., 1]
    ez, et = pb.edge["z"], pb.edge["theta"]
    Zt = (pb.tips["z"][1] - pb.tips["z"][0]) * pb.mpp
    Tt = (pb.tips["theta"][1] - pb.tips["theta"][0]) * pb.mpp
    nc = max(len(panels), 3)
    fig, axes = plt.subplots(2, nc, figsize=(3.8 * nc, 8), constrained_layout=True, squeeze=False)
    for ax, (title, m) in zip(axes[0], panels):
        mm_ = m - m[gin].mean()
        pm = ax.pcolormesh(Ze, Te, em(mm_), cmap=D.diverging_cmap(), norm=norm, shading="flat")
        ax.add_patch(plt.Rectangle((ez, et), Zt - 2 * ez, Tt - 2 * et, fill=False, ls="--", lw=0.8))
        fig.colorbar(pm, ax=ax, shrink=0.8, label="α (interior mean of log α removed)")
        ax.set_aspect("equal"); ax.set_title(title, fontsize=9)
        ax.set_xlabel("z (mm)", fontsize=8); ax.set_ylabel("θ (mm)", fontsize=8); ax.tick_params(labelsize=7)
    for ax in axes[0][len(panels):]:
        ax.axis("off")
    for ax, (cn, _, _) in zip(axes[1], COMPONENTS):
        Mx = np.array([[cross[(gf, gs, cn)] for gs in names] for gf in names])
        im = ax.imshow(Mx, cmap="viridis", vmin=min(0, Mx.min()), vmax=1)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{Mx[i, j]:+.2f}\n(α=1 {base[(names[j], cn)]:+.2f})", ha="center", va="center",
                        fontsize=8, color="w" if Mx[i, j] < 0.6 else "k")
        ax.set_xticks(range(3), [SHORT[n] for n in names]); ax.set_yticks(range(3), [SHORT[n] for n in names])
        ax.set(xlabel="scored on", ylabel="fitted to", title=f"{cn}: interior variance explained")
        fig.colorbar(im, ax=ax, shrink=0.8)
    for ax in axes[1][3:]:
        ax.axis("off")
    what = (f"synthetic '{synthetic}'" if synthetic else "measured") if not str(stem).startswith("stage2_stage3b")\
        else f"3D synthetic '{stem.split('_groups')[0].replace('stage2_', '')}'"
    txt = "\n".join(f"{k.replace('_dominant', '-dominant')}: r {v['r']:.2f}, amplitude ratio {v['amplitude_ratio']:.2f}"
                    for k, v in agree.items())
    fig.suptitle(f"{T.name}: stage 2 {'direction' if '_groups_' in stem else 'sequence'} test, basis {spec} ({what})\nmap agreement (interior), "
                 f"second map on first:\n{txt}", fontsize=10)
    fig.savefig(os.path.join(out, stem + ".png"), dpi=110)
    plt.close(fig)


def run_groups(pool, pb, data, ref, spec, folders, s, t, gin, true_log, out, tag, T, a, t_start, set_name="loading"):
    """Direction test: alpha fitted to each loading-ratio group, scored on every group."""
    B, mn = basis(s, t, spec)
    gfold = {gname: [f for f in folders if f[:5] in defs] for gname, defs in GROUP_SETS[set_name].items()}
    maps, coefs, pred, fits = {}, {}, {}, {}
    for gname, fs in gfold.items():
        t0 = time.time()
        fit = Fit(pool, pb, fs, spec, data)
        sol = least_squares(fit.fun, np.zeros(B.shape[-1]), jac=fit.jac, max_nfev=a.max_nfev, xtol=1e-3, ftol=1e-4)
        coefs[gname] = sol.x
        maps[gname] = B @ sol.x
        pred[gname] = {f: E for f, _, E in pool.map(_solve, [(f, spec, sol.x, None) for f in folders])}
        fits[gname] = dict(nfev=int(sol.nfev), message=str(sol.message), fit_s=time.time() - t0,
                           train_cost=float(np.mean(sol.fun ** 2)))
        print(f"  fitted to {gname} ({len(sol.x)} coef, {sol.nfev} it, {time.time() - t0:.0f} s)", flush=True)

    def ex(f, E, k):
        return explained(data[f]["E"][k], E[k], np.ones(data[f]["E"].shape[1], bool))[0]

    # cross-prediction: median explained share of each group's protocols, per component, for each fitted map
    cross = {}
    for gfit in gfold:
        for gscore, fs in gfold.items():
            for k, (cn, _, _) in enumerate(COMPONENTS):
                cross[(gfit, gscore, cn)] = float(np.median([ex(f, pred[gfit][f], k) for f in fs]))
    base = {(gs, cn): float(np.median([ex(f, ref[f], k) for f in fs]))
            for gs, fs in gfold.items() for k, (cn, _, _) in enumerate(COMPONENTS)}
    # map agreement on the interior (mean of log alpha removed: the level is not identified)
    L = {gname: m[gin] - m[gin].mean() for gname, m in maps.items()}
    names = list(gfold)
    agree = {}
    for i, g1 in enumerate(names):
        for g2 in names[i + 1:]:
            agree[f"{g1} vs {g2}"] = dict(r=float(np.corrcoef(L[g1], L[g2])[0, 1]),
                                          rms_log=float(np.sqrt(np.mean((L[g1] - L[g2]) ** 2))),
                                          amplitude_ratio=float(L[g2] @ L[g1] / (L[g1] @ L[g1])))
    vs_truth = {}
    if true_log is not None:
        tr = true_log[gin] - true_log[gin].mean()
        vs_truth = {gname: dict(r=float(np.corrcoef(L[gname], tr)[0, 1]),
                                rms_log=float(np.sqrt(np.mean((L[gname] - tr) ** 2)))) for gname in names}
    sd = {gname: float(L[gname].std()) for gname in names}
    for k_, v in agree.items():
        print(f"  maps {k_}: r {v['r']:.2f}, amplitude ratio {v['amplitude_ratio']:.2f}, "
              f"RMS difference of log alpha {v['rms_log']:.3f}", flush=True)
    for gname in names:
        print(f"  {gname} map: SD of log alpha {sd[gname]:.3f}"
              + (f", vs truth r {vs_truth[gname]['r']:.2f} RMS {vs_truth[gname]['rms_log']:.3f}" if vs_truth else ""),
              flush=True)
    for cn, _, _ in COMPONENTS:
        print(f"  {cn} explained (rows: fitted to; columns: scored on {', '.join(names)}; alpha = 1 in brackets)")
        for gfit in names:
            print(f"    {gfit:15s} " + "  ".join(f"{cross[(gfit, gs, cn)]:+.2f} ({base[(gs, cn)]:+.2f})" for gs in names))

    stem = f"stage2_{tag}_{'groups' if set_name == 'loading' else set_name}_{spec}"
    with open(os.path.join(out, stem + ".json"), "w") as fh:
        json.dump(dict(test=T.name, source=tag, basis=spec, n_coef=int(B.shape[-1]), groups=gfold, fits=fits,
                       coef={k: v.tolist() for k, v in coefs.items()}, map_sd_log=sd, map_agreement=agree,
                       map_vs_truth=vs_truth,
                       cross_explained={f"{gf} -> {gs}": {cn: cross[(gf, gs, cn)] for cn, _, _ in COMPONENTS}
                                        for gf in names for gs in names},
                       alpha1_explained={gs: {cn: base[(gs, cn)] for cn, _, _ in COMPONENTS} for gs in names},
                       runtime_s=time.time() - t_start), fh, indent=2)
    np.savez_compressed(os.path.join(out, stem + ".npz"), interior=gin,
                        truth_log_alpha=true_log if true_log is not None else np.zeros(0),
                        **{f"log_alpha_{k}": v for k, v in maps.items()})

    plot_groups(pb, gin, true_log, maps, names, cross, base, agree, spec, stem, out, T, a.synthetic)
    print(f"wrote {os.path.join(out, stem)}.* ({time.time() - t_start:.0f} s)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--synthetic", choices=("uniform", "smooth", "bands", "bands_theta", "bands_shape"))
    src.add_argument("--measured", action="store_true")
    src.add_argument("--data", help="synthetic positions from another model (.npz with pos_<folder> in image px, "
                                    "and optionally truth_log_alpha at the membrane Gauss points), e.g. from "
                                    "stage3b_thickness_synthetic.py; used as given (noise already added)")
    ap.add_argument("--bases", nargs="+", default=["1", "2", "3", "4"],
                    help='total degree ("4") or maximum degree in z and theta ("2x6")')
    ap.add_argument("--noise-px", type=float, default=0.06)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-nfev", type=int, default=15)
    ap.add_argument("--groups", nargs="?", const="loading", choices=tuple(GROUP_SETS),
                    help="direction test (loading, the default) or sequence test, with the first basis of --bases")
    a = ap.parse_args()
    T = testcfg.get()
    test = os.path.relpath(T.root, testcfg.ROOT)
    np.seterr(over="ignore", invalid="ignore")
    out = T.local_dir
    os.makedirs(out, exist_ok=True)
    tag = (f"synthetic_{a.synthetic}" if a.synthetic else
           os.path.basename(a.data).replace("_data.npz", "").replace(".npz", "") if a.data else "measured")
    df = pd.read_csv(os.path.join(T.fit_dir, "stretch_inputs_data.csv"))
    P = pd.read_csv(os.path.join(T.fit_dir, "params_by_stretch_input.csv")).set_index("input")
    params = P.loc[INPUT, list(M.PARAM_NAMES)].to_numpy(float)
    t_start = time.time()
    pb = Problem(T, df)
    folders = list(pb.prot)
    train = [f for f in folders if f[:5] in TRAIN]
    held = [f for f in folders if f not in train]
    g, s, t = mesh_coords(pb)
    gin = interior_gp(pb, g)

    # ---- data: measured interior strains, or FE data for a known alpha plus position noise
    xb_override, true_log = None, None
    data = {}

    def use_positions(f, x_px):
        """Synthetic data: positions x_px (image px) replace the measured ones, boundary included."""
        xm = pb.to_mesh(x_px)
        xb = pb.prot[f]["xb"].copy()
        bnd = ~pb.mesh.free
        xb[bnd] = xm[bnd]
        xb_override[f] = xb
        pb.prot[f]["xb"] = xb
        Em, _ = strain(pb.X0, x_px, np.ones(len(pb.X0), bool), pb.shape, pb.window_mm, pb.mpp)
        pb.prot[f]["Em"], pb.prot[f]["xp"] = Em, x_px

    if a.data:
        src = np.load(a.data)
        true_log = src["truth_log_alpha"] if "truth_log_alpha" in src and src["truth_log_alpha"].size else None
        xb_override = {}
        for f in folders:
            use_positions(f, src[f"pos_{f}"])
    if a.synthetic:
        true_log = truth(a.synthetic, g, s, t)
        rng = np.random.default_rng(1)
        xb_override = {}
        if a.synthetic == "bands_theta":                       # alpha on the SMC term only
            p_smc = np.zeros(8); p_smc[3], p_smc[4], p_smc[7] = params[3], params[4], params[7]
            p_smc[2], p_smc[6] = params[2], params[6]          # exponents kept positive; their C1 are 0
            p_rest = params.copy(); p_rest[3] = 0.0
            gen_params = [p_rest, p_smc]
            gen_alpha = np.stack([np.ones_like(true_log), np.exp(true_log)], -1)
        elif a.synthetic == "bands_shape":                     # another bulk fit's response shape
            gen_params = P.loc[SHAPE_INPUT, list(M.PARAM_NAMES)].to_numpy(float)
            gen_alpha = np.exp(true_log)
        else:
            gen_params, gen_alpha = params, np.exp(true_log)
        for f in folders:
            xf, _ = pb.solve(f, gen_params, warm=False, alpha=gen_alpha)
            use_positions(f, pb.to_px(xf) + rng.normal(0, a.noise_px, (len(pb.X0), 2)))
            print(f"  synthetic {a.synthetic}: {f} ({time.time() - t_start:.0f} s)", flush=True)
    for f in folders:
        p = pb.prot[f]
        E = np.stack([p["Em"][..., i, j][p["mask"]] for _, i, j in COMPONENTS])
        data[f] = dict(E=np.nan_to_num(E), sd=E.std(1))

    results = []
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(test, xb_override, params)) as pool:
        # alpha = 1 reference (stage 1) on every protocol
        ref = {f: E for f, _, E in pool.map(_solve, [(f, "1", np.zeros(0), None) for f in folders])}
        if a.groups:
            run_groups(pool, pb, data, ref, a.bases[0], folders, s, t, gin, true_log, out, tag, T, a, t_start, a.groups)
            return
        for deg in a.bases:
            B, mn = basis(s, t, deg)
            fit = Fit(pool, pb, train, deg, data)
            t0 = time.time()
            sol = least_squares(fit.fun, np.zeros(B.shape[-1]), jac=fit.jac, max_nfev=a.max_nfev, xtol=1e-3,
                                ftol=1e-4, x_scale=1.0)
            c = sol.x
            log_a = B @ c
            pred = {f: E for f, _, E in pool.map(_solve, [(f, deg, c, None) for f in folders])}
            row = dict(basis=deg, n_coef=len(c), coef=c.tolist(), modes=mn, nfev=int(sol.nfev),
                       message=str(sol.message), fit_s=time.time() - t0,
                       train_cost=float(np.mean(sol.fun ** 2)),
                       alpha_interior_range=[float(np.exp(log_a[gin].min())), float(np.exp(log_a[gin].max()))],
                       alpha_interior_sd_log=float(log_a[gin].std()))
            for role, fs in (("train", train), ("held_out", held)):
                for k, (cname, _, _) in enumerate(COMPONENTS):
                    ex = [explained(data[f]["E"][k], pred[f][k], np.ones(data[f]["E"].shape[1], bool))[0] for f in fs]
                    ex1 = [explained(data[f]["E"][k], ref[f][k], np.ones(data[f]["E"].shape[1], bool))[0] for f in fs]
                    row[f"{role}_{cname}_explained"] = float(np.median(ex))
                    row[f"{role}_{cname}_explained_alpha1"] = float(np.median(ex1))
            if true_log is not None:
                tr, fl = true_log[gin] - true_log[gin].mean(), log_a[gin] - log_a[gin].mean()
                Bc = np.c_[np.ones(gin.sum()), B[gin]]                 # constant: the level is not compared
                coef_p, *_ = np.linalg.lstsq(Bc, tr, rcond=None)
                proj = Bc @ coef_p
                proj = proj - proj.mean()
                row.update(truth_sd_log=float(tr.std()), recovered_corr=float(np.corrcoef(tr, fl)[0, 1])
                           if fl.std() > 0 and tr.std() > 0 else float("nan"),
                           recovered_rms_log=float(np.sqrt(np.mean((fl - tr) ** 2))),
                           projection_rms_log=float(np.sqrt(np.mean((proj - tr) ** 2))),
                           fit_vs_projection_rms_log=float(np.sqrt(np.mean((fl - proj) ** 2))))
            results.append(dict(row, log_alpha=log_a))
            print(f"basis {deg} ({len(c)} coef, {sol.nfev} it, {row['fit_s']:.0f} s): alpha {row['alpha_interior_range'][0]:.2f}"
                  f"-{row['alpha_interior_range'][1]:.2f} | held-out explained "
                  + " ".join(f"{cn} {row[f'held_out_{cn}_explained']:+.2f} (alpha=1 {row[f'held_out_{cn}_explained_alpha1']:+.2f})"
                             for cn, _, _ in COMPONENTS)
                  + (f" | truth SD {row['truth_sd_log']:.3f}, recovered r {row['recovered_corr']:.2f}, "
                     f"RMS {row['recovered_rms_log']:.3f} (projection {row['projection_rms_log']:.3f}, "
                     f"fit vs projection {row['fit_vs_projection_rms_log']:.3f})" if true_log is not None else ""),
                  flush=True)

    with open(os.path.join(out, f"stage2_{tag}.json"), "w") as fh:
        json.dump(dict(test=T.name, source=tag, noise_px=a.noise_px if a.synthetic else None, training=train,
                       held_out=held, params=dict(zip(M.PARAM_NAMES, map(float, params))),
                       fits=[{k: v for k, v in r.items() if k != "log_alpha"} for r in results],
                       runtime_s=time.time() - t_start), fh, indent=2)
    np.savez_compressed(os.path.join(out, f"stage2_{tag}.npz"), gauss_points=g, interior=gin,
                        truth_log_alpha=true_log if true_log is not None else np.zeros(0),
                        **{f"log_alpha_basis{r['basis']}": r["log_alpha"] for r in results})

    # ---- figure: alpha maps (element means), truth first when synthetic; one color scale for all panels
    em = lambda L: np.exp(L.mean(2))  # noqa: E731
    panels = ([("truth", em(true_log))] if true_log is not None else []) + \
             [(f"basis {r['basis']} ({r['n_coef']} coef)", em(r["log_alpha"])) for r in results]
    allv = np.concatenate([np.log(v[gin.all(2)]) for _, v in panels])
    span = max(np.exp(np.max(np.abs(allv))) - 1, 0.05)
    Xg = pb.mesh.Xg
    Ze, Te = Xg[..., 0], Xg[..., 1]
    fig, axes = plt.subplots(1, len(panels) + 1, figsize=(3.6 * (len(panels) + 1), 4.2), constrained_layout=True)
    from matplotlib.colors import TwoSlopeNorm
    norm = TwoSlopeNorm(1.0, 1 - span, 1 + span)
    for ax, (title, V) in zip(axes, panels):
        pm = ax.pcolormesh(Ze, Te, V, cmap=D.diverging_cmap(), norm=norm, shading="flat")
        ez, et = pb.edge["z"], pb.edge["theta"]
        Zt = (pb.tips["z"][1] - pb.tips["z"][0]) * pb.mpp
        Tt = (pb.tips["theta"][1] - pb.tips["theta"][0]) * pb.mpp
        ax.add_patch(plt.Rectangle((ez, et), Zt - 2 * ez, Tt - 2 * et, fill=False, ls="--", lw=0.8))
        fig.colorbar(pm, ax=ax, shrink=0.8, label="α")
        ax.set_aspect("equal"); ax.set_title(title, fontsize=9)
        ax.set_xlabel("z (mm)", fontsize=8); ax.set_ylabel("θ (mm)", fontsize=8); ax.tick_params(labelsize=7)
    ax = axes[-1]
    degs = np.arange(len(results))
    for k, (cn, _, _) in enumerate(COMPONENTS):
        ax.plot(degs, [r[f"held_out_{cn}_explained"] for r in results], "o-", color=f"C{k}", label=f"{cn} held out")
        ax.axhline(results[0][f"held_out_{cn}_explained_alpha1"], color=f"C{k}", ls=":", lw=1)
    if true_log is not None:
        ax2 = ax.twinx()
        ax2.plot(degs, [r["recovered_rms_log"] for r in results], "ks--", label="RMS error of log α")
        ax2.plot(degs, [r["projection_rms_log"] for r in results], "k^:", label="best possible (projection)")
        ax2.set_ylabel("RMS error of log α (interior)", fontsize=8)
        ax2.legend(fontsize=7, loc="lower right")
    ax.set(xlabel="basis (total degree, or max degree z x θ)", ylabel="interior variance explained (median)",
           title="held out (dotted: α = 1)")
    ax.set_xticks(degs, [f"{r['basis']}\n({r['n_coef']})" for r in results])
    ax.legend(fontsize=7, loc="upper left"); ax.grid(alpha=0.3)
    what = f"synthetic '{a.synthetic}', noise {a.noise_px:g} px" if a.synthetic else "measured"
    fig.suptitle(f"{T.name}: stage 2, smooth stiffness factor α(x) from the interior strains ({what}); "
                 f"training Def11, Def14, Def18; dashed: interior", fontsize=10)
    fig.savefig(os.path.join(out, f"stage2_{tag}.png"), dpi=110)
    print(f"wrote {out} ({time.time() - t_start:.0f} s)")


if __name__ == "__main__":
    main()
