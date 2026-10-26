"""
Stage 1 check: does any uniform material explain the interior strain field, or only the bulk fit fails?

Stage 1 (stage1_uniform_fe.py) loads a uniform material, the bulk fit, by the measured boundary positions and finds that
it leaves most of the interior strain variance unexplained. A uniform material with a different response shape (more
anisotropic, other exponents) carries the uneven tine loading further into the interior and could, in principle,
produce bands. Here the uniform material is fitted to the interior strain fields instead of the forces.

Model: the fixed eight-parameter model, the same parameters at every point. With displacements prescribed on the whole
boundary the FE solution depends only on the response shape, not on the stiffness level, so C1_smc is held at its bulk
value (it sets the level) and the other seven are free: Cgr, C1_el, C1_col (relative to C1_smc), C2_el, C2_smc,
C2_col, gamma. Bounds as in the bulk fit (C1_el >= 0.1 kPa, exponents 1e-3 to 80, gamma 0-90 deg).

Objective: at the peak of the training protocols, every interior point (edge_mm) and strain component (E_zz,
E_thetatheta, E_ztheta), measured minus FE, each component and protocol divided by the measured interior SD so all
nine blocks weigh equally. FE and measured strains go through the same Gaussian window. Solved with scipy
least_squares (finite differences; each FE solve starts from the previous solution of the same protocol).

Training: Def11, Def14, Def18 (equibiaxial, theta-dominant, z-dominant). Held out: the other six peaks.
Reported for every protocol, bulk vs best uniform: the explained share as in stage 1. Also: the bulk force misfit
of the best uniform material (closed-form biaxial with the DIC affine stretch, level refitted, which is the only
freedom left), against the bulk fit's RMS.

    .venv/bin/python constitutive_fit/stage1b_best_uniform.py --test N1739/test1
Outputs: <test>/constitutive_fit/local/stage1b_best_uniform.json / .csv, stage1b_best_uniform.png (explained shares),
stage1b_peak_maps_{zz,thetatheta}.png (measured, best uniform FE, residual at every protocol's peak)
"""
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import constitutive as M  # noqa: E402
import dic2d as D  # noqa: E402
from membrane_fe import Mesh  # noqa: E402
from stage1_uniform_fe import COMPONENTS, INPUT, explained, fill_edge, strain  # noqa: E402

TRAIN = ("Def11", "Def14", "Def18")
FREE = (0, 1, 2, 5, 4, 6, 7)          # Cgr, C1_el, C2_el, C1_col, C2_smc, C2_col, gamma (C1_smc fixed)
LO = np.array([0, 0.1, 1e-3, 0, 1e-3, 0, 1e-3, 0])
HI = np.array([300, 3000, 80, 3000, 80, 3000, 80, 90])
MAXITER = 5000                        # membrane solver iterations per solve (a solve that stalls at an extreme field is cut off)
CGR_FLOOR = 1e-3                      # Cgr = exp(x) - CGR_FLOOR >= 0 (lets Cgr reach 0 on a log scale)


class Problem:
    """Peak frames of every protocol on the common FE mesh, with the measured interior strains."""

    def __init__(self, T, df):
        self.tips = T.config["tine_tips_px"]
        self.window_mm = testcfg.analysis()["window_mm"]
        self.edge = testcfg.edge_mm()
        self.H = T.config["thickness_mm"]["H"]
        self.prot = {}
        for s, g in df.groupby("set", sort=True):
            folder = g.folder.iloc[0]
            R, _ = D.load_track(T.results, folder)
            X0, shape = R["X0"], tuple(R["shape"])
            mpp = 1.0 / float(R["px_per_mm"])
            if not self.prot:
                self.X0, self.shape, self.mpp = X0, shape, mpp
                z, th = D.specimen_mm(X0[:, 0], X0[:, 1], self.tips, mpp)
                self.mesh = Mesh(np.stack([z.reshape(shape), th.reshape(shape)], -1)[::-1])
                self.alpha = np.ones(self.mesh.gauss_points.shape[:3])
                self.interior = D.interior_mask(X0, shape, self.tips, mpp, self.edge)
            times = np.round(R["times"], 1)
            k = int(np.nonzero(times == round(g.t.iloc[-1], 1))[0][0])
            xp, valid = R["pos"][k], R["valid"][k]
            xm = self.to_mesh(np.where(valid[:, None], xp, np.nan))
            vm = valid.reshape(shape)[::-1]
            xb = xm.copy()
            for sl in ((0, slice(None)), (-1, slice(None)), (slice(None), 0), (slice(None), -1)):
                xb[sl] = fill_edge(xm[sl], vm[sl])
            bnd = ~self.mesh.free
            Fa = D.affine_fit(self.mesh.Xg[bnd], xb[bnd])
            x0 = self.mesh.Xg @ Fa.T + xb[bnd].mean(0) - self.mesh.Xg[bnd].mean(0) @ Fa.T
            Em, okm = strain(X0, xp, valid, shape, self.window_mm, mpp)
            m = self.interior & okm
            sd = {c: float(np.nanstd(Em[..., i, j][m])) for c, i, j in COMPONENTS}
            self.prot[folder] = dict(xb=xb, x0=x0, xp=xp, Em=Em, ok=okm, mask=m, sd=sd, warm=None)

    def to_mesh(self, x_px):
        z, th = D.specimen_mm(x_px[:, 0], x_px[:, 1], self.tips, self.mpp)
        return np.stack([z.reshape(self.shape), th.reshape(self.shape)], -1)[::-1]

    def to_px(self, xm):
        xm = xm[::-1].reshape(-1, 2)
        return np.stack([xm[:, 0] / self.mpp + self.tips["z"][0], self.tips["theta"][1] - xm[:, 1] / self.mpp], 1)

    def solve(self, folder, params, warm=True, alpha=None, start=None):
        """FE solution for one protocol's peak. alpha: stiffness factor at the Gauss points (default 1); start:
        initial node positions (default: the last warm solution if warm, else the affine map of the boundary)."""
        p = self.prot[folder]
        alpha = self.alpha if alpha is None else alpha
        if start is None:
            start = p["warm"] if (warm and p["warm"] is not None) else p["x0"]
        xf, info = self.mesh.solve(p["xb"], alpha, params, self.H, start, gtol=1e-8, maxiter=MAXITER)
        if not info["success"] or info["nit"] == 0:
            xs = p["x0"]
            for b in np.linspace(0.25, 1.0, 4):
                xs, info = self.mesh.solve(p["x0"] + b * (p["xb"] - p["x0"]), alpha, params, self.H, xs,
                                           gtol=1e-8, maxiter=MAXITER)
            xf = xs
        if warm:
            p["warm"] = xf
        Ef, _ = strain(self.X0, self.to_px(xf), np.ones(len(self.X0), bool), self.shape, self.window_mm, self.mpp)
        return xf, Ef

    def residual(self, folder, Ef):
        p = self.prot[folder]
        return np.concatenate([np.nan_to_num((p["Em"][..., i, j] - Ef[..., i, j])[p["mask"]]) / p["sd"][c]
                               for c, i, j in COMPONENTS])


def unpack(x, base):
    p = base.copy()
    for k, idx in enumerate(FREE):
        if idx == 7:
            p[idx] = x[k]
        elif idx == 0:
            p[idx] = np.exp(x[k]) - CGR_FLOOR
        else:
            p[idx] = np.exp(x[k])
    return np.clip(p, LO, HI)


def pack(p):
    return np.array([p[i] if i == 7 else np.log(p[i] + CGR_FLOOR) if i == 0 else np.log(max(p[i], LO[i], 1e-6))
                     for i in FREE])


def bounds():
    lo, hi = [], []
    for i in FREE:
        if i == 7:
            lo.append(LO[i]); hi.append(HI[i])
        elif i == 0:
            lo.append(np.log(CGR_FLOOR)); hi.append(np.log(HI[i] + CGR_FLOOR))
        else:
            lo.append(np.log(max(LO[i], 1e-3))); hi.append(np.log(HI[i]))
    return np.array(lo), np.array(hi)


def bulk_rms(df, params, w_t, w_z, H):
    """Bulk force misfit of a response shape with its level refitted (stress is linear in Cgr and the C1s)."""
    lz, lt = df[f"{INPUT}_z"].to_numpy(), df[f"{INPUT}_t"].to_numpy()
    t = np.r_[df.fz.to_numpy() * lz / (w_t * H), df.ft.to_numpy() * lt / (w_z * H)]
    m = np.r_[M.cauchy_biaxial(lz, lt, params)]
    s = float(m @ t / (m @ m))
    return float(np.sqrt(np.mean((s * m - t) ** 2))), s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", default="tine_span")
    ap.add_argument("--max-nfev", type=int, default=40, help="least-squares iterations (each: 3 FE solves)")
    a = ap.parse_args()
    T = testcfg.get()
    np.seterr(over="ignore", invalid="ignore")
    out = T.local_dir
    os.makedirs(out, exist_ok=True)
    df = pd.read_csv(os.path.join(T.fit_dir, "stretch_inputs_data.csv"))
    P = pd.read_csv(os.path.join(T.fit_dir, "params_by_stretch_input.csv")).set_index("input")
    bulk = P.loc[INPUT, list(M.PARAM_NAMES)].to_numpy(float)
    w_t, w_z = T.config["widths_mm"][a.width]
    t_start = time.time()
    pb = Problem(T, df)
    train = [f for f in pb.prot if f[:5] in TRAIN]
    print(f"{T.name}: training on {', '.join(train)}; {sum(pb.prot[f]['mask'].sum() for f in train)} interior "
          f"points x 3 components", flush=True)

    calls = dict(n=0, best=np.inf)

    def fun(x):
        p = unpack(x, bulk)
        r = np.concatenate([pb.residual(f, pb.solve(f, p)[1]) for f in train])
        calls["n"] += 1
        cost = float(np.mean(r ** 2))
        if cost < calls["best"]:
            calls["best"] = cost
            print(f"  eval {calls['n']:4d} ({time.time() - t_start:6.0f} s): mean squared normalized residual "
                  f"{cost:.4f} | " + " ".join(f"{n}={v:.3g}" for n, v in zip(M.PARAM_NAMES, p)), flush=True)
        return r

    lo, hi = bounds()
    x0 = np.clip(pack(bulk), lo + 1e-6, hi - 1e-6)
    sol = least_squares(fun, x0, bounds=(lo, hi), diff_step=1e-2, x_scale=np.where(np.array(FREE) == 7, 10.0, 1.0),
                        max_nfev=a.max_nfev, xtol=1e-4, ftol=1e-4)
    best = unpack(sol.x, bulk)
    print(f"least squares: {sol.message} ({sol.nfev} evaluations, {time.time() - t_start:.0f} s)", flush=True)

    # every protocol, bulk vs best uniform (cold starts, so both are solved the same way)
    rows, maps = [], {}
    for f in pb.prot:
        p = pb.prot[f]
        for name, prm in (("bulk", bulk), ("best_uniform", best)):
            xf, Ef = pb.solve(f, prm, warm=False)
            row = dict(folder=f, material=name, role="train" if f in train else "held out")
            for c, i, j in COMPONENTS:
                ex, pat, _, _ = explained(p["Em"][..., i, j], Ef[..., i, j], p["mask"])
                row[f"{c}_explained"], row[f"{c}_pattern_explained"] = ex, pat
            rows.append(row)
            if name == "best_uniform":
                maps[f] = dict(Ef=Ef, xf=pb.to_px(xf))
        print(f"  {f}: " + " | ".join(
            f"{r['material']} " + " ".join(f"{c} {r[f'{c}_explained']:+.2f}" for c, _, _ in COMPONENTS)
            for r in rows[-2:]), flush=True)
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(out, "stage1b_best_uniform.csv"), index=False)
    H = T.config["thickness_mm"]["H"]
    rms_bulk, _ = bulk_rms(df, bulk, w_t, w_z, H)
    rms_best, scale = bulk_rms(df, best, w_t, w_z, H)
    summ = {}
    for name in ("bulk", "best_uniform"):
        for role in ("train", "held out"):
            d = S[(S.material == name) & (S.role == role)]
            summ[f"{name}, {role}"] = {c: float(d[f"{c}_explained"].median()) for c, _, _ in COMPONENTS}
    res = dict(test=T.name, training=train, bulk=dict(zip(M.PARAM_NAMES, map(float, bulk))),
               best_uniform=dict(zip(M.PARAM_NAMES, map(float, best))), level_scale_for_forces=scale,
               bulk_force_rms_kPa=dict(bulk=rms_bulk, best_uniform_level_refitted=rms_best),
               explained_median=summ, optimizer=dict(message=str(sol.message), nfev=int(sol.nfev),
                                                     cost=float(np.mean(sol.fun ** 2))),
               runtime_s=time.time() - t_start)
    with open(os.path.join(out, "stage1b_best_uniform.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print(json.dumps({k: res[k] for k in ("best_uniform", "bulk_force_rms_kPa", "explained_median")}, indent=1))

    # ---- figure: explained share per protocol, bulk vs best uniform
    folders = list(pb.prot)
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.5), constrained_layout=True)
    for ax, (c, _, _) in zip(axes, COMPONENTS):
        for name, mk in (("bulk", "o-"), ("best_uniform", "s-")):
            d = S[S.material == name].set_index("folder").loc[folders]
            ax.plot(np.arange(len(folders)), d[f"{c}_explained"], mk, label=name.replace("_", " "))
        for k, f in enumerate(folders):
            if f in train:
                ax.axvspan(k - 0.4, k + 0.4, color="0.85", zorder=0)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(np.arange(len(folders)), [f[:5] for f in folders], fontsize=8)
        ax.set_ylim(max(-2, ax.get_ylim()[0]), 1.0)
        ax.set(title=f"{c}: interior variance explained", ylabel="explained share")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, title="gray: training protocols", title_fontsize=8)
    fig.suptitle(f"{T.name}: uniform material fitted to the interior strains (training Def11, Def14, Def18) "
                 f"vs the bulk fit")
    fig.savefig(os.path.join(out, "stage1b_best_uniform.png"), dpi=110)
    plt.close(fig)

    # ---- peak maps with the best uniform material (as stage1_peak_maps_*.png)
    for c, i, j in COMPONENTS[:2]:
        V = {f: dict(Em=pb.prot[f]["Em"][..., i, j], Ef=maps[f]["Ef"][..., i, j]) for f in folders}
        for f in folders:
            V[f]["res"] = V[f]["Em"] - V[f]["Ef"]
        both = np.concatenate([np.r_[v["Em"][pb.interior], v["Ef"][pb.interior]] for v in V.values()])
        lo_, hi_ = np.nanpercentile(both, 1), np.nanpercentile(both, 99)
        r = np.nanpercentile(np.abs(np.concatenate([v["res"][pb.interior] for v in V.values()])), 99)
        fig, axes = plt.subplots(3, len(folders), figsize=(2.9 * len(folders), 9.3), squeeze=False,
                                 constrained_layout=True)
        for col, f in enumerate(folders):
            for row, (kind, pos) in enumerate((("Em", pb.prot[f]["xp"]), ("Ef", maps[f]["xf"]),
                                               ("res", pb.prot[f]["xp"]))):
                ax = axes[row, col]
                zd, td = D.specimen_mm(pos[:, 0], pos[:, 1], pb.tips, pb.mpp)
                Zd, Td = zd.reshape(pb.shape), td.reshape(pb.shape)
                pm = ax.pcolormesh(D.cell_corners(Zd), D.cell_corners(Td), V[f][kind],
                                   cmap=D.diverging_cmap() if kind == "res" else "viridis",
                                   vmin=-r if kind == "res" else lo_, vmax=r if kind == "res" else hi_,
                                   shading="flat", edgecolors="face")
                ax.contour(Zd, Td, pb.interior.astype(float), levels=[0.5], colors="k" if kind == "res" else "w",
                           linewidths=0.7, linestyles="--")
                fig.colorbar(pm, ax=ax, shrink=0.75, pad=0.02).ax.tick_params(labelsize=6)
                ax.set_aspect("equal")
                ax.tick_params(labelsize=6)
                if row == 0:
                    ax.set_title(f.replace("_", " ", 1) + (" (train)" if f in train else ""), fontsize=8)
                if col == 0:
                    ax.set_ylabel(f"{['measured', 'best uniform FE', 'measured − FE'][row]}\nθ (mm)", fontsize=8)
                if row == 2:
                    ax.set_xlabel("z (mm)", fontsize=7)
        fig.suptitle(f"{T.name}: {c} at each protocol's peak, deformed configuration; uniform material fitted to "
                     f"the interior strains of Def11, Def14, Def18", fontsize=10)
        fig.savefig(os.path.join(out, f"stage1b_peak_maps_{c.split('_')[1]}.png"), dpi=110)
        plt.close(fig)
    print(f"wrote {out} ({time.time() - t_start:.0f} s)")


if __name__ == "__main__":
    main()
