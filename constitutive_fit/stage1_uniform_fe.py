"""
Stage 1 of the local-property stages: how much of the measured strain field a uniform material already explains.

A finite element model of the tracked domain (membrane_fe.py; nodes = the DIC grid, Q1 elements) made of one uniform
material, the test's bulk fit (DIC affine F, params_by_stretch_input.csv), is loaded by the measured positions of the
outermost ring of grid points (invalid boundary points filled along the edge). With displacements prescribed on the
whole boundary, the solution does not depend on the stiffness level or on H, only on the shape of the response
(exponents, gamma, ratios between the constants). The interior is solved from a smooth start (affine map of the
boundary), never from the measured interior positions.

The FE and the measured positions go through the same strain pipeline (Gaussian window_mm of analysis_settings.json).
Per frame and strain component (E_zz, E_thetatheta, E_ztheta), in the interior (edge_mm from the tine-tip lines) and
in the edge band:
    explained share = 1 - sum (E_meas - E_FE)^2 / sum (E_meas - mean E_meas)^2
The strain field at the first frame of each protocol (rest, after bridging) gives an upper bound on the noise
variance. Whatever variance the uniform material leaves unexplained, above that bound, is the evidence for material
heterogeneity (stiffness or thickness) that the later stages try to map.

Also reported, as a consistency check that involves H and the widths: the FE force resultants across the middle
element column (z) and the middle element row (theta), integral of the first PK membrane tension along the cut, against
the measured forces times the cut's share of the tine span. Mid-domain cuts are used because the edge reactions are
dominated by point loads where tines locally overstretch the tissue.

Frames: per protocol Def11-Def19, the loading-branch frames at 50%, 75% and 100% of the branch.

    .venv/bin/python constitutive_fit/stage1_uniform_fe.py --test N1739/test1
Reads the full-domain tracks in <test>/DIC_results/ and <test>/constitutive_fit/stretch_inputs_data.csv and
params_by_stretch_input.csv. Outputs: <test>/constitutive_fit/local/stage1_uniform_fe.csv / .json, stage1_peak_maps.npz,
stage1_peak_maps_{zz,thetatheta,ztheta}.png (every protocol's peak, deformed), stage1_explained.png
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import constitutive as M  # noqa: E402
import dic2d as D  # noqa: E402
from membrane_fe import Mesh  # noqa: E402

INPUT = "dic_affine"
FRACTIONS = (0.5, 0.75, 1.0)
COMPONENTS = (("E_zz", 0, 0), ("E_thetatheta", 1, 1), ("E_ztheta", 0, 1))


def fill_edge(v, ok):
    """Fill invalid entries of a 1D edge array (n, 2) by linear interpolation along the edge."""
    v = v.copy()
    i = np.arange(len(v))
    for c in range(2):
        v[~ok, c] = np.interp(i[~ok], i[ok], v[ok, c])
    return v


def strain(X0, x_px, valid, shape, window_mm, mpp):
    F, ok = D.deformation_gradient_field_mm(X0, x_px, valid, shape, window_mm, mpp)
    F = D.to_specimen(F)
    return D.green_lagrange(F), ok


def explained(Em, Ef, mask):
    m = mask & np.isfinite(Em) & np.isfinite(Ef)
    d = Em[m] - Em[m].mean()
    r = Em[m] - Ef[m]
    return (1 - np.sum(r ** 2) / np.sum(d ** 2), 1 - np.sum((r - r.mean()) ** 2) / np.sum(d ** 2),
            np.var(Em[m]), np.mean(r ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", default="tine_span")
    a = ap.parse_args()
    T = testcfg.get()
    np.seterr(over="ignore", invalid="ignore")
    out = T.local_dir
    os.makedirs(out, exist_ok=True)
    df = pd.read_csv(os.path.join(T.fit_dir, "stretch_inputs_data.csv"))
    P = pd.read_csv(os.path.join(T.fit_dir, "params_by_stretch_input.csv")).set_index("input")
    params = P.loc[INPUT, list(M.PARAM_NAMES)].to_numpy(float)
    H = T.config["thickness_mm"]["H"]
    w_t, w_z = T.config["widths_mm"][a.width]
    tips = T.config["tine_tips_px"]
    window_mm = testcfg.analysis()["window_mm"]
    edge = testcfg.edge_mm()

    rows, maps = [], {}
    mesh = None
    t_start = time.time()
    for s, g in df.groupby("set", sort=True):
        folder = g.folder.iloc[0]
        R, _ = D.load_track(T.results, folder)
        X0, shape = R["X0"], tuple(R["shape"])
        mpp = 1.0 / float(R["px_per_mm"])
        times = np.round(R["times"], 1)
        ny, nx = shape
        if mesh is None:
            z, th = D.specimen_mm(X0[:, 0], X0[:, 1], tips, mpp)
            Xg = np.stack([z.reshape(shape), th.reshape(shape)], -1)[::-1]   # rows: theta increasing
            mesh = Mesh(Xg)
            alpha = np.ones(mesh.gauss_points.shape[:3])
            interior = D.interior_mask(X0, shape, tips, mpp, edge)
            Lz, Lt = Xg[0, -1, 0] - Xg[0, 0, 0], Xg[-1, 0, 1] - Xg[0, 0, 1]

        def to_mesh(x_px):
            z, th = D.specimen_mm(x_px[:, 0], x_px[:, 1], tips, mpp)
            return np.stack([z.reshape(shape), th.reshape(shape)], -1)[::-1]

        def to_px(xm):
            xm = xm[::-1].reshape(-1, 2)
            return np.stack([xm[:, 0] / mpp + tips["z"][0], tips["theta"][1] - xm[:, 1] / mpp], 1)

        tl = g.t.to_numpy()
        # rest frame of the protocol: noise bound
        k_rest = int(np.nonzero(times == round(tl[0], 1))[0][0])
        Er, okr = strain(X0, R["pos"][k_rest], R["valid"][k_rest], shape, window_mm, mpp)
        noise = {c: float(np.nanvar(Er[..., i, j][interior & okr])) for c, i, j in COMPONENTS}
        for frac in FRACTIONS:
            r_ = g.iloc[int(round(frac * (len(g) - 1)))]
            k = int(np.nonzero(times == round(r_.t, 1))[0][0])
            xp, valid = R["pos"][k], R["valid"][k]
            xm = to_mesh(np.where(valid[:, None], xp, np.nan))
            vm = valid.reshape(shape)[::-1]
            xb = xm.copy()
            for sl in ((0, slice(None)), (-1, slice(None)), (slice(None), 0), (slice(None), -1)):
                xb[sl] = fill_edge(xm[sl], vm[sl])
            bnd = ~mesh.free
            Fa = D.affine_fit(mesh.Xg[bnd], xb[bnd])
            c0 = xb[bnd].mean(0) - mesh.Xg[bnd].mean(0) @ Fa.T
            x0 = mesh.Xg @ Fa.T + c0
            t0 = time.time()
            xf, info = mesh.solve(xb, alpha, params, H, x0, gtol=1e-7)
            if not info["success"] or info["nit"] == 0:
                # a kinked boundary can invert an element of the smooth start: step the boundary from its affine
                # map to the measured positions, each step started from the previous solution
                xs = x0
                for b in np.linspace(0.25, 1.0, 4):
                    xs, info = mesh.solve(x0 + b * (xb - x0), alpha, params, H, xs, gtol=1e-7)
                xf = xs
            N = M.first_pk_membrane(mesh.F(xf), params, H).mean(2)                  # (ne_r, ne_c, 2, 2), mN/mm
            fz_fe = N[:, N.shape[1] // 2, 0, 0].sum() * abs(mesh.hy)
            ft_fe = N[N.shape[0] // 2, :, 1, 1].sum() * abs(mesh.hx)
            fz_meas, ft_meas = r_.fz * Lt / w_t, r_.ft * Lz / w_z
            Em, okm = strain(X0, xp, valid, shape, window_mm, mpp)
            Ef, okf = strain(X0, to_px(xf), np.ones(len(X0), bool), shape, window_mm, mpp)
            dist = np.linalg.norm(to_px(xf) - xp, axis=1).reshape(shape) * mpp * 1e3          # um
            Fa_m = D.affine_fit(X0[valid], xp[valid])
            het = np.linalg.norm(xp - (X0 - X0[valid].mean(0)) @ Fa_m.T - xp[valid].mean(0), axis=1)
            het = het.reshape(shape) * mpp * 1e3
            row = dict(set=int(s), folder=folder, frac=frac, t=float(r_.t), fz_mN=float(r_.fz), ft_mN=float(r_.ft),
                       lambda_z=float(r_[f"{INPUT}_z"]), lambda_theta=float(r_[f"{INPUT}_t"]),
                       fe_converged=info["success"], fe_iterations=info["nit"],
                       fe_max_free_force_mN=info["max_free_force"],
                       cut_ratio_z=float(fz_fe / fz_meas) if fz_meas > 0 else np.nan,
                       cut_ratio_theta=float(ft_fe / ft_meas) if ft_meas > 0 else np.nan,
                       disp_resid_um_median=float(np.nanmedian(dist[interior & valid.reshape(shape)])),
                       disp_nonaffine_um_median=float(np.nanmedian(het[interior & valid.reshape(shape)])))
            for c, i, j in COMPONENTS:
                for region, mk in (("interior", interior), ("edge", ~interior)):
                    ex, exp_pat, var, mse = explained(Em[..., i, j], Ef[..., i, j], mk & okm)
                    row[f"{c}_{region}_explained"] = ex
                    row[f"{c}_{region}_pattern_explained"] = exp_pat
                    row[f"{c}_{region}_var"] = var
                    row[f"{c}_{region}_resid_var"] = mse
                row[f"{c}_noise_var"] = noise[c]
            rows.append(row)
            if frac == 1.0:
                maps[folder] = dict(Em=Em, Ef=Ef, xp=xp, xf=to_px(xf))
            print(f"  {folder} {int(frac * 100):3d}%: FE {time.time() - t0:5.1f} s ({info['nit']} it) | interior explained "
                  + " ".join(f"{c} {row[f'{c}_interior_explained']:+.2f} ({row[f'{c}_interior_pattern_explained']:+.2f})"
                             for c, _, _ in COMPONENTS)
                  + f" | mid-cut force / measured z {row['cut_ratio_z']:.2f} θ {row['cut_ratio_theta']:.2f}",
                  flush=True)

    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(out, "stage1_uniform_fe.csv"), index=False)
    peak = S[S.frac == 1.0]
    summ = {}
    for c, _, _ in COMPONENTS:
        for region in ("interior", "edge"):
            v, rv = peak[f"{c}_{region}_var"], peak[f"{c}_{region}_resid_var"]
            nv = peak[f"{c}_noise_var"]
            summ[f"{c}_{region}"] = dict(
                explained_median=float(peak[f"{c}_{region}_explained"].median()),
                pattern_explained_median=float(peak[f"{c}_{region}_pattern_explained"].median()),
                explained_range=[float(peak[f"{c}_{region}_explained"].min()),
                                 float(peak[f"{c}_{region}_explained"].max())],
                unexplained_above_noise_share_median=float(np.median((rv - nv).clip(0) / v)))
    with open(os.path.join(out, "stage1_uniform_fe.json"), "w") as f:
        json.dump(dict(test=T.name, params=dict(zip(M.PARAM_NAMES, map(float, params))), input=INPUT,
                       window_mm=window_mm, edge_mm=edge, fractions=FRACTIONS, peak_summary=summ,
                       cut_ratio_peak_median=dict(z=float(peak.cut_ratio_z.median()),
                                                       theta=float(peak.cut_ratio_theta.median())),
                       runtime_s=time.time() - t_start), f, indent=2)
    np.savez_compressed(os.path.join(out, "stage1_peak_maps.npz"),
                        **{f"{d}_{k}": v for d, m in maps.items() for k, v in m.items()})

    # ---- peak maps, one figure per strain component: rows measured / uniform FE / measured - FE, one column per
    # protocol, each drawn on its deformed configuration (measured and residual on the measured positions, FE on the
    # FE positions); one color scale for measured and FE across all protocols, one symmetric scale for the residual
    tag = {"E_zz": "zz", "E_thetatheta": "thetatheta", "E_ztheta": "ztheta"}
    for c, i, j in COMPONENTS:
        V = {f: dict(Em=m["Em"][..., i, j], Ef=m["Ef"][..., i, j], res=m["Em"][..., i, j] - m["Ef"][..., i, j])
             for f, m in maps.items()}
        both = np.concatenate([np.r_[v["Em"][interior], v["Ef"][interior]] for v in V.values()])
        lo, hi = np.nanpercentile(both, 1), np.nanpercentile(both, 99)
        r = np.nanpercentile(np.abs(np.concatenate([v["res"][interior] for v in V.values()])), 99)
        fig, axes = plt.subplots(3, len(maps), figsize=(2.9 * len(maps), 9.3), squeeze=False, constrained_layout=True)
        for col, (folder, m) in enumerate(maps.items()):
            for row, (kind, pos) in enumerate((("Em", m["xp"]), ("Ef", m["xf"]), ("res", m["xp"]))):
                ax = axes[row, col]
                zd, td = D.specimen_mm(pos[:, 0], pos[:, 1], tips, mpp)
                Zd, Td = zd.reshape(shape), td.reshape(shape)
                pm = ax.pcolormesh(D.cell_corners(Zd), D.cell_corners(Td), V[folder][kind],
                                   cmap=D.diverging_cmap() if kind == "res" else "viridis",
                                   vmin=-r if kind == "res" else lo, vmax=r if kind == "res" else hi,
                                   shading="flat", edgecolors="face")
                ax.contour(Zd, Td, interior.astype(float), levels=[0.5], colors="k" if kind == "res" else "w",
                           linewidths=0.7, linestyles="--")
                fig.colorbar(pm, ax=ax, shrink=0.75, pad=0.02).ax.tick_params(labelsize=6)
                ax.set_aspect("equal")
                ax.tick_params(labelsize=6)
                if row == 0:
                    ax.set_title(folder.replace("_", " ", 1), fontsize=8)
                if col == 0:
                    ax.set_ylabel(f"{['measured', 'uniform FE', 'measured − FE'][row]}\nθ (mm)", fontsize=8)
                if row == 2:
                    ax.set_xlabel("z (mm)", fontsize=7)
        fig.suptitle(f"{T.name}: {c} at each protocol's peak, deformed configuration; uniform bulk material, "
                     f"measured boundary positions (dashed: interior, edge_mm {edge['z']:g} / {edge['theta']:g} mm)",
                     fontsize=10)
        fig.savefig(os.path.join(out, f"stage1_peak_maps_{tag[c]}.png"), dpi=110)
        plt.close(fig)

    # ---- figure 2: explained share by protocol and load level; mid-cut forces
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5), constrained_layout=True)
    lab = [f"{f[:5]}" for f in peak.folder]
    for ax, (c, _, _) in zip(axes[:3], COMPONENTS):
        for frac, mk in zip(FRACTIONS, ("o", "s", "^")):
            d = S[S.frac == frac]
            ax.plot(np.arange(len(d)), d[f"{c}_interior_explained"], mk + "-", label=f"interior, {int(frac * 100)}%")
        ax.plot(np.arange(len(peak)), peak[f"{c}_interior_pattern_explained"], "C3*:", label="interior pattern only, 100%")
        ax.plot(np.arange(len(peak)), peak[f"{c}_edge_explained"], "kx--", label="edge band, 100%")
        ceil = 1 - peak[f"{c}_noise_var"] / peak[f"{c}_interior_var"]
        ax.plot(np.arange(len(peak)), ceil, color="0.5", lw=3, alpha=0.5, label="ceiling (rest-frame noise bound)")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(np.arange(len(peak)), lab, fontsize=8)
        ax.set_ylim(-0.5, 1.05)
        ax.set(title=f"{c}: variance explained", ylabel="explained share")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    ax = axes[3]
    ax.plot(np.arange(len(peak)), peak.cut_ratio_z, "o-", label="z")
    ax.plot(np.arange(len(peak)), peak.cut_ratio_theta, "s-", label="θ")
    ax.axhline(1, color="k", lw=0.6)
    ax.set_xticks(np.arange(len(peak)), lab, fontsize=8)
    ax.set(title="FE mid-cut force / measured share (peak)", ylabel="ratio")
    ax.legend(); ax.grid(alpha=0.3)
    fig.suptitle(f"{T.name}: stage 1, uniform bulk material ({INPUT}) with measured boundary positions")
    fig.savefig(os.path.join(out, "stage1_explained.png"), dpi=110)
    print(json.dumps(summ, indent=1))
    print(f"wrote {out} ({time.time() - t_start:.0f} s)")


if __name__ == "__main__":
    main()
