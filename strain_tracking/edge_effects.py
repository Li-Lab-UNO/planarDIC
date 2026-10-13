"""
How far tine effects reach into the specimen: signatures only the tines produce, measured against
distance from each rest-state tine-tip line, to choose the interior edge distance (edge_mm).

At each protocol's peak (frame of maximum mean axial stretch), F from the analysis_settings.json
window. For each edge (left / right rakes pull axially; top / bottom rakes circumferentially),
on every grid line parallel to that edge, at distance d from its tine-tip line:

  tine-pitch amplitude   the part of the stretch normal to the edge (lambda_z at the left/right
                         edges, lambda_theta at the top/bottom edges) that follows the tine
                         positions along the edge: regression of the detrended profile on a
                         template that is +1 at each measured tine and -1 midway between tines.
                         Discrete grips produce this pattern; tissue bands do not.
  shear                  mean |E_ztheta| along the line (the rakes' biaxial pull is shear-free).

Reference levels come from the specimen's core (lines at least 3 mm from both of that edge's tine
lines). The strain-window floor is the distance within which the window is cut off by the grid
edge (2.5 sigma of the Gaussian, sigma = window_mm / sqrt(12), plus the grid's inset from the tips).

The interior edge itself (edge_mm of analysis_settings.json) is the user's choice: this script reports the
distances measured here next to the current setting (dashed green on the plots, columns current_edge_mm and
within_current in the .csv) and changes nothing. To change the interior, edit edge_mm, or pass
--edge-mm Z THETA to peak_j_figure.py / compare_peak_j.py for one run.

    ../.venv/bin/python edge_effects.py --tests N1739/test1 N1690/test1
Output: comparisons/edge_effects_<tests>.png / .csv / _profiles.csv
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import testcfg  # noqa: E402

DEFS = [f"Def{i}" for i in range(11, 20)]
FAMILY = {**{d: "equibiaxial" for d in ("Def11", "Def15", "Def19")},
          **{d: "axial reduced" for d in ("Def12", "Def13", "Def14")},
          **{d: "circumferential reduced" for d in ("Def16", "Def17", "Def18")}}
FAM_COLOR = {"equibiaxial": "#1f1f1e", "axial reduced": "#2a78d6", "circumferential reduced": "#eb6834"}
EDGES = {"left": ("axial", 0), "right": ("axial", 0), "bottom": ("circumferential", 1), "top": ("circumferential", 1)}
MUTED, GRID = "#6b6b67", "#d9d9d6"


def peak_frame(R):
    X0, best = R["X0"], (-1.0, 0)
    for k in range(len(R["times"])):
        v = R["valid"][k]
        if v.sum() < 10:
            continue
        A = D.affine_fit(X0[v], R["pos"][k][v])
        lz = np.sqrt((A.T @ A)[0, 0])
        if lz > best[0]:
            best = (lz, k)
    return best[1]


def fields(T, defn, window_mm):
    R, _ = D.load_track(T.images, T.folder(defn), "tine_region")
    X0, sh = R["X0"], tuple(R["shape"])
    mpp = D.load_scale(T.scale_file)["mm_per_px"]
    k = peak_frame(R)
    F, ok = D.deformation_gradient_field_mm(X0, R["pos"][k], R["valid"][k], sh, window_mm, mpp)
    F = D.to_specimen(F)
    C = np.einsum("...ki,...kj->...ij", F, F)
    ok = ok & (R["valid"].mean(0).reshape(sh) >= 0.5)
    nan = lambda a: np.where(ok, a, np.nan)  # noqa: E731
    z, th = D.specimen_mm(X0[:, 0].reshape(sh), X0[:, 1].reshape(sh), T.setting("tine_tips_px"), mpp)
    return dict(lz=nan(np.sqrt(C[..., 0, 0])), lt=nan(np.sqrt(C[..., 1, 1])),
                ezt=nan(0.5 * C[..., 0, 1]), z=z, th=th, step_mm=float(R["step"]) * mpp, mpp=mpp)


def tine_positions(T, side, mpp):
    """Tine tip positions along an edge (mm, in that edge's coordinate: theta for left/right, z for top/bottom)."""
    tips = T.setting("tine_tips_px")
    c = np.array([p[0] for p in T.setting("calibration_tines_px", side, "tips")], float)
    return np.sort((tips["theta"][1] - c) * mpp if side in ("left", "right") else (c - tips["z"][0]) * mpp)


def template(s, tines):
    """+1 at each tine, -1 midway between neighbors (piecewise cosine in the local spacing)."""
    g = np.full(s.shape, np.nan)
    for a, b in zip(tines[:-1], tines[1:]):
        m = (s >= a) & (s <= b)
        g[m] = np.cos(2 * np.pi * (s[m] - a) / (b - a))
    return g


def edge_profiles(f, T, spans):
    """Per edge and grid line: distance d (mm), tine-pitch amplitude of the normal stretch, mean |E_zt|."""
    rows = []
    for side, (kind, _) in EDGES.items():
        normal = f["lz"] if kind == "axial" else f["lt"]
        tines = tine_positions(T, side, f["mpp"])
        if side in ("left", "right"):
            lines = [(normal[:, j], f["th"][:, j], f["ezt"][:, j], np.nanmean(f["z"][:, j])) for j in range(normal.shape[1])]
            L = spans[0]
        else:
            lines = [(normal[i, :], f["z"][i, :], f["ezt"][i, :], np.nanmean(f["th"][i, :])) for i in range(normal.shape[0])]
            L = spans[1]
        for val, s, ezt, pos in lines:
            d = pos if side in ("left", "bottom") else L - pos
            g = template(s, tines)
            m = np.isfinite(val) & np.isfinite(g)
            amp = np.nan
            if m.sum() >= 8:
                y = val[m] / np.nanmean(val[m]) - 1
                Xd = np.c_[np.ones(m.sum()), s[m] - s[m].mean(), g[m]]   # offset + linear trend + tine template
                amp = float(np.linalg.lstsq(Xd, y, rcond=None)[0][2])
            rows.append(dict(edge=side, kind=kind, d_mm=float(d), far_mm=float(L - d),
                             tine_amp=amp, abs_ezt=float(np.nanmean(np.abs(ezt)))))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tests", nargs="+", default=["N1739/test1", "N1690/test1"])
    ap.add_argument("--dmax", type=float, default=4.0)
    a = ap.parse_args()
    S = testcfg.analysis()
    edge = testcfg.edge_mm()
    sigma = S["window_mm"] / np.sqrt(12)
    allrows, info = [], {}
    for name in a.tests:
        testcfg.select(name)
        T = testcfg.get()
        tips = T.setting("tine_tips_px")
        mpp = D.load_scale(T.scale_file)["mm_per_px"]
        spans = ((tips["z"][1] - tips["z"][0]) * mpp, (tips["theta"][1] - tips["theta"][0]) * mpp)
        R0, _ = D.load_track(T.images, T.folder("Def11"), "tine_region")
        inset = float(np.min([(R0["X0"][:, 0].min() - tips["z"][0]), (tips["z"][1] - R0["X0"][:, 0].max()),
                              (R0["X0"][:, 1].min() - tips["theta"][0]), (tips["theta"][1] - R0["X0"][:, 1].max())])) * mpp
        info[name] = dict(floor=inset + 2.5 * sigma, subset=float(R0.get("subset_mm", np.nan)),
                          step=float(R0.get("step_mm", np.nan)), sample=T.sample, test=T.name.split("_")[0])
        for d in DEFS:
            for r in edge_profiles(fields(T, d, S["window_mm"]), T, spans):
                allrows.append(dict(test=name, protocol=d, family=FAMILY[d], **r))
        print(name, "done", flush=True)
    df = pd.DataFrame(allrows)

    # reference levels from the core of each edge's direction, and the distance where each measure settles
    summary = []
    for (name, kind), g in df.groupby(["test", "kind"]):
        core = g[(g.d_mm >= 3) & (g.far_mm >= 3)]
        amp_ref = float(np.sqrt(np.nanmean(core.tine_amp ** 2)))
        ezt_ref = float(np.nanmedian(core.abs_ezt))
        near = g[g.d_mm <= a.dmax].copy()
        near["bin"] = (near.d_mm / 0.2).round() * 0.2
        prof = near.groupby("bin").agg(amp=("tine_amp", lambda v: np.sqrt(np.nanmean(v ** 2))),
                                       ezt=("abs_ezt", "median"))
        def settle(col, ref):
            """First distance where the measure is at the core level and stays there for the next 0.4 mm."""
            thr = 2 * ref if col == "amp" else 1.5 * ref
            v = prof[col].dropna()
            for i in range(len(v)):
                if (v.iloc[i:i + 3] <= thr).all():
                    return float(v.index[i])
            return np.nan
        summary.append(dict(test=name, edges=kind, tine_amp_core_rms=amp_ref, tine_amp_at_edge=float(prof.amp.dropna().iloc[0]),
                            tine_settles_mm=settle("amp", amp_ref), ezt_core=ezt_ref,
                            ezt_at_edge=float(prof.ezt.dropna().iloc[0]), ezt_settles_mm=settle("ezt", ezt_ref),
                            window_floor_mm=info[name]["floor"]))
    sm = pd.DataFrame(summary)
    sm["current_edge_mm"] = [edge["z"] if k == "axial" else edge["theta"] for k in sm.edges]
    sm["within_current"] = sm.tine_settles_mm <= sm.current_edge_mm + 1e-9

    # figure: rows = tests; columns = tine-pitch amplitude (axial edges, circumferential edges), |E_zt|
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.labelcolor": MUTED,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig, ax = plt.subplots(len(a.tests), 3, figsize=(16, 4.6 * len(a.tests)), constrained_layout=True, squeeze=False)
    for r, name in enumerate(a.tests):
        g = df[(df.test == name) & (df.d_mm <= a.dmax)]
        for c, (col, kind, ylabel) in enumerate([("tine_amp", "axial", "tine-pitch amplitude of λ_z (fraction)"),
                                                 ("tine_amp", "circumferential", "tine-pitch amplitude of λ_θ (fraction)"),
                                                 ("abs_ezt", None, "mean |E_zθ| along the line")]):
            a_ = ax[r, c]
            gg = g if kind is None else g[g.kind == kind]
            for fam, color in FAM_COLOR.items():
                h = gg[gg.family == fam].copy()
                h["bin"] = (h.d_mm / 0.2).round() * 0.2
                val = h.groupby("bin")[col].apply(lambda v: np.sqrt(np.nanmean(v ** 2)) if col == "tine_amp"
                                                  else np.nanmedian(v))
                a_.plot(val.index, val.values, "-o", ms=3, lw=1.5, color=color, label=fam)
            ks = ["axial", "circumferential"] if kind is None else [kind]
            for k_ in ks:
                row = sm[(sm.test == name) & (sm.edges == k_)].iloc[0]
                ref = row.tine_amp_core_rms if col == "tine_amp" else row.ezt_core
                a_.axhline(2 * ref if col == "tine_amp" else 1.5 * ref, color=MUTED, ls=":", lw=1)
            a_.axvline(info[name]["floor"], color="#e34948", ls="--", lw=1)
            if kind is not None:
                cur = edge["z"] if kind == "axial" else edge["theta"]
                a_.axvline(cur, color="#2b8a3e", ls="--", lw=1.3)
                a_.text(cur + 0.05, a_.get_ylim()[1] * 0.85, f"current edge_mm {cur:g}", color="#2b8a3e",
                        fontsize=8.5, va="top")
            a_.text(info[name]["floor"] + 0.05, a_.get_ylim()[1] * 0.95, "window floor", color="#e34948",
                    fontsize=8.5, va="top")
            a_.set_xlabel("distance from tine-tip line (mm)"); a_.set_ylabel(ylabel); a_.grid(alpha=.3)
            title = {0: "axial edges (left, right rakes)", 1: "circumferential edges (top, bottom rakes)",
                     2: "all edges, shear"}[c]
            a_.set_title(f"{info[name]['sample'].split('_')[0]} {info[name]['test']} · {title}", fontsize=11, loc="left")
            if r == 0 and c == 0:
                a_.legend(fontsize=9)
    fig.suptitle(f"Tine effects vs distance from the tine-tip lines · {S['subset_mm']:.1f} mm subset, "
                 f"{S['window_mm']:g} mm strain window · dotted: core reference (2 × RMS amplitude, 1.5 × median shear)"
                 f"\ngreen: current interior edge (edge_mm z {edge['z']:g}, θ {edge['theta']:g} mm, analysis_settings.json)",
                 fontsize=12.5)
    tag = "_".join(n.replace("/", "-") for n in a.tests)
    out = os.path.join(testcfg.ROOT, "comparisons", f"edge_effects_{tag}")
    fig.savefig(out + ".png", dpi=110)
    df.to_csv(out + "_profiles.csv", index=False, float_format="%.5f")
    sm.to_csv(out + ".csv", index=False, float_format="%.4f")
    with pd.option_context("display.width", 200):
        print(sm.round(4).to_string(index=False))
    for n, v in info.items():
        print(f"{n}: subset {v['subset']:.3f} mm, step {v['step']:.3f} mm, window floor {v['floor']:.2f} mm")
    print("\nTine-pitch modulation vs the current interior edge (edge_mm in analysis_settings.json):")
    for kind, key, where in (("axial", "z", "left/right"), ("circumferential", "theta", "top/bottom")):
        g = sm[sm.edges == kind]
        over = ", ".join(f"{t} ({d:g} mm)" for t, d in zip(g.test, g.tine_settles_mm) if d > edge[key] + 1e-9)
        print(f"  {where} edges: settles at {g.tine_settles_mm.min():g}-{g.tine_settles_mm.max():g} mm; "
              f"edge_mm.{key} = {edge[key]:g} mm; beyond it: {over or 'none'}")
    print("wrote", out + ".png / .csv / _profiles.csv")


if __name__ == "__main__":
    main()
