"""
Parameter estimation for the bulk four-fiber family model (constitutive.py;
Jadidi et al. 2019, Eq. 2) from the CellScale 4-point tracking (Coords1.csv) over
several protocols. Homogeneous (single-material) fit.

Kinematics: F at the center of the bilinear element through the 4 CellScale points,
relative to ONE reference configuration (start of the first protocol listed), then
lambda_z = sqrt(C_zz), lambda_theta = sqrt(C_thetatheta).
Stress: t_zz = (F_z - F_z0) lambda_z / (w_theta H),  t_tt = (F_t - F_t0) lambda_t / (w_z H)
  F_z0, F_t0: force at the start of each protocol (load-cell offset at rest).
  w_theta, w_z: tine-span widths (default) or full specimen widths (--width specimen).
Data: loading frames of each protocol, from the Stretch start to the rake peak
(compare_stretch_inputs.loading_end).
Fit: multi-start least squares for the global minimum, then nonparametric bootstrap
(resample frames, refit) as in Jadidi et al. 2019. Also a gamma profile and the
linearized parameter correlations.

    .venv/bin/python constitutive_fit/fit_4pt.py --test N1739/test1 --sets 2-10
Thickness and widths default to the test's test_config.json.
Outputs: <test>/constitutive_fit/
"""
import argparse
import glob
import json
import os
import sys

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
from compare_stretch_inputs import loading_end  # noqa: E402

TEST_ROOT = testcfg.get().root
WIDTHS = {k: tuple(v) for k, v in testcfg.get().config.get("widths_mm", {}).items()}   # (w_theta for F_z, w_z for F_theta), mm
THICKNESS = testcfg.get().config.get("thickness_mm", {})
LO = np.array([0, 0, 1e-3, 0, 1e-3, 0, 1e-3, 0])
HI = np.array([300, 3000, 80, 3000, 80, 3000, 80, 90])
STRESS_PARAMS = [0, 1, 3, 5]
# Jadidi et al. 2019 Table 1 (flattening-adjusted, age-group averages), for reference
TABLE1 = {"16.3": [3.41, 40.78, 0.04, 0.13, 4.18, 10.46, 0.97, 61.66],
          "25.0": [12.55, 29.80, 0.32, 6.62, 3.37, 3.66, 3.68, 56.55],
          "36.0": [9.38, 34.43, 0.75, 6.88, 4.03, 4.00, 4.60, 55.11],
          "46.5": [9.63, 39.41, 1.82, 10.93, 6.31, 5.40, 7.21, 50.58],
          "56.1": [5.05, 41.09, 2.91, 14.75, 8.24, 7.72, 10.16, 50.23],
          "65.0": [3.96, 45.24, 6.30, 15.91, 10.51, 8.19, 18.84, 48.35],
          "75.5": [3.05, 38.53, 10.34, 26.55, 10.42, 7.56, 24.82, 50.96]}


def parse_sets(s):
    out = []
    for part in s.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return out


def load_data(sets, width):
    C = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "Coords*.csv"))[0], skiprows=1)
    tC = C["Time(S)"].round(1).to_numpy()
    P = C[["X1", "Y1", "X3", "Y3", "X4", "Y4", "X2", "Y2"]].to_numpy(float).reshape(-1, 4, 2)  # TL,TR,BR,BL
    Dd = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "*Data.csv"))[0], encoding="latin1")
    Dd = Dd.drop_duplicates("Time_S").set_index("Time_S")
    ph = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "Logs", "*PhaseTiming.csv"))[0])
    w_t, w_z = WIDTHS[width]
    rows, ref = [], None
    for s in sets:
        st = ph[(ph.Set == s) & (ph.Phase == "Stretch")].iloc[0]
        rc = ph[(ph.Set == s) & (ph.Phase == "Recover")].iloc[0]
        t0 = float(st["StartTime(s)"])
        t1 = loading_end(Dd, t0, float(rc["StartTime(s)"]))          # rake peak, not the logged boundary
        idx = np.nonzero((tC >= t0 - 1e-6) & (tC <= t1 + 1e-6))[0]
        if ref is None:
            ref = P[idx[0]]
        d = Dd.reindex(tC[idx])
        fz0, ft0 = d.X_Force_mN.iloc[0], d.Y_Force_mN.iloc[0]
        for j, i in enumerate(idx):
            E = D.green_lagrange(D.bilinear_F_center(ref, P[i]))
            lz, lt = np.sqrt(1 + 2 * E[0, 0]), np.sqrt(1 + 2 * E[1, 1])
            fz, ft = d.X_Force_mN.iloc[j] - fz0, d.Y_Force_mN.iloc[j] - ft0
            rows.append(dict(set=s, name=st.SetName.replace(" save cycle", "").replace(" save cyle", ""),
                             t=tC[i], lz=lz, lt=lt, Ezt=E[0, 1], fz=fz, ft=ft, fz0=fz0, ft0=ft0,
                             tzz_per_H=fz * lz / w_t, ttt_per_H=ft * lt / w_z))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="2-10", help="protocol Set numbers from the PhaseTiming log (2-10 = Def11-19)")
    ap.add_argument("--H", type=float, default=THICKNESS.get("H"), help="mm (default: the test's thickness_mm.H)")
    ap.add_argument("--H-sd", type=float, default=THICKNESS.get("H_sd"))
    ap.add_argument("--width", choices=tuple(WIDTHS), default="tine_span")
    ap.add_argument("--starts", type=int, default=400)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--out", default=testcfg.get().fit_dir)
    a = ap.parse_args()
    if a.H is None or a.width not in WIDTHS:
        sys.exit("set thickness_mm.H and widths_mm in the test's test_config.json (or pass --H)")
    os.makedirs(a.out, exist_ok=True)
    np.seterr(over="ignore", invalid="ignore")

    df = load_data(parse_sets(a.sets), a.width)
    df["tzz"], df["ttt"] = df.tzz_per_H / a.H, df.ttt_per_H / a.H
    lz, lt, tzz, ttt = (df[c].to_numpy() for c in ("lz", "lt", "tzz", "ttt"))
    allidx = np.arange(len(df))
    y = np.r_[tzz, ttt]

    def resid(p, idx):
        mz, mt = M.cauchy_biaxial(lz[idx], lt[idx], p)
        return np.r_[mz - tzz[idx], mt - ttt[idx]]

    rng = np.random.default_rng(0)
    sols = []
    for s in range(a.starts):
        p0 = np.r_[10 ** rng.uniform(-0.5, 1.5), 10 ** rng.uniform(0.3, 2.3), 10 ** rng.uniform(-1.5, 1.3),
                   10 ** rng.uniform(-0.5, 2.3), 10 ** rng.uniform(-1, 1.5), 10 ** rng.uniform(-1, 2),
                   10 ** rng.uniform(-1, 1.5), rng.uniform(5, 85)]
        r = least_squares(resid, np.clip(p0, LO + 1e-6, HI - 1e-6), args=(allidx,), bounds=(LO, HI),
                          x_scale="jac", max_nfev=5000)
        if np.isfinite(r.cost):
            sols.append((2 * r.cost, r.x))
    sols.sort(key=lambda s: s[0])
    best_rss, best = sols[0]
    r2 = 1 - best_rss / np.sum((y - y.mean()) ** 2)
    near = np.array([x for rss, x in sols if rss <= 1.05 * best_rss])

    boots = []
    for b in range(a.boot):
        idx = rng.choice(allidx, size=allidx.size, replace=True)
        boots.append(least_squares(resid, best, args=(idx,), bounds=(LO, HI), x_scale="jac", max_nfev=5000).x)
    boots = np.array(boots)

    # gamma profile
    prof = []
    for g in np.arange(0, 91, 5):
        bestr, bx = np.inf, None
        for s in range(15):
            p0 = best.copy() if s == 0 else sols[min(s * 7, len(sols) - 1)][1].copy()
            f = lambda q: resid(np.r_[q, g], allidx)
            r = least_squares(f, np.clip(p0[:7], LO[:7] + 1e-6, HI[:7] - 1e-6), bounds=(LO[:7], HI[:7]),
                              x_scale="jac", max_nfev=4000)
            if 2 * r.cost < bestr:
                bestr, bx = 2 * r.cost, r.x
        prof.append((g, bestr, bx))

    # linearized correlations (free parameters only)
    r = least_squares(resid, best, args=(allidx,), bounds=(LO, HI), x_scale="jac")
    free = [i for i in range(8) if LO[i] + 1e-6 < best[i] < HI[i] - 1e-6]
    Jf = r.jac[:, free]
    cov = np.linalg.pinv(Jf.T @ Jf) * best_rss / (y.size - len(free))
    sd = np.sqrt(np.diag(cov)); corr = cov / np.outer(sd, sd)

    # ---- report
    n_per = df.groupby("set").size()
    print(f"sets {a.sets}: {len(df)} loading frames ({', '.join(f'{k}:{v}' for k, v in n_per.items())}); "
          f"reference = start of set {df.set.iloc[0]} (t = {df.t.iloc[0]:.0f} s); widths '{a.width}'; H = {a.H} mm")
    print(f"force offsets at protocol start (mN): z {df.groupby('set').fz0.first().round(1).tolist()}, "
          f"theta {df.groupby('set').ft0.first().round(1).tolist()}")
    print(f"max |E_ztheta| in data {np.abs(df.Ezt).max():.4f}")
    print(f"\nbest fit: R^2 = {r2:.4f}, RMS {np.sqrt(best_rss / y.size):.2f} kPa "
          f"(peak stresses {tzz.max():.1f} / {ttt.max():.1f} kPa); {len(near)}/{len(sols)} starts within 5% of best")
    print(f"\n{'param':10s} {'best':>9s} {'boot median':>11s} {'boot 2.5%':>10s} {'97.5%':>9s} {'near-opt range':>20s}  thickness ±")
    res = {}
    for i, (n, u) in enumerate(zip(M.PARAM_NAMES, M.PARAM_UNITS)):
        lo_b, md, hi_b = np.percentile(boots[:, i], [2.5, 50, 97.5])
        hsd = best[i] * a.H_sd / a.H if i in STRESS_PARAMS else 0.0
        print(f"{n:10s} {best[i]:9.3f} {md:11.3f} {lo_b:10.3f} {hi_b:9.3f} {near[:, i].min():9.3f}-{near[:, i].max():<9.3f}  "
              f"{hsd:6.3f} {u}")
        res[n] = dict(value=float(best[i]), unit=u, bootstrap_median=float(md), bootstrap_95=[float(lo_b), float(hi_b)],
                      near_opt_range=[float(near[:, i].min()), float(near[:, i].max())], thickness_sd=float(hsd))
    print("\nper-protocol RMS residual (kPa) [z, theta]:")
    per = {}
    for s, g in df.groupby("set"):
        ii = g.index.to_numpy(); rr = resid(best, ii).reshape(2, -1)
        rz, rt = np.sqrt((rr ** 2).mean(1))
        per[int(s)] = dict(name=g['name'].iloc[0], rms_z=float(rz), rms_theta=float(rt))
        print(f"  set {s:2d} {g['name'].iloc[0]:24s} {rz:6.2f} {rt:6.2f}")
    print("\ngamma profile (RMS kPa):", " ".join(f"{g:.0f}°:{np.sqrt(rss / y.size):.2f}" for g, rss, _ in prof))
    names = [M.PARAM_NAMES[i] for i in free]
    strong = [(names[i], names[j], corr[i, j]) for i in range(len(free)) for j in range(i + 1, len(free))
              if abs(corr[i, j]) > 0.9]
    print("strong parameter correlations (|r| > 0.9):", ", ".join(f"{p}-{q} {c:+.2f}" for p, q, c in strong) or "none")

    json.dump(dict(model="Jadidi 2019 bulk four-fiber family (Eq. 2)", data="CellScale 4-point (Coords1.csv)",
                   sets=a.sets, reference_time_s=float(df.t.iloc[0]), H_mm=a.H, H_sd_mm=a.H_sd, width_choice=a.width,
                   widths_mm=dict(w_theta=WIDTHS[a.width][0], w_z=WIDTHS[a.width][1]), n_frames=int(len(df)),
                   R2=float(r2), rms_kPa=float(np.sqrt(best_rss / y.size)), params=res, per_protocol=per,
                   gamma_profile=[dict(gamma=float(g), rms_kPa=float(np.sqrt(rss / y.size))) for g, rss, _ in prof],
                   strong_correlations=[dict(a=p, b=q, r=float(c)) for p, q, c in strong]),
              open(os.path.join(a.out, "fit_4pt_params.json"), "w"), indent=2)
    df.to_csv(os.path.join(a.out, "fit_4pt_data.csv"), index=False)

    # ---- plots
    sets = list(df.set.unique())
    fig, ax = plt.subplots(figsize=(12, 8.5), constrained_layout=True)
    cmap = plt.get_cmap("tab10")
    from matplotlib.lines import Line2D
    handles = []
    for c, s in enumerate(sets):
        g = df[df.set == s]; ii = g.index.to_numpy()
        # model drawn along a lightly smoothed stretch path: the CellScale 4-point tracking has
        # occasional backward jumps (~0.5% stretch) that would otherwise break the curves
        from scipy.signal import savgol_filter
        win = min(11, len(ii) - (1 - len(ii) % 2))
        sz, sl = savgol_filter(lz[ii], win, 2), savgol_filter(lt[ii], win, 2)
        mz, mt = M.cauchy_biaxial(sz, sl, best)
        col = cmap(c % 10)
        ax.plot(g["lz"], g["tzz"], "o", ms=3, color=col, alpha=.35, mec="none")
        ax.plot(g["lt"], g["ttt"], "^", ms=3, color=col, alpha=.35, mec="none")
        ax.plot(sz, mz, "-", color=col, lw=1.8)
        ax.plot(sl, mt, "--", color=col, lw=1.8)
        handles.append(Line2D([], [], color=col, lw=2,
                              label=f"set {s}: {g['name'].iloc[0]}  (RMS {per[int(s)]['rms_z']:.1f} / {per[int(s)]['rms_theta']:.1f} kPa)"))
    handles += [Line2D([], [], color="k", marker="o", ls="-", ms=4, label=r"axial: data $t_{zz}$ vs $\lambda_z$ (o), model (solid)"),
                Line2D([], [], color="k", marker="^", ls="--", ms=4,
                       label=r"circumferential: data $t_{\theta\theta}$ vs $\lambda_\theta$ ($\triangle$), model (dashed)")]
    ax.legend(handles=handles, fontsize=9, loc="upper left")
    ax.set_xlabel(r"stretch ($\lambda_z$ or $\lambda_\theta$)"); ax.set_ylabel("Cauchy stress (kPa)"); ax.grid(alpha=.3)
    ax.set_title(f"Bulk four-fiber model, single homogeneous fit to CellScale 4-point data, sets {a.sets}\n"
                 f"R² = {r2:.3f}, RMS {np.sqrt(best_rss / y.size):.2f} kPa; H = {a.H} mm; widths: {a.width}; "
                 f"all 8 parameters optimized (γ = {best[7]:.1f}°)")
    fig.savefig(os.path.join(a.out, "fit_4pt_stress_stretch.png"), dpi=90)

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6), constrained_layout=True)
    g_ = [p[0] for p in prof]; rms_ = [np.sqrt(p[1] / y.size) for p in prof]
    ax[0].plot(g_, rms_, "o-"); ax[0].axvline(best[7], color="k", lw=.6)
    ax[0].set_xlabel(r"$\gamma$ (deg), fixed for this diagnostic only"); ax[0].set_ylabel("RMS residual (kPa)")
    ax[0].set_title(r"$\gamma$ profile: other 7 refitted at each fixed $\gamma$" "\n"
                    rf"line = fitted $\gamma$ = {best[7]:.1f}° (all 8 free)", fontsize=10); ax[0].grid(alpha=.3)
    at_bound = (best <= LO + 1e-6) | (best >= HI - 1e-6)
    lab = [n.replace("_", "\n") + ("\n(at bound)" if b else "") for n, b in zip(M.PARAM_NAMES, at_bound)]
    rel = boots / np.where(at_bound, np.nan, best)
    ax[1].boxplot([rel[:, i] if not at_bound[i] else [np.nan] for i in range(8)], labels=lab, showfliers=False)
    ax[1].axhline(1, color="k", lw=.6); ax[1].set_title(f"Bootstrap spread ({a.boot} resamples), value / best")
    ax[1].grid(alpha=.3)
    im = ax[2].imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    ax[2].set_xticks(range(len(names))); ax[2].set_xticklabels(names, rotation=60, fontsize=8)
    ax[2].set_yticks(range(len(names))); ax[2].set_yticklabels(names, fontsize=8)
    for i in range(len(names)):
        for j in range(len(names)):
            ax[2].text(j, i, f"{corr[i, j]:+.2f}", ha="center", va="center", fontsize=7)
    ax[2].set_title("Linearized parameter correlation"); fig.colorbar(im, ax=ax[2], shrink=.8)
    fig.savefig(os.path.join(a.out, "fit_4pt_identifiability.png"), dpi=90)
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
