"""
Stage 0 of the local-property stages: what the test's own bulk data determine, for the fixed eight-parameter model.

The model form never changes (no parameter is dropped; a parameter may go to its lower bound, 0). This script only
reports, per parameter, how far it can move before the fit gets worse: a profile of the bulk misfit. For each
parameter and each value on a grid spanning its bounds, that parameter is held fixed and the other seven are refitted
(continuation outward from the best fit, plus the best fit itself as a seed). The fit is the bulk fit of
compare_stretch_inputs.py (DIC affine F, tared forces, tine-span widths, the test's H, Def11-Def19 loading branches,
elastin-active family: C1_el >= EL_MIN).

Reported per parameter: the range where the RMS stays within +1% and +5% of the best, whether that range reaches a
bound (0 or the upper bound), the parameter that moves most along the profile at the +5% ends (what it trades off
with), and the protocol-bootstrap 95% interval from params_by_stretch_input.csv for comparison. A profile point below
the best RMS means the saved fit was not the global minimum and is flagged.

    .venv/bin/python constitutive_fit/stage0_identifiability.py --test N1739/test1
Reads <test>/constitutive_fit/stretch_inputs_data.csv and params_by_stretch_input.csv (run compare_stretch_inputs.py
first). Outputs: <test>/constitutive_fit/local/stage0_identifiability.csv / .json / .png
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
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import constitutive as M  # noqa: E402

INPUT = "dic_affine"
LO = np.array([0, 0.1, 1e-3, 0, 1e-3, 0, 1e-3, 0])        # compare_stretch_inputs.LO with C1_el >= EL_MIN
HI = np.array([300, 3000, 80, 3000, 80, 3000, 80, 90])
LEVELS = (0.01, 0.05)     # RMS increase thresholds
N_GRID = 31


def grid_for(i, best):
    """Coarse grid over the bounds plus a fine one around the best value (x0.33 to x3, or +-20 deg for gamma)."""
    if i == 7:
        g = np.r_[np.linspace(0, 90, N_GRID), np.linspace(best - 20, best + 20, 41)]
    else:
        lo = LO[i] if LO[i] > 0 else 1e-2
        g = np.geomspace(lo, HI[i], N_GRID - (LO[i] == 0))
        if best > lo:
            g = np.r_[g, best * np.geomspace(1 / 3, 3, 41)]
        if LO[i] == 0:
            g = np.r_[0.0, g]
    return np.unique(np.clip(g, LO[i], HI[i]))


def make_resid(lz, lt, tzz, ttt):
    def resid(p):
        mz, mt = M.cauchy_biaxial(lz, lt, p)
        return np.r_[mz - tzz, mt - ttt]
    return resid


def fit_fixed(resid, i, val, seeds):
    """Best fit with parameter i held at val, from each seed. Returns (rms, params)."""
    free = np.arange(8) != i
    best = (np.inf, None)
    for s in seeds:
        p = np.array(s, float)
        p[i] = val

        def r_(q):
            p[free] = q
            return resid(p)
        x0 = np.clip(p[free], LO[free] + 1e-9, HI[free] - 1e-9)
        r = least_squares(r_, x0, bounds=(LO[free], HI[free]), x_scale="jac", max_nfev=3000)
        if 2 * r.cost < best[0]:
            q = p.copy()
            q[free] = r.x
            best = (2 * r.cost, q)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", default="tine_span")
    a = ap.parse_args()
    T = testcfg.get()
    np.seterr(over="ignore", invalid="ignore")
    df = pd.read_csv(os.path.join(T.fit_dir, "stretch_inputs_data.csv"))
    P = pd.read_csv(os.path.join(T.fit_dir, "params_by_stretch_input.csv")).set_index("input")
    H = T.config["thickness_mm"]["H"]
    w_t, w_z = T.config["widths_mm"][a.width]
    lz, lt = df[f"{INPUT}_z"].to_numpy(), df[f"{INPUT}_t"].to_numpy()
    tzz, ttt = df.fz.to_numpy() * lz / (w_t * H), df.ft.to_numpy() * lt / (w_z * H)
    resid = make_resid(lz, lt, tzz, ttt)
    n = 2 * len(lz)
    rms = lambda sse: np.sqrt(sse / n)  # noqa: E731

    # best fit: the saved one, polished, and every input's saved fit as a seed
    seeds = [np.clip(P.loc[k, list(M.PARAM_NAMES)].to_numpy(float), LO, HI) for k in P.index]
    p_saved = np.clip(P.loc[INPUT, list(M.PARAM_NAMES)].to_numpy(float), LO, HI)
    best = (np.inf, None)
    for s in [p_saved] + seeds:
        r = least_squares(resid, np.clip(s, LO + 1e-9, HI - 1e-9), bounds=(LO, HI), x_scale="jac", max_nfev=5000)
        if 2 * r.cost < best[0]:
            best = (2 * r.cost, r.x)
    sse0, p0 = best
    print(f"{T.name}: best RMS {rms(sse0):.4f} kPa (saved {P.loc[INPUT, 'rms_kPa']:.4f}); "
          + " ".join(f"{k}={v:.3g}" for k, v in zip(M.PARAM_NAMES, p0)), flush=True)

    profiles, rows, improved = {}, [], []
    for i, name in enumerate(M.PARAM_NAMES):
        g = np.unique(np.r_[grid_for(i, p0[i]), p0[i]])
        k0 = int(np.argmin(np.abs(g - p0[i])))
        R = np.full(len(g), np.nan)
        Q = np.full((len(g), 8), np.nan)
        R[k0], Q[k0] = sse0, p0
        for direction in (1, -1):
            prev = p0
            k = k0 + direction
            while 0 <= k < len(g):
                sse, q = fit_fixed(resid, i, g[k], [prev, p0])
                R[k], Q[k] = sse, q
                prev = q
                k += direction
        if np.nanmin(R) < sse0 * (1 - 1e-4):
            improved.append(dict(param=name, value=float(g[np.nanargmin(R)]), rms=float(rms(np.nanmin(R)))))
        rel = np.sqrt(R / sse0) - 1                       # relative RMS increase
        profiles[name] = dict(grid=g, rel=rel, params=Q)
        row = dict(param=name, unit=M.PARAM_UNITS[i], best=p0[i], lower_bound=LO[i], upper_bound=HI[i],
                   boot_lo95=P.loc[INPUT, f"{name}_lo95"], boot_hi95=P.loc[INPUT, f"{name}_hi95"])
        for lev in LEVELS:
            inside = rel <= lev
            # contiguous range around the best
            lo_k = k0
            while lo_k - 1 >= 0 and inside[lo_k - 1]:
                lo_k -= 1
            hi_k = k0
            while hi_k + 1 < len(g) and inside[hi_k + 1]:
                hi_k += 1
            tag = f"{int(lev * 100)}pct"
            row[f"lo_{tag}"], row[f"hi_{tag}"] = g[lo_k], g[hi_k]
            row[f"at_lower_bound_{tag}"] = bool(lo_k == 0)
            row[f"at_upper_bound_{tag}"] = bool(hi_k == len(g) - 1)
            if lev == LEVELS[-1]:
                # trade-off partner: largest log change of another parameter across the +5% range
                ends = Q[[lo_k, hi_k]]
                with np.errstate(divide="ignore", invalid="ignore"):
                    ch = np.abs(np.log((np.abs(ends[1]) + 1e-3) / (np.abs(ends[0]) + 1e-3)))
                ch[i] = -1
                j = int(np.nanargmax(ch))
                row["trades_off_with"] = M.PARAM_NAMES[j]
                row["partner_range"] = f"{ends[0, j]:.3g} to {ends[1, j]:.3g}"
                lb, ub = row[f"at_lower_bound_{tag}"], row[f"at_upper_bound_{tag}"]
                row["status"] = ("unconstrained" if lb and ub else "upper limit only" if lb
                                 else "lower limit only" if ub else "determined")
        rows.append(row)
        print(f"  {name:7s} best {p0[i]:8.3g} | +5% range {row['lo_5pct']:.3g} to {row['hi_5pct']:.3g} "
              f"({row['status']}; trades off with {row['trades_off_with']} {row['partner_range']}) | "
              f"bootstrap {row['boot_lo95']:.3g} to {row['boot_hi95']:.3g}", flush=True)

    S = pd.DataFrame(rows)
    out = T.local_dir
    os.makedirs(out, exist_ok=True)
    S.to_csv(os.path.join(out, "stage0_identifiability.csv"), index=False)
    np.savez_compressed(os.path.join(out, "stage0_profiles.npz"),
                        **{f"{k}_{f}": v[f] for k, v in profiles.items() for f in ("grid", "rel", "params")})
    with open(os.path.join(out, "stage0_identifiability.json"), "w") as f:
        json.dump(dict(test=T.name, input=INPUT, H_mm=H, widths_mm=[w_t, w_z], frames=int(len(lz)),
                       best_rms_kPa=float(rms(sse0)), saved_rms_kPa=float(P.loc[INPUT, "rms_kPa"]),
                       best=dict(zip(M.PARAM_NAMES, map(float, p0))), levels=LEVELS,
                       profile_found_lower_rms=improved,
                       status={r["param"]: r["status"] for r in rows}), f, indent=2)

    fig, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)
    for i, (ax, name) in enumerate(zip(axes.flat, M.PARAM_NAMES)):
        pr, r = profiles[name], rows[i]
        x = pr["grid"]
        ax.axvspan(r["boot_lo95"], r["boot_hi95"], color="C1", alpha=0.2, label="bootstrap 95%")
        ax.plot(x, 100 * pr["rel"], "o-", ms=3, color="C0")
        for lev, ls in zip(LEVELS, (":", "--")):
            ax.axhline(100 * lev, color="k", ls=ls, lw=0.8, label=f"+{int(lev * 100)}% RMS")
        ax.axvline(p0[i], color="C3", lw=1, label="best")
        if i != 7:
            ax.set_xscale("symlog", linthresh=max(LO[i], 1e-2))
        top = max(30, 100 * np.nanmin(pr["rel"][pr["rel"] > LEVELS[-1]], initial=0.3) * 1.5)
        ax.set_ylim(-1, min(top, 60))
        ax.set(xlabel=f"{name} ({M.PARAM_UNITS[i]})", ylabel="RMS increase (%)",
               title=f"{name}: {r['status']}\ntrades off with {r['trades_off_with']}")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8)
    fig.suptitle(f"{T.name}: profile of the bulk misfit per parameter (others refitted), "
                 f"best RMS {rms(sse0):.3f} kPa, DIC affine F, Def11–Def19")
    fig.savefig(os.path.join(out, "stage0_identifiability.png"), dpi=110)
    if improved:
        print("  profile found a lower RMS than the saved fit:", improved)
    print("wrote", out)


if __name__ == "__main__":
    main()
