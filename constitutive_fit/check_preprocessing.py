"""
Sensitivity of the frozen bulk baseline (baseline.py) to force taring and light
per-protocol smoothing before fitting.

Variants (applied per protocol, loading branch, on the frozen data set):
    baseline        tare = force at protocol start; no smoothing (as frozen)
    no_tare         raw load-cell forces
    tare_first5     tare = mean force over the first 5 frames
    sg11_force      Savitzky-Golay (11 frames, 2nd order) on forces
    sg11_both       same filter on forces and stretches
    sg21_both       21-frame filter on forces and stretches
    no_tare_sg11    raw forces + 11-frame filter on both
Each variant is refit with the elastin-active multi-start (seeded from the baseline) and
compared with the baseline value and its protocol-bootstrap 95% interval.

    .venv/bin/python constitutive_fit/check_preprocessing.py --test N1739/test1 [--inputs dic_affine cs4] [--tracks tine_region|standard]
Outputs: preprocessing_sensitivity.csv / .png in the chosen baseline's folder
(constitutive_fit/baseline/ by default, the primary full-domain baseline)
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.signal import savgol_filter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import baseline as B  # noqa: E402
import compare_stretch_inputs as CSI  # noqa: E402
import constitutive as M  # noqa: E402

TEST_ROOT = testcfg.get().root
VARIANTS = ["baseline", "no_tare", "tare_first5", "sg11_force", "sg11_both", "sg21_both", "no_tare_sg11"]
LABELS = {"baseline": "baseline (tare at start)", "no_tare": "no tare", "tare_first5": "tare: mean of first 5 frames",
          "sg11_force": "SG-11 on force", "sg11_both": "SG-11 on force + stretch",
          "sg21_both": "SG-21 on force + stretch", "no_tare_sg11": "no tare + SG-11 on both"}


def prepare(df, raw, inp, variant):
    """Per-protocol force/stretch series for a variant -> (lz, lt, fz, ft)."""
    out = {k: [] for k in ("lz", "lt", "fz", "ft")}
    for s, g in df.groupby("set", sort=False):
        lz, lt = g[f"{inp}_z"].to_numpy(), g[f"{inp}_t"].to_numpy()
        rz, rt = raw.loc[g.index, "fz_raw"].to_numpy(), raw.loc[g.index, "ft_raw"].to_numpy()
        if variant.startswith("no_tare"):
            fz, ft = rz, rt
        elif variant == "tare_first5":
            fz, ft = rz - rz[:5].mean(), rt - rt[:5].mean()
        else:
            fz, ft = g["fz"].to_numpy(), g["ft"].to_numpy()       # tare at start (frozen)
        win = 21 if "sg21" in variant else 11
        if "sg" in variant:
            fz, ft = savgol_filter(fz, win, 2), savgol_filter(ft, win, 2)
            if "both" in variant or variant == "no_tare_sg11":
                lz, lt = savgol_filter(lz, win, 2), savgol_filter(lt, win, 2)
        for k, v in zip(("lz", "lt", "fz", "ft"), (lz, lt, fz, ft)):
            out[k].append(v)
    return [np.concatenate(out[k]) for k in ("lz", "lt", "fz", "ft")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", default=["dic_affine", "cs4"])
    ap.add_argument("--starts", type=int, default=40)
    ap.add_argument("--tracks", choices=tuple(B.TRACKS), default="tine_region",
                    help="which frozen baseline to perturb (default: full domain)")
    a = ap.parse_args()
    np.seterr(over="ignore", invalid="ignore")
    base, df = B.load(a.tracks)
    out = B.TRACKS[a.tracks][1]
    st = base["settings"]
    H, w_t, w_z = st["H_mm"], st["w_theta_mm"], st["w_z_mm"]
    Dd = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "*Data.csv"))[0], encoding="latin1")
    Dd = Dd.drop_duplicates("Time_S").set_index("Time_S")
    raw = pd.DataFrame({"fz_raw": Dd.X_Force_mN.reindex(df.t.round(1)).to_numpy(),
                        "ft_raw": Dd.Y_Force_mN.reindex(df.t.round(1)).to_numpy()}, index=df.index)
    rng = np.random.default_rng(0)
    rows = []
    for inp in a.inputs:
        ref = base["inputs"][inp]
        seed = np.array([ref["params"][n] for n in M.PARAM_NAMES])
        for v in VARIANTS:
            lz, lt, fz, ft = prepare(df, raw, inp, v)
            tzz, ttt = fz * lz / (w_t * H), ft * lt / (w_z * H)
            p, r2, rms = CSI.fit(lz, lt, tzz, ttt, a.starts, rng, seeds=[seed])
            row = dict(input=inp, variant=v, label=LABELS[v], R2=r2, rms_kPa=rms)
            outside = []
            for i, n in enumerate(M.PARAM_NAMES):
                lo, hi = ref["bootstrap_95"][n]
                row[n] = p[i]
                row[n + "_pct_change"] = 100 * (p[i] / seed[i] - 1) if seed[i] else np.nan
                if not lo <= p[i] <= hi:
                    outside.append(n)
            row["outside_baseline_95"] = ", ".join(outside)
            rows.append(row)
            print(f"{inp:10s} {LABELS[v]:30s} R2 {r2:.4f} RMS {rms:.3f} | "
                  + " ".join(f"{n} {p[i]:.3g} ({row[n + '_pct_change']:+.0f}%)" for i, n in enumerate(M.PARAM_NAMES))
                  + (f" | outside 95%: {row['outside_baseline_95']}" if outside else ""), flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(out, "preprocessing_sensitivity.csv"), index=False)

    fig, ax = plt.subplots(len(a.inputs), 1, figsize=(15, 4.2 * len(a.inputs)), constrained_layout=True, squeeze=False)
    x = np.arange(len(M.PARAM_NAMES))
    for k, inp in enumerate(a.inputs):
        a_ = ax[k, 0]; ref = base["inputs"][inp]
        seed = np.array([ref["params"][n] for n in M.PARAM_NAMES])
        lo = np.array([100 * (ref["bootstrap_95"][n][0] / seed[i] - 1) for i, n in enumerate(M.PARAM_NAMES)])
        hi = np.array([100 * (ref["bootstrap_95"][n][1] / seed[i] - 1) for i, n in enumerate(M.PARAM_NAMES)])
        a_.bar(x, hi - lo, bottom=lo, width=.8, color="0.88", label="baseline protocol-bootstrap 95%")
        for j, v in enumerate(VARIANTS[1:]):
            r = T[(T.input == inp) & (T.variant == v)].iloc[0]
            a_.plot(x + (j - 2.5) * 0.11, [r[n + "_pct_change"] for n in M.PARAM_NAMES], "o", ms=6,
                    color=plt.get_cmap("tab10")(j), label=f"{LABELS[v]} (RMS {r.rms_kPa:.2f} kPa)")
        a_.axhline(0, color="k", lw=.8)
        a_.set_xticks(x); a_.set_xticklabels(M.PARAM_NAMES); a_.set_ylabel("change vs baseline (%)")
        a_.set_ylim(-110, 160); a_.grid(alpha=.3); a_.legend(fontsize=8, ncol=2, loc="upper left")
        a_.set_title(f"{ref['label']}, {a.tracks} baseline: preprocessing sensitivity "
                     f"(baseline RMS {ref['rms_kPa']:.2f} kPa)")
    fig.savefig(os.path.join(out, "preprocessing_sensitivity.png"), dpi=90)
    print("wrote", os.path.join(out, "preprocessing_sensitivity.csv"), "and .png")


if __name__ == "__main__":
    main()
