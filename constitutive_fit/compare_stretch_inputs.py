"""
Sensitivity of the bulk four-fiber parameters to the global stretch input.

Same model, stress definition, protocols (default Set 2-10 = Def11-19, loading branch),
force tare and reference configuration (start of the first protocol) for every input;
only the stretch measure changes. The measure's lambda is used both in the model and
in the current cross-section area: t_zz = F_z lambda_z / (w_theta H).

Stretch inputs (all relative to the common reference):
    cs4         CellScale 4-point tracking (Coords1.csv), bilinear F at element center
    dic4        DIC displacement at the same 4 material points, bilinear F
    dic_mean    mean of the local DIC stretches over the tracked region
    dic_affine  least-squares affine F of all valid DIC points
    dic_edges   change in tracked-region axial length / circumferential width (edge rows/cols)
    rake        reported rake separation (X/Y_Size_um)

Needs track_protocols.py output in Images/DefNN_*/strain_results/ (--tracks standard, the
test's roi_px.standard) or track_tine_region.py output (--tracks tine_region, grid out to the
tine tips; outputs go to stretch_inputs_tine_region/). --table-only writes the stretch table
and stops before any fitting.
Loading branch: from the logged Stretch start to the rake peak (loading_end), which is the logged
Recover start in every protocol except Def12, where the motion continues 0.5 s past it.
Fit: multi-start least squares restricted to the physiological solution family with the
elastin term active (C1_el > EL_MIN). The model also has a near-equal minimum with
elastin switched off (C1_el ~ 0, gamma ~ 34 deg), which is excluded as non-physiological.
Each input is also re-seeded from every input's best elastin-active solution so the
family is found reliably. Uncertainty: protocol-level bootstrap (resample whole
protocols with replacement, refit from the best solution); resamples that fall into the
elastin-off minimum are discarded and counted.

    ../.venv/bin/python compare_stretch_inputs.py --test N1739/test1
Outputs: Images/constitutive_fit/stretch_inputs/
"""
import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import constitutive as M  # noqa: E402
import dic2d as D  # noqa: E402

TEST_ROOT = testcfg.get().root
IMAGES = os.path.join(TEST_ROOT, "Images")
WIDTHS = {k: tuple(v) for k, v in testcfg.get().config.get("widths_mm", {}).items()}   # (w_theta, w_z), mm
THICKNESS = testcfg.get().config.get("thickness_mm", {})
LO = np.array([0, 0, 1e-3, 0, 1e-3, 0, 1e-3, 0])
HI = np.array([300, 3000, 80, 3000, 80, 3000, 80, 90])
INPUTS = ["cs4", "dic4", "dic_mean", "dic_affine", "dic_edges", "rake"]
LABELS = {"cs4": "CellScale 4-point", "dic4": "DIC 4-point", "dic_mean": "DIC region-mean λ",
          "dic_affine": "DIC affine F", "dic_edges": "DIC region edges", "rake": "rake (X/Y_Size)"}
EL_MIN = 0.1           # kPa; C1_el above this = elastin active (physiological solution family)


def lam_from_F(F):
    C = F.T @ F
    return np.sqrt(C[0, 0]), np.sqrt(C[1, 1])


TRACK_SUFFIX = {"standard": "", "tine_region": "_tine_region"}
LOADING_SLACK_S = D.LOADING_SLACK_S   # s past the logged Recover start searched for the rake peak


def loading_end(Dd, t0, t_recover, slack=LOADING_SLACK_S):
    """End of the loading branch: the time the rakes peak (largest summed rake stretch), searched from
    the Stretch start to `slack` s past the logged Recover start. The logged boundary can lag the motion:
    in Def12 the rakes and forces keep rising until 864.5 s, 0.5 s after the logged 864.0 s."""
    w = Dd.loc[t0 - 1e-6:t_recover + slack]
    rake = w.X_Size_um / w.X_Size_um.iloc[0] + w.Y_Size_um / w.Y_Size_um.iloc[0]
    return round(float(rake.idxmax()), 1)


def build_table(sets, tracks="tine_region"):
    C = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "Coords*.csv"))[0], skiprows=1)
    C["t"] = C["Time(S)"].round(1)
    C = C.drop_duplicates("t").set_index("t")
    P_all = C[["X1", "Y1", "X3", "Y3", "X4", "Y4", "X2", "Y2"]]          # TL, TR, BR, BL
    Dd = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "*Data.csv"))[0], encoding="latin1")
    Dd = Dd.drop_duplicates("Time_S").set_index("Time_S")
    ph = pd.read_csv(glob.glob(os.path.join(TEST_ROOT, "Logs", "*PhaseTiming.csv"))[0])
    dfolders = sorted(d for d in os.listdir(IMAGES) if re.match(r"Def\d+_", d))
    # set number -> deformation folder: Def folders are numbered by Stretch phase order
    stretch = ph[ph.Phase == "Stretch"].reset_index(drop=True)
    set_to_folder = {int(r.Set): dfolders[i] for i, r in stretch.iterrows() if i < len(dfolders)}
    # (sets with several cycles map to their last cycle; only single-cycle sets are used here)
    rows, ref = [], {}
    for s in sets:
        fd = set_to_folder[s]
        st = ph[(ph.Set == s) & (ph.Phase == "Stretch")].iloc[0]
        rc = ph[(ph.Set == s) & (ph.Phase == "Recover")].iloc[0]
        t0 = round(float(st["StartTime(s)"]), 1)
        t1 = loading_end(Dd, t0, round(float(rc["StartTime(s)"]), 1))
        R, _ = D.load_track(IMAGES, fd, tracks)      # standard: sliced from the full domain if not tracked
        X0, shape, times = R["X0"], tuple(R["shape"]), np.round(R["times"], 1)
        pos, valid = R["pos"], R["valid"]
        win = D.window_points(float(R["step"]))
        xs, ys = X0[:, 0].reshape(shape)[0], X0[:, 1].reshape(shape)[:, 0]
        if not ref:
            ref["t"] = times[0]
            ref["P"] = P_all.loc[times[0]].to_numpy(float).reshape(4, 2)
            ref["size"] = (Dd.X_Size_um[times[0]], Dd.Y_Size_um[times[0]])
            ref["X0"] = X0
            ref["W"] = X0[:, 0].max() - X0[:, 0].min()                   # region edges at reference
            ref["Hh"] = X0[:, 1].max() - X0[:, 1].min()
        k_load = np.nonzero((times >= t0 - 1e-6) & (times <= t1 + 1e-6))[0]
        Pg = pos.reshape((len(times),) + shape + (2,))
        fz0, ft0 = Dd.X_Force_mN[times[k_load[0]]], Dd.Y_Force_mN[times[k_load[0]]]
        for k in k_load:
            t = times[k]
            r = dict(set=s, folder=fd, name=st.SetName.replace(" save cycle", "").replace(" save cyle", ""), t=t,
                     fz=Dd.X_Force_mN[t] - fz0, ft=Dd.Y_Force_mN[t] - ft0)
            # cs4
            r["cs4_z"], r["cs4_t"] = lam_from_F(D.bilinear_F_center(ref["P"], P_all.loc[t].to_numpy(float).reshape(4, 2)))
            # dic4: DIC positions at the CellScale reference points
            Pd = np.stack([RegularGridInterpolator((ys, xs), pos[k, :, c].reshape(shape))(ref["P"][:, ::-1])
                           for c in (0, 1)], 1)
            r["dic4_z"], r["dic4_t"] = lam_from_F(D.bilinear_F_center(ref["P"], Pd))
            # dic_mean
            F, ok = D.deformation_gradient_field(X0, pos[k], valid[k], shape, win)
            Cm = np.einsum("...ki,...kj->...ij", F, F)
            r["dic_mean_z"] = np.nanmean(np.sqrt(Cm[..., 0, 0])[ok])
            r["dic_mean_t"] = np.nanmean(np.sqrt(Cm[..., 1, 1])[ok])
            # dic_affine
            v = valid[k]
            r["dic_affine_z"], r["dic_affine_t"] = lam_from_F(D.affine_fit(X0[v], pos[k][v]))
            # dic_edges
            r["dic_edges_z"] = (Pg[k, :, -1, 0].mean() - Pg[k, :, 0, 0].mean()) / ref["W"]
            r["dic_edges_t"] = (Pg[k, -1, :, 1].mean() - Pg[k, 0, :, 1].mean()) / ref["Hh"]
            # rake
            r["rake_z"] = Dd.X_Size_um[t] / ref["size"][0]
            r["rake_t"] = Dd.Y_Size_um[t] / ref["size"][1]
            rows.append(r)
        print(f"  {fd}: {len(k_load)} loading frames", flush=True)
    return pd.DataFrame(rows), ref


def fit(lz, lt, tzz, ttt, starts, rng, seeds=()):
    """Best elastin-active (C1_el > EL_MIN) least-squares solution from seeds + random starts.
    Returns (params, R2, RMS) or (None, nan, nan) if no elastin-active minimum was found."""
    def resid(p):
        mz, mt = M.cauchy_biaxial(lz, lt, p)
        return np.r_[mz - tzz, mt - ttt]
    cands = list(seeds)
    for _ in range(starts):
        cands.append(np.r_[10 ** rng.uniform(-0.5, 1.5), 10 ** rng.uniform(0.3, 2.3), 10 ** rng.uniform(-1.5, 1.3),
                           10 ** rng.uniform(-0.5, 2.3), 10 ** rng.uniform(-1, 1.5), 10 ** rng.uniform(-1, 2),
                           10 ** rng.uniform(-1, 1.5), rng.uniform(5, 85)])
    best = (np.inf, None)
    for p0 in cands:
        r = least_squares(resid, np.clip(p0, LO + 1e-6, HI - 1e-6), bounds=(LO, HI), x_scale="jac", max_nfev=5000)
        if r.x[1] > EL_MIN and 2 * r.cost < best[0]:
            best = (2 * r.cost, r.x)
    if best[1] is None:
        return None, np.nan, np.nan
    y = np.r_[tzz, ttt]
    return best[1], 1 - best[0] / np.sum((y - y.mean()) ** 2), np.sqrt(best[0] / y.size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="2-10")
    ap.add_argument("--H", type=float, default=THICKNESS.get("H"), help="mm (default: the test's thickness_mm.H)")
    ap.add_argument("--width", choices=tuple(WIDTHS), default="tine_span")
    ap.add_argument("--starts", type=int, default=150)
    ap.add_argument("--boot", type=int, default=200, help="protocol-level bootstrap resamples")
    ap.add_argument("--tracks", choices=tuple(TRACK_SUFFIX), default="tine_region",
                    help="DIC tracks: standard ROI, or tine_region (grid out to the tine tips)")
    ap.add_argument("--table-only", action="store_true", help="write the stretch table, no fitting")
    ap.add_argument("--out", default=None, help="default constitutive_fit/stretch_inputs[_tine_region]")
    a = ap.parse_args()
    if a.H is None or a.width not in WIDTHS:
        sys.exit("set thickness_mm.H and widths_mm in the test's test_config.json (or pass --H)")
    a.out = a.out or os.path.join(IMAGES, "constitutive_fit", "stretch_inputs" + TRACK_SUFFIX[a.tracks])
    np.seterr(over="ignore", invalid="ignore")
    os.makedirs(a.out, exist_ok=True)
    lo_, _, hi_ = a.sets.partition("-")
    sets = list(range(int(lo_), int(hi_ or lo_) + 1))

    cache = os.path.join(a.out, "stretch_inputs_data.csv")
    print("building stretch table (DIC, CellScale, rake)")
    df, ref = build_table(sets, a.tracks)
    df.to_csv(cache, index=False)
    if a.table_only:
        print("wrote", cache, "(table only, no fitting)")
        return
    w_t, w_z = WIDTHS[a.width]
    rng = np.random.default_rng(0)
    sets_arr = df.set.to_numpy()
    data = {}
    for inp in INPUTS:
        lz, lt = df[f"{inp}_z"].to_numpy(), df[f"{inp}_t"].to_numpy()
        data[inp] = (lz, lt, df.fz.to_numpy() * lz / (w_t * a.H), df.ft.to_numpy() * lt / (w_z * a.H))
    # pass 1: random starts; pass 2: re-seed every input from all inputs' elastin-active bests
    first = {inp: fit(*data[inp], a.starts, rng)[0] for inp in INPUTS}
    seeds = [p for p in first.values() if p is not None]
    results = {}
    for inp in INPUTS:
        lz, lt, tzz, ttt = data[inp]
        best, r2, rms = fit(lz, lt, tzz, ttt, 0, rng, seeds=seeds)
        boots, flipped = [], 0
        for b in range(a.boot):
            pick = rng.choice(sets, size=len(sets), replace=True)
            idx = np.concatenate([np.nonzero(sets_arr == s)[0] for s in pick])
            pb, _, _ = fit(lz[idx], lt[idx], tzz[idx], ttt[idx], 0, rng, seeds=[best])
            if pb is None:                                   # refit fell to the elastin-off minimum
                pb, _, _ = fit(lz[idx], lt[idx], tzz[idx], ttt[idx], 0, rng, seeds=seeds)
            if pb is None:
                flipped += 1
            else:
                boots.append(pb)
        boots = np.array(boots)
        results[inp] = dict(best=best, R2=r2, rms=rms, lo=np.percentile(boots, 2.5, 0), hi=np.percentile(boots, 97.5, 0),
                            peak_z=float(lz.max()), peak_t=float(lt.max()), n_boot=len(boots), flipped=flipped)
        print(f"{LABELS[inp]:20s} R2 {r2:.4f} RMS {rms:.2f} kPa | peak λz {lz.max():.3f} λθ {lt.max():.3f} | "
              f"bootstrap {len(boots)} kept, {flipped} elastin-off discarded | "
              + " ".join(f"{n}={v:.3g}" for n, v in zip(M.PARAM_NAMES, best)), flush=True)

    # ---- tables
    rows = []
    for inp, r in results.items():
        row = dict(input=inp, label=LABELS[inp], R2=r["R2"], rms_kPa=r["rms"], peak_lambda_z=r["peak_z"],
                   peak_lambda_theta=r["peak_t"], bootstrap_kept=r["n_boot"], bootstrap_elastin_off=r["flipped"])
        for i, n in enumerate(M.PARAM_NAMES):
            row[n], row[n + "_lo95"], row[n + "_hi95"] = r["best"][i], r["lo"][i], r["hi"][i]
        rows.append(row)
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(a.out, "params_by_stretch_input.csv"), index=False)
    pd.set_option("display.width", 200)
    print("\n", T[["label", "R2", "rms_kPa", "peak_lambda_z", "peak_lambda_theta", *M.PARAM_NAMES]].round(3).to_string(index=False))

    # ---- plots
    cols = {k: plt.get_cmap("tab10")(i) for i, k in enumerate(INPUTS)}
    # 1) stretch histories (loading branches concatenated)
    fig, ax = plt.subplots(2, 1, figsize=(15, 7.5), constrained_layout=True, sharex=True)
    xi = np.arange(len(df))
    for inp in INPUTS:
        ax[0].plot(xi, df[f"{inp}_z"], color=cols[inp], lw=1.2, label=LABELS[inp])
        ax[1].plot(xi, df[f"{inp}_t"], color=cols[inp], lw=1.2)
    for s in sets:
        i0 = np.nonzero(sets_arr == s)[0][0]
        for a_ in ax:
            a_.axvline(i0, color="k", lw=.4, alpha=.4)
        ax[0].text(i0 + 2, 0.97, f"set {s}", transform=ax[0].get_xaxis_transform(), fontsize=8, va="top")
    ax[0].set_ylabel(r"$\lambda_z$"); ax[1].set_ylabel(r"$\lambda_\theta$"); ax[1].set_xlabel("loading frame (protocols concatenated)")
    ax[0].legend(ncol=4, fontsize=9); [a_.grid(alpha=.3) for a_ in ax]
    ax[0].set_title("Global stretch inputs, loading branches of each protocol, relative to the common reference")
    fig.savefig(os.path.join(a.out, "stretch_inputs_histories.png"), dpi=90)

    # 2) parameters by input (elastin-active solution) with protocol-level bootstrap 95% intervals
    fig, ax = plt.subplots(2, 4, figsize=(17, 8.5), constrained_layout=True)
    for i, (n, u) in enumerate(zip(M.PARAM_NAMES, M.PARAM_UNITS)):
        a_ = ax.ravel()[i]
        for j, inp in enumerate(INPUTS):
            r = results[inp]
            a_.plot([j, j], [r["lo"][i], r["hi"][i]], color=cols[inp], lw=2, alpha=.4)
            a_.plot(j, r["best"][i], "o", ms=7, color=cols[inp])
        a_.set_xticks(range(len(INPUTS))); a_.set_xticklabels([LABELS[k] for k in INPUTS], rotation=60, ha="right", fontsize=8)
        a_.set_title(f"{n} ({u})"); a_.grid(alpha=.3)
    fig.suptitle("Fitted parameters by stretch input, elastin-active solution "
                 "(dots: best fit; bars: protocol-level bootstrap 95%)")
    fig.savefig(os.path.join(a.out, "params_by_stretch_input.png"), dpi=90)

    # 3) model responses along three reference paths, each input's parameters
    fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True, sharey=True)
    lz_ = np.linspace(1, 1.35, 100); lt_ = np.linspace(1, 1.175, 100)
    for inp in INPUTS:
        p = results[inp]["best"]
        mz, mt = M.cauchy_biaxial(lz_, lt_, p)                       # proportional (equibiaxial-type) path
        ax[0].plot(lz_, mz, color=cols[inp], lw=1.6, label=LABELS[inp])
        ax[0].plot(lt_, mt, color=cols[inp], lw=1.6, ls="--")
        mz, _ = M.cauchy_biaxial(lz_, np.ones_like(lz_), p)           # axial only
        ax[1].plot(lz_, mz, color=cols[inp], lw=1.6)
        _, mt = M.cauchy_biaxial(np.ones_like(lt_), lt_, p)           # circumferential only
        ax[2].plot(lt_, mt, color=cols[inp], lw=1.6, ls="--")
    ax[0].set_title(r"proportional path to $\lambda_z$ = 1.35, $\lambda_\theta$ = 1.175" "\n"
                    r"(solid $t_{zz}$ vs $\lambda_z$, dashed $t_{\theta\theta}$ vs $\lambda_\theta$)", fontsize=10)
    ax[1].set_title(r"axial only ($\lambda_\theta$ = 1): $t_{zz}$ vs $\lambda_z$", fontsize=10)
    ax[2].set_title(r"circumferential only ($\lambda_z$ = 1): $t_{\theta\theta}$ vs $\lambda_\theta$", fontsize=10)
    for a_ in ax:
        a_.set_xlabel("stretch"); a_.grid(alpha=.3)
    ax[0].set_ylabel("model Cauchy stress (kPa)"); ax[0].legend(fontsize=8)
    fig.suptitle("Model response from each input's parameters (elastin-active solution)")
    fig.savefig(os.path.join(a.out, "model_response_by_stretch_input.png"), dpi=90)
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
