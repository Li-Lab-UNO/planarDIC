"""
Frozen bulk (homogeneous) four-fiber baseline for comparing future parameterizations,
e.g. local strain-based or zone-based fits.

freeze   Copy the current compare_stretch_inputs.py results into
         <test>/constitutive_fit/baseline/ together with the exact data set
         (stretches, tared forces, stresses), the settings, per-protocol errors, and
         SHA-256 checksums of the data and the scripts that produced them.
         Primary baseline input: DIC affine F; the other inputs are kept for context.
evaluate Score any parameter set on the frozen data set with the same stress definition:
         overall R^2 / RMS and per-protocol RMS, side by side with the baseline.

    .venv/bin/python constitutive_fit/baseline.py --test N1739/test1 freeze
    .venv/bin/python constitutive_fit/baseline.py --test N1739/test1 freeze --tracks standard     (standard_roi/baseline/)
    .venv/bin/python constitutive_fit/baseline.py --test N1739/test1 evaluate --params "5.9,3.6,1.8,9.4,9.9,1.8,7.2,49"
    .venv/bin/python constitutive_fit/baseline.py --test N1739/test1 evaluate --params-json my_fit.json --input dic_affine

In Python:  import testcfg; testcfg.select("N1739/test1"); import baseline; baseline.evaluate([...])  -> dict
            baseline.load(tracks="tine_region") / evaluate(..., tracks="tine_region") for the
            baseline built on the tracks that reach the tine tips (track_tine_region.py)
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import constitutive as M  # noqa: E402

TC = testcfg.get()
SRC = TC.fit_dir                                   # compare_stretch_inputs.py output
OUT = os.path.join(TC.fit_dir, "baseline")
# DIC tracks -> (stretch_inputs folder, baseline folder, ROI setting of the tracked grid)
TRACKS = {"tine_region": (SRC, OUT, TC.config.get("roi_px", {}).get("tine_region")),
          "standard": (os.path.join(SRC, "standard_roi"), os.path.join(SRC, "standard_roi", "baseline"),
                       TC.config.get("roi_px", {}).get("standard"))}
PRIMARY = "dic_affine"
SCRIPTS = ["dic2d.py", "track_protocols.py", "compare_stretch_inputs.py", "constitutive.py", "baseline.py"]


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _stress(df, inp, H, w_t, w_z):
    lz, lt = df[f"{inp}_z"].to_numpy(), df[f"{inp}_t"].to_numpy()
    return lz, lt, df["fz"].to_numpy() * lz / (w_t * H), df["ft"].to_numpy() * lt / (w_z * H)


def _score(params, df, inp, H, w_t, w_z):
    lz, lt, tzz, ttt = _stress(df, inp, H, w_t, w_z)
    mz, mt = M.cauchy_biaxial(lz, lt, np.asarray(params, float))
    r = np.r_[mz - tzz, mt - ttt]
    y = np.r_[tzz, ttt]
    per = {}
    for s, idx in df.groupby("set").indices.items():
        per[int(s)] = dict(name=df["name"].iloc[idx[0]],
                           rms_z=float(np.sqrt(np.mean((mz[idx] - tzz[idx]) ** 2))),
                           rms_theta=float(np.sqrt(np.mean((mt[idx] - ttt[idx]) ** 2))))
    return dict(R2=float(1 - np.sum(r ** 2) / np.sum((y - y.mean()) ** 2)), rms_kPa=float(np.sqrt(np.mean(r ** 2))),
                per_protocol=per)


def freeze(H=None, H_sd=None, widths=None, tracks="tine_region"):
    """Defaults: the test's thickness_mm and widths_mm.tine_span (test_config.json)."""
    H = TC.setting("thickness_mm", "H") if H is None else H
    H_sd = TC.setting("thickness_mm", "H_sd") if H_sd is None else H_sd
    widths = tuple(TC.setting("widths_mm", "tine_span")) if widths is None else widths
    SRC, OUT, roi = TRACKS[tracks]
    os.makedirs(OUT, exist_ok=True)
    data_csv = os.path.join(OUT, "baseline_data.csv")
    shutil.copyfile(os.path.join(SRC, "stretch_inputs_data.csv"), data_csv)
    T = pd.read_csv(os.path.join(SRC, "params_by_stretch_input.csv"))
    df = pd.read_csv(data_csv)
    w_t, w_z = widths
    inputs = {}
    for _, row in T.iterrows():
        p = [float(row[n]) for n in M.PARAM_NAMES]
        sc = _score(p, df, row["input"], H, w_t, w_z)
        inputs[row["input"]] = dict(
            label=row["label"], params=dict(zip(M.PARAM_NAMES, p)),
            bootstrap_95={n: [float(row[n + "_lo95"]), float(row[n + "_hi95"])] for n in M.PARAM_NAMES},
            peak_lambda=dict(z=float(row["peak_lambda_z"]), theta=float(row["peak_lambda_theta"])), **sc)
    base = dict(
        description=f"Homogeneous bulk four-fiber family fit (Jadidi et al. 2019, Eq. 2), {TC.config.get('sample', TC.name)}",
        frozen=datetime.date.today().isoformat(), primary_input=PRIMARY,
        settings=dict(protocols="PhaseTiming Sets 2-10 = Def11-Def19",
                      branch="loading only: Stretch start to the rake peak (compare_stretch_inputs.loading_end)",
                      reference=f"first frame of {df['folder'].iloc[0].split('_')[0]} (t = {df['t'].iloc[0]:.1f} s), "
                                "common to all protocols",
                      force_tare="force at the start of each protocol subtracted",
                      stress="t_zz = F_z*lambda_z/(w_theta*H), t_thetatheta = F_theta*lambda_theta/(w_z*H)",
                      H_mm=H, H_sd_mm=H_sd, w_theta_mm=w_t, w_z_mm=w_z, width_choice="tine_span",
                      solution_family="elastin active (C1_el > 0.1 kPa); elastin-off minimum excluded",
                      objective="unweighted sum of squared Cauchy-stress residuals, both directions",
                      uncertainty="protocol-level bootstrap, 200 resamples",
                      dic=dict(subset_px=61, grid_px=15, strain_window_px=75, roi_px=roi, tracks=tracks),
                      param_units=dict(zip(M.PARAM_NAMES, M.PARAM_UNITS))),
        inputs=inputs,
        checksums=dict(data={"baseline_data.csv": _sha(data_csv)},
                       scripts={s: _sha(testcfg.script_path(s)) for s in SCRIPTS}))
    with open(os.path.join(OUT, "baseline.json"), "w") as f:
        json.dump(base, f, indent=2)
    # human-readable table
    rows = [dict(input=k, label=v["label"], R2=v["R2"], rms_kPa=v["rms_kPa"], **v["params"]) for k, v in inputs.items()]
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "baseline_params.csv"), index=False)
    for fn in ("params_by_stretch_input.png", "fit_stress_stretch_dic_affine.png"):
        if os.path.exists(os.path.join(SRC, fn)):
            shutil.copyfile(os.path.join(SRC, fn), os.path.join(OUT, fn))
    print(f"baseline frozen in {OUT} (primary input: {PRIMARY})")
    print(pd.DataFrame(rows).round(3).to_string(index=False))


def load(tracks="tine_region"):
    OUT = TRACKS[tracks][1]
    with open(os.path.join(OUT, "baseline.json")) as f:
        base = json.load(f)
    df = pd.read_csv(os.path.join(OUT, "baseline_data.csv"))
    if _sha(os.path.join(OUT, "baseline_data.csv")) != base["checksums"]["data"]["baseline_data.csv"]:
        raise RuntimeError("baseline_data.csv does not match the frozen checksum")
    return base, df


def evaluate(params, input=PRIMARY, verbose=True, tracks="tine_region"):
    """Score a parameter vector (order: constitutive.PARAM_NAMES) on the frozen data set."""
    base, df = load(tracks)
    s = base["settings"]
    sc = _score(params, df, input, s["H_mm"], s["w_theta_mm"], s["w_z_mm"])
    ref = base["inputs"][input]
    if verbose:
        print(f"evaluated on the frozen {input} data set ({len(df)} loading frames, Def11-19, {tracks} tracks)")
        print(f"  candidate: R2 {sc['R2']:.4f}, RMS {sc['rms_kPa']:.3f} kPa")
        print(f"  baseline:  R2 {ref['R2']:.4f}, RMS {ref['rms_kPa']:.3f} kPa")
        print(f"  {'param':10s} {'candidate':>10s} {'baseline':>10s} {'baseline 95%':>18s}")
        for n, v in zip(M.PARAM_NAMES, params):
            lo, hi = ref["bootstrap_95"][n]
            flag = "" if lo <= v <= hi else "  outside"
            print(f"  {n:10s} {v:10.3f} {ref['params'][n]:10.3f}   [{lo:7.3f}, {hi:7.3f}]{flag}")
        print("  per-protocol RMS z / theta (kPa), candidate vs baseline:")
        for k, v in sc["per_protocol"].items():
            b = ref["per_protocol"][str(k)]
            print(f"    set {k:2d} {v['name']:22s} {v['rms_z']:5.2f} / {v['rms_theta']:5.2f}   vs   "
                  f"{b['rms_z']:5.2f} / {b['rms_theta']:5.2f}")
    return sc


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f_ = sub.add_parser("freeze")
    f_.add_argument("--tracks", choices=tuple(TRACKS), default="tine_region")
    e = sub.add_parser("evaluate")
    e.add_argument("--tracks", choices=tuple(TRACKS), default="tine_region")
    e.add_argument("--params", help="comma-separated, order: " + ", ".join(M.PARAM_NAMES))
    e.add_argument("--params-json", help="JSON file with a 'params' dict keyed by parameter name")
    e.add_argument("--input", default=PRIMARY)
    a = ap.parse_args()
    if a.cmd == "freeze":
        freeze(tracks=a.tracks)
    else:
        if a.params_json:
            pj = json.load(open(a.params_json))["params"]
            p = [pj[n]["value"] if isinstance(pj[n], dict) else pj[n] for n in M.PARAM_NAMES]
        else:
            p = [float(v) for v in a.params.split(",")]
        evaluate(p, a.input, tracks=a.tracks)
