"""
Compare strain bands across tests at the peak of one deformation.

For each test, the exported fields of the full-domain track (export_fields.py on
DefNN_*.npz) are taken at the frame where the mean axial stretch peaks. Maps are
drawn on the reference grid in mm from the left (z) and bottom (theta) tine tips, theta up.

Band measures, per test (valid field points only):
  concentration    E / spatial mean of E, so specimens with different peak strains compare
  vertical bands   E_zz averaged along theta, vs z   (bands running in theta, varying along z)
  horizontal bands E_thetatheta averaged along z, vs theta
  contrast         peak-to-peak and SD of each profile, as a fraction of the mean
  control          the same profile taken in the other direction (E_zz vs theta, E_thetatheta vs z)
  wavelength       spacing of the strongest periodic component of each profile (FFT of the
                   detrended profile; only meaningful if it is well below the profile length)
Band positions are not comparable between specimens unless they were mounted in the same
orientation; amplitude, direction and spacing are.

    .venv/bin/python strain_tracking/compare_bands.py --sample N1739 --tests test1 test2 test3 --def Def11
Outputs: comparisons/bands_<Def>_<tests>.png / .csv
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
import testcfg  # noqa: E402

OUT = os.path.join(testcfg.ROOT, "comparisons")


def load(test, defn, suffix):
    T = testcfg.Test(testcfg._resolve(test))
    fd = T.folder(defn)
    p = os.path.join(T.results, f"{fd}{suffix}_fields.npz")
    if not os.path.exists(p):
        sys.exit(f"missing {p}; run export_fields.py on {fd}{suffix}.npz first")
    d = np.load(p, allow_pickle=True)
    lz = np.where(d["valid_field"], d["lambda_z"], np.nan)
    k = int(np.nanargmax(np.nanmean(lz.reshape(len(lz), -1), 1)))       # peak mean axial stretch
    mpp = float(d["mm_per_px"])
    tips = T.setting("tine_tips_px")
    z = (d["X"] - tips["z"][0]) * mpp                                     # mm from the left tips
    th = (tips["theta"][1] - d["Y"]) * mpp                                # mm from the bottom tips, up
    ok = d["valid_field"][k]
    E = {c: np.where(ok, d[c][k], np.nan) for c in ("E_zz", "E_thetatheta")}
    return dict(name=T.name, short=T.name.split("_")[0], z=z, th=th, E=E, t=float(d["times"][k]),
                Fz=float(d["force_z_mN"][k]), Ft=float(d["force_theta_mN"][k]), ok=ok,
                lz=float(np.nanmean(d["lambda_z"][k][ok])), lt=float(np.nanmean(d["lambda_theta"][k][ok])))


def wavelength(prof, step):
    p = prof[np.isfinite(prof)]
    if len(p) < 8:
        return np.nan
    x = np.arange(len(p))
    p = p - np.polyval(np.polyfit(x, p, 1), x)                            # remove the linear trend
    A = np.abs(np.fft.rfft(p * np.hanning(len(p))))
    f = np.fft.rfftfreq(len(p), d=step)
    j = 2 + int(np.argmax(A[2:])) if len(A) > 3 else 1                    # skip DC and the longest mode
    return 1 / f[j]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", default="N1739", help="sample folder (unique prefix) for bare test names")
    ap.add_argument("--tests", nargs="+", default=["test1", "test2", "test3"],
                    help="test names in --sample, or SAMPLE/TEST")
    ap.add_argument("--def", dest="defn", default="Def11")
    ap.add_argument("--suffix", default="", help="track suffix ('_standard' for the standard grid)")
    a = ap.parse_args()
    R = [load(t if "/" in t else f"{a.sample}/{t}", a.defn, a.suffix) for t in a.tests]
    n = len(R)
    fig, ax = plt.subplots(n, 6, figsize=(27, 4.6 * n), constrained_layout=True, squeeze=False,
                           gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 1.15, 1.15]))
    rows = []
    for i, r in enumerate(R):
        z, th = r["z"], r["th"]
        step = abs(z[0, 1] - z[0, 0])
        ext = (z.min(), z.max(), th.min(), th.max())
        for j, (c, lab) in enumerate((("E_zz", r"$E_{zz}$"), ("E_thetatheta", r"$E_{\theta\theta}$"))):
            E = r["E"][c]
            m = np.nanmean(E)
            im = ax[i, j].pcolormesh(z, th, E, cmap="viridis", shading="nearest")
            fig.colorbar(im, ax=ax[i, j], shrink=.85)
            ax[i, j].set_title(f"{r['short']}: {lab} (mean {m:.3f})")
            im2 = ax[i, j + 2].pcolormesh(z, th, E / m, cmap="RdBu_r", vmin=0.6, vmax=1.4, shading="nearest")
            fig.colorbar(im2, ax=ax[i, j + 2], shrink=.85)
            ax[i, j + 2].set_title(f"{r['short']}: {lab} / mean")
            for q in (ax[i, j], ax[i, j + 2]):
                q.set_aspect("equal"); q.set_xlim(ext[0], ext[1]); q.set_ylim(ext[2], ext[3])
                q.set_xlabel("z, axial (mm)"); q.set_ylabel(r"$\theta$, circumferential (mm)")
        Ezz, Ett = r["E"]["E_zz"], r["E"]["E_thetatheta"]
        zc, tc = np.nanmean(z, 0), np.nanmean(th, 1)
        pz = np.nanmean(Ezz, 0) / np.nanmean(Ezz)          # vertical bands: E_zz vs z
        pt = np.nanmean(Ett, 1) / np.nanmean(Ett)          # horizontal bands: E_tt vs theta
        cz = np.nanmean(Ezz, 1) / np.nanmean(Ezz)          # control: E_zz vs theta
        ct = np.nanmean(Ett, 0) / np.nanmean(Ett)          # control: E_tt vs z
        a_ = ax[i, 4]
        a_.plot(zc, pz, "-o", ms=3, label=r"$E_{zz}$ vs z (vertical bands)")
        a_.plot(zc, ct, "--", color="0.5", label=r"control: $E_{\theta\theta}$ vs z")
        a_.axhline(1, color="k", lw=.6); a_.set_ylim(0.6, 1.4); a_.grid(alpha=.3)
        a_.set_xlabel("z, axial (mm)"); a_.set_ylabel("profile / mean"); a_.legend(fontsize=8)
        a_.set_title(f"{r['short']}: profiles along z")
        a_ = ax[i, 5]
        a_.plot(tc, pt, "-o", ms=3, color="C3", label=r"$E_{\theta\theta}$ vs $\theta$ (horizontal bands)")
        a_.plot(tc, cz, "--", color="0.5", label=r"control: $E_{zz}$ vs $\theta$")
        a_.axhline(1, color="k", lw=.6); a_.set_ylim(0.6, 1.4); a_.grid(alpha=.3)
        a_.set_xlabel(r"$\theta$, circumferential (mm)"); a_.set_ylabel("profile / mean"); a_.legend(fontsize=8)
        a_.set_title(f"{r['short']}: profiles along $\\theta$")
        rows.append(dict(test=r["name"], deformation=a.defn, t_s=r["t"], force_z_mN=r["Fz"], force_theta_mN=r["Ft"],
                         mean_lambda_z=r["lz"], mean_lambda_theta=r["lt"],
                         mean_E_zz=float(np.nanmean(Ezz)), mean_E_thetatheta=float(np.nanmean(Ett)),
                         Ezz_vs_z_p2p=float(np.nanmax(pz) - np.nanmin(pz)), Ezz_vs_z_sd=float(np.nanstd(pz)),
                         Ezz_vs_theta_p2p_control=float(np.nanmax(cz) - np.nanmin(cz)),
                         Ett_vs_theta_p2p=float(np.nanmax(pt) - np.nanmin(pt)), Ett_vs_theta_sd=float(np.nanstd(pt)),
                         Ett_vs_z_p2p_control=float(np.nanmax(ct) - np.nanmin(ct)),
                         Ezz_vs_z_wavelength_mm=wavelength(pz, step), Ett_vs_theta_wavelength_mm=wavelength(pt, step),
                         domain_z_mm=float(z.max() - z.min()), domain_theta_mm=float(th.max() - th.min()),
                         map_sd_Ezz_over_mean=float(np.nanstd(Ezz) / np.nanmean(Ezz)),
                         map_sd_Ett_over_mean=float(np.nanstd(Ett) / np.nanmean(Ett))))
    tag = f"bands_{a.defn}_{'_'.join(r['short'] for r in R)}"
    fig.suptitle(f"{a.defn} at peak mean axial stretch: Green-Lagrange strain maps and band profiles "
                 f"(full-domain tracks; z = axial, θ = circumferential, up)", fontsize=13)
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, tag + ".png"), dpi=80)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, tag + ".csv"), index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(df.drop(columns=["test"]).round(3).T.rename(columns={i: R[i]["short"] for i in range(n)}))
    print("wrote", os.path.join(OUT, tag + ".png"), "and .csv")


if __name__ == "__main__":
    main()
