"""
Compare measure_thickness.py results with the manual readings in each sample's
Thickness, …/Thickess measurements.xlsx (5 caliper lines on each of 2 images per test).

Each image is matched by the photo name in its label (e.g. "Image 1 (DSC_8630)"). The
N1739 workbook's labels were copied from N1690's and name no N1739 photo, so there
image 1 / image 2 are the first / second photo (by name) in the test's folder. The automatic value per image
pools the stations of its selected strips (see measure_thickness.py); unselected strips are listed in the CSV.

    .venv/bin/python thickness/compare_manual.py

Run measure_thickness.py on the thickness folders first. Writes
comparisons/thickness_auto_vs_manual.csv (per image), _by_test.csv and .png.
"""
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import openpyxl  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure_thickness as MT  # noqa: E402
OUT = os.path.join(ROOT, "comparisons")


def manual_readings(xlsx):
    """{test: [(label, [5 mm readings of image 1]), (label, [image 2])]} and {test: [px, px]}."""
    ws = openpyxl.load_workbook(xlsx, data_only=True).active
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    readings, scale = {}, {}
    for r, row in enumerate(rows):
        for c, v in enumerate(row):
            if isinstance(v, str) and re.match(r"Test\s*\d", v.strip()):
                t = int(re.findall(r"\d", v)[0])
                imgs, cur = [], None
                for row2 in rows[r + 1:]:
                    lab, mm = row2[c], row2[c + 2] if c + 2 < len(row2) else None
                    if isinstance(row2[c + 1], str) and row2[c + 1].startswith(("Average", "S.D")):
                        break
                    if isinstance(lab, str) and lab.startswith("Image"):
                        cur = []
                        imgs.append((lab, cur))
                    if cur is not None and isinstance(mm, (int, float)):
                        cur.append(float(mm))
                readings[t] = imgs
            if isinstance(v, str) and v.startswith("Distance in pixels"):
                # this test's column block starts one column left of the label
                px = [rows[r][c + 1], rows[r + 1][c + 1]]
                scale[c] = [float(p) for p in px if isinstance(p, (int, float))]
    # scale blocks in column order = tests 1, 2, 3
    scale = {i + 1: scale[c] for i, c in enumerate(sorted(scale))}
    return readings, scale


def main():
    per_image, per_test = [], []
    for xlsx in sorted(glob.glob(os.path.join(ROOT, "N*", "Thickness*", "*.xlsx"))):
        sample = os.path.basename(os.path.dirname(os.path.dirname(xlsx)))
        readings, _ = manual_readings(xlsx)
        for t, imgs in sorted(readings.items()):
            tdir = [d for d in glob.glob(os.path.join(os.path.dirname(xlsx), "*"))
                    if os.path.isdir(d) and re.sub(r"\D", "", os.path.basename(d)) == str(t)][0]
            odir = MT.output_location(tdir)
            summ = os.path.join(odir, "summary.csv")
            auto = pd.read_csv(summ)
            photos = sorted(f for f in os.listdir(tdir) if f.lower().endswith((".jpg", ".jpeg")))
            pooled_auto = []
            for k, (label, man) in enumerate(imgs):
                # the label names the photo (e.g. "Image 1 (DSC_8630)") unless it is not in
                # this folder (N1739's labels were copied from N1690): then go by order
                named = [f for f in photos if os.path.splitext(f)[0] in label]
                img = named[0] if named else photos[k]
                a = auto[auto.image == img]
                c = a[a.selected]
                v = np.concatenate([pd.read_csv(os.path.join(odir, f"{os.path.splitext(img)[0]}_strip{int(k)}_stations.csv"))
                                    .query("ok").thickness_mm.values for k in c.strip])
                pooled_auto.append(v)
                per_image.append(dict(
                    sample=sample, test=t, image=img, by="label" if named else "order",
                    manual_mean=np.mean(man),
                    manual_sd=np.std(man, ddof=1), manual_n=len(man),
                    auto_mean=v.mean(), auto_sd=v.std(ddof=1), auto_median=np.median(v), auto_n=len(v),
                    auto_strips=", ".join(str(int(k)) for k in c.strip),
                    auto_px_per_mm=c.px_per_mm.iloc[0], selected_by=c.selected_by.iloc[0],
                    diff_mm=v.mean() - np.mean(man),
                    other_strips="; ".join(f"{int(r.strip)}: {r.mean_mm:.2f}±{r.sd_mm:.2f}"
                                           for r in a[~a.selected].itertuples())))
            man_all = np.concatenate([m for _, m in imgs])
            auto_all = np.concatenate(pooled_auto)
            per_test.append(dict(sample=sample, test=t, manual_mean=man_all.mean(),
                                 manual_sd=man_all.std(ddof=1), manual_n=len(man_all),
                                 auto_mean=auto_all.mean(), auto_sd=auto_all.std(ddof=1),
                                 auto_n=len(auto_all), diff_mm=auto_all.mean() - man_all.mean(),
                                 diff_pct=100 * (auto_all.mean() / man_all.mean() - 1)))
    pi, pt = pd.DataFrame(per_image), pd.DataFrame(per_test)
    os.makedirs(OUT, exist_ok=True)
    pi.to_csv(os.path.join(OUT, "thickness_auto_vs_manual.csv"), index=False)
    pt.to_csv(os.path.join(OUT, "thickness_auto_vs_manual_by_test.csv"), index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    for s, mk in zip(sorted(pi["sample"].unique()), "os"):
        d = pi[pi["sample"] == s]
        ax.errorbar(d.manual_mean, d.auto_mean, xerr=d.manual_sd, yerr=d.auto_sd, fmt=mk,
                    capsize=3, label=s)
    lim = [0.8, 3.6]
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.set(xlim=lim, ylim=lim, xlabel="manual (mm, 5 lines)", ylabel="automatic (mm, selected strips)",
           title="per image, mean ± SD")
    ax.legend()
    ax = axes[1]
    x = np.arange(len(pt))
    ax.bar(x - 0.2, pt.manual_mean, 0.4, yerr=pt.manual_sd, capsize=3, label="manual (n=10)")
    ax.bar(x + 0.2, pt.auto_mean, 0.4, yerr=pt.auto_sd, capsize=3, label="automatic (same 2 images)")
    ax.set_xticks(x, [f"{s[:5]}\ntest{t}" for s, t in zip(pt["sample"], pt.test)])
    ax.set(ylabel="thickness (mm)", title="per test, mean ± SD")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "thickness_auto_vs_manual.png"), dpi=130)

    pd.set_option("display.width", 200)
    print(pi.drop(columns=["other_strips"]).round(3).to_string(index=False))
    print()
    print(pt.round(3).to_string(index=False))
    d = pi.diff_mm
    print(f"\nper image: bias {d.mean():+.3f} mm, SD of differences {d.std(ddof=1):.3f} mm, "
          f"limits of agreement {d.mean() - 1.96 * d.std(ddof=1):+.2f} to {d.mean() + 1.96 * d.std(ddof=1):+.2f} mm")


if __name__ == "__main__":
    main()
