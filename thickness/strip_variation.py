"""
How much, and over what distance, wall thickness varies along the thickness strips, as a reference for local
stiffness maps (a membrane stiffness map cannot be told apart from thickness variation of the same size and scale).

The strips are cut from tissue next to the tested specimen, so they give statistics of H, not H on the tested region.
Some strips were photographed more than once. Per test (the strips that measure_thickness.py pools: selected, scale
reliable):

  1. Duplicates. Every pair of strips is aligned along the centerline (shift, and both directions in case the strip
     was flipped) to the best correlation of their thickness profiles. A pair is one physical strip photographed twice
     when the correlation is >= DUP_R, the means differ by <= DUP_MEAN and the lengths by <= DUP_LEN. Groups are
     the connected pairs.
  2. Measurement noise, from the duplicates: d(s) = H_1(s) - H_2(s) on the overlap. The SD of d about its mean,
     divided by sqrt(2), is the per-photo noise; the mean of d over the strip mean is the photo-to-photo offset
     (mostly scale).
  3. Per physical strip (the longest photo stands for it): mean, SD, the SD with the noise removed and its CV, the
     share of variance in a linear trend along the strip, and the correlation length (1/e point of the
     autocorrelation of the mean-removed profile, linearly interpolated; ACF_MAX_FRAC of the profile at most).
  4. Between strips of a test: the SD of the strip means.

    .venv/bin/python thickness/strip_variation.py

Run measure_thickness.py on the thickness folders first. Writes comparisons/thickness_strip_variation.csv (per
physical strip), _pairs.csv (aligned pairs with their match statistics), _by_test.csv and .png.
"""
import glob
import itertools
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "comparisons")

STEP_MM = 0.1         # station spacing (measure_thickness.STEP_MM)
MIN_OVERLAP = 0.7     # aligned profiles must overlap by at least this fraction of the shorter one
DUP_R = 0.8           # duplicate: profile correlation at least this ...
DUP_MEAN = 0.10       # ... means within this fraction ...
DUP_LEN = 0.10        # ... and lengths within this fraction
ACF_MAX_FRAC = 0.5    # autocorrelation evaluated to this fraction of the profile length


def load_strips():
    """One row per strip that measure_thickness.py pools, with its profile (thickness on a STEP_MM grid)."""
    rows = []
    for d in sorted(glob.glob(os.path.join(ROOT, "N*", "Thickness*", "thickness_auto", "test*"))):
        sample = d.split(os.sep)[-4]
        test = int(re.sub(r"\D", "", os.path.basename(d)))
        summ = pd.read_csv(os.path.join(d, "summary.csv"))
        for r in summ[summ.selected & summ.scale_ok].itertuples():
            st = pd.read_csv(os.path.join(d, f"{os.path.splitext(r.image)[0]}_strip{int(r.strip)}_stations.csv"))
            s = st.s_mm.to_numpy()
            k = np.round((s - s.min()) / STEP_MM).astype(int)   # from the first station: s_mm can sit on x.x5
            h = np.full(k.max() + 1, np.nan)
            h[k] = np.where(st.ok, st.thickness_mm, np.nan)
            rows.append(dict(sample=sample, test=test, image=r.image, strip=int(r.strip), length_mm=r.length_mm,
                             px_per_mm=r.px_per_mm, h=h))
    return rows


def align(a, b):
    """Best alignment of profile b to a (shift and direction): (r, shift, flipped, overlap mask pair)."""
    best = (-np.inf, 0, False)
    nmin = int(MIN_OVERLAP * min(len(a), len(b)))
    for flip in (False, True):
        bb = b[::-1] if flip else b
        for sh in range(-(len(bb) - nmin), len(a) - nmin + 1):
            i0, i1 = max(0, sh), min(len(a), sh + len(bb))
            x, y = a[i0:i1], bb[i0 - sh:i1 - sh]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < nmin:
                continue
            r = np.corrcoef(x[ok], y[ok])[0, 1]
            if r > best[0]:
                best = (r, sh, flip)
    return best


def overlap(a, b, sh, flip):
    bb = b[::-1] if flip else b
    i0, i1 = max(0, sh), min(len(a), sh + len(bb))
    x, y = a[i0:i1], bb[i0 - sh:i1 - sh]
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok]


def corr_length(h, noise_var=0.0):
    """1/e length (mm) of the autocorrelation of the mean-removed profile. Lag 0 uses the variance with the
    measurement noise removed, so uncorrelated noise does not shorten it. NaN if it never falls to 1/e."""
    x = h - np.nanmean(h)
    var = max(np.nanvar(x) - noise_var, 1e-12)
    n = len(x)
    lags = np.arange(1, int(ACF_MAX_FRAC * n))
    acf = np.array([np.nanmean(x[:-L] * x[L:]) for L in lags]) / var
    acf = np.concatenate([[1.0], acf])
    below = np.nonzero(acf < np.exp(-1))[0]
    if not len(below):
        return np.nan, acf
    j = below[0]
    l0, l1, a0, a1 = j - 1, j, acf[j - 1], acf[j]
    return STEP_MM * (l0 + (a0 - np.exp(-1)) / (a0 - a1) * (l1 - l0)), acf


def main():
    strips = load_strips()
    pairs, groups = [], []
    for (sample, test), idx in itertools.groupby(range(len(strips)),
                                                 key=lambda i: (strips[i]["sample"], strips[i]["test"])):
        idx = list(idx)
        parent = {i: i for i in idx}

        def root(i):
            while parent[i] != i:
                i = parent[i]
            return i

        for i, j in itertools.combinations(idx, 2):
            a, b = strips[i], strips[j]
            r, sh, flip = align(a["h"], b["h"])
            x, y = overlap(a["h"], b["h"], sh, flip)
            ma, mb = np.nanmean(a["h"]), np.nanmean(b["h"])
            dmean = abs(ma - mb) / np.mean([ma, mb])
            dlen = abs(a["length_mm"] - b["length_mm"]) / np.mean([a["length_mm"], b["length_mm"]])
            dup = r >= DUP_R and dmean <= DUP_MEAN and dlen <= DUP_LEN
            d = x - y
            pairs.append(dict(sample=sample, test=test, a=f"{a['image']}#{a['strip']}", b=f"{b['image']}#{b['strip']}",
                              r=r, shift_mm=sh * STEP_MM, flipped=flip, mean_a=ma, mean_b=mb, dmean_frac=dmean,
                              dlen_frac=dlen, duplicate=dup, n_overlap=len(d),
                              offset_frac=d.mean() / np.mean([ma, mb]) if dup else np.nan,
                              noise_sd_mm=d.std(ddof=1) / np.sqrt(2) if dup else np.nan))
            if dup:
                parent[root(j)] = root(i)
        for g in sorted({root(i) for i in idx}):
            groups.append((sample, test, [i for i in idx if root(i) == g]))
    P = pd.DataFrame(pairs)

    dups = P[P.duplicate]
    noise_sd = float(np.sqrt(np.mean(dups.noise_sd_mm ** 2))) if len(dups) else 0.0
    # noise autocorrelation length from the duplicate differences
    ncl = []
    for p in dups.itertuples():
        a = next(s for s in strips if f"{s['image']}#{s['strip']}" == p.a and s["sample"] == p.sample and s["test"] == p.test)
        b = next(s for s in strips if f"{s['image']}#{s['strip']}" == p.b and s["sample"] == p.sample and s["test"] == p.test)
        x, y = overlap(a["h"], b["h"], int(round(p.shift_mm / STEP_MM)), p.flipped)
        ncl.append(corr_length(x - y)[0])
    noise_cl = float(np.nanmedian(ncl)) if ncl else np.nan

    rows, acfs = [], []
    for sample, test, members in groups:
        rep = max(members, key=lambda i: np.isfinite(strips[i]["h"]).sum())
        h = strips[rep]["h"]
        v = h[np.isfinite(h)]
        s = np.arange(len(h)) * STEP_MM
        ok = np.isfinite(h)
        coef = np.polyfit(s[ok], h[ok], 1)
        trend_share = np.var(np.polyval(coef, s[ok])) / np.var(v)
        sd_true = np.sqrt(max(v.var(ddof=1) - noise_sd ** 2, 0.0))
        cl, acf = corr_length(h, noise_sd ** 2)
        acfs.append((sample, test, len(rows), acf))
        rows.append(dict(sample=sample, test=test, strip_id=len([r for r in rows if r["sample"] == sample and r["test"] == test]) + 1,
                         photos="; ".join(f"{strips[i]['image']}#{strips[i]['strip']}" for i in members),
                         n_photos=len(members), length_mm=strips[rep]["length_mm"], mean_mm=v.mean(),
                         sd_mm=v.std(ddof=1), sd_true_mm=sd_true, cv_true=sd_true / v.mean(),
                         slope_mm_per_mm=coef[0], trend_share=trend_share, corr_length_mm=cl))
    S = pd.DataFrame(rows)

    by = []
    for (sample, test), g in S.groupby(["sample", "test"]):
        # within-strip variance (noise removed) and between-strip variance of the means, per physical strip
        within = np.mean(g.sd_true_mm ** 2)
        between = g.mean_mm.var(ddof=1) if len(g) > 1 else np.nan
        m = g.mean_mm.mean()
        by.append(dict(sample=sample, test=test, n_strips=len(g), n_photos=g.n_photos.sum(), mean_mm=m,
                       cv_within=np.sqrt(within) / m, cv_between=np.sqrt(between) / m if len(g) > 1 else np.nan,
                       cv_total=np.sqrt(within + (between if len(g) > 1 else 0.0)) / m,
                       corr_length_mm_median=g.corr_length_mm.median(),
                       strip_means="; ".join(f"{x:.2f}" for x in g.mean_mm)))
    B = pd.DataFrame(by)

    os.makedirs(OUT, exist_ok=True)
    S.to_csv(os.path.join(OUT, "thickness_strip_variation.csv"), index=False)
    P.to_csv(os.path.join(OUT, "thickness_strip_variation_pairs.csv"), index=False)
    B.to_csv(os.path.join(OUT, "thickness_strip_variation_by_test.csv"), index=False)

    # figure: profiles per test (one line per physical strip, duplicates overlaid), CV and correlation length
    keys = list(B[["sample", "test"]].itertuples(index=False, name=None))
    nc = 4
    nr = int(np.ceil(len(keys) / nc)) + 1
    fig, axes = plt.subplots(nr, nc, figsize=(4 * nc, 2.8 * nr), squeeze=False)
    for ax, (sample, test) in zip(axes.flat, keys):
        for k, (smp, tst, members) in enumerate([g for g in groups if g[0] == sample and g[1] == test]):
            ref = members[0]
            for i in members:
                h = strips[i]["h"]
                sh, flip = 0, False
                if i != ref:
                    _, sh, flip = align(strips[ref]["h"], h)
                hh = h[::-1] if flip else h
                ax.plot((np.arange(len(hh)) + sh) * STEP_MM, hh, color=f"C{k}", lw=1, alpha=0.8)
        ax.set(title=f"{sample[:5]} test{test}", xlabel="along strip (mm)", ylabel="H (mm)")
    for ax in axes.flat[len(keys):]:
        ax.axis("off")
    ax = axes[-1, 0]
    ax.axis("on")
    x = np.arange(len(B))
    ax.bar(x - 0.2, B.cv_within, 0.4, label="within strip")
    ax.bar(x + 0.2, B.cv_between.fillna(0), 0.4, label="between strips")
    ax.set_xticks(x, [f"{s[:5]}\nt{t}" for s, t in zip(B["sample"], B.test)], fontsize=7)
    ax.set(ylabel="CV of H (noise removed)", title="thickness variation")
    ax.legend(fontsize=8)
    ax = axes[-1, 1]
    ax.axis("on")
    ax.plot(S.mean_mm, S.corr_length_mm, "o")
    ax.axhline(noise_cl, color="k", ls=":", lw=1, label=f"noise ({noise_cl:.2f} mm)")
    ax.set(xlabel="strip mean H (mm)", ylabel="correlation length (mm)", title="per physical strip")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "thickness_strip_variation.png"), dpi=130)

    pd.set_option("display.width", 220)
    print(P.drop(columns=["n_overlap"]).round(3).to_string(index=False))
    print(f"\nmeasurement noise (per photo, from {len(dups)} duplicate pairs): SD {noise_sd:.3f} mm, "
          f"correlation length {noise_cl:.2f} mm; photo-to-photo offset "
          f"{100 * dups.offset_frac.abs().median():.1f}% median, {100 * dups.offset_frac.abs().max():.1f}% max")
    print()
    print(S.round(3).to_string(index=False))
    print()
    print(B.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
