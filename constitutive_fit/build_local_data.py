"""
Stage 0 of the local estimation (local_estimation_handoff.md): assemble the local data set
and the band diagnostics. No constitutive fitting.

Frames, force tare and the global dic_affine stretches come from a frozen baseline data set
(baseline.load(), checksum-verified; default the tine_region baseline, whose DIC domain is
the same as the local tracks; --baseline standard for the original ROI), so every frame here is a baseline frame
(Def11-Def19 loading branches, 1489 frames). Local F comes from the tine-region tracks
(track_tine_region.py), which share the baseline's common reference (first frame of the test's reference_deformation)
and reach 5-15 px from the tine tips.

Zones fill the region between the tine tips: each grid point owns the cell bounded by the
midpoints to its neighbors, and the outer cells extend to the tip lines. The zones are a grid
of patches (default 10 x 10 in z x theta, 0.85 x 1.00 mm; set with --patches N_Z N_THETA, e.g. 6 x 7 for
1.42 mm squares, slightly larger than the 75 px strain
window). Cells that straddle a patch edge are split by exact area, so the patches tile the
region with no gaps or overlap.

Axes: theta increases upward in the image (dic2d.THETA_SIGN). Stored grids are [row, col] =
[theta, z] with row 0 at the bottom tine tips and col 0 at the left tine tips; F is in the
specimen frame (dic2d.to_specimen).

    ../.venv/bin/python build_local_data.py --test N1739/test1
Outputs: Images/constitutive_fit/local/
    local_data.npz          F per frame on the grid, masks, tared forces, zones, checksums
    patches.json            patch edges (mm from the left / bottom tips), coverage, strain concentrations
    stage0_peak_maps.png    lambda_z, lambda_theta, shear angle, J_2 at each protocol's peak
    stage0_band_profiles.png  means along z vs theta (and along theta vs z, the control)
    stage0_patches.png      patch grid, patch-mean strain concentration, coverage, within-patch spread
    stage0_patch_layout.png   patch configuration on the reference image, and a zoom on the grid cells
"""
import argparse
import hashlib
import json
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import baseline as B  # noqa: E402
import dic2d as D  # noqa: E402
from track_tine_region import SUFFIX, TINE_TIPS  # noqa: E402

IMAGES = testcfg.get().images
OUT = os.path.join(IMAGES, "constitutive_fit", "local")
MAP_SETS = (2, 3, 4, 5, 7, 8, 9)            # peak equibiaxial (Def11) + each ratio protocol
PATCHES = (10, 10)                          # default (z, theta) patch counts: 0.85 x 1.00 mm; --patches to change


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cell_edges(centers, lo, hi):
    """Cell boundaries for lattice points: midpoints, outer cells extended to [lo, hi]."""
    mid = 0.5 * (centers[1:] + centers[:-1])
    return np.r_[lo, mid, hi]


def kin(F):
    """lambda_z, lambda_theta, shear angle (deg, change of the initially right angle z-theta), J_2."""
    C = np.einsum("...ki,...kj->...ij", F, F)
    lz, lt = np.sqrt(C[..., 0, 0]), np.sqrt(C[..., 1, 1])
    shear = np.degrees(np.arcsin(np.clip(C[..., 0, 1] / (lz * lt), -1, 1)))
    return dict(lambda_z=lz, lambda_theta=lt, shear_deg=shear, J2=np.linalg.det(F))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--baseline", choices=tuple(B.TRACKS), default="tine_region",
                    help="frozen baseline supplying frames, tare and global stretches")
    ap.add_argument("--patches", nargs=2, type=int, default=PATCHES, metavar=("N_Z", "N_THETA"),
                    help="patch grid over the tine region, counts along z and theta (default 10 x 10, "
                         "0.85 x 1.00 mm; 6 x 7 gives 1.42 mm squares)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    np.seterr(invalid="ignore", divide="ignore")
    warnings.filterwarnings("ignore", "Mean of empty slice")
    mpp = D.load_scale()["mm_per_px"]

    base, df = B.load(a.baseline)
    print(f"baseline: {a.baseline} ({B.TRACKS[a.baseline][1]})")
    df = df.reset_index(drop=True)
    df["t"] = df.t.round(1)
    folders = list(dict.fromkeys(df.folder))
    ref_label = f"reference image {folders[0].split('_')[0]}, t = {df.t.iloc[0]:.1f} s"
    lo_z, hi_z = TINE_TIPS["z"]
    lo_t, hi_t = TINE_TIPS["theta"]

    # ---------------------------------------------------------------- local F per baseline frame
    F_all, ok_all, vt_all, src = [], [], [], {}
    for fd in folders:
        path = os.path.join(IMAGES, fd, "strain_results", fd + SUFFIX + ".npz")
        R = np.load(path, allow_pickle=True)
        src[os.path.relpath(path, IMAGES)] = _sha(path)
        X0, shape, times = R["X0"], tuple(R["shape"]), np.round(R["times"], 1)
        win = D.window_points(float(R["step"]))
        kmap = {t: k for k, t in enumerate(times)}
        for t in df.t[df.folder == fd]:
            k = kmap[t]
            F, ok = D.deformation_gradient_field(X0, R["pos"][k], R["valid"][k], shape, win)
            # specimen frame (theta up); rows reordered so row 0 is at the bottom tine tips
            F_all.append(D.to_specimen(F)[::-1].astype(np.float32)); ok_all.append(ok[::-1])
            vt_all.append(R["valid"][k].reshape(shape)[::-1])
        print(f"  {fd}: {int((df.folder == fd).sum())} loading frames", flush=True)
    F_all, ok_all, vt_all = np.array(F_all), np.array(ok_all), np.array(vt_all)
    ny, nx = shape
    Ximg, Yimg = X0[:, 0].reshape(shape)[::-1], X0[:, 1].reshape(shape)[::-1]   # image px, rows bottom-up
    xs, ys = Ximg[0], Yimg[:, 0]
    zc, tc = xs - lo_z, hi_t - ys          # px from the left / bottom tine tips (theta increases upward)
    ex, ey = cell_edges(zc, 0, hi_z - lo_z), cell_edges(tc, 0, hi_t - lo_t)
    area = np.outer(np.diff(ey), np.diff(ex))                    # tributary area (px^2)
    full_win = np.zeros(shape, bool)
    full_win[win:-win, win:-win] = True                          # strain window not truncated
    K = kin(F_all.astype(np.float64))
    print(f"local F: {F_all.shape[0]} frames, grid {nx} x {ny}; F available in "
          f"{ok_all.mean()*100:.1f}% of point-frames (min per frame {ok_all.mean((1, 2)).min()*100:.1f}%)")

    # ---------------------------------------------------------------- band profiles
    lz_g, lt_g = df.dic_affine_z.to_numpy(), df.dic_affine_t.to_numpy()
    Ezz, Ett = 0.5 * (K["lambda_z"] ** 2 - 1), 0.5 * (K["lambda_theta"] ** 2 - 1)
    Ezz_g, Ett_g = 0.5 * (lz_g ** 2 - 1), 0.5 * (lt_g ** 2 - 1)
    # strain concentration: local / global Green strain, high-load frames of each protocol
    hi_t_ = np.zeros(len(df), bool); hi_z_ = np.zeros(len(df), bool)
    for s, idx in df.groupby("set").indices.items():
        hi_t_[idx] = Ett_g[idx] > 0.5 * Ett_g[idx].max()
        hi_z_[idx] = Ezz_g[idx] > 0.5 * Ezz_g[idx].max()
    hi_t_ &= Ett_g > 0.02                                       # skip near-zero theta strain
    hi_z_ &= Ezz_g > 0.02
    conc_t = np.nanmean(Ett[hi_t_] / Ett_g[hi_t_, None, None], 0)
    conc_z = np.nanmean(Ezz[hi_z_] / Ezz_g[hi_z_, None, None], 0)
    ok_any = ok_all.mean(0) > 0.5
    conc_t[~ok_any] = np.nan; conc_z[~ok_any] = np.nan

    def zmean(Q):                       # area-weighted mean along z (per theta row)
        w = np.where(np.isfinite(Q), area, 0)
        return np.nansum(Q * w, -1) / w.sum(-1)

    def tmean(Q):
        w = np.where(np.isfinite(Q), area, 0)
        return np.nansum(Q * w, -2) / w.sum(-2)

    # ---------------------------------------------------------------- patches
    Wz_px, Wt_px = hi_z - lo_z, hi_t - lo_t
    nz_p, nt_p = a.patches
    pz, pt = np.linspace(0, Wz_px, nz_p + 1), np.linspace(0, Wt_px, nt_p + 1)     # patch edges (px from tips)

    def overlap(edges, pe):            # (n_patches, n_cells): overlap length of each cell with each patch (px)
        return np.clip(np.minimum(edges[None, 1:], pe[1:, None]) - np.maximum(edges[None, :-1], pe[:-1, None]),
                       0, None)
    # PW[i, j, row, col]: area of grid cell (row, col) inside patch (i = theta, j = z), px^2. Cells that
    # straddle a patch edge are split exactly, so the patches tile the tine region with no gaps or overlap.
    PW = np.einsum("ar,bc->abrc", overlap(ey, pt), overlap(ex, pz))
    assert np.allclose(PW.sum((0, 1)), area)

    def patch_mean(Q):                 # Q (..., ny, nx) -> (..., nt_p, nz_p), area-weighted over finite values
        fin = np.isfinite(Q)
        num = np.einsum("abrc,...rc->...ab", PW, np.where(fin, Q, 0))
        return num / np.einsum("abrc,...rc->...ab", PW, fin.astype(float))
    coverage = np.einsum("abrc,rc->ab", PW, ok_all.mean(0)) / PW.sum((2, 3))
    pc_t, pc_z = patch_mean(conc_t), patch_mean(conc_z)
    sd_t = np.sqrt(np.maximum(patch_mean(conc_t ** 2) - pc_t ** 2, 0))
    sd_z = np.sqrt(np.maximum(patch_mean(conc_z ** 2) - pc_z ** 2, 0))
    side_z, side_t = Wz_px / nz_p * mpp, Wt_px / nt_p * mpp
    print(f"patches: {nz_p} (z) x {nt_p} (theta) = {nz_p * nt_p}, {side_z:.3f} x {side_t:.3f} mm "
          f"(aspect {side_z / side_t:.3f}); F coverage per patch min {coverage.min() * 100:.1f}%")
    win_mm = (2 * win + 1) * float(R["step"]) * mpp
    if min(side_z, side_t) < win_mm:
        print(f"  note: patches are smaller than the {win_mm:.2f} mm strain window, so neighboring patch "
              "means share displacement data")
    print(f"  patch E_thetatheta concentration {np.nanmin(pc_t):.2f}-{np.nanmax(pc_t):.2f}, "
          f"E_zz concentration {np.nanmin(pc_z):.2f}-{np.nanmax(pc_z):.2f}")

    # ---------------------------------------------------------------- save data set
    script_sha = {s: _sha(testcfg.script_path(s)) for s in ("build_local_data.py", "track_tine_region.py", "dic2d.py")}
    np.savez_compressed(
        os.path.join(a.out, "local_data.npz"),
        F=F_all, ok=ok_all, valid_track=vt_all, full_window=full_win,
        X_image_px=Ximg, Y_image_px=Yimg, shape=np.array(shape),
        z_mm=np.broadcast_to(zc * mpp, shape), theta_mm=np.broadcast_to((tc * mpp)[:, None], shape),
        cell_edges_z_mm=ex * mpp, cell_edges_theta_mm=ey * mpp, area_px2=area, area_mm2=area * mpp ** 2,
        grid_layout="[row, col] = [theta, z]; row 0 at the bottom tine tips, col 0 at the left tine tips; "
                    "F in the specimen frame (index 0 = z, 1 = theta, theta increasing upward)",
        patch_weights_mm2=(PW * mpp ** 2).astype(np.float32),
        patch_edges_z_mm=pz * mpp, patch_edges_theta_mm=pt * mpp,
        patch_layout="patch (i, j): i = theta row from the bottom tine tips, j = z column from the left tine tips; "
                     "patch_weights_mm2[i, j, row, col] = area of grid cell (row, col) inside patch (i, j)",
        patch_coverage=coverage, patch_conc_theta=pc_t, patch_conc_z=pc_z,
        tine_tips_z_px=np.array(TINE_TIPS["z"]), tine_tips_theta_px=np.array(TINE_TIPS["theta"]),
        mm_per_px=mpp, window_px=(2 * win + 1) * float(R["step"]),
        t=df.t.to_numpy(), set=df.set.to_numpy(), folder=df.folder.to_numpy().astype(str),
        force_z_mN=df.fz.to_numpy(), force_theta_mN=df.ft.to_numpy(),
        lambda_z_dic_affine=lz_g, lambda_theta_dic_affine=lt_g,
        conc_theta=conc_t, conc_z=conc_z,
        axis_convention=D.AXIS_CONVENTION,
        baseline=a.baseline, source_checksums=json.dumps(src), baseline_data_sha=base["checksums"]["data"]["baseline_data.csv"],
        script_checksums=json.dumps(script_sha))
    with open(os.path.join(a.out, "patches.json"), "w") as f:
        json.dump(dict(description="Rectangular patches tiling the region between the tine tips",
                       layout="patch (i, j): i = theta row from the bottom tine tips (theta increasing upward), "
                              "j = z column from the left tine tips",
                       n_z=nz_p, n_theta=nt_p, side_z_mm=side_z, side_theta_mm=side_t,
                       edges_z_mm=[round(v, 4) for v in pz * mpp], edges_theta_mm=[round(v, 4) for v in pt * mpp],
                       tine_tips_px=TINE_TIPS, mm_per_px=mpp, strain_window_mm=(2 * win + 1) * float(R["step"]) * mpp,
                       coverage=np.round(coverage, 4).tolist(),
                       conc_theta=np.round(pc_t, 3).tolist(), conc_z=np.round(pc_z, 3).tolist()), f, indent=2)

    # ---------------------------------------------------------------- figures
    Xmm, Ymm = ex * mpp, ey * mpp                             # cell edges, mm from the left / bottom tips
    ycen, xcen = tc * mpp, zc * mpp
    Wz, Wt = (hi_z - lo_z) * mpp, (hi_t - lo_t) * mpp

    def frame_axes(ax_):
        ax_.set_xlim(0, Wz); ax_.set_ylim(0, Wt); ax_.set_aspect("equal")

    peak = {s: idx[np.argmax((lz_g[idx] - 1) + (lt_g[idx] - 1))] for s, idx in df.groupby("set").indices.items()}
    qs = [("lambda_z", r"$\lambda_z$", "viridis"), ("lambda_theta", r"$\lambda_\theta$", "viridis"),
          ("shear_deg", r"shear angle $\gamma_{z\theta}$ (°)", "RdBu_r"), ("J2", r"$J_2$ = det F", "viridis")]
    fig, ax = plt.subplots(len(MAP_SETS), 4, figsize=(17, 3.6 * len(MAP_SETS)), constrained_layout=True,
                           squeeze=False)
    lims = {}                                                    # common color limits per column
    for q, _, _ in qs:
        allv = np.concatenate([np.where(ok_all[peak[s]], K[q][peak[s]], np.nan).ravel() for s in MAP_SETS])
        lo, hi = np.nanpercentile(allv, [1, 99])
        if q == "shear_deg":
            hi = max(abs(lo), abs(hi)); lo = -hi
        lims[q] = (lo, hi)
    for r_, s in enumerate(MAP_SETS):
        k = peak[s]
        for c_, (q, lab, cmap) in enumerate(qs):
            v = np.where(ok_all[k], K[q][k], np.nan)
            m = ax[r_, c_].pcolormesh(Xmm, Ymm, v, cmap=cmap, vmin=lims[q][0], vmax=lims[q][1])
            frame_axes(ax[r_, c_]); fig.colorbar(m, ax=ax[r_, c_], shrink=.85)
            ax[r_, c_].set_title(f"{df.folder[k]}  {lab}\nglobal $\\lambda_z$ {lz_g[k]:.3f}, "
                                 f"$\\lambda_\\theta$ {lt_g[k]:.3f}", fontsize=9)
            if r_ == len(MAP_SETS) - 1:
                ax[r_, c_].set_xlabel("z from left tine tips (mm)")
            if c_ == 0:
                ax[r_, c_].set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    fig.suptitle("Local kinematics at each protocol's peak (reference configuration, region between the tine tips; "
                 "color limits common to each column)")
    fig.savefig(os.path.join(a.out, "stage0_peak_maps.png"), dpi=80)

    fig, ax = plt.subplots(2, 4, figsize=(18, 8.5), constrained_layout=True)
    cols = {s: plt.get_cmap("tab10")(i) for i, s in enumerate(MAP_SETS)}
    for c_, (q, lab, _) in enumerate(qs):
        for s in MAP_SETS:
            k = peak[s]
            v = np.where(ok_all[k], K[q][k], np.nan)
            ax[0, c_].plot(ycen, zmean(v), color=cols[s], lw=1.4, label=df.folder[k].split("_", 1)[0])
            ax[1, c_].plot(xcen, tmean(v), color=cols[s], lw=1.4)
            if q in ("lambda_z", "lambda_theta"):                # machine (rake) stretch behind the profiles
                g = df["rake_z" if q == "lambda_z" else "rake_t"][k]
                for r_ in (0, 1):
                    ax[r_, c_].axhline(g, color="0.6", ls="--", lw=1, zorder=0)
        ax[0, c_].set_title(f"{lab}: mean along z vs $\\theta$")
        ax[1, c_].set_title(f"{lab}: mean along $\\theta$ vs z (control)")
        ax[0, c_].set_xlabel(r"$\theta$ from bottom tine tips (mm)"); ax[1, c_].set_xlabel("z from left tine tips (mm)")
        ax[0, c_].set_xlim(0, Wt); ax[1, c_].set_xlim(0, Wz)
        ax[0, c_].grid(alpha=.3); ax[1, c_].grid(alpha=.3)
    h, _ = ax[0, 0].get_legend_handles_labels()
    fig.legend(h + [plt.Line2D([], [], color="0.6", ls="--", lw=1)],
               [l_.get_label() for l_ in h] + ["protocol machine stretch\n(rake separation)"],
               loc="outside right upper", fontsize=9)
    fig.suptitle("Band profiles at each protocol's peak, region between the tine tips "
                 "(dashed gray: each protocol's machine stretch at the same frame)")
    fig.savefig(os.path.join(a.out, "stage0_band_profiles.png"), dpi=85)

    # ---------------------------------------------------------------- patches
    import cv2
    ref_img = cv2.cvtColor(cv2.imread(D.find_frames([os.path.join(IMAGES, folders[0])])[1][0]), cv2.COLOR_BGR2RGB)
    PZ, PT = pz * mpp, pt * mpp
    pzc, ptc = 0.5 * (PZ[1:] + PZ[:-1]), 0.5 * (PT[1:] + PT[:-1])
    hist_t, hist_z = patch_mean(K["lambda_theta"]), patch_mean(K["lambda_z"])     # (frames, nt_p, nz_p)
    fs = 8 if max(nz_p, nt_p) <= 7 else 6                         # label size in the patch maps
    shape_word = "square patches" if abs(side_z / side_t - 1) < 0.02 else "patches"

    def grid_lines(ax_, color="k"):
        for v in PZ:
            ax_.axvline(v, color=color, lw=.8)
        for v in PT:
            ax_.axhline(v, color=color, lw=.8)

    def patch_map(ax_, V, title, cmap, lim, fmt):
        m = ax_.pcolormesh(PZ, PT, V, cmap=cmap, vmin=lim[0], vmax=lim[1])
        frame_axes(ax_); grid_lines(ax_); fig.colorbar(m, ax=ax_, shrink=.8)
        for i in range(nt_p):
            for j in range(nz_p):
                r, g, b_, _ = m.cmap(m.norm(V[i, j]))
                ax_.text(pzc[j], ptc[i], fmt.format(V[i, j]), ha="center", va="center", fontsize=fs,
                         color="w" if 0.299 * r + 0.587 * g + 0.114 * b_ < 0.45 else "k")
        ax_.set_title(title, fontsize=10); ax_.set_xlabel("z from left tine tips (mm)")

    fig, ax = plt.subplots(2, 4, figsize=(22, 11), constrained_layout=True)
    a_ = ax[0, 0]
    a_.imshow(ref_img[lo_t:hi_t, lo_z:hi_z], extent=(0, Wz, 0, Wt))
    frame_axes(a_); grid_lines(a_, "w")
    for i in range(nt_p):
        for j in range(nz_p):
            a_.text(pzc[j], ptc[i], f"{i},{j}", ha="center", va="center", fontsize=fs - 1, color="w", weight="bold")
    a_.set_title(f"{nz_p} × {nt_p} {shape_word}, {side_z:.2f} × {side_t:.2f} mm\n"
                 f"(label: θ row, z column; {ref_label})", fontsize=10)
    a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    for a_, V, lab in ((ax[0, 1], pc_t, r"$E_{\theta\theta}$ / global $E_{\theta\theta}$"),
                       (ax[0, 2], pc_z, r"$E_{zz}$ / global $E_{zz}$")):
        d = np.nanmax(np.abs(V - 1))
        patch_map(a_, V, f"patch-mean strain concentration {lab}\n(high-load frames, Def11–19)", "RdBu_r",
                  (1 - d, 1 + d), "{:.2f}")
    patch_map(ax[0, 3], 100 * coverage, "F available (% of patch area, all frames)", "Greens",
              (min(90, np.nanmin(100 * coverage)), 100), "{:.0f}")
    for a_, V, lab in ((ax[1, 0], sd_t, r"$E_{\theta\theta}$"), (ax[1, 1], sd_z, r"$E_{zz}$")):
        patch_map(a_, V, f"spread inside each patch: sd of the {lab} concentration", "magma_r",
                  (0, np.nanmax(V)), "{:.2f}")
    ax[1, 0].set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    for a_, H, g, lab, cm, idx in ((ax[1, 2], hist_t, lt_g, r"\lambda_\theta", "viridis", "row"),
                                   (ax[1, 3], hist_z, lz_g, r"\lambda_z", "plasma", "col")):
        n = nt_p if idx == "row" else nz_p
        for i in range(nt_p):
            for j in range(nz_p):
                c = plt.get_cmap(cm)((i if idx == "row" else j) / max(n - 1, 1))
                a_.plot(g, H[:, i, j], ".", ms=.6, color=c)
        a_.plot([1, g.max()], [1, g.max()], "k", lw=.8)
        sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, n - 1))
        fig.colorbar(sm, ax=a_, shrink=.8, label=r"patch $\theta$ row" if idx == "row" else "patch z column")
        a_.set_xlabel(f"global ${lab}$ (dic_affine)"); a_.set_ylabel(f"patch-mean ${lab}$")
        a_.set_title(f"patch-mean ${lab}$, all {nz_p * nt_p} patches, 1489 frames", fontsize=10); a_.grid(alpha=.3)
    fig.suptitle(f"{shape_word.capitalize()} between the tine tips ({nz_p} × {nt_p}); "
                 f"strain window {(2 * win + 1) * float(R['step']) * mpp:.2f} mm")
    fig.savefig(os.path.join(a.out, "stage0_patches.png"), dpi=85)

    # ---------------------------------------------------------------- patch layout
    pad = 110                                                    # px of image shown beyond the tine tips
    fig, ax = plt.subplots(1, 2, figsize=(19, 10), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.35, 1]))
    a_ = ax[0]
    a_.imshow(ref_img[lo_t - pad:hi_t + pad, lo_z - pad:hi_z + pad],
              extent=(-pad * mpp, Wz + pad * mpp, -pad * mpp, Wt + pad * mpp))
    for i in range(nt_p):
        for j in range(nz_p):
            a_.add_patch(Rectangle((PZ[j], PT[i]), PZ[j + 1] - PZ[j], PT[i + 1] - PT[i],
                                   fc="none", ec="w", lw=1.6))
            a_.text(pzc[j], ptc[i], f"({i},{j})", ha="center", va="center", fontsize=fs + 2, color="w", weight="bold",
                    bbox=dict(fc="k", ec="none", alpha=.35, pad=1.5))
    a_.add_patch(Rectangle((0, 0), Wz, Wt, fc="none", ec="r", lw=2.2, label="tine tips (region boundary)"))
    a_.plot([], [], color="w", lw=1.6, label=f"patch, {side_z:.2f} × {side_t:.2f} mm")
    a_.set_xlim(-pad * mpp, Wz + pad * mpp); a_.set_ylim(-pad * mpp, Wt + pad * mpp); a_.set_aspect("equal")
    a_.annotate("", xy=(Wz * .22, -pad * mpp * .55), xytext=(0, -pad * mpp * .55),
                arrowprops=dict(arrowstyle="->", color="w", lw=2))
    a_.text(Wz * .23, -pad * mpp * .55, "z (axial)", color="w", va="center", fontsize=11, weight="bold")
    a_.annotate("", xy=(-pad * mpp * .55, Wt * .2), xytext=(-pad * mpp * .55, 0),
                arrowprops=dict(arrowstyle="->", color="w", lw=2))
    a_.text(-pad * mpp * .55, Wt * .21, r"$\theta$ (circ.)", color="w", ha="center", va="bottom", fontsize=11,
            weight="bold", rotation=90)
    a_.legend(loc="upper right", fontsize=9, facecolor="0.25", labelcolor="w", framealpha=.85)
    a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    a_.set_title(f"{nz_p} × {nt_p} = {nz_p * nt_p} {shape_word} tiling the region between the tine tips "
                 f"({Wz:.2f} × {Wt:.2f} mm)\nlabel (i, j): i = θ row from the bottom, j = z column from the left; "
                 + ref_label, fontsize=11)
    # zoom: bottom-left 2 x 2 patches with the DIC grid points and their cells
    a_ = ax[1]
    zx, zy = PZ[2], PT[2]
    a_.imshow(ref_img[lo_t:hi_t, lo_z:hi_z], extent=(0, Wz, 0, Wt), alpha=.55)
    for e in Xmm:
        a_.plot([e, e], [0, Wt], color="0.35", lw=.5)
    for e in Ymm:
        a_.plot([0, Wz], [e, e], color="0.35", lw=.5)
    a_.plot(*np.meshgrid(xcen, ycen), ".", color="C0", ms=5)
    for i in range(nt_p):
        for j in range(nz_p):
            a_.add_patch(Rectangle((PZ[j], PT[i]), PZ[j + 1] - PZ[j], PT[i + 1] - PT[i], fc="none", ec="r", lw=2))
            if PZ[j] < zx and PT[i] < zy:
                a_.text(pzc[j], ptc[i], f"({i},{j})", ha="center", va="center", fontsize=12, color="r",
                        weight="bold", bbox=dict(fc="w", ec="none", alpha=.7))
    a_.set_xlim(-.15, zx + .15); a_.set_ylim(-.15, zy + .15); a_.set_aspect("equal")
    a_.plot([], [], ".", color="C0", ms=6, label=f"DIC grid point (15 px = {15 * mpp:.3f} mm spacing)")
    a_.plot([], [], color="0.35", lw=.8, label="grid cell (tributary area of a point)")
    a_.plot([], [], color="r", lw=2, label="patch boundary")
    a_.legend(loc="upper right", fontsize=9, framealpha=.9)
    a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    a_.set_title(f"Zoom on the bottom-left corner: a patch spans {side_z / (15 * mpp):.1f} × "
                 f"{side_t / (15 * mpp):.1f} grid spacings;\n"
                 "cells cut by a patch boundary are split by area; outer cells extend to the tine tips",
                 fontsize=11)
    fig.savefig(os.path.join(a.out, "stage0_patch_layout.png"), dpi=90)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
