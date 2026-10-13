"""
Stage 2 setup (local_estimation_handoff.md): assemble the virtual-field equations for the interior
patches and check how well they constrain the patch stiffness factors. No fitting.

Model: membrane equilibrium in the reference configuration, Div N = 0, with N = H t F^-T the first
Piola-Kirchhoff membrane tension (kPa*mm = mN/mm; constitutive.first_pk_membrane). Ladder level (a):
one factor alpha_k per patch scales Cgr, C1_el, C1_smc and C1_col of the full-domain bulk baseline
(dic_affine), with the exponents and gamma fixed, so in grid cell c

    N_c = sum_k phi_kc alpha_k N_base(F_c),   phi_kc = share of cell c's area inside patch k

and every equation below is linear in alpha.

Domain: the interior patches only (--exclude-ring outer rings removed), i.e. the model boundary is
moved inward past the tine-adjacent margin. Cells are integrated with their area share inside the
domain. Three blocks per frame:

  cut_z      for each full column of cells: sum_rows N_zZ dTheta = f_z F_z,  f_z = domain height / w_theta
  cut_theta  for each full row of cells:    sum_cols N_thetaTheta dZ = f_theta F_theta,  f_theta = width / w_z
             (the baseline's uniform-flux assumption, applied to the cuts through the domain)
  vf         interior virtual fields u* = sin(m pi xi) sin(n pi eta) e_i (m, n = 1..modes, i = z, theta),
             zero on the domain boundary: sum_cells A_c N_iJ du*_i/dX_J = 0. No width assumption.

F comes from the full-domain tracks with a (2 win + 1)-point window (--win 1 = 45 px by default; the synthetic
test showed the 75 px window of Stage 0 roughly doubles the spurious contrast of a uniform material).

Checks:
  * conditioning: singular values of each block and of all blocks together (rows scaled to unit norm,
    then each block to equal total weight; provisional until noise-based weights exist), and maps of the
    patch patterns the equations constrain least (smallest right singular vectors);
  * consistency: residuals of each block at alpha = 1 (bulk) and at the Stage 1 interior alpha map,
    relative to the size of the terms. Evaluations only; nothing is fitted.

    ../.venv/bin/python stage2_setup.py --test N1739/test1 [--exclude-ring 1] [--every 4] [--modes 8] [--win 1]
Outputs: Images/constitutive_fit/local/stage2_setup/
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import baseline as B  # noqa: E402
import constitutive as M  # noqa: E402
import dic2d as D  # noqa: E402
from track_tine_region import SUFFIX, TINE_TIPS  # noqa: E402

IMAGES = testcfg.get().images
LOCAL = testcfg.get().local_dir


def measured_positions(df, frames, mpp, with_valid=False):
    """Tracked positions of the full-domain tracks for the selected baseline frames, in the local layout
    (row 0 at the bottom tine tips) and specimen mm (z from the left tips, theta up from the bottom tips)."""
    lo_z, _ = TINE_TIPS["z"]
    _, hi_t = TINE_TIPS["theta"]
    cache, out, val = {}, [], []
    for k in frames:
        fd, t = df.folder.iloc[k], round(float(df.t.iloc[k]), 1)
        if fd not in cache:
            R = np.load(os.path.join(IMAGES, fd, "strain_results", fd + SUFFIX + ".npz"), allow_pickle=True)
            cache[fd] = (R, {tt: i for i, tt in enumerate(np.round(R["times"], 1))}, tuple(R["shape"]))
        R, kmap, shape = cache[fd]
        p = R["pos"][kmap[t]].reshape(shape + (2,))[::-1]
        out.append(np.stack([(p[..., 0] - lo_z) * mpp, (hi_t - p[..., 1]) * mpp], -1))
        val.append(R["valid"][kmap[t]].reshape(shape)[::-1])
    return (np.array(out), np.array(val)) if with_valid else np.array(out)


def window_F(Xg, x, valid, win):
    """F (specimen frame) on the grid from positions x (frames, ny, nx, 2) with a (2 win + 1)-point window."""
    F = np.full(x.shape[:3] + (2, 2), np.nan)
    for n in range(len(x)):
        Fi, ok = D.deformation_gradient_field(Xg.reshape(-1, 2), x[n].reshape(-1, 2), valid[n].ravel(),
                                              x.shape[1:3], win)
        Fi[~ok] = np.nan
        F[n] = Fi
    return F


def geometry(d, k_ring=1, modes=8):
    """Interior domain, patch area shares per cell, cuts and virtual-field gradients (see module docstring)."""
    W = d["patch_weights_mm2"].astype(float)                      # (nt, nz, ny, nx)
    nt, nz = W.shape[:2]
    ii, jj = np.meshgrid(np.arange(nt), np.arange(nz), indexing="ij")
    inner = (ii >= k_ring) & (ii < nt - k_ring) & (jj >= k_ring) & (jj < nz - k_ring)
    area = d["area_mm2"]
    phi = W[inner] / area                                          # (npat, ny, nx) area share per patch
    ez, et = d["cell_edges_z_mm"], d["cell_edges_theta_mm"]
    pz, pt = d["patch_edges_z_mm"], d["patch_edges_theta_mm"]
    Z0, Z1, T0, T1 = pz[k_ring], pz[nz - k_ring], pt[k_ring], pt[nt - k_ring]
    Lz, Lt = Z1 - Z0, T1 - T0
    zc, tc = d["z_mm"], d["theta_mm"]
    xi, eta = (zc - Z0) / Lz, (tc - T0) / Lt
    grads, names = [], []
    for m in range(1, modes + 1):
        for n in range(1, modes + 1):
            dZ = m * np.pi / Lz * np.cos(m * np.pi * xi) * np.sin(n * np.pi * eta)
            dT = n * np.pi / Lt * np.sin(m * np.pi * xi) * np.cos(n * np.pi * eta)
            for i, comp in enumerate(("z", "theta")):
                g = np.zeros((2, 2) + zc.shape)
                g[i, 0], g[i, 1] = dZ, dT
                grads.append(g); names.append(f"u*_{comp} m={m} n={n}")
    return dict(inner=inner, npat=int(inner.sum()), phi=phi, frac=phi.sum(0), area=area,
                dz=np.diff(ez), dt=np.diff(et), pz=pz, pt=pt, nt=nt, nz=nz, domain=(Z0, Z1, T0, T1),
                full_col=(ez[:-1] >= Z0 - 1e-9) & (ez[1:] <= Z1 + 1e-9),
                full_row=(et[:-1] >= T0 - 1e-9) & (et[1:] <= T1 + 1e-9),
                grads=np.array(grads), vf_names=names)


def assemble(g, F, Fz, Ft, P, H, w_theta, w_z):
    """Virtual-field equations for deformation gradients F (frames, ny, nx, 2, 2; specimen frame, NaN where
    unavailable) and forces Fz, Ft (mN). Returns (Gz, Gt, Gv, bz, bt): coefficient arrays
    (frames, rows, npat) and right-hand sides (frames, rows); the virtual-field right-hand side is zero."""
    Z0, Z1, T0, T1 = g["domain"]
    f_z, f_t = (T1 - T0) / w_theta, (Z1 - Z0) / w_z
    good = np.isfinite(F[..., 0, 0])
    Nb = M.first_pk_membrane(np.where(good[..., None, None], F, np.eye(2)), P, H)
    Nb = np.where(good[..., None, None], Nb, 0.0)
    phi = g["phi"]
    Gz = np.einsum("krc,frc,r->fck", phi, Nb[..., 0, 0], g["dt"])[:, g["full_col"]]
    Gt = np.einsum("krc,frc,c->frk", phi, Nb[..., 1, 1], g["dz"])[:, g["full_row"]]
    Gv = np.einsum("krc,frcij,vijrc,rc->fvk", phi, Nb, g["grads"], g["area"])
    bz = np.repeat((f_z * np.asarray(Fz))[:, None], Gz.shape[1], 1)
    bt = np.repeat((f_t * np.asarray(Ft))[:, None], Gt.shape[1], 1)
    return Gz, Gt, Gv, bz, bt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exclude-ring", type=int, default=1)
    ap.add_argument("--every", type=int, default=4, help="use every n-th loading frame")
    ap.add_argument("--modes", type=int, default=8,
                    help="sine modes per direction for the virtual fields (4 is too few to resolve 8 x 8 patches; "
                         "see stage2_synthetic.py)")
    ap.add_argument("--win", type=int, default=1,
                    help="strain window half-width in grid points: 1 = 45 px (default for Stage 2; lower bias in the "
                         "synthetic test), 2 = 75 px (as in Stage 0)")
    ap.add_argument("--out", default=os.path.join(LOCAL, "stage2_setup"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    d = np.load(os.path.join(LOCAL, "local_data.npz"))
    base, df = B.load(str(d["baseline"]))
    st = base["settings"]
    H, w_t, w_z = st["H_mm"], st["w_theta_mm"], st["w_z_mm"]
    P = np.array([base["inputs"]["dic_affine"]["params"][n] for n in M.PARAM_NAMES])

    # ---------------------------------------------------------------- geometry and assembly
    g = geometry(d, a.exclude_ring, a.modes)
    k_, inner, npat, nt, nz = a.exclude_ring, g["inner"], g["npat"], g["nt"], g["nz"]
    pz, pt, vf_names = g["pz"], g["pt"], g["vf_names"]
    Z0, Z1, T0, T1 = g["domain"]
    Lz, Lt = Z1 - Z0, T1 - T0
    f_z, f_t = Lt / w_t, Lz / w_z
    print(f"domain z {Z0:.3f}-{Z1:.3f} mm, theta {T0:.3f}-{T1:.3f} mm ({Lz:.2f} x {Lt:.2f} mm); {npat} patches; "
          f"force fractions f_z = {f_z:.3f}, f_theta = {f_t:.3f}; {g['full_col'].sum()} z-cuts, "
          f"{g['full_row'].sum()} theta-cuts")
    sets_all = d["set"]
    fidx = np.arange(0, len(sets_all), a.every)
    Xg = np.stack([d["z_mm"], d["theta_mm"]], -1)
    xm, vm = measured_positions(df.reset_index(drop=True), fidx, float(d["mm_per_px"]), with_valid=True)
    F = window_F(Xg, xm, vm, a.win)                                # (nf, ny, nx, 2, 2), specimen frame
    print(f"F from the tracks with a {(2 * a.win + 1) * 15} px window")
    missing = (~np.isfinite(F[..., 0, 0]) & (g["frac"] > 0)[None]).mean()
    Fz, Ft = d["force_z_mN"][fidx], d["force_theta_mN"][fidx]
    sets = sets_all[fidx]
    print(f"{len(fidx)} frames (every {a.every}); interior cells without F: {100 * missing:.3f}% (left out)")
    Gz, Gt, Gv, bz, bt = assemble(g, F, Fz, Ft, P, H, w_t, w_z)
    blocks = dict(cut_z=(Gz.reshape(-1, npat), bz.ravel()), cut_theta=(Gt.reshape(-1, npat), bt.ravel()),
                  vf=(Gv.reshape(-1, npat), np.zeros(Gv.shape[0] * Gv.shape[1])))

    # ---------------------------------------------------------------- conditioning
    def scaled(G):
        n = np.linalg.norm(G, axis=1, keepdims=True)
        keep = n[:, 0] > 0
        return G[keep] / n[keep] / np.sqrt(keep.sum())            # unit rows, equal total weight per block
    sv, Vt_all = {}, None
    for name, (G, _) in blocks.items():
        sv[name] = np.linalg.svd(scaled(G), compute_uv=False)
    Gall = np.vstack([scaled(G) for G, _ in blocks.values()])
    _, s_all, Vt_all = np.linalg.svd(Gall, full_matrices=False)
    sv["all"] = s_all
    for name, s in sv.items():
        print(f"  {name:9s}: {len(blocks[name][1]) if name != 'all' else Gall.shape[0]:6d} rows; "
              f"sigma_max/sigma_min = {s[0] / s[-1]:.3g}; directions with sigma > 1e-2 sigma_max: "
              f"{int((s > 1e-2 * s[0]).sum())}/{npat}, > 1e-3: {int((s > 1e-3 * s[0]).sum())}/{npat}")

    # ---------------------------------------------------------------- consistency at alpha = 1 and Stage 1 alpha
    s1 = pd.read_csv(os.path.join(LOCAL, f"stage1_interior{k_}" if k_ else "stage1", "stage1_patches.csv"))
    a1 = np.full((nt, nz), np.nan)
    a1[s1.i_theta, s1.j_z] = s1.alpha
    cand = {"alpha = 1 (bulk)": np.ones(npat), "Stage 1 alpha (interior)": a1[inner]}
    cons = {}
    for cname, al in cand.items():
        cons[cname] = {}
        for name, (G, b) in blocks.items():
            r = G @ al - b
            if name == "vf":                                      # relative to the size of the terms
                scale = np.sqrt(np.mean((np.abs(G) @ np.abs(al)) ** 2))
            else:
                scale = np.sqrt(np.mean(b ** 2))
            cons[cname][name] = float(np.sqrt(np.mean(r ** 2)) / scale)
        print(f"  {cname:26s}: relative RMS residual  " + "  ".join(f"{k}: {v:.3f}" for k, v in cons[cname].items()))
    # per-protocol view of the cut blocks at alpha = 1
    per = {}
    for name, (G, b) in (("cut_z", blocks["cut_z"]), ("cut_theta", blocks["cut_theta"])):
        ncut = Gz.shape[1] if name == "cut_z" else Gt.shape[1]
        rs = np.repeat(sets, ncut)
        per[name] = {int(s): float(np.sqrt(np.mean((G[rs == s].sum(1) - b[rs == s]) ** 2)) /
                                   np.sqrt(np.mean(b[rs == s] ** 2))) for s in np.unique(sets)}

    # ---------------------------------------------------------------- save
    np.savez_compressed(os.path.join(a.out, "stage2_setup.npz"), G_cut_z=Gz, G_cut_theta=Gt, G_vf=Gv,
                        b_cut_z=bz, b_cut_theta=bt, frames=fidx, sets=sets, inner=inner,
                        singular_values_all=s_all, right_singular_vectors=Vt_all, vf_names=np.array(vf_names),
                        domain_mm=np.array([Z0, Z1, T0, T1]), force_fractions=np.array([f_z, f_t]),
                        params_baseline=P, H_mm=H)
    with open(os.path.join(a.out, "stage2_setup.json"), "w") as f:
        json.dump(dict(description="Stage 2 virtual-field equations for the interior patches (ladder (a), exponents "
                       "and gamma fixed at the full-domain baseline); conditioning and consistency checks, no fit",
                       domain_mm=dict(z=[Z0, Z1], theta=[T0, T1]), n_patches=npat, frames_used=int(len(fidx)),
                       every=a.every, modes=a.modes, window_px=(2 * a.win + 1) * 15, force_fractions=dict(f_z=f_z, f_theta=f_t),
                       rows=dict(cut_z=int(Gz.size // npat), cut_theta=int(Gt.size // npat), vf=int(Gv.size // npat)),
                       condition_numbers={k: float(v[0] / v[-1]) for k, v in sv.items()},
                       relative_residuals=cons, cut_residual_alpha1_by_protocol=per), f, indent=2)

    # ---------------------------------------------------------------- figure
    fig = plt.figure(figsize=(20, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 5)
    a_ = fig.add_subplot(gs[0, :2])
    for name, s in sv.items():
        a_.semilogy(np.arange(1, len(s) + 1), s / s[0], "o-", ms=3, label=f"{name} (κ = {s[0] / s[-1]:.2g})")
    a_.set_xlabel("singular value index"); a_.set_ylabel("σ / σ_max"); a_.grid(alpha=.3, which="both")
    a_.legend(fontsize=9); a_.set_title(f"Singular values of the scaled equations ({npat} patch factors)", fontsize=10)
    a_ = fig.add_subplot(gs[1, :2]); a_.axis("off")
    txt = ["Relative RMS residual (no fit):"]
    for cname, v in cons.items():
        txt.append(f"  {cname}:  " + ",  ".join(f"{k} {x:.3f}" for k, x in v.items()))
    txt.append("")
    txt.append("cut residual at α = 1 by protocol (z / θ):")
    for s_ in per["cut_z"]:
        txt.append(f"  set {s_:2d}  {per['cut_z'][s_]:.3f} / {per['cut_theta'][s_]:.3f}")
    a_.text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=9)
    show = [(-1, "weakest"), (-2, "2nd weakest"), (-3, "3rd weakest"), (0, "strongest"), (1, "2nd strongest"),
            (2, "3rd strongest")]
    for n_, (idx, lab) in enumerate(show):
        a_ = fig.add_subplot(gs[n_ // 3, 2 + n_ % 3])
        v = np.full((nt, nz), np.nan)
        v[inner] = Vt_all[idx]
        mval = np.nanmax(np.abs(v))
        a_.set_facecolor("0.85")
        m = a_.pcolormesh(pz, pt, v, cmap="PuOr_r", vmin=-mval, vmax=mval)
        a_.set_xlim(0, pz[-1]); a_.set_ylim(0, pt[-1]); a_.set_aspect("equal")
        a_.set_title(f"{lab} direction (σ/σ_max = {s_all[idx] / s_all[0]:.2g})", fontsize=10)
        a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
        fig.colorbar(m, ax=a_, shrink=.8)
    fig.suptitle(f"Stage 2 setup: virtual-field equations for the {npat} interior patches (no fit); "
                 f"{len(fidx)} frames, {a.modes}×{a.modes} sine modes per direction")
    fig.savefig(os.path.join(a.out, "stage2_setup.png"), dpi=85)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
