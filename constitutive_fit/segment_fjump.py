"""
Computational dissection of the tine-region domain by deformation-gradient jump
(Witzenburg, Dhume, Lake & Barocas 2016, IEEE Trans Med Imaging 35:29-41).

Mesh: bilinear quadrilateral elements on the tracked material points of the full-domain tracks
(track_tine_region.py). Each element spans 2 x 2 grid cells (30 px, ~0.48 mm); the last column of
elements spans the one remaining cell. F is evaluated at each element-edge midpoint from the
element's own bilinear interpolation of the tracked positions (their Eqs. 1-2), and the jump across
a shared edge is the double contraction of the difference with itself, ||F_1 - F_2||^2 (Eq. 3).
A bilinear field is continuous along a shared edge, so the jump comes from the normal derivative.

Jumps: averaged over the last 5 loading frames of each of Def11-Def19 (the same loading branches as
the bulk baselines), elements with any invalid node skipped.

Component scaling (default; departs from the paper's Eq. 3): before summing over the four components
of F, each component's squared jump is divided by its mean square over the interior links (two elements
in from the domain edge), per protocol. With equal weights, the larger axial stretches make F_zz
dominate the sum and the circumferential structure is missed; scaling makes every component count
equally. --unscaled restores the paper's ||F_1 - F_2||^2.

Each protocol's jump field is then divided by its 99th percentile (the paper normalizes per extension;
the percentile keeps one noisy edge from setting the scale) and the nine are summed.

Segmentation: the element mesh becomes an unweighted network (node = element, link = shared edge).
Edge betweenness (number of shortest paths through a link) is multiplied by the summed jump and the
link with the largest product is removed; betweenness is recomputed and the step repeated. Each time
the network splits, the partition and its modularity (Newman) on the original network are recorded.
As in the paper, modularity is reported to help choose the number of subdomains; it does not give a
consistent stopping value, so the choice is left to the user.

Image axes: the segmentation is invariant to the sign of theta. Plots use theta increasing upward,
measured from the bottom tine tips.

    ../.venv/bin/python segment_fjump.py --test N1739/test1 [--max-parts 10] [--frames 5] [--unscaled] [--uniform-jump]
Outputs: <test>/constitutive_fit/local/segmentation/
    segmentation.npz      element geometry, summed and per-protocol jumps, labels for every partition
    segmentation.json     partitions (sizes, modularity, removal step) and settings
    fjump_segmentation.png  summed jump map, modularity, partitions for 2..max-parts subdomains
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components, shortest_path  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import testcfg  # noqa: E402
import baseline as B  # noqa: E402
import dic2d as D  # noqa: E402
from track_tine_region import SUFFIX, TINE_TIPS  # noqa: E402

OUT = os.path.join(testcfg.get().local_dir, "segmentation")


def node_lines(n, stride):
    """Grid indices used as element corners: every `stride`-th line, plus the last one."""
    idx = list(range(0, n, stride))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return np.array(idx)


def edge_F(Xn, xn):
    """
    F at the four edge midpoints of bilinear elements with axis-aligned rectangular reference shape.
    Xn, xn: (..., ne_r, ne_c, 4, 2) reference / current corners ordered (bl, br, tr, tl) in grid index
    (row up = +1). Returns dict side -> F (..., ne_r, ne_c, 2, 2) for sides 'bottom', 'right', 'top', 'left'.
    """
    hx = Xn[..., 1, 0] - Xn[..., 0, 0]                 # element width (reference, image x)
    hy = Xn[..., 3, 1] - Xn[..., 0, 1]                 # element height (reference, image y; sign kept)
    x0, x1, x2, x3 = (xn[..., k, :] for k in range(4))
    out = {}
    # local coords s, t in [0, 1]; x = (1-s)(1-t) x0 + s(1-t) x1 + s t x2 + (1-s) t x3
    for side, (s, t) in dict(bottom=(.5, 0.), right=(1., .5), top=(.5, 1.), left=(0., .5)).items():
        dxds = (1 - t) * (x1 - x0) + t * (x2 - x3)
        dxdt = (1 - s) * (x3 - x0) + s * (x2 - x1)
        F = np.stack([dxds / hx[..., None], dxdt / hy[..., None]], -1)   # columns d x / d X, d x / d Y
        out[side] = F
    return out


def edge_betweenness(Dm, sigma, links):
    """Unweighted edge betweenness over ordered pairs from distance and path-count matrices."""
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(sigma > 0, 1.0 / sigma, 0.0)
    b = np.zeros(len(links))
    for e, (u, v) in enumerate(links):
        for p, q in ((u, v), (v, u)):
            m = (Dm[:, p][:, None] + 1 + Dm[q, :][None, :]) == Dm
            b[e] += np.sum(m * (sigma[:, p][:, None] * sigma[q, :][None, :]) * inv)
    return b


def path_counts(A, Dm):
    """Number of shortest paths between every pair (BFS layers as matrix products)."""
    n = A.shape[0]
    P = np.eye(n)
    S = P.copy()
    finite = np.isfinite(Dm)
    kmax = int(Dm[finite].max())
    with np.errstate(all="ignore"):                    # Accelerate BLAS can raise spurious FP flags
        for k in range(1, kmax + 1):
            P = (P @ A) * (Dm == k)
            S += P
    if not np.isfinite(S).all():
        raise FloatingPointError("non-finite shortest-path counts")
    return S


def modularity(A0, labels):
    m2 = A0.sum()
    deg = A0.sum(1)
    Q = 0.0
    for c in np.unique(labels):
        idx = labels == c
        Q += A0[np.ix_(idx, idx)].sum() / m2 - (deg[idx].sum() / m2) ** 2
    return Q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=2, help="grid cells per element side")
    ap.add_argument("--frames", type=int, default=5, help="last loading frames per protocol averaged")
    ap.add_argument("--max-parts", type=int, default=10,
                    help="stop once this many subdomains of at least --min-frac of the elements have formed")
    ap.add_argument("--min-frac", type=float, default=0.02,
                    help="pieces smaller than this fraction of the elements are recorded but not counted")
    ap.add_argument("--unscaled", action="store_true",
                    help="paper's Eq. 3: sum the squared jumps of the F components without scaling")
    ap.add_argument("--uniform-jump", action="store_true",
                    help="null control: give every link the same jump, so only geometry drives the cuts")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    mpp = D.load_scale()["mm_per_px"]
    _, df = B.load("tine_region")
    df["t"] = df.t.round(1)
    lo_z, hi_z = TINE_TIPS["z"]
    lo_t, hi_t = TINE_TIPS["theta"]

    # ---------------------------------------------------------------- jumps per protocol
    jumps, folders = [], list(dict.fromkeys(df.folder))
    for fd in folders:
        R = np.load(os.path.join(testcfg.get().results, fd + SUFFIX + ".npz"), allow_pickle=True)
        X0, shape, times = R["X0"], tuple(R["shape"]), np.round(R["times"], 1)
        rows, cols = node_lines(shape[0], a.stride), node_lines(shape[1], a.stride)
        kmap = {t: k for k, t in enumerate(times)}
        ks = [kmap[t] for t in df.t[df.folder == fd].to_numpy()[-a.frames:]]
        Xg = X0.reshape(shape + (2,))
        Xn = np.stack([Xg[np.ix_(rows[:-1], cols[:-1])], Xg[np.ix_(rows[:-1], cols[1:])],
                       Xg[np.ix_(rows[1:], cols[1:])], Xg[np.ix_(rows[1:], cols[:-1])]], 2)
        ne_r, ne_c = len(rows) - 1, len(cols) - 1
        jh, jv = np.zeros((ne_r, ne_c - 1, 2, 2)), np.zeros((ne_r - 1, ne_c, 2, 2))   # per F component
        nh, nv = np.zeros((ne_r, ne_c - 1)), np.zeros((ne_r - 1, ne_c))
        for k in ks:
            xg = R["pos"][k].reshape(shape + (2,))
            vg = R["valid"][k].reshape(shape)
            xn = np.stack([xg[np.ix_(rows[:-1], cols[:-1])], xg[np.ix_(rows[:-1], cols[1:])],
                           xg[np.ix_(rows[1:], cols[1:])], xg[np.ix_(rows[1:], cols[:-1])]], 2)
            ok = (vg[np.ix_(rows[:-1], cols[:-1])] & vg[np.ix_(rows[:-1], cols[1:])]
                  & vg[np.ix_(rows[1:], cols[1:])] & vg[np.ix_(rows[1:], cols[:-1])])
            Fe = edge_F(Xn, xn)
            # horizontal neighbors share element 1's right side and element 2's left side
            d = (Fe["right"][:, :-1] - Fe["left"][:, 1:]) ** 2
            m = ok[:, :-1] & ok[:, 1:]
            jh += np.where(m[..., None, None], d, 0); nh += m
            # vertical neighbors (grid row + 1) share element 1's top and element 2's bottom
            d = (Fe["top"][:-1] - Fe["bottom"][1:]) ** 2
            m = ok[:-1] & ok[1:]
            jv += np.where(m[..., None, None], d, 0); nv += m
        jh, jv = jh / np.maximum(nh, 1)[..., None, None], jv / np.maximum(nv, 1)[..., None, None]
        if a.unscaled:
            s2 = np.ones((2, 2))
        else:                                          # typical squared jump of each component, interior links
            s2 = np.concatenate([jh[2:-2, 2:-3].reshape(-1, 2, 2), jv[2:-3, 2:-2].reshape(-1, 2, 2)]).mean(0)
        jh, jv = (jh / s2).sum((-1, -2)), (jv / s2).sum((-1, -2))
        scale = np.percentile(np.r_[jh.ravel(), jv.ravel()], 99)
        jumps.append((jh / scale, jv / scale))
        print(f"  {fd}: frames {times[ks[0]]:.1f}-{times[ks[-1]]:.1f} s, jump 99th percentile {scale:.2e}",
              flush=True)
    Jh = sum(j[0] for j in jumps)
    Jv = sum(j[1] for j in jumps)

    # ---------------------------------------------------------------- network
    n = ne_r * ne_c
    eid = np.arange(n).reshape(ne_r, ne_c)
    links = [(eid[r, c], eid[r, c + 1]) for r in range(ne_r) for c in range(ne_c - 1)] + \
            [(eid[r, c], eid[r + 1, c]) for r in range(ne_r - 1) for c in range(ne_c)]
    wj = np.r_[Jh.ravel(), Jv.ravel()]
    if a.uniform_jump:
        wj = np.ones_like(wj)
        print("null control: uniform jump on every link")
    links = np.array(links)
    A0 = np.zeros((n, n))
    A0[links[:, 0], links[:, 1]] = A0[links[:, 1], links[:, 0]] = 1
    alive = np.ones(len(links), bool)
    parts, last_nc, last_big, step = [], 1, 1, 0
    min_size = max(1, int(np.ceil(a.min_frac * n)))
    print(f"network: {n} elements ({ne_c} x {ne_r}), {len(links)} links; removing links ...", flush=True)
    while True:
        A = np.zeros((n, n))
        L = links[alive]
        A[L[:, 0], L[:, 1]] = A[L[:, 1], L[:, 0]] = 1
        Dm = shortest_path(csr_matrix(A), unweighted=True, directed=False)
        sigma = path_counts(A, Dm)
        b = edge_betweenness(Dm, sigma, L)
        k = np.flatnonzero(alive)[np.argmax(b * wj[alive])]
        alive[k] = False
        step += 1
        A[links[k, 0], links[k, 1]] = A[links[k, 1], links[k, 0]] = 0
        nc, lab = connected_components(csr_matrix(A), directed=False)
        if nc > last_nc:
            sizes = np.bincount(lab)
            nbig = int((sizes >= min_size).sum())
            last_nc = nc
            if nbig > last_big:                        # a new subdomain of meaningful size formed
                Q = modularity(A0, lab)
                parts.append(dict(n_parts=nbig, n_pieces=int(nc), step=step, modularity=float(Q),
                                  labels=lab.copy(), sizes=sorted(sizes.tolist(), reverse=True)))
                print(f"  step {step:4d}: {nbig} subdomains >= {min_size} elements ({nc} pieces), sizes "
                      f"{sorted(sizes.tolist(), reverse=True)[:nbig]}, modularity {Q:.3f}", flush=True)
                last_big = nbig
                if nbig >= a.max_parts:
                    break

    # ---------------------------------------------------------------- save
    ecz = (Xn[..., 0, 0] + Xn[..., 1, 0]) / 2
    ect = (Xn[..., 0, 1] + Xn[..., 3, 1]) / 2
    ez_px = np.r_[X0.reshape(shape + (2,))[0, cols, 0]]
    ey_px = np.r_[X0.reshape(shape + (2,))[rows, 0, 1]]
    np.savez_compressed(
        os.path.join(a.out, "segmentation.npz"),
        node_cols=cols, node_rows=rows, element_corner_x_px=ez_px, element_corner_y_px=ey_px,
        jump_horizontal=Jh, jump_vertical=Jv,
        jump_per_protocol_h=np.array([j[0] for j in jumps]), jump_per_protocol_v=np.array([j[1] for j in jumps]),
        folders=np.array(folders), **{f"labels_{p['n_parts']}": p["labels"].reshape(ne_r, ne_c) for p in parts},
        layout="element (row, col) in image grid order (row 0 at the top); labels per element")
    with open(os.path.join(a.out, "segmentation.json"), "w") as f:
        json.dump(dict(method="Witzenburg et al. 2016 deformation-gradient-jump segmentation "
                              "(betweenness x jump, Girvan-Newman link removal)",
                       elements=f"{ne_c} x {ne_r} bilinear quads, {a.stride} x {a.stride} grid cells "
                                f"({a.stride * 15} px = {a.stride * 15 * mpp:.3f} mm)",
                       frames_per_protocol=a.frames,
                       components="unscaled (paper's Eq. 3)" if a.unscaled else
                                  "each F component scaled by its mean-square interior jump, per protocol",
                       normalization="per protocol, by the 99th percentile",
                       subdomains_counted=f"pieces of at least {min_size} elements ({a.min_frac:.0%})",
                       partitions=[{k: v for k, v in p.items() if k != "labels"} for p in parts]), f, indent=2)

    # ---------------------------------------------------------------- figure
    zmm = lambda x: (x - lo_z) * mpp                      # noqa: E731
    tmm = lambda y: (hi_t - y) * mpp                      # noqa: E731  theta up, from the bottom tips
    Wz, Wt = (hi_z - lo_z) * mpp, (hi_t - lo_t) * mpp
    segs, vals = [], []
    for r in range(ne_r):                                 # vertical segments between horizontal neighbors
        for c in range(ne_c - 1):
            x = ez_px[c + 1]
            segs.append([(zmm(x), tmm(ey_px[r])), (zmm(x), tmm(ey_px[r + 1]))]); vals.append(Jh[r, c])
    for r in range(ne_r - 1):                             # horizontal segments between vertical neighbors
        for c in range(ne_c):
            y = ey_px[r + 1]
            segs.append([(zmm(ez_px[c]), tmm(y)), (zmm(ez_px[c + 1]), tmm(y))]); vals.append(Jv[r, c])
    vals = np.array(vals)
    show = [p for p in parts if p["n_parts"] >= 2][:8]
    ncol = 4
    nrow = 1 + int(np.ceil(len(show) / ncol))
    fig, ax = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 5.0 * nrow), constrained_layout=True, squeeze=False)
    a_ = ax[0, 0]
    lc = LineCollection(segs, array=vals, cmap="magma_r", linewidths=2.2,
                        norm=plt.Normalize(0, np.percentile(vals, 99)))
    a_.add_collection(lc); fig.colorbar(lc, ax=a_, shrink=.85, label="summed normalized jump")
    a_.set_title(f"Sum of normalized F jumps, Def11–Def19\n(last {a.frames} loading frames each)", fontsize=10)
    a_ = ax[0, 1]
    a_.plot([p["n_parts"] for p in parts], [p["modularity"] for p in parts], "o-")
    a_.set_xlabel("number of subdomains"); a_.set_ylabel("modularity"); a_.grid(alpha=.3)
    a_.set_title("Modularity when each new subdomain forms", fontsize=10)
    a_ = ax[0, 2]
    a_.plot([p["step"] for p in parts], [p["n_parts"] for p in parts], "s-")
    a_.set_xlabel("links removed"); a_.set_ylabel("number of subdomains"); a_.grid(alpha=.3)
    a_.set_title("Links removed before each split", fontsize=10)
    ax[0, 3].axis("off")
    ax[0, 3].text(0, .95, "\n".join(f"{p['n_parts']:2d} subdomains: sizes {p['sizes'][:p['n_parts']][:6]}"
                                    + (" …" if p['n_parts'] > 6 else "") + f"  (Q {p['modularity']:.2f})"
                                    for p in parts), va="top", fontsize=8, family="monospace")
    edges_z = np.r_[[zmm(x) for x in ez_px]]
    edges_t = np.r_[[tmm(y) for y in ey_px]]
    for k_, p in enumerate(show):
        a_ = ax[1 + k_ // ncol, k_ % ncol]
        lab = p["labels"].reshape(ne_r, ne_c)
        order = np.argsort(-np.bincount(lab.ravel()))               # color by size, largest first
        rank = np.empty_like(order); rank[order] = np.arange(len(order))
        R_ = np.where(rank[lab] < p["n_parts"], rank[lab], np.nan).astype(float)   # small pieces -> gray
        a_.set_facecolor("0.8")
        a_.pcolormesh(edges_z, edges_t, R_, cmap="tab10", vmin=-.5, vmax=9.5)
        a_.set_title(f"{p['n_parts']} subdomains (modularity {p['modularity']:.2f})\n"
                     f"gray: {p['n_pieces'] - p['n_parts']} pieces < {min_size} elements", fontsize=10)
    for a_ in ax.ravel():
        if a_.has_data() and a_ not in (ax[0, 1], ax[0, 2]):
            a_.set_xlim(0, Wz); a_.set_ylim(0, Wt); a_.set_aspect("equal")
            a_.set_xlabel("z from left tine tips (mm)"); a_.set_ylabel(r"$\theta$ from bottom tine tips (mm)")
    for a_ in ax.ravel()[4 + len(show):]:
        a_.axis("off")
    fig.suptitle("Deformation-gradient-jump segmentation (Witzenburg et al. 2016) of the region between the tine tips"
                 + ("" if a.unscaled else "; each F component scaled by its typical jump"))
    fig.savefig(os.path.join(a.out, "fjump_segmentation.png"), dpi=85)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
