"""
Loading branch of every analyzed deformation: from the Stretch start to the rake peak (dic2d.rake_peak),
found automatically. Draws them shaded on the rake forces (DIC_results/loading_branches.png). With
--review (runDIC.sh: review = yes), the branches are shown in a window first, and wrong ones are set by
clicking their start and end. Branches set by hand are saved in the test's loading_branches.json, the
same file the MATLAB version writes; --loading-only tracking (track_protocols.py) and the fits
(compare_stretch_inputs.py, fit_4pt.py) use them in place of the automatic ones.

    .venv/bin/python strain_tracking/review_branches.py --test N1757/test1 --review
    .venv/bin/python strain_tracking/review_branches.py --test N1757/test1 --defs Def11 Def12
"""
import argparse
import datetime
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import review  # noqa: E402
import testcfg  # noqa: E402

FZ, FT = "#1f77b4", "#ff7f0e"        # F_z (axial, X) and F_theta (circumferential, Y)


def branches(T, defs):
    """One dict per deformation: key, name, automatic branch (a0, a1), Recover start, plot window, branch used."""
    ph = pd.read_csv(glob.glob(os.path.join(T.root, "Logs", "*_PhaseTiming.csv"))[0])
    Dd = pd.read_csv(T.first("*Data.csv"), encoding="latin1").drop_duplicates("Time_S").set_index("Time_S")
    stretch = ph[ph.Phase == "Stretch"].reset_index(drop=True)
    out = []
    for n in defs:
        st = stretch.iloc[n - 1]
        rc = ph[(ph.Phase == "Recover") & (ph.Set == st.Set) & (ph.Cycle == st.Cycle)].iloc[0]
        a0, t_rec = round(float(st["StartTime(s)"]), 1), round(float(rc["StartTime(s)"]), 1)
        a1 = D.rake_peak(Dd, a0, t_rec)
        t0, t1 = D.branch_by_hand(T.root, f"Def{n:02d}") or (a0, a1)
        dur = t_rec - a0
        out.append(dict(key=f"Def{n:02d}", name=re.sub(r"\s+save\s+cyc?le$", "", st.SetName, flags=re.I),
                        a0=a0, a1=a1, t_rec=t_rec, t0=t0, t1=t1,
                        w0=max(Dd.index[0], min(t0, a0) - max(2, 0.15 * dur)),
                        w1=min(Dd.index[-1], max(t1, t_rec) + max(4, 0.4 * dur))))
    return Dd, out


def by_hand(b):
    return abs(b["t0"] - b["a0"]) > 0.05 or abs(b["t1"] - b["a1"]) > 0.05


def draw(ax, Dd, b, span, shade):
    """Rake forces around one branch, with `span` shaded and the logged Recover start dotted."""
    w = Dd.loc[b["w0"]:b["w1"]]
    h1, = ax.plot(w.index, w.X_Force_mN, color=FZ, lw=1.2, label="F_z (axial)")
    h2, = ax.plot(w.index, w.Y_Force_mN, color=FT, lw=1.2, label="F_θ (circumferential)")
    hp = ax.axvspan(*span, color=shade, alpha=0.35, lw=0, label="loading branch")
    ax.axvline(b["t_rec"], ls=":", color="0.4")
    ax.set_xlim(b["w0"], b["w1"])
    ax.set_xlabel("time (s)")
    ax.set_ylabel("force (mN)")
    ax.grid(True, alpha=0.4)
    return [h1, h2, hp]


def plot_all(fig, Dd, info, T, top=0.9, bottom=0.06):
    """Every branch on the rake forces, as small multiples."""
    nc = int(np.ceil(np.sqrt(len(info))))
    nr = int(np.ceil(len(info) / nc))
    gs = fig.add_gridspec(nr, nc, top=top, bottom=bottom, left=0.05, right=0.98, hspace=0.5, wspace=0.25)
    for i, b in enumerate(info):
        ax = fig.add_subplot(gs[i // nc, i % nc])
        hand = by_hand(b)
        h = draw(ax, Dd, b, (b["t0"], b["t1"]), "#f28c8c" if hand else "#8cd98c")
        ax.set_title(f"{b['key']} {b['name']}{', set by hand' if hand else ''}", fontsize=10)
        if i == 0:
            ax.legend(handles=h, loc="upper left", fontsize=8)
    fig.suptitle(f"{T.name.replace('_', ' ')}\nshaded: loading branch (green automatic, red set by hand); "
                 "dotted: logged Recover start", fontsize=review.TEXT_FONT)


def click_branch(Dd, b):
    """Start and end of one loading branch, clicked; "auto" for the automatic one, None if closed."""
    fig = review.figure(f"runDIC: loading branch of {b['key']}")
    ax = fig.add_axes([0.06, 0.08, 0.92, 0.8])
    h = draw(ax, Dd, b, (b["a0"], b["a1"]), "0.8")
    h[2].set_label("automatic branch")
    ax.legend(handles=h, loc="upper left")
    pts = review.click_points(fig, ax, 2, f"{b['key']} {b['name']}: click where loading starts, then where it "
                              "ends (the force peak)", mark=lambda x, y: ax.axvline(x, color="r", lw=1.5),
                              keys={"a": ("auto", "A goes back to the automatic branch")})
    if pts is None or pts == "auto":
        return pts
    t = sorted(round(min(max(x, b["w0"]), b["w1"]), 1) for x, _ in pts)   # Data.csv rows are 0.1 s apart
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--defs", nargs="+", default=[f"Def{n}" for n in range(11, 20)],
                    help="deformations to review (default: Def11-Def19)")
    ap.add_argument("--review", action="store_true", help="show the branches and ask whether they are right")
    a = ap.parse_args()
    if a.review:
        review.enable()
    T = testcfg.get()
    defs = [int(d.lower().replace("def", "")) for d in a.defs]
    Dd, info = branches(T, defs)
    f = os.path.join(T.root, "loading_branches.json")
    old = json.load(open(f)) if os.path.exists(f) else {}

    reviewed = review.enabled()
    while reviewed:
        fig = review.figure("runDIC: loading branches", toolbar=False)
        plot_all(fig, Dd, info, T, top=0.88, bottom=0.2)
        answer = review.ask(fig, "Does each shaded band cover one loading curve, from where the forces start to "
                            "rise to where they peak?", ["Continue", "No, fix some"])
        if answer != "No, fix some":
            break
        items = [f"{b['key']}  {b['name']}{'  (set by hand)' if by_hand(b) else ''}" for b in info]
        for item in review.pick("runDIC", "Loading branches to click again:", items, "Click them"):
            b = info[items.index(item)]
            t = click_branch(Dd, b)
            if t == "auto":
                b["t0"], b["t1"] = b["a0"], b["a1"]
            elif t:
                b["t0"], b["t1"] = t

    # the branches set by hand, with any of other deformations kept
    hand = {k: v for k, v in (old.get("by_hand") or {}).items() if k not in {b["key"] for b in info}}
    for b in info:
        if by_hand(b):
            hand[b["key"]] = [b["t0"], b["t1"]]
            print(f"  {b['key']} {b['name']}: {b['t0']:.1f}-{b['t1']:.1f} s, set by hand "
                  f"(automatic: {b['a0']:.1f}-{b['a1']:.1f} s)")
    out = {"reviewed": datetime.datetime.now().strftime("%Y-%m-%d %H:%M") if reviewed else old.get("reviewed"),
           "by_hand": dict(sorted(hand.items()))}
    if out != old:
        with open(f, "w") as fh:
            json.dump({k: v for k, v in out.items() if v is not None}, fh, indent=2)
            fh.write("\n")

    for b in info:                         # a loading-only track that ends before the branch does
        npz = glob.glob(os.path.join(T.results, f"{b['key']}_*.npz"))
        npz = [q for q in npz if not q.endswith(("_fields.npz", "_strain.npz"))]
        if npz and by_hand(b):
            t_last = float(np.load(npz[0], allow_pickle=True)["times"].max())
            if t_last < b["t1"] - 0.05:
                print(f"  note: the {b['key']} track ends at {t_last:.1f} s, before the branch set by hand "
                      f"({b['t1']:.1f} s): track {b['key']} again (runDIC.sh --steps track --defs {b['key']})")

    import matplotlib
    if not reviewed:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14, 10))
    plot_all(fig, Dd, info, T)
    os.makedirs(T.results, exist_ok=True)
    png = os.path.join(T.results, "loading_branches.png")
    fig.savefig(png, dpi=90)
    plt.close(fig)
    print(f"{len(info)} loading branches, {sum(map(by_hand, info))} set by hand -> loading_branches.json, "
          f"DIC_results/loading_branches.png")


if __name__ == "__main__":
    main()
