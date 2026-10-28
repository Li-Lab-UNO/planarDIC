"""
Tracking report of a test: the quality of every tracked deformation, and a window that plays them back
side by side (the Python version of the MATLAB version's viewTracking).

Quality, over each deformation's loading branch (Stretch start to the rake peak, or as set by hand in
loading_branches.json): the share of grid points tracked per frame and the median match quality (ZNCC)
per frame, their worst frame, and a status: Good when the worst frame has at least GOOD = (93 %, 0.95),
Problem below PROBLEM = (80 %, 0.90), Check between (the MATLAB version's limits). Written to
DIC_results/tracking_report.csv, with a one-line verdict printed.

The window shows the deformations as tiles: the tracked points over the images (green: tracked, red:
failed in that frame and filled in from the neighbors), every 4th grid row and column, and the stretches
(affine fit), valid points and ZNCC of the frame shown. One slider moves all of them by percent of their
loading branch; Play runs them, paced by the clock (a slow computer shows fewer frames, not a slower
playback). It opens at the worst frame of the worst deformation. A click on a tile selects it; "Open
large" shows the selected one in a big viewer. runDIC.sh's report step, right after tracking and the
artifact check, writes the report and, with review = yes, opens the window (--ask): when a deformation is
marked Check or Problem, closing it asks whether to proceed to the test's figures, exports and fit.

    .venv/bin/python strain_tracking/view_tracking.py --test N1757/test1             report and window
    .venv/bin/python strain_tracking/view_tracking.py --test N1757/test1 --no-window report only
"""
import argparse
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dic2d as D  # noqa: E402
import review  # noqa: E402
import testcfg  # noqa: E402

GOOD = (93.0, 0.95)         # worst frame: valid points (%) and median ZNCC at least this -> Good
PROBLEM = (80.0, 0.90)      # below this -> Problem; between -> Check
LINE_EVERY = 4              # grid lines drawn: every 4th row and column
TILE_PX = 360               # tile size (image px): images are read at about this size
SPEEDS = [1, 2, 5, 10]
COL = {"good": "#268c26", "check": "#d98c00", "problem": "#d91a1a"}
STATUS = ["good", "check", "problem"]


def grade(x, good, problem):
    """0 good, 1 check, 2 problem (higher is better for both measures)."""
    return 0 if x >= good else 1 if x >= problem else 2


def loading_frames(T, folder, times):
    """Indices of the frames of the loading branch."""
    hand = D.branch_by_hand(T.root, folder)
    t0 = hand[0] if hand else times[0]
    t1 = D.loading_end_time(T.root, folder)
    k = np.nonzero((times >= t0 - 0.05) & (times <= t1 + 0.05))[0]
    return k if len(k) else np.arange(len(times))


def report(T, defs):
    """One dict per tracked deformation (the quality of its loading branch), and the verdict line."""
    rows = []
    for d in defs:
        try:
            folder = T.folder(d)
        except SystemExit:
            continue
        path = os.path.join(T.results, folder + ".npz")
        if not os.path.exists(path):
            continue
        R = np.load(path, allow_pickle=True)
        k = loading_frames(T, folder, np.round(R["times"], 1))
        v = 100 * R["valid"][k].mean(1)
        z = np.nanmedian(R["zncc"][k], 1)
        kv, kz = int(np.argmin(v)), int(np.argmin(z))
        sv, sz = grade(v[kv], GOOD[0], PROBLEM[0]), grade(z[kz], GOOD[1], PROBLEM[1])
        w = kz if sz > sv or (sz == sv and v[kv] >= GOOD[0]) else kv          # frame of the worse measure
        rows.append(dict(name=folder, frames=len(k), branch=str(R["branch"]) if "branch" in R.files else "",
                         valid_median=round(float(np.median(v)), 2), valid_worst=round(float(v[kv]), 2),
                         zncc_median=round(float(np.median(z)), 4), zncc_worst=round(float(z[kz]), 4),
                         worst_frame=int(k[w]), worst_fraction=round(w / max(len(k) - 1, 1), 4),
                         status=STATUS[max(sv, sz)]))
    if not rows:
        sys.exit(f"no tracks in {T.results}: run the track step first")
    pd.DataFrame(rows).to_csv(os.path.join(T.results, "tracking_report.csv"), index=False)
    bad = [r for r in rows if r["status"] != "good"]
    if not bad:
        v = (f"all {len(rows)} deformations good (worst frame {min(r['valid_worst'] for r in rows):.1f}% valid, "
             f"ZNCC {min(r['zncc_worst'] for r in rows):.3f})")
    else:
        v = "; ".join(f"{r['name'].split('_')[0]} {r['status'].upper()} (worst frame {r['valid_worst']:.1f}% valid, "
                      f"ZNCC {r['zncc_worst']:.3f})" for r in bad) + ". Check the tracking report (view_tracking.py)"
    return rows, "Tracking: " + v


def grid_lines(x, y, shape, every=LINE_EVERY):
    """Every `every`-th row and column of a point grid (stored row by row, shape = (ny, nx)) as one
    NaN-separated line."""
    X, Y = x.reshape(shape), y.reshape(shape)
    r, c = slice(0, None, every), slice(0, None, every)
    nan_r = np.full((X[r].shape[0], 1), np.nan)
    nan_c = np.full((X[:, c].T.shape[0], 1), np.nan)
    gx = np.r_[np.hstack([X[r], nan_r]).ravel(), np.hstack([X[:, c].T, nan_c]).ravel()]
    gy = np.r_[np.hstack([Y[r], nan_r]).ravel(), np.hstack([Y[:, c].T, nan_c]).ravel()]
    return gx, gy


class Tile:
    """One deformation in an axes: the image (the tine region, read at about `size` px), the tine-tip
    rectangle, grid lines, points and the readout of the frame shown."""

    def __init__(self, ax, T, row, size, small=True):
        self.ax, self.row, self.small = ax, row, small
        R = np.load(os.path.join(T.results, row["name"] + ".npz"), allow_pickle=True)
        self.pos, self.valid, self.zncc = R["pos"], R["valid"], R["zncc"]
        self.X0, self.shape = R["X0"], tuple(int(s) for s in R["shape"])
        times = np.round(R["times"], 1)
        self.k = loading_frames(T, row["name"], times)
        self.times = times
        self.paths = [os.path.join(T.images, row["name"], str(p)) for p in R["paths"]]
        (z0, z1), (y0, y1) = T.setting("tine_tips_px", "z"), T.setting("tine_tips_px", "theta")
        m = 0.35 * max(z1 - z0, y1 - y0)
        h, w = cv2.imread(self.paths[0], cv2.IMREAD_GRAYSCALE).shape
        self.c0, self.c1 = max(0, int(z0 - m)), min(w, int(z1 + m))
        self.r0, self.r1 = max(0, int(y0 - m)), min(h, int(y1 + m))
        ds = max(1, (max(self.c1 - self.c0, self.r1 - self.r0)) // size)
        self.red = 4 if ds >= 4 else 2 if ds >= 2 else 1          # decode the jpg at 1/2 or 1/4 size (faster)
        self.ds = self.red * max(1, ds // self.red)
        self.im = ax.imshow(self.read(0), cmap="gray", vmin=0, vmax=255, interpolation="nearest")
        ax.plot(self.tx(np.array([z0, z1, z1, z0, z0])), self.ty(np.array([y0, y0, y1, y1, y0])), "--",
                color="#ffbf1a", lw=1)
        self.hgrid, = ax.plot([], [], "-", color=(0.3, 0.9, 1, 0.6), lw=0.6)
        self.hgood, = ax.plot([], [], ".", color="#33e633", ms=2.5 if small else 4)
        self.hbad, = ax.plot([], [], "x", color="#ff3333", ms=4 if small else 6, mew=1.2)
        ax.set_xticks([])
        ax.set_yticks([])
        yy, xx = np.indices(self.shape)
        self.sub = ((yy % 2 == 0) & (xx % 2 == 0)).ravel() if small else np.ones(xx.size, bool)
        self.shown = None

    def tx(self, x):
        return (np.asarray(x) - self.c0) / self.ds

    def ty(self, y):
        return (np.asarray(y) - self.r0) / self.ds

    def read(self, j):
        flag = {1: cv2.IMREAD_GRAYSCALE, 2: cv2.IMREAD_REDUCED_GRAYSCALE_2, 4: cv2.IMREAD_REDUCED_GRAYSCALE_4}
        img = cv2.imread(self.paths[j], flag[self.red])
        s = self.ds // self.red
        return img[self.r0 // self.red:self.r1 // self.red:s, self.c0 // self.red:self.c1 // self.red:s]

    def frame(self, f):
        """Frame index of fraction f of the loading branch."""
        return int(self.k[int(round(f * (len(self.k) - 1)))])

    def draw(self, f, grid, points, selected=False, force=False):
        j = self.frame(f)
        if j == self.shown and not force:
            return
        self.shown = j
        self.im.set_data(self.read(j))
        x, y, v = self.pos[j, :, 0], self.pos[j, :, 1], self.valid[j]
        g = v & self.sub
        self.hgood.set_data(self.tx(x[g]), self.ty(y[g])) if points else self.hgood.set_data([], [])
        self.hbad.set_data(self.tx(x[~v]), self.ty(y[~v])) if points else self.hbad.set_data([], [])
        if grid:
            gx, gy = grid_lines(x, y, self.shape)
            self.hgrid.set_data(self.tx(gx), self.ty(gy))
        else:
            self.hgrid.set_data([], [])
        F = D.affine_fit(self.X0[v], self.pos[j, v])
        C = F.T @ F
        n = list(self.k).index(j) + 1
        self.ax.set_xlabel(f"frame {n}/{len(self.k)}   valid {100 * v.mean():.1f}%   ZNCC "
                           f"{np.nanmedian(self.zncc[j]):.3f}\nλ_z {np.sqrt(C[0, 0]):.3f}   λ_θ {np.sqrt(C[1, 1]):.3f}"
                           f"   t = {self.times[j]:.1f} s", fontsize=9 if self.small else 11, family="monospace",
                           color="0.2")
        st = self.row["status"]
        self.ax.set_title(f"{'> ' if selected else ''}{self.row['name'].replace('_', ' ')}   {st.upper()}",
                          color=COL[st], fontsize=11 if self.small else 13, weight="bold" if selected else "normal")
        for s in self.ax.spines.values():
            s.set_color(COL[st])
            s.set_linewidth(0.5 if st == "good" else 3)


class Report:
    def __init__(self, T, rows):
        from matplotlib.widgets import CheckButtons, Slider
        self.plt = review.pyplot()
        self.T, self.rows = T, rows
        sev = [STATUS.index(r["status"]) for r in rows]
        self.sel = sorted(range(len(rows)), key=lambda i: (-sev[i], rows[i]["valid_worst"]))[0]
        self.f = 1.0 if sev[self.sel] == 0 else rows[self.sel]["worst_fraction"]
        self.playing, self.speed = False, 1
        fig = self.fig = review.figure(f"runDIC: tracking report, {T.sample}/{T.name.split('_')[0]}", toolbar=False)
        msg, sub, c = self.verdict(sev)
        fig.text(0.01, 0.975, msg, fontsize=11, weight="bold", color=c, va="top")
        fig.text(0.01, 0.945, sub, fontsize=review.TEXT_FONT, color="0.35", va="top")
        self.checks = CheckButtons(fig.add_axes([0.86, 0.935, 0.13, 0.055], frameon=False), ["Grid", "Points"],
                                   [True, True])
        self.checks.on_clicked(lambda _: self.draw(force=True))
        n = len(rows)
        nc = int(np.ceil(np.sqrt(n)))
        nr = int(np.ceil(n / nc))
        gs = fig.add_gridspec(nr, nc, left=0.01, right=0.99, bottom=0.12, top=0.9, hspace=0.32, wspace=0.04)
        self.tiles = [Tile(fig.add_subplot(gs[i // nc, i % nc]), T, r, TILE_PX) for i, r in enumerate(rows)]
        W, _ = fig.get_size_inches()                    # buttons as in MATLAB's own dialogs (review.button)
        self.bplay = review.button(fig, "Play", 0.15, 0.2, review.button_width("Pause"))
        self.bplay.on_clicked(lambda _: self.play())
        self.bspeed = review.button(fig, f"{SPEEDS[self.speed]}x", 0.25 + review.button_width("Pause"), 0.2)
        self.bspeed.on_clicked(lambda _: self.next_speed())
        self.slider = Slider(fig.add_axes([0.17, 0.035, 0.58, 0.025]), "", 0, 1, valinit=self.f, valfmt="%.2f")
        self.slider.valtext.set_visible(False)
        self.slider.on_changed(self.on_slider)
        ax_mark = fig.add_axes([0.17, 0.008, 0.58, 0.022])
        ax_mark.set_xlim(0, 1)
        ax_mark.set_ylim(0, 1)
        ax_mark.axis("off")
        for r in rows:                                   # the worst frame of each marked deformation
            if r["status"] != "good":
                ax_mark.plot(r["worst_fraction"], 0.5, "^", ms=8, color=COL[r["status"]])
        self.pct = fig.text(0.76, 0.04, "", fontsize=review.TEXT_FONT)
        self.blarge = review.button_row(fig, ["Open large"], W - 0.15, 0.2)[0]
        self.blarge.on_clicked(lambda _: self.large())
        fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.timer = fig.canvas.new_timer(interval=40)
        self.timer.add_callback(self.tick)
        self.draw(force=True)

    def verdict(self, sev):
        nb, nck = sev.count(2), sev.count(1)
        excl = len(self.T.config.get("dic", {}).get("exclude_px", []))
        art = ("Artifact check: nothing excluded." if not excl
               else f"Artifact check: {excl} region(s) excluded and tracked again.")
        if nb + nck == 0:
            return (f"All {len(sev)} deformations tracked well. {art}",
                    "Press Play to watch them; the slider moves all of them through their loading branch.", COL["good"])
        p = []
        if nb:
            p.append(f"{nb} deformation{' has' if nb == 1 else 's have'} a tracking problem")
        if nck:
            who = f"{nck}" if nb else f"{nck} deformation{'s' if nck > 1 else ''}"
            p.append(f"{who} need{'s' if nck == 1 else ''} a look")
        return (f"{', '.join(p)}. {art}", f"Opened at the worst frame of {self.rows[self.sel]['name'].split('_')[0]}. "
                "Play it before using its results.", COL["problem" if nb else "check"])

    def draw(self, force=False):
        grid, points = self.checks.get_status()
        for i, t in enumerate(self.tiles):
            t.draw(self.f, grid, points, i == self.sel, force)
        self.pct.set_text(f"{100 * self.f:.0f}% of loading")
        self.fig.canvas.draw_idle()

    def on_slider(self, val):
        self.f = float(val)
        self.draw()

    def on_click(self, e):
        for i, t in enumerate(self.tiles):
            if e.inaxes is t.ax and i != self.sel:
                self.sel = i
                self.draw(force=True)

    def next_speed(self):
        self.speed = (self.speed + 1) % len(SPEEDS)
        self.bspeed.label.set_text(f"{SPEEDS[self.speed]}x")
        self.fig.canvas.draw_idle()

    def play(self):
        self.playing = not self.playing
        self.bplay.label.set_text("Pause" if self.playing else "Play")
        if self.playing:
            if self.f >= 1:
                self.f = 0.0
            self.dur = max(t.times[t.k[-1]] - t.times[t.k[0]] for t in self.tiles)   # longest branch (s)
            self.f0, self.t0 = self.f, time.monotonic()
            self.timer.start()
        else:
            self.timer.stop()
        self.fig.canvas.draw_idle()

    def tick(self):
        f = min(1.0, self.f0 + (time.monotonic() - self.t0) * SPEEDS[self.speed] / self.dur)   # paced by the clock
        self.slider.set_val(f)
        if f >= 1:
            self.play()

    def large(self):
        Large(self.T, self.rows[self.sel], self.f)

    def show(self):
        review.wait(self.fig)


class Large:
    """One deformation in a big window, with its own slider."""

    def __init__(self, T, row, f):
        from matplotlib.widgets import Slider
        fig = self.fig = review.figure(f"runDIC: {row['name']}")
        self.tile = Tile(fig.add_axes([0.02, 0.12, 0.96, 0.82]), T, row, 1400, small=False)
        self.slider = Slider(fig.add_axes([0.1, 0.03, 0.8, 0.03]), "", 0, 1, valinit=f)
        self.slider.valtext.set_visible(False)
        self.slider.on_changed(lambda v: self.draw(float(v)))
        fig.canvas.mpl_connect("scroll_event", lambda e: review.zoom(self.tile.ax, e))
        fig._large = self                                  # keep the widgets alive with the window
        self.draw(f)
        fig.show()

    def draw(self, f):
        self.tile.draw(f, True, True, True, force=True)
        self.fig.canvas.draw_idle()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--defs", nargs="+", default=[f"Def{n}" for n in range(11, 20)],
                    help="deformations in the report (default: Def11-Def19)")
    ap.add_argument("--no-window", action="store_true", help="write the report only")
    ap.add_argument("--ask", action="store_true",
                    help="after the window, when a deformation is marked Check or Problem, ask whether to proceed "
                         "with this test; stopping exits with an error (runDIC.sh then skips its later steps)")
    a = ap.parse_args()
    T = testcfg.get()
    rows, verdict = report(T, a.defs)
    print(f"{verdict} -> DIC_results/tracking_report.csv")
    if a.no_window:
        return
    review.enable()
    if not review.enabled():
        return
    Report(T, rows).show()
    bad = [r for r in rows if r["status"] != "good"]
    if a.ask and bad:
        who = f"{T.sample.split('_')[0]}/{T.name.split('_')[0]}"
        marks = ", ".join(f"{r['name'].split('_')[0]} {r['status'].upper()}" for r in bad)
        fig = review.figure("runDIC: tracking", size=(10, 1.9), toolbar=False)
        if review.ask(fig, f"{who}: {marks}.\n\nProceed with its figures, exports and fit, or stop this test here?",
                      ["Proceed", "Stop this test"]) == "Stop this test":
            sys.exit(f"stopped after the tracking report ({marks})")


if __name__ == "__main__":
    main()
