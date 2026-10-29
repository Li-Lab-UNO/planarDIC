"""
Measure a test's wall thickness from its photos (<SAMPLE>/Thickness, <ID>/test <N>/), interactively: the
Python version of the MATLAB version's measureThickness, step for step. For each photo the user sets the
scale (a line across one grid spacing or ruler tick, then the scale check: a longer line between grid lines,
with the grid spacings in between counted), crops the strip (two corners; the view stays on it) and traces
its centerline; the edges are found automatically (strip_edges.py) and can be refined (remove a section,
redraw the centerline, click the edges by hand). One strip per crop; a photo can hold more than one.

The thickness is the mean of every position used, all strips of all photos pooled, saved in the test's
test_config.json as thickness_mm (H, H_sd, source). When the test already has a thickness, the window
first asks whether to keep it or measure again. When nothing is measured (no photos, or the window closed)
and the test has no thickness, a box asks to type one (source "typed"). Results go to <test>/thickness/:
stations.csv (every position, used or why not), summary.csv (per strip, with its scale and crop) and an
overlay per photo.

    .venv/bin/python thickness/measure_interactive.py --test N1757/test1
"""
import argparse
import os
import re
import sys
import textwrap

import cv2
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), "strain_tracking"))
import strip_edges as E  # noqa: E402
import review  # noqa: E402
import testcfg  # noqa: E402

DISPLAY_PX = 2000                       # displayed image size (px) for the whole photo
UNITS = {"mm": 1.0, "cm": 10.0, "in": 25.4}
RED, ORANGE, GRAY, GREEN = "#ff3333", "#ff9900", "#8c8c8c", "#33d94d"
CLICK = "#ff0000"                # the clicks and the line following the cursor


def photos_of(root):
    """The thickness photo folder of a test, <sample>/Thickness*/test <N>, and its .jpg photos. Names are matched
    ignoring case, with or without the space ("Test 1", "test1"), since Linux file names are case-sensitive."""
    m = re.match(r"(?i)test(\d+)", os.path.basename(root))
    if not m:
        return None, []

    def subfolders(parent, pattern):
        return sorted(x for x in os.listdir(parent)
                      if re.search(pattern, x, re.I) and os.path.isdir(os.path.join(parent, x)))

    sample = os.path.dirname(root)
    for d in subfolders(sample, r"^thickness"):
        for t in subfolders(os.path.join(sample, d), r"^test\s*(\d+)$"):
            if int(re.search(r"\d+", t)[0]) != int(m[1]):
                continue
            f = os.path.join(sample, d, t)
            p = sorted(x for x in os.listdir(f) if re.search(r"\.jpe?g$", x, re.I))
            if p:
                return f, p
    return None, []


def used(S):
    return S["ok"] & ~S["outlier"] & ~S["removed"]


def chords(S, sel):
    """Chords (edge to edge) at the selected positions, as one NaN-separated line."""
    k = np.nonzero(sel)[0]
    if not len(k):
        return [], []
    xs = np.c_[S["x"][k] - S["rm"][k] * S["nx"][k], S["x"][k] + S["rp"][k] * S["nx"][k], np.full(len(k), np.nan)]
    ys = np.c_[S["y"][k] - S["rm"][k] * S["ny"][k], S["y"][k] + S["rp"][k] * S["ny"][k], np.full(len(k), np.nan)]
    return xs.ravel(), ys.ravel()


class Measure:
    def __init__(self, T, photo_dir, photos):
        from matplotlib.widgets import Button, RadioButtons, TextBox
        self.plt = review.pyplot()
        self.T, self.photo_dir, self.photos = T, photo_dir, photos
        p = review.prefs().get("thickness_square", {"value": "", "unit": "in"})
        self.i_photo, self.state, self.pts, self.title, self.ins = -1, "", [], "", []
        self.L1 = self.size_mm = self.pxmm = np.nan
        self.n_sq, self.scale_ok = 0, False
        self.S, self.edge_a, self.strip_no, self.msg, self.frame = None, [], 1, "", None
        self.strips, self.photo_strips, self.finished = [], [], False

        fig = self.fig = review.figure(f"runDIC: thickness, {T.name}")
        ax = self.ax = fig.add_axes([0.01, 0.02, 0.65, 0.95])
        ax.axis("off")
        self.img = ax.imshow(np.zeros((2, 2, 3), np.uint8), interpolation="nearest")
        line = lambda *a, **k: ax.plot([], [], *a, **k)[0]  # noqa: E731
        self.h_line = line("-o", color=CLICK, mfc=CLICK, mec="w", lw=2, ms=5)
        self.h_ok = line("-", color=RED, lw=0.8)
        self.h_out = line("-", color=ORANGE, lw=1.2)
        self.h_ex = line("-", color=GRAY, lw=0.8)
        self.h_edge = line(".", color=RED, ms=3)
        self.h_done = line("-", color=GREEN, lw=0.6)
        self.h_frame = line("--", color="w", lw=1)                    # the crop, once accepted
        self.h_prev = line("--", color=CLICK, lw=1.5)                 # follows the cursor
        x0, w = self.x0, self.w = 0.68, 0.30                        # the panel on the right (layout() places it)
        self.t_title = fig.text(x0, 0.975, "", fontsize=11, weight="bold", va="top")
        self.t_ins = fig.text(x0, 0.885, "", fontsize=review.TEXT_FONT, va="top")
        self.t_size = fig.text(x0, 0.632, "Size", fontsize=review.TEXT_FONT, va="center")
        self.box = TextBox(fig.add_axes([x0 + 0.04, 0.615, 0.08, 0.035]), "", initial=p.get("value", ""))
        self.box.on_text_change(lambda _: self.t_err.set_text(""))
        self.radio = RadioButtons(fig.add_axes([x0 + 0.14, 0.585, 0.06, 0.075], frameon=False), list(UNITS),
                                  active=list(UNITS).index(p.get("unit", "in")) if p.get("unit") in UNITS else 2)
        self.t_err = fig.text(x0, 0.565, "", fontsize=review.TEXT_FONT, color="#cc1a1a", va="top")
        self.buttons, self.actions = [], []
        for k in range(6):                                  # sizes in inches, set in layout()
            b = Button(fig.add_axes([x0, 0.48 - k * 0.052, w, 0.045]), "")
            b.label.set_fontsize(review.BUTTON_FONT)
            b.on_clicked(lambda _e, k=k: self.actions[k]() if k < len(self.actions) else None)
            self.buttons.append(b)
        self.t_stat = fig.text(x0, 0.2, "", fontsize=review.TEXT_FONT, va="top", family="monospace")
        self.b_view = Button(fig.add_axes([x0, 0.02, 0.12, 0.04]), "Whole photo")
        self.b_view.label.set_fontsize(review.BUTTON_FONT)
        self.b_view.on_clicked(lambda _: self.whole_photo())
        fig.text(x0, 0.075, "Scroll to zoom where the pointer is.", fontsize=review.TEXT_FONT, color="0.4", va="center")
        fig.canvas.mpl_connect("button_press_event", self.on_press)
        fig.canvas.mpl_connect("key_press_event", self.on_key)
        fig.canvas.mpl_connect("motion_notify_event", self.on_move)
        fig.canvas.mpl_connect("scroll_event", self.on_scroll)
        fig.canvas.mpl_connect("resize_event", lambda _e: self.layout())
        self.next_photo()

    def run(self):
        """Returns the strips measured (a list per photo), or None when the window was closed first."""
        review.wait(self.fig)
        return self.photo_strips if self.finished else None

    # ================================================================== photos and states
    def next_photo(self):
        if self.i_photo >= 0:
            self.photo_strips.append(self.strips)
        self.i_photo += 1
        if self.i_photo >= len(self.photos):
            self.finished = True
            self.plt.close(self.fig)
            return
        self.strips, self.strip_no, self.S, self.pts, self.edge_a, self.msg, self.frame = [], 1, None, [], [], "", None
        self.t_title.set_text(f"Loading {self.photos[self.i_photo]} ...")
        self.fig.canvas.draw()
        bgr = cv2.imread(os.path.join(self.photo_dir, self.photos[self.i_photo]))
        self.rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
        self.H, self.W = self.rgb.shape[:2]
        self.ds = max(1, int(np.ceil(max(self.H, self.W) / DISPLAY_PX)))
        self.small = self.rgb[::self.ds, ::self.ds]
        self.whole_photo()
        self.set_state("scale1")

    def set_state(self, st):
        self.state, self.pts = st, []
        if st in ("seed", "edgeA"):
            self.edge_a = []
        self.render()

    def render(self):
        st = self.state
        self.title = f"{self.photos[self.i_photo]} (photo {self.i_photo + 1} of {len(self.photos)})"
        show_size = st == "scale1"
        for a in (self.box.ax, self.radio.ax):
            a.set_visible(show_size)
        self.box.set_active(show_size)
        self.radio.active = show_size                         # (set_active picks a button here)
        self.t_size.set_visible(show_size)
        self.t_err.set_text("")
        stat, B = "", []
        n = len(self.pts)
        if st == "scale1":
            ins = ["Step 1 of 5: Scale, 1 grid spacing", "",
                   "1. Click on a grid line, then on the next one, straight across one grid spacing (or one ruler "
                   "tick).",
                   "2. Type the size of that spacing and choose its unit, then press Next.", "",
                   "Right-click to remove a marker."]
            if n == 2:
                stat = f"one grid spacing: {self.L1:.1f} px"
            B = [("Next", self.scale_next), ("Skip this photo", self.next_photo)]
        elif st == "scale2":
            ins = ["Step 2 of 5: Scale check", "",
                   "1. Click on a grid line near the strip, then on another grid line as far away as the photo "
                   "allows.",
                   "2. Check the count of grid spacings between your clicks (below) and press Confirm.", "",
                   "Right-click to remove a marker."]
            B = [("Back", lambda: self.set_state("scale1"))]
            self.scale_ok = False
            if n == 2:
                L2 = np.hypot(*np.subtract(self.pts[1], self.pts[0]))
                p, N, ratio = E.scale_from_lines(self.L1, L2, self.size_mm)
                if abs(ratio - N) > 0.25:
                    ins += ["", f"{ratio:.2f} grid spacings between your clicks: not a whole number, so a click is off "
                            "a grid line (here or in step 1). Draw again."]
                    B = [("Draw again", lambda: self.set_state("scale2"))] + B
                else:
                    self.n_sq, self.pxmm, self.scale_ok = N, p, True
                    u = self.unit()
                    ins += ["", f"{N} grid spacings between your clicks = {round(N * self.size_mm / UNITS[u], 4):g} {u} "
                            f"({N * self.size_mm:.2f} mm): scale {self.pxmm:.2f} px/mm."]
                    B = [("Confirm", lambda: self.set_state("frame")),
                         ("Draw again", lambda: self.set_state("scale2"))] + B
        elif st == "frame":
            ins = ["Step 3 of 5: Crop", "",
                   "1. Click two opposite corners of an ROI around the strip.",
                   "2. Press Use this crop. The view then stays on it.", "",
                   "Crop whole strips only; leave out pieces cut off at the edge of the photo. If this photo has a "
                   "second strip of this test, you can crop it after this one."]
            B = [("No crop (whole photo)", lambda: self.use_frame(None))]
            if n == 2:
                B = [("Use this crop", lambda: self.use_frame(self.pts)),
                     ("Draw again", lambda: self.set_state("frame"))] + B
        elif st == "seed":
            ins = ["Step 4 of 5: Trace the centerline", "",
                   "1. Click along the middle of the strip, from one end to the other.",
                   "2. Press Done (or Enter, or double-click the last point). The edges are then found for you.", "",
                   "Right-click to remove a marker."]
            stat = f"{n} click(s)"
            B = [("Done", self.seed_done), ("Skip this photo", self.next_photo)]
        elif st == "check":
            ins = ["Step 5 of 5: Check the thickness measurements", "",
                   "Each red line should cross the strip and end on both of its edges.", "",
                   "Orange lines were left out automatically (a loose fiber or a cut end). Gray lines are ones you "
                   "removed."]
            stat = self.strip_stats()
            B = [("Yes, use this strip", lambda: self.use_strip(False)),
                 ("Use it, and crop another strip", lambda: self.use_strip(True)),
                 ("Remove a section", lambda: self.set_state("cut")),
                 ("Redraw the centerline", lambda: self.set_state("seed")),
                 ("Click the edges myself", lambda: self.set_state("edgeA")),
                 ("Skip this strip", self.skip_strip)]
        elif st == "cut":
            ins = ["Remove a section", "", "Click the two ends of the stretch to leave out, along the strip. Its "
                   "lines turn gray."]
            stat = self.strip_stats()
            B = [("Cancel", lambda: self.set_state("check"))]
        elif st == "edgeA":
            ins = ["Click the edges yourself, 1 of 2", "",
                   "Click along one edge of the strip from end to end, then press Done.", "",
                   "Right-click to remove a marker."]
            B = [("Done", self.edge_a_done), ("Cancel", lambda: self.set_state("seed"))]
        else:                                                # edgeB
            ins = ["Click the edges yourself, 2 of 2", "",
                   "Now click along the other edge from end to end, then press Done.", "",
                   "Right-click to remove a marker."]
            B = [("Done", self.edge_b_done), ("Cancel", lambda: self.set_state("seed"))]
        if self.msg:
            stat = f"{stat}\n\n{self.msg}"
        self.ins = ins
        self.t_stat.set_text(stat)
        self.actions = [f for _, f in B]
        for k, b in enumerate(self.buttons):
            on = k < len(B)
            b.ax.set_visible(on)
            b.set_active(on)
            if on:
                b.label.set_text(B[k][0])
        self.layout()
        self.draw_overlays()
        self.fig.canvas.draw_idle()

    def wrap(self, lines, fontsize):
        """Lines wrapped to the panel's width (characters of about 0.55 font sizes)."""
        n = max(20, int(self.w * self.fig.get_figwidth() * 72 / (0.55 * fontsize)))
        return "\n".join(textwrap.fill(t, n) if t else "" for t in lines)

    def layout(self):
        """Place the panel from the top: title, instructions, the size row (scale, part 1), the buttons shown
        and the readout, each below the one before, at the window's current size."""
        fig, x0, w = self.fig, self.x0, self.w
        r = fig.canvas.get_renderer()
        bottom = lambda t: t.get_window_extent(r).transformed(fig.transFigure.inverted()).y0  # noqa: E731
        W, H = fig.get_size_inches()
        bh, step = review.BUTTON_H / H, (review.BUTTON_H + 0.08) / H      # buttons as in MATLAB's own dialogs
        self.t_title.set_text(self.wrap([self.title], 11))
        self.t_ins.set_y(bottom(self.t_title) - 0.03)
        self.t_ins.set_text(self.wrap(self.ins, review.TEXT_FONT))
        y = bottom(self.t_ins) - 0.025
        if self.box.ax.get_visible():
            self.t_size.set_y(y - 0.0375)
            self.box.ax.set_position([x0 + 0.04, y - 0.0375 - bh / 2, 0.08, bh])
            self.radio.ax.set_position([x0 + 0.14, y - 0.075, 0.06, 0.075])
            y -= 0.085
        self.t_err.set_y(y)
        y -= 0.035
        shown = [b for b in self.buttons if b.ax.get_visible()]
        bw = min(w, max([review.button_width(b.label.get_text()) for b in shown], default=1) / W)
        for b in shown:                               # one width, the widest label's, text on the left: a tidy column
            b.ax.set_position([x0, y - bh, bw, bh])
            b.label.set_horizontalalignment("left")
            b.label.set_x(0.12 / (bw * W))
            y -= step
        self.t_stat.set_y(y - 0.015)
        self.b_view.ax.set_position([x0, 0.15 / H, review.button_width("Back to the crop") / W, bh])

    # ================================================================== actions
    def unit(self):
        return self.radio.value_selected

    def scale_next(self):
        try:
            v = float(self.box.text.strip().replace(",", "."))
        except ValueError:
            v = np.nan
        if len(self.pts) < 2:
            self.err("Click across one grid spacing first.")
            return
        if not (np.isfinite(v) and v > 0):
            self.err("Enter the size of the grid spacing, for example 0.5.")
            return
        self.size_mm = v * UNITS[self.unit()]
        self.L1 = np.hypot(*np.subtract(self.pts[1], self.pts[0]))
        review.prefs(thickness_square={"value": self.box.text.strip(), "unit": self.unit()})
        self.set_state("scale2")

    def err(self, text):
        self.t_err.set_text(text)
        self.fig.canvas.draw_idle()

    def use_frame(self, P):
        """P: two opposite corners, or None for no frame."""
        if P is None:
            self.frame = None
        else:
            P = np.asarray(P)
            self.frame = [max(-0.5, P[:, 0].min()), max(-0.5, P[:, 1].min()),
                          min(self.W - 0.5, P[:, 0].max()), min(self.H - 0.5, P[:, 1].max())]
        self.set_state("seed")
        self.whole_photo()

    def seed_done(self):
        if len(self.pts) < 2:
            self.err("Click at least twice along the strip.")
            return
        self.t_stat.set_text("Finding the edges ...")
        self.fig.canvas.draw()
        self.new_strip(E.measure_strip(self.rgb, np.array(self.pts), self.pxmm))
        if used(self.S).sum() < 10:
            self.msg = "Fewer than 10 positions found: redraw the\ncenterline, or click the edges yourself."

    def new_strip(self, S):
        S["removed"] = np.zeros(len(S["ok"]), bool)
        self.S, self.msg, self.state, self.pts = S, "", "check", []
        self.zoom_to_strip()
        self.render()

    def zoom_to_strip(self):
        """Show the strip just measured, with 2 mm around it."""
        m = 2 * self.pxmm
        S = self.S
        xs = [S["x"].min() - m, S["x"].max() + m]
        ys = [S["y"].min() - m, S["y"].max() + m]
        side = max(xs[1] - xs[0], (ys[1] - ys[0]) * self.W / self.H)
        cx, cy = np.mean(xs), np.mean(ys)
        self.set_view(cx - side / 2, cx + side / 2, cy - side * self.H / self.W / 2, cy + side * self.H / self.W / 2)

    def edge_a_done(self):
        if len(self.pts) < 2:
            self.err("Click at least twice along the edge.")
            return
        self.edge_a, self.state, self.pts = list(self.pts), "edgeB", []
        self.render()

    def edge_b_done(self):
        if len(self.pts) < 2:
            self.err("Click at least twice along the edge.")
            return
        S = E.measure_manual(np.array(self.edge_a), np.array(self.pts), self.pxmm)
        self.edge_a = []
        self.new_strip(S)

    def use_strip(self, another):
        """Keep the strip; then the next photo, or (another) frame a further strip of this photo, with its scale
        kept."""
        if used(self.S).sum() < 10:
            self.err("A strip needs at least 10 positions.")
            return
        self.strips.append(dict(photo=self.photos[self.i_photo], strip=self.strip_no, S=self.S, pxmm=self.pxmm,
                                scale=f"{self.n_sq} squares of {self.size_mm:.4g} mm", frame=self.frame))
        if not another:
            self.next_photo()
            return
        self.msg = (f"Strip saved ({np.mean(self.S['thickness_mm'][used(self.S)]):.3f} mm). "
                    "Crop the next strip.")
        self.strip_no, self.S, self.frame = self.strip_no + 1, None, None
        self.whole_photo()
        self.set_state("frame")

    def skip_strip(self):
        self.S = None
        self.next_photo()

    def strip_stats(self):
        S = self.S
        u = used(S)
        t = S["thickness_mm"][u]
        return (f"thickness {np.mean(t):.3f} +- {np.std(t, ddof=1) if len(t) > 1 else 0:.3f} mm\n"
                f"{u.sum()} positions used, {(S['outlier'] & ~S['removed']).sum()} flagged,\n"
                f"{S['removed'].sum()} removed, {(~S['ok']).sum()} without edges")

    # ================================================================== mouse, keys, view
    def on_press(self, e):
        tb = self.fig.canvas.toolbar
        if e.inaxes is not self.ax or (tb is not None and tb.mode) or e.xdata is None:
            return
        if e.button == 3:
            self.undo()
        elif e.button == 1 and e.dblclick:
            self.primary()                                    # its first click placed the point
        elif e.button == 1:
            self.add_click((e.xdata, e.ydata))

    def add_click(self, p):
        st = self.state
        if st in ("scale1", "scale2", "frame"):
            if len(self.pts) >= 2:
                self.pts = []
            self.pts.append(p)
            if st == "scale1" and len(self.pts) == 2:
                self.L1 = np.hypot(*np.subtract(self.pts[1], self.pts[0]))
            self.render()
        elif st in ("seed", "edgeA", "edgeB"):
            self.pts.append(p)
            self.render()
        elif st == "cut":
            self.pts.append(p)
            if len(self.pts) == 2:
                S = self.S
                sa, sb = (S["s"][np.argmin(np.hypot(S["x"] - q[0], S["y"] - q[1]))] for q in self.pts)
                S["removed"] |= (S["s"] >= min(sa, sb)) & (S["s"] <= max(sa, sb))
                self.state, self.pts = "check", []
            self.render()

    def undo(self):
        if self.pts:
            self.pts.pop()
            self.render()

    def primary(self):
        """Double-click or Enter: proceed, as the step's first button does, once its clicks are placed."""
        st = self.state
        if st == "scale1":
            self.scale_next()
        elif st == "scale2" and self.scale_ok:
            self.set_state("frame")
        elif st == "frame" and len(self.pts) == 2:
            self.use_frame(self.pts)
        elif st == "seed":
            self.seed_done()
        elif st == "edgeA":
            self.edge_a_done()
        elif st == "edgeB":
            self.edge_b_done()

    def on_key(self, e):
        if self.box.capturekeystrokes:                       # typing the square's size
            return
        if e.key in ("enter", "return"):
            self.primary()
        elif e.key in ("backspace", "delete"):
            self.undo()

    def on_move(self, e):
        """Crosshair over the photo; preview from the last click to the cursor: the scale line (with the live
        square count on the long line), the frame, or the next segment of a centerline or edge."""
        from matplotlib.backend_tools import Cursors
        inside = e.inaxes is self.ax and e.xdata is not None
        self.fig.canvas.set_cursor(Cursors.SELECT_REGION if inside else Cursors.POINTER)
        if not self.pts:
            return
        if not inside:
            self.h_prev.set_data([], [])
            self.fig.canvas.draw_idle()
            return
        q, p, st = self.pts[-1], (e.xdata, e.ydata), self.state
        if st in ("scale1", "scale2") and len(self.pts) == 1:
            self.h_prev.set_data([q[0], p[0]], [q[1], p[1]])
            d = np.hypot(p[0] - q[0], p[1] - q[1])
            self.t_stat.set_text(f"{d / self.L1:.2f} grid spacings" if st == "scale2" else f"{d:.1f} px")
        elif st == "frame" and len(self.pts) == 1:
            self.h_prev.set_data([q[0], p[0], p[0], q[0], q[0]], [q[1], q[1], p[1], p[1], q[1]])
        elif st in ("seed", "edgeA", "edgeB"):
            self.h_prev.set_data([q[0], p[0]], [q[1], p[1]])
        else:
            return
        self.fig.canvas.draw_idle()

    def on_scroll(self, e):
        if e.inaxes is not self.ax or e.xdata is None:
            return
        (x0, x1), (y1, y0) = self.ax.get_xlim(), self.ax.get_ylim()      # y axis points down
        f = 1 / 1.25 if e.button == "up" else 1.25
        w, h = min(self.W, (x1 - x0) * f), min(self.H, (y1 - y0) * f)
        if w < 50:
            return
        a, b = (e.xdata - x0) / (x1 - x0), (e.ydata - y0) / (y1 - y0)
        nx = min(max(e.xdata - a * w, -0.5), self.W - 0.5 - w)
        ny = min(max(e.ydata - b * h, -0.5), self.H - 0.5 - h)
        self.set_view(nx, nx + w, ny, ny + h)

    def set_view(self, x0, x1, y0, y1):
        x0, x1 = max(x0, -0.5), min(x1, self.W - 0.5)
        y0, y1 = max(y0, -0.5), min(y1, self.H - 0.5)
        self.ax.set_xlim(x0, x1)
        self.ax.set_ylim(y1, y0)
        self.refresh_image()

    def whole_photo(self):
        """The frame around the strips when one is set (with 5% around it), else the whole photo."""
        if self.frame is None:
            self.set_view(-0.5, self.W - 0.5, -0.5, self.H - 0.5)
        else:
            f = self.frame
            c = ((f[0] + f[2]) / 2, (f[1] + f[3]) / 2)
            side = 1.05 * max(f[2] - f[0], (f[3] - f[1]) * self.W / self.H)
            self.set_view(c[0] - side / 2, c[0] + side / 2, c[1] - side * self.H / self.W / 2,
                          c[1] + side * self.H / self.W / 2)
        self.b_view.label.set_text("Whole photo" if self.frame is None else "Back to the crop")

    def refresh_image(self):
        """The whole photo subsampled, or, zoomed in, the visible part at up to full resolution."""
        (x0, x1), (y1, y0) = self.ax.get_xlim(), self.ax.get_ylim()
        if x1 - x0 > 0.6 * self.W:
            d = self.ds
            h, w = self.small.shape[:2]
            self.img.set_data(self.small)
            self.img.set_extent((-0.5 * d, (w - 0.5) * d, (h - 0.5) * d, -0.5 * d))
        else:
            cx = (max(0, int(x0)), min(self.W - 1, int(np.ceil(x1))))
            cy = (max(0, int(y0)), min(self.H - 1, int(np.ceil(y1))))
            st = max(1, int(np.ceil((cx[1] - cx[0]) / DISPLAY_PX)))
            crop = self.rgb[cy[0]:cy[1] + 1:st, cx[0]:cx[1] + 1:st]
            h, w = crop.shape[:2]
            self.img.set_data(crop)
            self.img.set_extent((cx[0] - 0.5 * st, cx[0] + (w - 0.5) * st, cy[0] + (h - 0.5) * st, cy[0] - 0.5 * st))
        self.ax.set_xlim(x0, x1)
        self.ax.set_ylim(y1, y0)
        self.fig.canvas.draw_idle()

    def draw_overlays(self):
        for h in (self.h_line, self.h_ok, self.h_out, self.h_ex, self.h_edge, self.h_prev):
            h.set_data([], [])
        st, P = self.state, np.array(self.pts).reshape(-1, 2)
        if st in ("scale1", "scale2", "frame", "seed", "edgeA", "edgeB", "cut") and len(P):
            self.h_line.set_data(P[:, 0], P[:, 1])
        if st == "frame" and len(P) == 2:
            self.h_line.set_data(P[[0, 1, 1, 0, 0], 0], P[[0, 0, 1, 1, 0], 1])
        f = self.frame
        self.h_frame.set_data(([f[0], f[2], f[2], f[0], f[0]], [f[1], f[1], f[3], f[3], f[1]])
                              if f is not None and st != "frame" else ([], []))
        if st == "edgeB" and self.edge_a:
            A = np.array(self.edge_a)
            self.h_ex.set_data(A[:, 0], A[:, 1])
        if st in ("check", "cut") and self.S is not None:
            S = self.S
            self.h_ok.set_data(*chords(S, used(S)))
            self.h_out.set_data(*chords(S, S["outlier"] & ~S["removed"]))
            self.h_ex.set_data(*chords(S, S["removed"] & S["ok"]))
            k = S["ok"] & ~S["removed"]
            self.h_edge.set_data(np.r_[S["x"][k] + S["rp"][k] * S["nx"][k], S["x"][k] - S["rm"][k] * S["nx"][k]],
                                 np.r_[S["y"][k] + S["rp"][k] * S["ny"][k], S["y"][k] - S["rm"][k] * S["ny"][k]])
        dx, dy = [], []                                      # strips already measured on this photo
        for a in self.strips:
            x, y = chords(a["S"], used(a["S"]))
            dx += list(x)
            dy += list(y)
        self.h_done.set_data(dx, dy)


def results(T, photo_dir, photo_strips):
    """Pool the strips; write <test>/thickness/ (stations.csv, summary.csv, an overlay per photo). Returns
    (H, sd, source) or None when no strip was measured."""
    strips = [a for ps in photo_strips for a in ps]
    if not strips:
        return None
    t, rows, summ = [], [], []
    for a in strips:
        S = a["S"]
        u = used(S)
        t += list(S["thickness_mm"][u])
        flag = np.array(S["reason"], dtype=object)
        flag[S["outlier"] & ~S["removed"]] = "outlier"
        flag[S["removed"]] = "removed"
        for i in range(len(S["ok"])):
            rows.append(dict(photo=a["photo"], strip=a["strip"], method=S["method"], px_per_mm=a["pxmm"],
                             x_px=S["x"][i], y_px=S["y"][i], thickness_mm=S["thickness_mm"][i], used=int(u[i]),
                             flag=flag[i]))
        v = S["thickness_mm"][u]
        summ.append(dict(photo=a["photo"], strip=a["strip"], method=S["method"], px_per_mm=a["pxmm"], scale=a["scale"],
                         frame_px=str([round(x) for x in a["frame"]]) if a["frame"] is not None else "",
                         n_used=int(u.sum()), mean_mm=v.mean(), sd_mm=v.std(ddof=1), median_mm=np.median(v)))
    t = np.array(t)
    n_ph = len({a["photo"] for a in strips})
    source = f"measured from {n_ph} photo(s), {len(strips)} strip(s), {len(t)} positions"
    out = os.path.join(T.root, "thickness")
    os.makedirs(out, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(out, "stations.csv"), index=False)
    pd.DataFrame(summ).to_csv(os.path.join(out, "summary.csv"), index=False)
    for ps in photo_strips:
        if ps:
            overlay(photo_dir, ps, os.path.join(out, os.path.splitext(ps[0]["photo"])[0] + "_overlay.png"))
    return float(t.mean()), float(t.std(ddof=1)), source


def overlay(photo_dir, ps, path):
    """The photo around its strips (3 mm margin), with the chords used."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    im = cv2.cvtColor(cv2.imread(os.path.join(photo_dir, ps[0]["photo"])), cv2.COLOR_BGR2RGB)
    xy = np.vstack([np.c_[a["S"]["x"], a["S"]["y"]] for a in ps])
    fr = ps[0]["frame"]
    if fr is not None:
        xy = np.vstack([xy, [fr[0], fr[1]], [fr[2], fr[3]]])
    m = 3 * ps[0]["pxmm"]
    cx = (max(0, int(xy[:, 0].min() - m)), min(im.shape[1] - 1, int(xy[:, 0].max() + m)))
    cy = (max(0, int(xy[:, 1].min() - m)), min(im.shape[0] - 1, int(xy[:, 1].max() + m)))
    st = max(1, int(np.ceil(max(cx[1] - cx[0], cy[1] - cy[0]) / 1600)))
    fig = Figure(figsize=(12, 10))
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 0.93])
    ax.imshow(im[cy[0]:cy[1] + 1:st, cx[0]:cx[1] + 1:st], extent=(cx[0] - 0.5, cx[1] + 0.5, cy[1] + 0.5, cy[0] - 0.5))
    ax.axis("off")
    if fr is not None:
        ax.plot([fr[0], fr[2], fr[2], fr[0], fr[0]], [fr[1], fr[1], fr[3], fr[3], fr[1]], "--", color="#ffd933", lw=1.5)
    lab = []
    for a in ps:
        S = a["S"]
        u = used(S)
        ax.plot(*chords(S, u), "-", color=RED, lw=0.5)
        ax.text(S["x"][0], S["y"][0], f" {a['strip']}", color="y", fontsize=14, weight="bold")
        v = S["thickness_mm"][u]
        lab.append(f"strip {a['strip']}: {v.mean():.3f} +- {v.std(ddof=1):.3f} mm (n = {u.sum()})")
    ax.set_title(f"{ps[0]['photo']}, {ps[0]['pxmm']:.2f} px/mm: {'; '.join(lab)}")
    fig.savefig(path, dpi=90)


def save(T, H, source, sd=None):
    """thickness_mm in the test's test_config.json: H, H_sd (when given) and source; other entries kept."""
    cfg = dict(T.config)
    t = dict(cfg.get("thickness_mm") or {})
    t["H"] = round(H, 5)
    if sd is not None and np.isfinite(sd):
        t["H_sd"] = round(sd, 4)
    else:
        t.pop("H_sd", None)
    t["source"] = source
    cfg["thickness_mm"] = t
    T.write_config(cfg)


def typed(T, why):
    """Ask for the thickness (mm) in a box; saved with source "typed". Returns it, or None."""
    while True:
        a = review.form("runDIC: thickness", f"{why} Type the wall thickness of {T.name.split('_')[0]} (mm), for "
                        "example from the thickness workbook, or leave it blank to skip.",
                        [("Thickness (mm)", "")], ok="Save", cancel="Skip")
        if not a or not a[0]:
            return None
        try:
            H = float(a[0].replace(",", "."))
        except ValueError:
            H = np.nan
        if 0 < H <= 20:
            save(T, H, "typed")
            return H


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    T = testcfg.get()
    test = T.name.split("_")[0]
    review.enable()
    if not review.enabled():
        sys.exit("measuring the thickness needs a display")
    old = (T.config.get("thickness_mm") or {})
    H_old = old.get("H")
    photo_dir, photos = photos_of(T.root)
    if H_old is not None:
        keep = True
        if photos:
            fig = review.figure("runDIC: thickness", size=(10, 1.9), toolbar=False)
            keep = review.ask(fig, f"{test} already has a thickness: {H_old:.3f} mm\n({old.get('source', '')}).\n\n"
                              "Keep it, or measure it again from the photos?",
                              ["Keep it and proceed", "Measure again"]) != "Measure again"
        if keep:
            print(f"{test}: kept {H_old:.3f} mm ({old.get('source', '')})")
            return
    if not photos:
        H = typed(T, f"No thickness photos for {test} (Thickness*/test {test[4:]}).")
        print(f"{test}: thickness {H:.3f} mm (typed)" if H else f"{test}: no thickness photos, no thickness")
        return
    strips = Measure(T, photo_dir, photos).run()
    r = results(T, photo_dir, strips) if strips is not None else None
    if r is None:
        if H_old is not None:
            print(f"{test}: nothing measured, kept {H_old:.3f} mm ({old.get('source', '')})")
            return
        H = typed(T, "No strip was measured." if strips is not None else "The window was closed.")
        print(f"{test}: thickness {H:.3f} mm (typed)" if H else f"{test}: no thickness (stresses NaN)")
        return
    H, sd, source = r
    save(T, H, source, sd)
    print(f"{test}: thickness {H:.3f} +- {sd:.3f} mm ({source}) -> test_config.json, {test}/thickness/")


if __name__ == "__main__":
    main()
