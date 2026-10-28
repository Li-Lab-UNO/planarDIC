"""
Check windows shown during a run: the sample label and notes (sample_info.py), the sorted-image table
(sort_images.py), the tine tips (make_config.py), the loading branches (review_branches.py), the thickness
measurement (thickness/measure_interactive.py) and the tracking report (view_tracking.py). They open when runDIC.sh runs with review = yes (it sets DIC_REVIEW=1) or a script is
given --review, and only where there is a display; otherwise the scripts run without asking.

Each window answers with buttons along its bottom edge; closing a window counts as its first answer.
In the click windows, the scroll wheel zooms where the pointer is, a right-click (or Backspace) removes
the last click, and the toolbar's zoom and pan work as usual (no clicks are taken while they are on).
"""
import os
import sys

ENV = "DIC_REVIEW"
BUTTON_H = 0.3                   # button height (in) and text size: as MATLAB's own dialogs (its step list)
BUTTON_FONT = 9
TEXT_FONT = 10                   # questions, labels and instructions


def enable():
    os.environ[ENV] = "1"


def enabled():
    """True when check windows are wanted and can be shown."""
    if os.environ.get(ENV) != "1":
        return False
    if sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        print("no display: the check windows are skipped")
        return False
    return True


def pyplot():
    """matplotlib.pyplot with a window backend: the native one on macOS (the Tk bundled with the system
    Python is too old there), Tk elsewhere (Ubuntu: sudo apt install python3-tk)."""
    import matplotlib
    matplotlib.use("macosx" if sys.platform == "darwin" else "TkAgg")
    import matplotlib.pyplot as plt
    return plt


def figure(title, size=(14, 8.5), toolbar=True):
    """A large window, maximized where the backend allows it; the toolbar (zoom, pan) only where there is
    something to zoom."""
    plt = pyplot()
    fig = plt.figure(title, figsize=size, facecolor="w")
    tb = fig.canvas.toolbar
    if not toolbar and hasattr(tb, "pack_forget"):
        tb.pack_forget()                              # Tk; the macOS window keeps it (made without one, it
    win = getattr(fig.canvas.manager, "window", None)  # can keep its event loop running once closed)
    try:
        if sys.platform.startswith("linux"):
            win.attributes("-zoomed", True)
        elif sys.platform == "win32":
            win.state("zoomed")
    except Exception:
        pass                                          # macOS, or a window manager without maximize
    return fig


def wait(fig):
    """Show fig, in front, and return once it is closed. plt.show() brings the window to the front; on
    macOS it notices a closed window only at the next mouse or keyboard event, so closing fig also stops
    its event loop at once (_macosx.stop, as Ctrl-C does)."""
    plt = pyplot()
    if sys.platform == "darwin":
        from matplotlib.backends import _macosx
        fig.canvas.mpl_connect("close_event", lambda _e: _macosx.stop())
    plt.show()


def button_width(label):
    """Width (in) of a button that fits `label` at BUTTON_FONT."""
    return max(0.95, 0.3 + 0.075 * len(label))


def button(fig, label, x, y, width=None):
    """A Button with its lower left at (x, y) inches from the figure's, BUTTON_H tall and as wide as its
    label (or `width`, in)."""
    from matplotlib.widgets import Button
    W, H = fig.get_size_inches()
    w = width or button_width(label)
    b = Button(fig.add_axes([x / W, y / H, w / W, BUTTON_H / H]), label)
    b.label.set_fontsize(BUTTON_FONT)
    return b


def button_row(fig, labels, x_right, y):
    """Buttons side by side, 0.1 in apart, the last one ending at x_right (in); in the order of labels."""
    out, x = [], x_right
    for label in reversed(labels):
        w = button_width(label)
        x -= w
        out.append(button(fig, label, x, y, w))
        x -= 0.1
    return out[::-1]


def ask(fig, question, answers):
    """Show `question` and one button per answer at the bottom of `fig`; return the answer clicked (the
    first one if the window is closed)."""
    plt = pyplot()
    W, H = fig.get_size_inches()                         # sizes in inches, so any window size looks the same
    fig.text(0.5, (0.15 + BUTTON_H + 0.12) / H, question, ha="center", va="bottom", fontsize=TEXT_FONT, wrap=True)
    total = sum(button_width(a) for a in answers) + 0.1 * (len(answers) - 1)
    chosen, buttons = [answers[0]], []
    for label, b in zip(answers, button_row(fig, answers, (W + total) / 2, 0.15)):

        def on_click(_event, label=label):
            chosen[0] = label
            plt.close(fig)
        b.on_clicked(on_click)
        buttons.append(b)                               # widgets stop working once garbage-collected
    wait(fig)
    return chosen[0]


def pick(title, question, items, ok):
    """Window with a check box per item and an `ok` button; returns the items checked ([] if closed)."""
    from matplotlib.widgets import CheckButtons
    plt = pyplot()
    fig = plt.figure(title, figsize=(3.6, 1.2 + 0.3 * len(items)), facecolor="w")
    fig.text(0.05, 0.95, question, va="top", fontsize=TEXT_FONT)
    boxes = CheckButtons(fig.add_axes([0.05, 0.6 / (1.2 + 0.3 * len(items)), 0.9, 0.3 * len(items) / (1.2 + 0.3 * len(items))],
                                      frameon=False), items,
                         label_props={"fontsize": [TEXT_FONT] * len(items)}, frame_props={"s": [60] * len(items)},
                         check_props={"s": [60] * len(items)})
    done = [False]
    b = button(fig, ok, (3.6 - button_width(ok)) / 2, 0.12)

    def on_ok(_event):
        done[0] = True
        plt.close(fig)
    b.on_clicked(on_ok)
    wait(fig)
    return [x for x, on in zip(items, boxes.get_status()) if on] if done[0] else []


def form(title, intro, fields, ok="Save", cancel="Cancel", width=11):
    """A window with one text box per field ((label, initial text) pairs) and two buttons; returns the texts
    when `ok` is pressed (or Enter in the last box), None when cancelled or closed."""
    from matplotlib.widgets import TextBox
    plt = pyplot()
    rows = len(fields)
    fig = plt.figure(title, figsize=(width, 1.3 + 0.42 * rows), facecolor="w")
    H = 1.3 + 0.42 * rows
    fig.text(0.02, 1 - 0.15 / H, intro, va="top", fontsize=TEXT_FONT, wrap=True)
    boxes = []
    for i, (label, initial) in enumerate(fields):
        y = 1 - (0.6 + 0.42 * (i + 1)) / H
        fig.text(0.02, y + 0.15 / H, label, va="center", fontsize=TEXT_FONT)
        tb = TextBox(fig.add_axes([0.3, y, 0.68, BUTTON_H / H]), "", initial=initial or "")
        tb.text_disp.set_fontsize(TEXT_FONT)
        boxes.append(tb)
    answer = [None]

    def finish(_=None):
        answer[0] = [b.text.strip() for b in boxes]
        plt.close(fig)
    b_ok, b_cancel = button_row(fig, [ok, cancel], width - 0.2, 0.15)
    b_ok.on_clicked(finish)
    b_cancel.on_clicked(lambda _: plt.close(fig))
    boxes[-1].on_submit(lambda _: finish())
    fig._widgets = (boxes, b_ok, b_cancel)                 # widgets stop working once garbage-collected
    wait(fig)
    return answer[0]


def prefs(**update):
    """Settings remembered between runs (~/.config/planarDIC/review.json); with keywords, update them."""
    import json
    f = os.path.join(os.path.expanduser("~"), ".config", "planarDIC", "review.json")
    try:
        p = json.load(open(f))
    except (OSError, ValueError):
        p = {}
    if update:
        p.update(update)
        os.makedirs(os.path.dirname(f), exist_ok=True)
        json.dump(p, open(f, "w"), indent=2)
    return p


def click_points(fig, ax, n, head, mark=None, keys=None):
    """Collect n clicks (x, y in data coordinates) on ax; None if the window is closed first. `mark(x, y)`
    draws a click and returns its artist (default: a red circle); `keys` maps a key to (value, hint): the
    key closes the window and returns value, and the hint is added to the title, e.g.
    {"a": ("auto", "A goes back to the automatic branch")}."""
    plt = pyplot()
    pts, marks, pressed = [], [], []
    keys = keys or {}
    mark = mark or (lambda x, y: ax.plot(x, y, "o", color="r", ms=10, mfc="none", mew=1.5)[0])
    extra = "".join(f"; {v[1]}" for v in keys.values())

    def status():
        ax.set_title(f"{head}\n{len(pts)} of {n} clicked. Scroll to zoom; right-click or Backspace removes "
                     f"the last one{extra}", fontsize=13)
        fig.canvas.draw_idle()

    def undo():
        if pts:
            pts.pop()
            marks.pop().remove()
            status()

    def on_press(e):
        tb = fig.canvas.toolbar
        if e.inaxes is not ax or (tb is not None and tb.mode) or len(pts) >= n:
            return
        if e.button == 1:
            pts.append((e.xdata, e.ydata))
            marks.append(mark(e.xdata, e.ydata))
            status()
            if len(pts) == n:
                ax.set_title(f"{head}\n{n} of {n} clicked", fontsize=13)
                t = fig.canvas.new_timer(interval=600)
                t.single_shot = True
                t.add_callback(plt.close, fig)
                t.start()
                fig._close_timer = t                    # keep a reference until it fires
        elif e.button == 3:
            undo()

    def on_key(e):
        if e.key in ("backspace", "delete"):
            undo()
        elif e.key and e.key.lower() in keys:
            pressed.append(keys[e.key.lower()][0])
            plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_press)
    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("scroll_event", lambda e: zoom(ax, e))
    status()
    wait(fig)
    if pressed:
        return pressed[0]
    return pts if len(pts) == n else None


def zoom(ax, e, factor=1.25):
    """Scroll-wheel zoom about the pointer."""
    if e.inaxes is not ax or e.xdata is None:
        return
    f = 1 / factor if e.button == "up" else factor
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    ax.set_xlim(e.xdata + (x0 - e.xdata) * f, e.xdata + (x1 - e.xdata) * f)
    ax.set_ylim(e.ydata + (y0 - e.ydata) * f, e.ydata + (y1 - e.ydata) * f)
    ax.figure.canvas.draw_idle()
