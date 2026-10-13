"""
Test selection and per-test settings, shared by strain_tracking/ and constitutive_fit/.

Every script works on one test folder. Sample folders (N1739_LFPA, N1690_RFPA, ...) sit next to
strain_tracking/ and constitutive_fit/, each holding test folders test1_..., test2_.... Select a
test with --test on the command line or with the environment variable DIC_TEST:
    --test N1739_LFPA/test2       sample/test, each as any unique prefix (e.g. N1739/test2)
    --test test2                  a bare test name works only if exactly one sample has it
    --test /path/to/testN_...     any test folder by path
Scripts that take a file inside a test folder (a tracking .npz) use the test that file is in.
There is no default, so outputs cannot land in the wrong test.
The --test option is read and removed from sys.argv before the script's own options are parsed.

Per-test settings live in <test folder>/test_config.json (make_config.py drafts it):
    sample                 label for figure titles
    roi_px                 {"standard": [x0, y0, x1, y1], "tine_region": [...]}, image px (track.py)
    tine_tips_px           rest-state tip lines {"z": [left, right], "theta": [top, bottom]}, image px
    calibration_tines_px   each tine's [centerline, tip], per rake (calibrate_scale.py, measure_tines.py)
    calibration_exclude    tine pairs left out of the scale calibration, e.g. ["Y3", "Y5"]
    widths_mm              {"tine_span": [w_theta, w_z], "specimen": [w_theta, w_z], "tine_span_by_rake": {...}};
                           w_theta divides F_z, w_z divides F_theta; tine_span = means of opposing
                           rakes (measure_tines.py)
    thickness_mm           {"H": ..., "H_sd": ...}
    reference_deformation  deformation whose first frame is the common reference, e.g. "Def11"
    dic                    optional tracker settings, e.g. {"mask_saturated": true, "exclude_px": [...]}
    notes                  caveats (dropped frames, scale uncertainty, ...), printed by status.py --notes
The image scale is in <test folder>/scale_calibration.json (calibrate_scale.py).

In Python:  import testcfg; testcfg.select("N1739/test1"); T = testcfg.get()
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                     # holds the code folders, the sample folders and comparisons/
TEST_RE = re.compile(r"test\d+_")
CODE_DIRS = [os.path.join(ROOT, "strain_tracking"), os.path.join(ROOT, "constitutive_fit")]
ENV = "DIC_TEST"
SETTINGS_FILE = os.path.join(ROOT, "analysis_settings.json")
_selected = None
_test = None


def tests(sample_dir):
    return sorted(d for d in os.listdir(sample_dir) if TEST_RE.match(d) and os.path.isdir(os.path.join(sample_dir, d)))


def samples():
    """Sample folders under ROOT: those holding at least one test folder."""
    return sorted(d for d in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, d)) and not d.startswith(".")
                  and tests(os.path.join(ROOT, d)))


def all_tests():
    return [f"{s}/{t}" for s in samples() for t in tests(os.path.join(ROOT, s))]


def _match(name, names, what, where):
    m = [d for d in names if d == name or d.startswith(name)]
    m = [d for d in m if d == name or d.startswith(name + "_")] or m
    if len(m) != 1:
        sys.exit(f"--test: {name!r} matches {m or 'no ' + what} in {where}; choose one of: {', '.join(names)}")
    return m[0]


def _resolve(name):
    """Absolute path of the test folder named by `name` (see the module docstring)."""
    if os.path.isdir(name) and TEST_RE.match(os.path.basename(os.path.abspath(name))):
        return os.path.abspath(name)
    if "/" in name:
        sample, test = name.rstrip("/").split("/", 1)
        sdir = os.path.join(ROOT, _match(sample, samples(), "sample folder", ROOT))
        return os.path.join(sdir, _match(test, tests(sdir), "test folder", sdir))
    hits = []
    for s in samples():
        m = [t for t in tests(os.path.join(ROOT, s)) if t == name or t.startswith(name + "_")]
        hits += [f"{s}/{t}" for t in m]
    if len(hits) != 1:
        sys.exit(f"--test: {name!r} matches {', '.join(hits) or 'no test folder'}; give it as SAMPLE/TEST, "
                 f"one of: {', '.join(all_tests())}")
    return os.path.join(ROOT, hits[0])


def _from_argv():
    for i, arg in enumerate(sys.argv):
        if arg == "--test" and i + 1 < len(sys.argv):
            name = sys.argv[i + 1]
            del sys.argv[i:i + 2]
            return name
        if arg.startswith("--test="):
            del sys.argv[i]
            return arg.split("=", 1)[1]
    return None


def select(name):
    """Choose the test folder for this process (before calling get())."""
    global _selected, _test
    _selected, _test = _resolve(name), None


class Test:
    def __init__(self, root):
        root = root if os.path.isabs(root) else os.path.join(ROOT, root)
        self.root = os.path.normpath(root)
        self.name = os.path.basename(self.root)                     # e.g. test1_08252026_PLC25_38-20
        self.sample = os.path.basename(os.path.dirname(self.root))  # e.g. N1739_LFPA
        self.images = os.path.join(self.root, "Images")
        self.fit_dir = os.path.join(self.images, "constitutive_fit")
        self.local_dir = os.path.join(self.fit_dir, "local")
        self.scale_file = os.path.join(self.root, "scale_calibration.json")
        self.config_file = os.path.join(self.root, "test_config.json")
        self.config = json.load(open(self.config_file)) if os.path.exists(self.config_file) else {}

    def setting(self, *keys):
        """Nested config value, e.g. setting("roi_px", "standard"); exits with a clear message if missing."""
        v = self.config
        for k in keys:
            if not isinstance(v, dict) or k not in v:
                sys.exit(f"{self.config_file}: missing setting {'.'.join(keys)}")
            v = v[k]
        return v

    def folder(self, defn):
        """Deformation folder name for 'Def11', 'Def1' (= Def01) or any unique folder prefix."""
        num = re.fullmatch(r"(?i)def(\d+)", defn)
        prefix = f"Def{int(num.group(1)):02d}_" if num else defn
        m = sorted(d for d in os.listdir(self.images)
                   if re.match(r"Def\d+_", d) and d.startswith(prefix) and os.path.isdir(os.path.join(self.images, d)))
        if len(m) != 1:
            sys.exit(f"{defn!r} matches {m or 'no folder'} in {self.images}")
        return m[0]

    def write_config(self, cfg):
        """Write test_config.json: objects and lists that fit within ~110 characters stay on one line."""
        def fmt(v, ind, lead=0):
            one = json.dumps(v, ensure_ascii=False)
            if lead + len(one) + ind <= 110 or not isinstance(v, (dict, list)):
                return one
            pad = " " * (ind + 2)
            if isinstance(v, dict):
                items = [f'{pad}{json.dumps(k)}: {fmt(x, ind + 2, len(k) + 4)}' for k, x in v.items()]
                return "{\n" + ",\n".join(items) + "\n" + " " * ind + "}"
            return "[\n" + ",\n".join(pad + fmt(x, ind + 2) for x in v) + "\n" + " " * ind + "]"
        with open(self.config_file, "w") as f:
            f.write(fmt(cfg, 0) + "\n")
        self.config = cfg

    def first(self, pattern):
        """First file in the test folder matching a glob pattern (e.g. '*Data.csv')."""
        import glob
        m = sorted(glob.glob(os.path.join(self.root, pattern)))
        if not m:
            sys.exit(f"no {pattern} in {self.root}")
        return m[0]


_ARGV_NAME = _from_argv()     # read at import, so the script's own argparse never sees --test


def _from_path(path):
    """Test folder (absolute path) containing a file, e.g. a tracking .npz, or None."""
    d = os.path.dirname(os.path.abspath(path))
    while d != os.path.dirname(d):
        if TEST_RE.match(os.path.basename(d)) and os.path.isdir(os.path.join(d, "Images")):
            return d
        d = os.path.dirname(d)
    return None


def get(path=None):
    """The selected test. Scripts that take a file inside a test folder pass it as `path`,
    and the test is then taken from where the file is (it must agree with --test if both are given)."""
    global _selected, _test
    if path is not None:
        p_name = _from_path(path)
        if p_name:
            other = _selected or (_resolve(_ARGV_NAME) if _ARGV_NAME else None)
            if other and other != p_name:
                sys.exit(f"{path} is in {p_name}, but the selected test is {other}")
            if _test is not None and _test.root != p_name:
                _test = None
            _selected = p_name
    if _test is None:
        name = _selected or _ARGV_NAME or os.environ.get(ENV) or os.environ.get("N1739_TEST")
        if not name:
            sys.exit(f"choose a test with --test (or {ENV}): {', '.join(all_tests())}")
        _selected = _resolve(name)
        _test = Test(_selected)
    return _test


def analysis():
    """Project-wide physical sizes (mm) from analysis_settings.json: subset_mm, step_mm, window_mm,
    edge_mm, low_spot_box_mm. Each test converts them to px with its own image scale."""
    with open(SETTINGS_FILE) as f:
        return json.load(f)


def edge_mm(override=None):
    """Interior edge distances {"z": mm from the left/right tine lines, "theta": mm from the top/bottom
    ones}: override (one number, a (z, theta) pair or a dict) if given, else edge_mm of analysis_settings.json."""
    e = analysis()["edge_mm"] if override is None else override
    if isinstance(e, dict):
        return {"z": float(e["z"]), "theta": float(e["theta"])}
    if isinstance(e, (list, tuple)):
        return {"z": float(e[0]), "theta": float(e[1])}
    return {"z": float(e), "theta": float(e)}


def script_path(name):
    """Path of a code file in strain_tracking/ or constitutive_fit/."""
    for d in CODE_DIRS:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(name)
