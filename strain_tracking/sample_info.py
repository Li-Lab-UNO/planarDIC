"""
The label (e.g. healthy, or slightly calcified) and notes (tears, re-mounting, anything unusual) of a test,
in a window filled in with what is saved: the Python version of the MATLAB version's info step. The label
goes into test_config.json's sample name ("N1757 LFPA, test1 (10/06/2026), labeled healthy"), the notes into
its notes list (status.py --notes prints them). Saving with an empty note box removes that note.

    .venv/bin/python strain_tracking/sample_info.py --test N1757/test1
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import review  # noqa: E402
import testcfg  # noqa: E402

EMPTY_NOTES = 2               # empty note boxes after the saved notes (at least 4 boxes)


def base_name(T):
    """'N1757 LFPA, test1 (10/06/2026)'."""
    m = re.match(r"(test\d+)_(\d\d)(\d\d)(\d{4})", T.name)
    return f"{T.sample.replace('_', ' ')}, " + (f"{m[1]} ({m[2]}/{m[3]}/{m[4]})" if m else T.name)


def label_of(T):
    """The label in the sample name ('' if none), without 'labeled'."""
    s, base = T.config.get("sample", ""), base_name(T)
    label = s[len(base):].lstrip(", ") if s.startswith(base) else ""
    return re.sub(r"^labeled\s+", "", label, flags=re.I)


def sample_name(T, label):
    label = label.strip()
    if label and not label.lower().startswith("labeled"):
        label = "labeled " + label
    return base_name(T) + (f", {label}" if label else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    T = testcfg.get()
    review.enable()
    if not review.enabled():
        sys.exit("the sample label and notes window needs a display")
    notes = list(T.config.get("notes") or [])
    fields = [("Label (e.g. healthy, or slightly calcified)", label_of(T))]
    fields += [(f"Note {i + 1}", n) for i, n in enumerate(notes)]
    fields += [(f"Note {i + 1}", "") for i in range(len(notes), max(4, len(notes) + EMPTY_NOTES))]
    a = review.form(f"runDIC: label and notes, {T.sample}/{T.name.split('_')[0]}",
                    f"{base_name(T)}: label, and notes, one per box (tears, re-mounting, anything unusual).", fields)
    if a is None:
        print("cancelled, nothing changed")
        return
    cfg = dict(T.config)
    cfg["sample"] = sample_name(T, a[0])
    new = [n for n in a[1:] if n]
    if new:
        cfg["notes"] = new
    else:
        cfg.pop("notes", None)
    T.write_config(cfg)
    print(f'label "{label_of(testcfg.Test(T.root))}", {len(new)} note(s) -> test_config.json')


if __name__ == "__main__":
    main()
