# Running tests on Ubuntu, from scratch

Every command below is run from the `DIC` folder in Ubuntu's Terminal:

```bash
cd ~/DIC
```

The examples use a sample **N1900_LFPA** and its **test1**; replace them with your own names.
`./run_pipeline.sh` runs steps 3–6 in one command; the steps are listed separately so a single one
can be rerun.

---

## 0. One-time setup on a new computer

Copy the `DIC` folder across without its `.venv` (that one belongs to the other computer). For a
clean run from scratch, also leave out every result, so the pipeline recomputes them. From the
computer that has the project (here to an Ubuntu machine over SSH; a USB drive path works the same
way), with `N1739_LFPA` as the sample:

```bash
rsync -av --exclude .venv --exclude pipeline_logs --exclude comparisons --exclude __pycache__ --exclude .DS_Store --exclude strain_results --exclude scale_calibration.json --exclude 'peak_J_over_mean_*' --exclude 'Images/constitutive_fit' --exclude thickness_auto --exclude '*.pptx' --include 'N1739_LFPA/***' --exclude 'N1[0-9][0-9][0-9]_*' ~/Documents/DIC/ user@ubuntu-machine:~/DIC/
```

This copies the code, the settings, and the sample's raw CellScale files, images and
`test_config.json` files. For a run that derives every setting from the raw data, add
`--exclude test_config.json` (the configs then hold no settings chosen on the other computer). Keep the order of the options: rsync uses the first rule that matches, so
the result folders are excluded before the sample is included, and the other samples after. Add
`-n` to see the file list without copying anything. Then, on the Ubuntu machine:

```bash
sudo apt install python3-venv
```
```bash
./run_pipeline.sh --setup
```

This builds `.venv` and installs `requirements.txt`. If `./run_pipeline.sh` says "permission
denied", run `chmod +x run_pipeline.sh` once.

---

## 1. Put the CellScale files in place

Make a sample folder next to the existing ones, and a test folder inside it. **The test folder name must
start with `test` + a number + `_`** (for example `test1_10142026_PLC25_40-12`), or the code won't find it.

```
DIC/
└── N1900_LFPA/
    ├── Thickness, 1900/             thickness photos (and workbook, if any)
    └── test1_10142026_PLC25_40-12/
        ├── ID_…Data.csv             from CellScale
        ├── ID_….tst                 from CellScale (the protocol; gives the programmed stretches)
        ├── Coords1.csv              from CellScale
        ├── Logs/                    from CellScale; must contain …_PhaseTiming.csv
        └── Images/                  all the CellScale .jpg files, in any subfolders
```

From here on, `--test N1900/test1` (or just `N1900` for all its tests) selects this test. Any unique
start of the folder names works.

Nothing else is needed: the pipeline derives every setting from these files. The tine tips and glint
masking come from the images, the thickness from the sample's workbook (if there is one), excluded tine
pairs from the calibration's outlier rule, and excluded regions from the artifact check. Only the label
is yours to give (`--label`). A `test_config.json` that is already there is kept and updated, so
settings entered by hand override the rules.

---

## 2. Check the plan

```bash
./run_pipeline.sh --dry-run N1900/test1
```

This prints every command without running anything.

---

## 3. Sort the images and draft `test_config.json`

```bash
./run_pipeline.sh --steps sort,config --label "labeled healthy" N1900/test1
```

- **sort** moves the images into 19 folders, `Def01_Preconditioning_c01` … `Def19_Equibiaxial_Save`
  (`Images/sort_manifest.csv` records every move). Safe to rerun.
- **config** finds the 20 tine tips in the first Def11 image, refines them, and writes `test_config.json`
  with the tip lines, the label and a glint check. With an existing config it updates only the tine
  settings and keeps everything else.

Open `Images/Def11_*/strain_results/config_tine_tips.png` and check that every red circle sits on a
tine tip. If the step failed with `SELF-CHECK FAILED`, the message names the rake. Copy a config from
another test, enter that rake's tines by hand (open the image in GIMP, `sudo apt install gimp`, which
shows the pixel position under the mouse; within about 20 px is enough), and rerun the step: it
refines them.

| Rake | `horizontal` | `out` | Each tine is written as |
|---|---|---|---|
| `top` | `false` | `-1` | `[x, y]`: tine x, then tip y |
| `bottom` | `false` | `1` | `[x, y]`: tine x, then tip y |
| `left` | `true` | `-1` | `[y, x]`: tine y, then tip x |
| `right` | `true` | `1` | `[y, x]`: tine y, then tip x |

---

## 4. Scale and tine spans

```bash
./run_pipeline.sh --steps scale,widths N1900/test1
```

Tine pairs that disagree with the rest are left out automatically and listed in the log. Open
`Images/Def01_…/strain_results/scale_calibration.png` to check: the pairs used should agree within a few
percent (right panel). To choose the pairs yourself, add them to the config, for example
`"calibration_exclude": ["Y3"]`, and rerun this step.

---

## 5. Track Def11–Def19

```bash
./run_pipeline.sh --steps track,artifacts N1900/test1
```

This is the long step: one job per deformation, `--jobs N` at a time (default: CPU count − 1).
Each deformation takes roughly 5–30 min. If you're connected over SSH, run it inside `tmux` so it
keeps going if the connection drops. Each `Images/DefNN_…/strain_results/` folder gets a
`DefNN_…_tine_region.npz` track.

**artifacts** then looks for local tracking errors, such as a speckle particle sliding over the tissue or
a glare patch whose edges shift with load. When it confirms one, it excludes that region in
`test_config.json` and tracks the test again. Its map is `Images/Def11_…/strain_results/artifact_check.png`.

---

## 6. Figures and exports

```bash
./run_pipeline.sh --steps jfig,exports,edge,compare N1900
```

- **jfig**: `Images/peak_J_over_mean_Def11-19.png / .csv`, peak areal stretch in all nine protocols.
- **exports**: for Def11, the CellScale comparison, field exports (`_fields.npz`, `_fields.mat`) and
  animations (`.mp4`), next to the track. `--export-defs Def11,Def15` exports more deformations.
- **edge**: `comparisons/edge_effects_<tests>.png / .csv`, how far tine effects reach.
- **compare**: `comparisons/peak_J_<samples>_w1.5mm.png / .csv`, all the tests side by side.

---

## Everything at once

```bash
./run_pipeline.sh --label "labeled healthy" N1900
```

It asks before starting (`--yes` skips the question), writes logs to `pipeline_logs/run_<date>_<time>/`
in each test folder, and ends with a summary. Check `config_tine_tips.png` and `scale_calibration.png` afterwards.

To see what each test has:

```bash
.venv/bin/python strain_tracking/status.py --notes
```

---

## Thickness (needed only for stress and fitting)

```bash
.venv/bin/python thickness/measure_thickness.py "N1900_LFPA/Thickness, 1900/test 1"
```

Review the overlays in `thickness_auto/`, then add the result to `test_config.json`, for example
`"thickness_mm": {"H": 1.85, "H_sd": 0.30, "source": "measure_thickness.py, 3 photos"}`. See
"Thickness from photos" in `README.md`. Constitutive fitting is separate: see `constitutive_fit/README.md`.
