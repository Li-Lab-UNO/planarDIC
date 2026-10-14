# Running tests on Ubuntu, from scratch

Every command below is run from the project folder (the git clone, `~/Documents/planarDIC` here) in Ubuntu's
Terminal.

The examples use a sample **N1900_LFPA** and its **test1**; replace them with your own names.
`./runDIC.sh` runs steps 3–8 in one command; the steps are listed separately so a single one
can be rerun. Its defaults on this computer (tests, jobs, steps, ...) are in `runDIC.ubuntu.ini`.

---

## 0. One-time setup on a new computer

The code and the settings files come from the git repository; the data does not (git ignores the
sample folders). Install git, the GitHub tool and Python's venv module, sign in, and clone:

```bash
sudo apt install git gh python3-venv
```
```bash
gh auth login
```
```bash
gh repo clone davidshuenli/planarDIC ~/Documents/planarDIC
```

Then copy each sample folder into the project folder, in any way that suits,
for example a USB drive. The raw CellScale files, the `Images/` folders, the `Thickness, <ID>/` photos and
workbook, and the `test_config.json` files are what the pipeline needs. For a clean run from scratch,
leave out (or delete afterwards) every result, so the pipeline recomputes them: in each test folder
`DIC_results/`, `DIC_logs/`, `constitutive_fit/` and `scale_calibration.json`, and
`Thickness, <ID>/thickness_auto/`. For a run that derives every setting from the raw data, leave out the
`test_config.json` files too (they then hold no settings chosen on the other computer). Then, on the
Ubuntu machine, build the Python environment:

```bash
./runDIC.sh --setup
```

This builds `.venv` and installs `requirements.txt`. If `./runDIC.sh` says "permission
denied", run `chmod +x runDIC.sh` once. Later code changes arrive with `git pull`; if this
computer's copy of the code was edited, `git fetch origin && git reset --hard origin/main` replaces it
with the repository's (the data is ignored by git and stays as it is).

---

## 1. Put the CellScale files in place

Make a sample folder next to the existing ones, and a test folder inside it. **The test folder name must
start with `test` + a number + `_`** (for example `test1_10142026_PLC25_40-12`), or the code won't find it.

```
planarDIC/ (the project folder)
└── N1900_LFPA/
    ├── Thickness, 1900/             thickness photos in test 1/, test 2/, … (and workbook, if any)
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
masking come from the images, the thickness from the photos in `Thickness, 1900/test 1` (the sample's
workbook where the photos give none, or with `thickness = workbook` in `runDIC.ubuntu.ini`), excluded tine
pairs from the calibration's outlier rule, and excluded regions from the artifact check. Only the label
is yours to give (`--label`). A `test_config.json` that is already there is kept and updated, so
settings entered by hand override the rules.

---

## 2. Check the plan

```bash
./runDIC.sh --dry-run N1900/test1
```

This prints every command without running anything.

---

## 3. Sort the images and draft `test_config.json`

```bash
./runDIC.sh --steps sort,config --label "labeled healthy" N1900/test1
```

- **sort** moves the images into 19 folders, `Def01_Preconditioning_c01` … `Def19_Equibiaxial_Save`
  (`Images/sort_manifest.csv` records every move). Safe to rerun.
- **config** finds the 20 tine tips in the first Def11 image, refines them, and writes `test_config.json`
  with the tip lines, the label and a glint check. With an existing config it updates only the tine
  settings and keeps everything else.

Open `DIC_results/config_tine_tips.png` and check that every red circle sits on a
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
./runDIC.sh --steps scale,widths N1900/test1
```

Tine pairs that disagree with the rest are left out automatically and listed in the log. Open
`DIC_results/scale_calibration.png` to check: the pairs used should agree within a few
percent (right panel). To choose the pairs yourself, add them to the config, for example
`"calibration_exclude": ["Y3"]`, and rerun this step.

---

## 5. Thickness

The thickness is needed only for the stresses (step 8), but it runs here, before tracking, so its overlays
can be checked while tracking runs; it needs only the photos and the test's config. The pipeline's
`thickness` step measures the photos in `Thickness, 1900/test <N>` and writes `thickness_mm` into each
test's config, in place of the workbook value (`thickness = auto` in `runDIC.ubuntu.ini`), or only where
the workbook gives none (`thickness = workbook`, or `--thickness workbook` for one run). Its log lists the photos used and those left out, with the reason.
To run it on its own:

```bash
./runDIC.sh --steps thickness N1900
```

Check the overlays in `Thickness, 1900/thickness_auto/test<N>/`: the yellow outline should follow the tissue
around every strip, and every strip marked `*` (all, by default) should be a piece of the test's specimen
(`strip_selection.json` in the photo folder can use one strip of a photo, or leave a photo out). See
"Thickness from photos" in `README.md`.

---

## 6. Track Def11–Def19

```bash
./runDIC.sh --steps track,artifacts N1900/test1
```

This is the long step: one job per deformation, `jobs` of them at a time (`runDIC.ubuntu.ini`;
`--jobs N` for one run), each matching on several threads. With all nine deformations of a test running
at once, a 380-frame deformation takes about 4 min on the Threadripper (about 7 min on the M4 Mac). If
you're connected over SSH, run it inside `tmux` so it keeps going if the connection drops. The fits use
only the loading branches: `--loading-only` (or `loading_only = yes` in `runDIC.ubuntu.ini`) stops each
deformation at the rake peak and tracks about half the frames. Each
deformation's track goes to the test's `DIC_results/` folder as `DefNN_….npz`.

**artifacts** then looks for local tracking errors, such as a speckle particle sliding over the tissue or
a glare patch whose edges shift with load. When it confirms one, it excludes that region in
`test_config.json` and tracks the test again. Its map is `DIC_results/artifact_check.png`.

---

## 7. Figures and exports

```bash
./runDIC.sh --steps jfig,exports,edge,compare N1900
```

- **jfig**: `DIC_results/peak_J_over_mean_Def11-19.png / .csv`, peak areal stretch in all nine protocols.
- **exports**: for Def11, the CellScale comparison, field exports (`_fields.npz`, `_fields.mat`) and
  animations (`.mp4`), next to the track. `--export-defs Def11,Def15` exports more deformations.
- **edge**: `comparisons/edge_effects_<tests>.png / .csv`, how far tine effects reach.
- **compare**: `comparisons/peak_J_<samples>_w1.5mm.png / .csv`, all the tests side by side.

---

## 8. Stresses and the constitutive fit

```bash
./runDIC.sh --steps fit N1900
```

This fits every test of the sample at once (about 0.5–2 min each) and skips, with the reason, a test that
doesn't have all of Def11–Def19 tracked or has no thickness. To fit one test directly:

```bash
.venv/bin/python constitutive_fit/compare_stretch_inputs.py --test N1900/test1
```

The Cauchy stresses come from the tared rake forces, the stretches, the tine spans (`widths_mm`) and the
thickness (`thickness_mm`): t_zz = F_z·lambda_z / (w_theta·H) and t_thetatheta = F_theta·lambda_theta / (w_z·H),
over the loading branch of Def11–Def19. The four-fiber model is fitted for each stretch input (CellScale
4-point, DIC 4-point, DIC region mean, DIC affine F, DIC region edges, rake), with a protocol-level
bootstrap. Results go to the test's `constitutive_fit/`: `params_by_stretch_input.csv / .png` (parameters
and 95% intervals) and `fit_stress_stretch_<input>.png` (data and model, every protocol on one plot). See
`constitutive_fit/README.md`.

---

## Everything at once (steps 3–8)

```bash
./runDIC.sh --label "labeled healthy" N1900
```

It asks before starting (`--yes` skips the question), writes logs to `DIC_logs/run_<date>_<time>/`
in each test folder, and ends with a summary that includes each test's fit (R² and RMSE). Check
`config_tine_tips.png` and `scale_calibration.png`, the thickness overlays (ready about two minutes in,
while tracking runs) and the fit plots.

To see what each test has:

```bash
.venv/bin/python strain_tracking/status.py --notes
```
