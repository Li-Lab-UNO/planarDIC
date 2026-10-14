# CellScale biaxial DIC analysis

Analysis code for CellScale planar biaxial tests of soft-tissue samples (pulmonary artery so far).
One copy of the code serves every sample and test folder; a new sample is a new folder next to the
others, and the code finds it automatically.

- `strain_tracking/`: DIC tracking, kinematics, figures and exports.
- `constitutive_fit/`: bulk and local constitutive fits.
- `thickness/`: wall thickness from the strip photos.
- `runDIC.sh`: runs the tracking pipeline for one or more tests (see "Running the pipeline").

## Layout

```
planarDIC/
├── runDIC.sh         the pipeline in one command
├── runDIC.mac.ini, runDIC.ubuntu.ini   runDIC.sh defaults for each computer
├── analysis_settings.json  physical sizes shared by every test (subset, grid step, strain window, ...)
├── requirements.txt        Python packages (.venv)
├── NEW_TEST_UBUNTU.md      a new test, step by step, on Ubuntu
├── strain_tracking/  constitutive_fit/  thickness/
├── comparisons/            cross-test and cross-sample figures and tables (DIC_logs/: their run logs)
└── <SAMPLE>/               e.g. N1739_LFPA: one folder per sample
    ├── Thickness, <ID>/    thickness photos in test <N>/, their results in thickness_auto/test<N>/
    │                       (and the manual workbook, if any)
    └── test<N>_<date>_<protocol>/       e.g. test1_08252026_PLC25_38-20
        ├── test_config.json             this test's settings (below)
        ├── scale_calibration.json       image scale (calibrate_scale.py)
        ├── Coords1.csv, *Data.csv, *.tst, Logs/   CellScale output
        ├── DIC_logs/run_<date>_<time>/     runDIC.sh logs for this test
        ├── DIC_results/              tracks (DefNN_*.npz), exports and figures, named by deformation;
        │                                config_tine_tips, scale_calibration, artifact_check, peak_J_over_mean
        ├── constitutive_fit/            fit outputs for this test (see constitutive_fit/README.md)
        └── Images/
            ├── Def01_… / … Def19_…/     images, one folder per load–unload cycle (sort_images.py)
            └── sort_manifest.csv        old → new path for every sorted image
```

A test folder name must start with `test`, a number and `_`. Outputs always go inside the test
folder they belong to; figures that span tests go to `comparisons/`. Only Def11–Def19 are
analyzed; Def01–Def10 are preconditioning.

## Choosing a test

Every script works on one test. Select it with `--test SAMPLE/TEST`, where each part is any
unique prefix of the folder name, or with the environment variable `DIC_TEST`:

```bash
cd strain_tracking
../.venv/bin/python track_protocols.py --test N1739/test2
```

```bash
DIC_TEST=N1739/test1 ../.venv/bin/python ../constitutive_fit/baseline.py evaluate --params "…"
```

`--test N1739/test2` selects `N1739_LFPA/test2_…`. A bare test name (`--test test2`) works only
when a single sample has that test, and a path to a test folder also works. There is no default
test, so outputs can't land in the wrong folder. Scripts that take a tracking `.npz`
(`compare_cellscale.py`, `export_fields.py`, `animate_fields.py`) use the test that the file is in.
In Python, use `import testcfg; testcfg.select("N1739/test1")` before importing the other modules.

## Physical analysis settings: `analysis_settings.json`

Subset size, grid step, strain window, interior edge, low-spot box and the resampling target are
set once, in mm, in `analysis_settings.json`. Each test converts them to px with its own image scale,
so every specimen is tracked and smoothed at the same physical scale regardless of camera zoom.
See `strain_tracking/README.md` (DIC settings).

The interior edge (`edge_mm`) is a user choice. `strain_tracking/edge_effects.py` reports how far
the tine effects reach in each test against it and changes nothing; edit `edge_mm` to change it,
or pass `--edge-mm Z TH` to `peak_j_figure.py` / `compare_peak_j.py` for one run. Every J/mean
figure and .csv states the interior it used.

## Per-test settings: `test_config.json`

| Key | Meaning | Set by |
|---|---|---|
| `sample` | label used in figure titles, e.g. `"N1739 LFPA, test1 (08/25/2026), labeled healthy"` | `make_config.py --label`, or by hand |
| `tine_tips_px` | rest-state tip lines: `z` = [left, right], `theta` = [top, bottom] (px); the tracking grid, the mm origin and the interior are measured from them | `make_config.py` |
| `calibration_tines_px` | each tine's [centerline, tip] (px), per rake: starting points for `calibrate_scale.py` and `measure_tines.py`, tine positions for `edge_effects.py` | `make_config.py` |
| `calibration_exclude` | tine pairs left out of the scale calibration, e.g. `["Y3", "Y5"]`. Optional: without it, `calibrate_scale.py`'s outlier rule chooses them and records them in `scale_calibration.json` | by hand, to override the rule |
| `roi_px.standard`, `roi_px.tine_region` | grid extents [x0, y0, x1, y1] for the single-deformation scripts (`track.py`) | `make_config.py` |
| `widths_mm.tine_span` | [w_θ, w_z] in mm, the load-bearing widths in the stress: w_θ divides F_z, w_z divides F_θ; means of the opposing rakes' spans | `measure_tines.py --write` |
| `widths_mm.tine_span_by_rake` | each rake's span (mm); opposing rakes differ because the tines are bent and applied unevenly | `measure_tines.py --write` |
| `widths_mm.specimen` | [w_θ, w_z] of the whole specimen (mm), for the specimen-width sensitivity check | by hand |
| `thickness_mm` | `H`, `H_sd` (mm) and their `source`; `workbook`, the workbook value the photo measurement replaced | `measure_thickness.py --write`, from the photos (the pipeline's `thickness` step); `make_config.py` first fills it from the sample's thickness workbook (mean and sample SD of the test's readings), which stays where the photos give none or with `thickness = workbook`; or by hand |
| `reference_deformation` | protocol whose first frame is the common reference (`Def11`) | `make_config.py` |
| `dic.mask_saturated` | `true` to leave saturated glint pixels out of DIC matching | `make_config.py`: on when the tissue in the rest frame has saturated glints, off when it has none or is overexposed (≥ 10% saturated) |
| `dic.exclude_px` | rectangles [x0, y0, x1, y1] (reference-frame px) left out of DIC matching, e.g. a speckle particle that slides over the tissue or a glare patch whose edges move with load. The rectangles move with the tracked tissue; neighboring subsets track on the remaining pixels | `artifact_check.py`, for confirmed tracking artifacts; or by hand |
| `dic.note`, `notes` | why the settings above were chosen; caveats such as dropped frames or scale uncertainty | by hand |

## Status

`status.py` reads the files on disk and prints, for every test: sorted folders, config, scale (and
excluded tine pairs), tine spans, thickness, tracked deformations, J/mean figure and exports.
`--notes` adds each test's label and notes.

```bash
.venv/bin/python strain_tracking/status.py
```
```bash
.venv/bin/python strain_tracking/status.py N1739 --notes
```

## Running the pipeline: `runDIC.sh`

Runs every tracking step for the tests given, in order, with one log per step:

| Step | Script | Does |
|---|---|---|
| `sort` | `sort_images.py` | sorts the images into Def01–Def19 folders (safe to re-run) |
| `config` | `make_config.py` | finds and refines the tine tips; with an existing config (`--update`) every other setting is kept |
| `scale` | `calibrate_scale.py` | image scale from Def01 and Def19, leaving out `calibration_exclude` |
| `widths` | `measure_tines.py --write` | tine spans → `widths_mm` |
| `thickness` | `measure_thickness.py --write` | thickness from the photos in `<SAMPLE>/Thickness, <ID>/test <N>` → `thickness_mm`, in place of the workbook value (`thickness = auto`, the default), or only where there is no workbook value (`thickness = workbook`, `--thickness workbook`) |
| `track` | `track_protocols.py` | Def11–Def19 against the common reference, one job per deformation, several at a time |
| `artifacts` | `artifact_check.py` | finds local tracking errors (a sliding speckle particle, a shifting glare patch) at the Def11 peak, confirms each against independent patch matching, excludes it in `test_config.json` and tracks that test again |
| `jfig` | `peak_j_figure.py` | peak J/mean figure per test |
| `exports` | `compare_cellscale.py`, `export_fields.py`, `animate_fields.py` | CellScale comparison, field exports and animations for Def11 |
| `edge` | `edge_effects.py` | tine-effect check across the tests run |
| `compare` | `compare_peak_j.py` | cross-test J/mean comparison of the samples run |

```bash
./runDIC.sh N1739
```
```bash
./runDIC.sh --steps track,jfig --jobs 6 N1739/test1 N1739/test2
```
```bash
./runDIC.sh --dry-run N1739
```

Each computer's defaults live in its settings file: `runDIC.mac.ini` on macOS and
`runDIC.ubuntu.ini` on Linux (tests, steps, jobs, deformations, status interval, loading-only,
confirmation, thickness source). `runDIC.sh` picks the file from the operating system; both are tracked by git, so
they can be edited on either computer and committed like code. Options on the command line override it
for one run (`--full-cycle` and `--ask` undo `loading_only = yes` and `yes = yes`); without tests on the
command line, its `tests` are run. `--config FILE` reads another settings file.

A bare sample name (`N1739`) runs all its tests. The script prints the plan and asks before
starting (`--yes` skips the question), writes each test's logs to `<test folder>/DIC_logs/run_<date>_<time>/` (one file per step, plus
`runDIC.log` with everything the script printed; the cross-test steps log to `comparisons/DIC_logs/`), skips the later
steps of a test whose step failed, and prints a summary. `--jobs` sets how many tracking or export
jobs run at once (default: `jobs` in the settings file; `auto` = CPU count − 1); `--defs` and
`--export-defs` choose the deformations;
`--setup` builds `.venv` first. `./runDIC.sh --help` lists every option. The parallel jobs share
the CPUs: each job gets (CPU count ÷ jobs) threads (`DIC_THREADS`), matches chunks of its subsets in
parallel on them and reads the next image while matching; OpenCV's and NumPy's own threading are limited
to match. A script run on its own uses every CPU. The subsets are matched independently, so the results
do not depend on the thread count (results from two computers agree to within rounding, far below
the DIC noise, but are not bit-identical). Tracking takes most of the time: with the nine deformations of
a test running at once, a 380-frame deformation takes about 7 min on an M4 Mac (10 cores) and about
4 min on a Threadripper 7970X (32 cores), so a whole test takes about as long as its longest deformation.

`--loading-only` tracks each deformation only to the end of loading (the frame where the rakes peak),
about half the frames. The fits use only the loading branch, and tracking runs forward, so the
loading frames come out the same as in a full-cycle track. Animations and CellScale comparisons then
cover loading only, and the unloading branch (closure back toward the rest state) is not tracked. The
track records `branch`, and `strain_tracking/status.py` marks loading-only tracks.

## Setting up a new test

`NEW_TEST_UBUNTU.md` walks through these steps command by command.

1. **Copy the CellScale output** into `<SAMPLE>/test<N>_…/`: `*Data.csv`, `Coords1.csv`, the `.tst` file,
   `Logs/` (with `*_PhaseTiming.csv`) and the images in `Images/`.
2. **Sort the images:** `sort_images.py --test SAMPLE/testN --dry-run`, then without `--dry-run`.
3. **Draft `test_config.json`:** `make_config.py --test SAMPLE/testN --label "labeled healthy"` finds the
   20 tine tips in the Def11 rest frame, refines them with the calibration's tip finder and writes the
   tip lines, ROIs and label. Check `DIC_results/config_tine_tips.png`. If its
   self-check fails it writes nothing: enter the flagged rake's tines in a config by hand (within about
   20 px) and run it again with `--update`, which refines them.
4. **Calibrate the scale:** `calibrate_scale.py --test SAMPLE/testN` (Def01 and Def19). Tine pairs that
   disagree with the rest are left out by its outlier rule and listed; check `scale_calibration.png`.
   `calibration_exclude` in the config overrides the rule.
5. **Measure the tine spans:** `measure_tines.py --test SAMPLE/testN --write`.
6. **Track:** `track_protocols.py --test SAMPLE/testN` (Def11–Def19), then `artifact_check.py --test
   SAMPLE/testN`: if it adds an excluded region (exit status 3), track the test again.
7. **Figures and exports:** `peak_j_figure.py --test SAMPLE/testN`, then `compare_cellscale.py`,
   `export_fields.py` and `animate_fields.py` on a track `.npz`.
8. **Thickness** (for stress and fitting only): `measure_thickness.py --recursive --write "<SAMPLE>/Thickness, <ID>"`
   measures the photos and writes each test's thickness, in place of the workbook value `make_config.py`
   entered (`--keep-workbook`: workbook values stay, the photos fill only tests without one). Check its overlays.

Steps 2–8 are what `runDIC.sh` runs (the thickness as its `thickness` step). Constitutive fitting is separate: see
`constitutive_fit/README.md`.

## Thickness from photos: `thickness/measure_thickness.py`

Measures the wall thickness of the strips in the `Thickness, …` photos, which is the width of the
tissue band, at stations every 0.1 mm along its centerline (10% trimmed at each end).
The scale comes from the background: the 0.5 in grid of the blue mat, or the mm edge of the steel
ruler (told from the inch edge automatically). When a Petri dish is in the frame, the grid spacing is
read from the mat outside it (the grid under a dish is often worn and its rim confuses the spacing).
Strips cut by the image border, and round blobs such as glare, are skipped. Loose fibers thinner
than about 0.5 mm are not counted as wall.

```bash
.venv/bin/python thickness/measure_thickness.py "N1739_LFPA/Thickness, 1739/test 1"
```

Writes `<SAMPLE>/Thickness, <ID>/thickness_auto/test<N>/` for the photos in `test <N>/`: an overlay per image (mask, centerline, chords),
per-station CSVs and `summary.csv` (mean, SD, median, IQR per strip). When an image holds several
strips, each is reported and `selected` marks the ones used: all of them, since the strips in a test's
photos are pieces of its specimen, unless `strip_selection.json` says otherwise
(`{"DSC_1234.JPG": "upper"}` uses one strip, picked by position: upper, lower, left, right or center;
`"exclude"` leaves a bad photo out). Strip-shaped regions with less than `STRIP_CONTRAST_FRAC` of the
strongest strip's contrast against the mat (stains, films) are dropped and listed. `scale_ok = False`
(printed as `(?)`) means the two grid directions disagree by more than 5%; check the overlay's scale
bar, and rerun that image with `--px-per-mm` if it is wrong.

`--write` puts the result into each test's `test_config.json` as `thickness_mm`: H and H_sd are the mean
and SD of all stations of the selected strips of the test's photos with a reliable scale (the folder
`test <N>` or `test<N>` is matched to the test folder `test<N>_…`). One such photo is enough; the log
lists the photos used and those left out with the reason (excluded, failed, no complete strip, scale
uncertain), or says that the test got no thickness (the config then keeps what it has, such as the
workbook value). A workbook value it replaces is kept in `thickness_mm.workbook`, and the log shows the
difference. With `--keep-workbook`, workbook values stay in use (one kept in `thickness_mm.workbook` is
put back) and the photos fill only tests without one. The pipeline's `thickness` step runs it for each
test (`thickness = auto`, the default, or `workbook` in the settings file).

On the steel ruler the tissue is told from the steel by brightness alone (wet films on the steel look
warm), with the edge `RULER_EDGE_FRAC` of the way from the steel's brightness to the tissue's.
On the blue mat, pixels both warmer and brighter than the mat are the confident tissue; each strip is
then grown (`grow_on_mat`) into the rest of its tissue, judged against the mat around it: pale layers
and blurred edges (the edge half-way through the blur), and a darker but warm shaded side. Growth
stays out of the mat's printed grid lines (found as straight warm runs along the grid directions),
joins pieces of one strip whose gap is tissue-colored all along (a faint middle in a blurred photo),
and never joins far-apart parts of a strip (the ends of a C stay apart). Isolated small pieces such
as printed numerals stay out. Blurred photos still give a less certain edge: check their overlays.

`thickness/compare_manual.py` compares the center strips with the workbooks' manual readings
(same two images per test) and writes `comparisons/thickness_auto_vs_manual*.csv` and `.png`.

## Environment

Python 3.9–3.12. Build `.venv` in the project folder (each computer needs its own; a copied `.venv`
from another operating system doesn't work):

```bash
./runDIC.sh --setup
```

which runs `python3 -m venv .venv` and `.venv/bin/python -m pip install -r requirements.txt`. On
Ubuntu (22.04 or 24.04), install the venv module first:

```bash
sudo apt install python3-venv
```

Use `.venv/bin/python -m pip`, not the `pip` script, which stops working if the folder is moved.
Figures are written to files (matplotlib Agg backend), so no display is needed, and the MP4 encoder
comes with the `imageio-ffmpeg` package. Linux file names are case sensitive; the code uses the
CellScale names as written (`Images/`, `Coords1.csv`, `*Data.csv`).

### Progress output

Tracking (`track.py`, `track_protocols.py`), `export_fields.py` and `animate_fields.py` show a progress
line for the current loop: frames done, percent, elapsed time, ETA and, for tracking, the valid fraction
and median ZNCC. `track_protocols.py` also numbers the deformations (`[3/9] Def13_…`). In a terminal the
line updates in place; when output goes to a file (`nohup`, `> log.txt`, `| tee`) a full line is
written every 30 s (`DIC_PROGRESS_EVERY` sets the interval). `DIC_PROGRESS=lines` forces full lines and
`DIC_PROGRESS=0` turns the line off.

`runDIC.sh` tracks each deformation as its own job with its own log, so the per-frame lines go to
the logs. While tracking runs it prints a status line every second (`--status N` for every N s; 0 turns
it off) with the deformations finished and the percent done of each running one; each job's log is
updated as often (at most every 10 s). It prints the total tracking time when tracking ends, and the total
runtime of the run at the end of the summary, both on screen and in each test's `runDIC.log`.
