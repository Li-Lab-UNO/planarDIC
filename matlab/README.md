# runDIC for MATLAB: stress–stretch data from CellScale biaxial tests

`runDIC` turns the raw output of a CellScale BioTester planar biaxial test into stress–stretch data,
ready for model fitting. It measures the wall thickness from photos of the strips, sorts the images,
finds the tine tips, calibrates the image scale, measures the specimen widths, tracks the tissue with
digital image correlation (DIC), checks the tracking for local errors and writes one row per image
frame of every loading curve: CSV files, and an Excel file in the layout of CellScale's own AllData
export. A tracking report shows the tracking of every protocol at the end.

It runs on Windows, macOS and Linux with **MATLAB R2021a or later and no toolboxes**.

For a first run, print `QUICKSTART.pdf`: every window in order, with screenshots and what to answer.

## Setting up (once)

1. Copy this `matlab` folder anywhere on your computer, for example `Documents\runDIC`.
2. Open MATLAB and go to that folder (the folder bar at the top, or `cd Documents\runDIC`).

## The data

Each sample has a folder holding one folder per test. A test folder's name must start with `test`, a
number and `_` (for example `test1_10142026_PLC25_40-12`), and it holds the CellScale files as the
tester wrote them:

```
N1900_LFPA\                          sample folder
├── Thickness, 1900\test 1\          photos of the test's strips for the thickness (see "Thickness")
└── test1_10142026_PLC25_40-12\      test folder
    ├── ID_…Data.csv                 forces and rake positions
    ├── ID_….tst                     the protocol
    ├── Coords1.csv                  CellScale's 4-point tracking
    ├── Logs\…_PhaseTiming.csv       start and end of every phase
    └── Images\                      all the .jpg images (in any subfolders)
```

## Running

In the MATLAB Command Window:

```matlab
runDIC
```

A window asks for a folder. Choose a **test folder** to analyze that test, or a **sample folder** to
analyze all its tests one after another. A second window then lists the steps, marks the ones already
done, and runs the ones you pick (all of them the first time). Then a window asks for the **label** of
each test (for example healthy) and **notes** (tears, re-mounting, anything unusual), filled in with
what was saved before; they go into `test_config.json` and the results' `README.txt`. Then the wall
thickness is measured from the test's photos (see "Thickness"). The folder window opens where the
folder of the last run is, also after MATLAB is restarted. You can also give the folder directly, which
runs the default steps without the step window:

```matlab
runDIC('C:\Data\N1900_LFPA\test1_10142026_PLC25_40-12')
```

A test takes about 15 minutes on a recent 10-core computer, most of it tracking, and longer on a
laptop with fewer cores. Progress is printed in the Command Window.

## What to check

After the images are sorted, `runDIC` shows a table of the folders Def01–Def19 with the number of
images in each and the frames missing from the tester's files: images, rows of `Data.csv` (forces,
rake positions) and rows of `Coords1.csv` (CellScale's 4 points), each a gap in their times. It then
**asks whether to proceed**. Folders missing anything are marked red. A few missing frames are normal
(the camera drops frames now and then); a folder with no images, or many missing, means the test or its
folders need a look before going on. "No, stop this test" stops that test only. The table is saved as
`Images\sort_summary.csv`.

When the tine tips have been found, `runDIC` shows them and asks **whether every red circle is on the
tip of a tine**. If not, answer "No, fix some rakes", pick the rakes that are wrong and click their tips
(see below).

Next, `runDIC` shows the rake forces of Def11–Def19 with each **loading branch** shaded: the part
exported, from the start of stretching to the force peak. Each band should start where the forces
begin to rise and end at their peak. If one doesn't, answer "No, fix some", pick the protocols and click
where loading starts and where it ends; a right-click (or Backspace) removes the last click, and A goes
back to the automatic branch. The run then proceeds by itself.

After a run, look at the pictures in the test's `DIC_results` folder:

- **`config_tine_tips.png`**: every red circle should sit on the tip of a tine, and the yellow
  rectangle should join the tips.
- **`loading_branches.png`**: the loading branches of all 19 deformations (Def01–Def10 are only shown,
  not asked about; red: set by hand).
- **`scale_calibration.png`**: in the right panel, the tine pairs should agree within a few percent.
  Pairs that disagree are left out automatically and listed in the Command Window.
- **`artifact_check.png`**: local tracking errors. A speckle particle sliding over the tissue or a
  glare patch whose edges shift with load is found automatically (green circle), left out of the
  tracking and the test tracked again; the Command Window says so. Blue circles are spots that were
  checked and found to be real.

If the tine tips can't be found automatically on one side, `runDIC` shows the image and asks you to
**click the 5 tine tips of that rake**. Click the end of each rod, where it meets the tissue; a red
circle marks each click, and a right-click (or Backspace) removes the last one. The clicks are refined
automatically. If the tips were found but some red circles are off, click all 20 yourself and redo the
later steps:

```matlab
runDIC(folder, 'tine_tips', 'click')
```

## The tracking report

At the end of a run, `runDIC` opens the **tracking report** of each test it tracked: the analyzed
deformations (Def11–Def19) side by side, each with the tracked points over the images (green: tracked,
red: failed in that frame and filled in from the neighbors), every 4th grid row and column, and the
stretches, the share of valid points and the match quality (ZNCC, 1 is a perfect match) of the frame
shown. One slider moves all of them through their loading branch together; Play runs them. The
preconditioning cycles (Def01–Def10) are tracked too, for the AllData file; they are not shown, but
their numbers are in `DIC_results\tracking_report.csv` and the Command Window says if any of them
tracked poorly.

Each deformation is marked **Good**, **Check** (amber) or **Problem** (red) from its worst frame. If any
is marked, the report opens at the worst frame of the worst one and says so at the top. A few red
crosses that come and go are normal; a patch of red that persists, or grid lines that kink, is a
tracking problem. Show a report with Check or Problem marks to your supervisor before using those
results. "Open large" shows the selected deformation in a big viewer.

The report reads the tracks already saved, so it can be opened again at any time without tracking:
run `runDIC`, pick the test (or sample) folder and choose only the step `report`, or type

```matlab
viewTracking
```

## The results

The results are in the test's `stress_stretch` folder:

| File | Contents |
|---|---|
| `all_protocols.csv` | every loading curve, one row per image frame |
| `Def11_Equibiaxial.csv`, … `Def19_….csv` | the same, one file per protocol |
| `stress_stretch.png` | every protocol on one plot, to check the data at a glance |
| `README.txt` | what the columns are; the label and notes, widths, thickness, scale, tracking verdict and any loading branches set by hand |
| `stress-stretch test<N>.AllData.xlsx` | the loading branches of Def01–Def19 in CellScale's AllData layout (see below) |

Each row is one image frame of a **loading curve**, from the start of stretching to the peak. The
stretches are relative to one reference state for every protocol: the first image of Def11. Def01–Def10
are preconditioning and are only in the AllData file (below).

| Column | Meaning |
|---|---|
| `protocol`, `protocol_name` | e.g. `Def12`, `1-0.75_ratio` (axial : circumferential) |
| `time_s` | test time (s) |
| `force_z_mN`, `force_theta_mN` | rake forces, set to zero at the start of each loading curve |
| `lambda_z`, `lambda_theta` | stretches from the DIC (axial z, circumferential theta) |
| `tension_z_N_per_m`, `tension_theta_N_per_m` | membrane tension: force × stretch / width |
| `stress_zz_kPa`, `stress_thetatheta_kPa` | Cauchy stress: tension / thickness |
| `lambda_z_cellscale`, `lambda_theta_cellscale` | stretches from CellScale's own 4-point tracking |
| `lambda_z_rake`, `lambda_theta_rake` | stretches from the rake separation |
| `valid_points_pct` | share of the DIC points tracked in that frame |

Axial (z) is the tester's X direction and circumferential (theta) its Y direction. The widths are the
tine spans measured from the images: the axial force is divided by the circumferential width and the
circumferential force by the axial width. Use `lambda_z` and `lambda_theta` (the DIC stretches) for
fitting; the CellScale and rake stretches are there for comparison.

### The AllData file

`stress-stretch test<N>.AllData.xlsx` (in `stress_stretch`) has the loading branch of every deformation,
Def01–Def19 (preconditioning included, for checking), in the layout of CellScale's AllData export, one
row per image frame. As in CellScale's export, each cycle is its own reference: the stretches are 1 and
the forces are set to zero at the first frame of the cycle.

| Columns | Meaning |
|---|---|
| `SetName` … `Temperature` | as in the tester's Data.csv (the last frame or two of a loading branch can carry the next phase's name, `1-Recover`: the log lags the motion) |
| `X1` … `Y4` | CellScale's 4 tracked points (px) |
| `lx`, `ly`, `kx`, `ky` | deformation gradient from the DIC: F11, F22, F12, F21 (x axial, y circumferential) |
| `txx`, `tyy`, `txy`, `tyx` | Cauchy stresses (kPa) with the DIC deformation, the tine spans and the thickness: force × lx / (width × thickness), … |
| `ShearAngleDeg` | angle between the deformed axes (90 undeformed) |
| `lx_4pt` … `tyx_4pt` | the same from CellScale's 4 points, with the rake separations at the start of the test as widths: CellScale's own values, for comparison |

The first 27 columns have the names and order of CellScale's export, so a script written for it reads
this file unchanged.

## Thickness

The stresses need the wall thickness. The step `thickness`, selected by default, measures it from the
test's photos (below); when the test already has a thickness, it first asks whether to keep it or
measure again. For a test without photos, `runDIC` opens a window at the start to type the thickness
(mm); leave it blank to skip that test (its stress columns are then `NaN`). The thickness is saved in
the test's `test_config.json`, with where it came from.

### Measuring from the photos

The photos of a test's strips, lying on their cut face on the blue cutting mat or the steel ruler, are
in the sample folder's `Thickness, <ID>\test <N>` folder (for example `Thickness, 1900\test 1`). Each
photo is opened in turn; scroll to zoom where the pointer is.

1. **Scale, part 1.** Click on one grid line and then on the next one, straight across one square (or
   across one ruler tick), and type the size of that square with its unit (a mat square is 0.5 in). The
   size is remembered for the next photo.
2. **Scale, part 2.** Click on a grid line and then on another one as far away as the photo allows,
   close to the strip. The squares in between are counted for you; confirm the count. If the ends are
   not on grid lines, you are asked to draw the line again.
3. **Frame the strip.** Click two opposite corners of a box around the strip to measure; the view then
   stays on it ("Back to the frame" returns to it). Measure complete strips only and leave out pieces
   cut off by the edge of the photo. "No frame" keeps the whole photo.
4. **Centerline.** Click along the middle of a strip, end to end (two clicks for a straight strip, more
   where it bends), and press Done. The edges are found automatically.
5. **Check.** The strip is shown with a red line across it at every position measured, edge to edge.
   If the red lines end on the edges of the strip, use the strip. Orange positions are far from the
   rest (a loose fiber, a cut end) and are left out automatically. Otherwise: **Remove a section** (click
   its two ends), **Redraw the centerline**, **Click the edges myself** (click along one edge, then the
   other), or **Skip this strip**.

The cursor is a crosshair over the photo. A double-click (or Enter) on the last point of a step proceeds,
as the step's first button does.

"Yes, use this strip" proceeds to the next photo. If a photo holds a second strip of the test, "Use it,
and frame another strip in this photo" frames the next one with the same scale. A strip needs at least
10 positions. The test's thickness is the mean of every position used, all strips of all photos
together, and is saved when the last photo is done; closing the window instead keeps the thickness
saved before (with none saved, a box asks to type one). The positions, a summary per strip (with its
frame) and a picture of each photo are saved in the test's `thickness` folder.

### Where the thickness is saved

The thickness is saved in the test's `test_config.json`, the same entry the Python version of this
pipeline uses:

```
"thickness_mm": {"H": 1.237, "H_sd": 0.251, "source": "measured from 3 photo(s), 3 strip(s), 349 positions"}
```

To use another thickness for one run, give it when running the test:

```matlab
runDIC('C:\Data\N1900_LFPA\test1_10142026_PLC25_40-12', 'thickness_mm', 1.25)
```

Without a thickness, the stress columns are `NaN`; everything else is written. To add or change a
thickness later without tracking again, run only the thickness and export steps:

```matlab
runDIC('C:\Data\N1900_LFPA', 'steps', {'thickness', 'export'})
```

## Running only some steps

The steps are, in order: `info`, `thickness`, `sort`, `tines`, `scale`, `widths`, `branches`, `track`,
`artifacts`, `export`, `report`. To run some of them:

```matlab
runDIC(folder, 'steps', {'track', 'export'})
```

Each step uses what the earlier ones wrote, so after changing something, run that step and the ones
after it.

## Settings

`DIC_settings.m` holds the settings, each explained in the file: the steps, the thickness, the
deformations, the DIC sizes (in mm), the limits of the tracking report and how often the tracking
progress line is updated. The defaults suit these tests. `track_deformations` (Def01–Def19) are
tracked; `deformations` (Def11–Def19) are analyzed: the CSVs, the loading-branch check and the
tracking report.

`review_sort = false`, `review_tines = false`, `review_branches = false`, `review_tracking = false` and
`ask_thickness = false` turn off the sorting check, the tine-tip question, the loading-branch question,
the tracking report window and the thickness window, for example when running many tests unattended;
leave `info` and `thickness` out of the steps then too.

The tracking needs about 1 GB of memory.

## Files written

| Where | What |
|---|---|
| `Images\Def01_…` … `Def19_…` | the sorted images (`sort_manifest.csv`: every move; `sort_summary.csv`: the completeness of every cycle) |
| `test_config.json` | sample label and notes, thickness, tine tips, tracking region, widths, excluded regions |
| `scale_calibration.json` | image scale (px per mm) |
| `loading_branches.json` | loading branches set by hand (start and end, s) |
| `DIC_results\DefNN_….mat` | the DIC tracks (point positions in every frame) |
| `DIC_results\*.png` | the pictures to check |
| `DIC_results\tracking_report.csv` | tracking quality of every deformation (its verdict is also in `stress_stretch\README.txt`) |
| `stress_stretch\` | the CSV files and the AllData file |
| `thickness\` | thickness measured from the photos: positions, summary per strip, a picture per photo |

`test_config.json` and `scale_calibration.json` are the same files the Python version of this pipeline
(`runDIC.sh`) writes, so either version can pick up a test the other one calibrated.
