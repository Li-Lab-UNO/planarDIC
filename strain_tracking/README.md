# Strain tracking for CellScale biaxial sequences

Incremental 2D subset DIC (IC-GN, affine subsets, ZNSSD, per-point reference updating) on the
5 Hz CellScale image stream. Results are compared against the CellScale 4-point tracking
(`Coords1.csv`) and the rake strain (`*Data.csv`).

Every script needs a test: `--test SAMPLE/testN` (e.g. `--test N1739/test2`), or the test that a given `.npz` is in (see the
top-level `README.md`). Test-specific settings come from that test's `test_config.json`.
`./runDIC.sh` (in the project folder) runs the scripts below in order for one or more tests.
The constitutive fits are in `../constitutive_fit/`.

## Axes and naming

| Axis | Image direction | Sign |
|---|---|---|
| Axial **z** | image x (tester X) | increases to the right |
| Circumferential **θ** | image y (tester Y) | increases **upward** (θ = −image y) |

Tracking works in image coordinates. `dic2d.to_specimen` converts F to the specimen frame,
and `dic2d.THETA_SIGN` holds the sign. Only the signs of F_zθ, F_θz, E_zθ, u_θ and φ_1
depend on this convention; stretches and J do not. Every export records its convention
in `axis_convention`.

Exports use ASCII names that are valid MATLAB identifiers:

- `u_z`, `u_theta` (and `_mm` versions)
- `F_zz`, `F_ztheta`, `F_thetaz`, `F_thetatheta`, where F[i][j] = ∂x_i/∂X_j
- `lambda_z`, `lambda_theta`
- `E_zz`, `E_thetatheta`, `E_ztheta`
- `force_z_mN`, `force_theta_mN`
- `phi_1`: angle of the major principal stretch from the axial direction

`X`, `Y`, `x` and `y` are image pixel positions. `calibrate_scale.py` keeps the tester's
X/Y names because it works in machine axes.

## Tracking results

Results are written to `<test>/DIC_results/`, named by deformation. Def11–Def19 are tracked
against one common reference, the first frame of the test's `reference_deformation` (Def11).

The default grid is the **full domain** between the tine tips (`roi_px.tine_region`, files
`DefNN_*.npz`, named after the deformation folder). The tracked region then covers the whole load path between the
rakes. That's what the local estimation and the primary bulk baseline need.

The **standard region** (`roi_px.standard`) is the full domain without its outer two rows and
columns. `dic2d.load_track(images, folder, "standard")` slices it from the full-domain track
when there's no separately tracked file (`DefNN_*.npz`). `compare_stretch_inputs.py --tracks
standard` and `validate_bands.py` use this. On N1739 test1, the sliced and the separately tracked
standard grids agree to a median of 0–0.1 px, and the stretch table differs by at most 6e-4.
`track_protocols.py --grid standard` tracks it separately when a frozen baseline needs it.

The full-domain grid is a lattice of the grid step centered between the rest-state tine tips
(`tine_tips_px`), at least 5 px inside them (`dic2d.tip_lattice_roi`), laid out on the resampled
image. `roi_px` in `test_config.json` is used only by the single-deformation scripts (`track.py`).

`track.py` tracks one deformation against its own first frame. It refuses to overwrite a
common-reference track; use `track_protocols.py` for Def11–Def19.

## Scripts

| File | Purpose |
|---|---|
| `testcfg.py` | Test selection (`--test SAMPLE/TEST`, `DIC_TEST`, or the test containing a given file) and per-test settings from `test_config.json`. Shared with `../constitutive_fit/` |
| `dic2d.py` | Core: IC-GN matcher, sequence tracker, deformation gradient field, axis convention and export names. `load_scale()` reads the selected test's `scale_calibration.json` |
| `track.py` | Tracks one deformation (or a `--t0`/`--t1` window) against its own first frame → `DIC_results/<name>[_tag].npz` (`--grid standard`: `<name>_standard[_tag].npz`) |
| `track_protocols.py` | Tracks Def11–19 (default) against the common reference: each protocol's first frame is matched to the reference frame, then the protocol is tracked from there. `--grid tine_region` (default) → `DefNN_*.npz`; `--grid standard` → `DefNN_*_standard.npz`. `--loading-only` stops each protocol at the end of loading, the rake peak (`branch` in the npz). `X0` = the common reference grid |
| `track_tine_region.py` | Same as `track_protocols.py` with the default grid; other scripts import its `ROI`, `TINE_TIPS` and `SUFFIX` |
| `compare_cellscale.py` | Strain vs time against CellScale and rake, point disagreement, peak strain maps, CSVs (written next to the .npz) |
| `export_fields.py` | Writes `*_fields.npz`, `*_fields.mat` (MATLAB) and `*_fields_peak.png`. Per frame, on the reference grid: displacement, F, λ_z, λ_θ, principal stretches and φ_1, J, Green–Lagrange strain and rake forces |
| `animate_fields.py` | MP4s of displacement, stretch and strain on the deformed image (`--quantity displacement stretch strain`, default all). The lower panels show ROI means against time and force. The displacement video instead shows the change in ROI length and width against the change in rake separation. Run `compare_cellscale.py` first for the CellScale and force traces |
| `calibrate_scale.py` | Image scale from tine-tip motion vs reported rake separation (default `--def Def01 Def19`), starting from `calibration_tines_px` → `<test>/scale_calibration.json`. Also writes a diagnostic `scale_calibration.png` to the first deformation's `DIC_results/`. Works from the images, with no DIC tracking. The result is the mean over the middle three tine pairs of both axes. A pair is left out when it is more than 6% from the median of all pairs, when Def01 and Def19 disagree by more than 5%, or when its tips weren't tracked throughout (`outliers()`); the camera's pixels are square, so such a pair has tip tracking errors or bent tines. The excluded pairs and the reason are printed and recorded (`excluded`, `excluded_by`). `calibration_exclude` in `test_config.json` or `--exclude Y3 Y5` replaces the rule; `--no-auto` turns it off |
| `artifact_check.py` | Finds local tracking errors at the reference protocol's peak: clusters of at least 6 grid points whose displacement departs from their 9 × 9 neighborhood median (peak ≥ 100 µm, at least 1 mm from every tine line), each confirmed by matching small speckle patches between the rest and peak frames independently of the tracker (median disagreement ≥ 3 px and ≥ 3 × that in a ring around it). A confirmed cluster's bounding box is added to `dic.exclude_px` (exit status 3: track the test again); `--report-only` changes nothing → `artifact_check.png` in the reference protocol's `DIC_results/`. It finds the sliding speckle particle of N1739 test3 and the glare patch of N1690 test3, and flags nothing in tests without one |
| `animate_calibration.py` | Video of the scale calibration for one tine pair (`--def Def01 --pair X3`): the tracked tips on the specimen, zoomed insets of the two tips against their rest positions, and the calibration plots (tip separation vs reported rake motion, residual vs force) drawn up to the current frame → `scale_calibration_<pair>.mp4` in the deformation's `DIC_results/` |
| `edge_effects.py` | How far tine effects reach in (`--tests N1739/test1 N1690/test1`): along every grid line parallel to an edge, the part of the normal stretch that follows the measured tine positions (tine-pitch amplitude) and mean \|E_zθ\|, against distance from that edge's tine-tip line, by protocol family, with core reference levels and the strain-window floor. Reports where each test's modulation settles next to the current `edge_mm` (green on the plots; `current_edge_mm`, `within_current` in the .csv; a printed summary naming any test beyond it). It never changes `edge_mm` → `comparisons/edge_effects_<tests>.png / .csv / _profiles.csv` |
| `make_config.py` | Drafts a test's `test_config.json` from its Def11 rest frame (`--label "labeled healthy"`): rods as straight dark lines in a background band (fallback bands when a large specimen covers the first, and wider smoothing for wide rods at high zoom), tips where each rod's dark core ends, refined with `calibrate_scale.find_tip` (tip coordinate only; its centerline sits on a rod edge); `tine_tips_px` = median refined tip per rake. Also sets `dic.mask_saturated` (on when the tissue in the rest frame has saturated glints, off when none or ≥ 10% saturated), unless the config already sets it. A tip off its rake's median is searched again (longer highlight gap, looser darkness limit, then the calibration's tip finder walked in from the median). Self-check: five rods per rake, tips within 35 px of their rake's median (tines bend), a clear step at every tip; on failure nothing is written (`--allow-warnings` overrides). `--update` replaces only the tine settings and ROIs of an existing config and starts a failing rake from its existing tines; `--review` shows the tips and asks whether they are right, and the rakes answered wrong (or failing the self-check without existing tines) are clicked (kept as clicked, not refined); `--out` writes elsewhere → `test_config.json`, `<test>/DIC_results/config_tine_tips.png` |
| `measure_tines.py` | Each rake's tine span in the reference rest frame (rod centerlines from both edges, 15 px behind the tips); w_θ and w_z = means of opposing rakes. `--write` updates `widths_mm` in `test_config.json` |
| `validate_bands.py` | Checks that the strain bands are not code artifacts, on the test's Def01 (first frame to the rake peak). Four tests: a uniform-stretch null test, recovery of imposed E_zz and E_θθ bands, a brute-force check of the strain calculation, and an independent template-matching tracker. Needs the Def01 track first (`track.py --def Def01`). Output: `bands_validation.png` in Def01's `DIC_results/` |
| `compare_bands.py` | Cross-test band comparison at the peak of one deformation (`--sample N1739 --tests test1 test2 test3 --def Def11`): E_zz and E_θθ maps, maps divided by their mean, and band profiles with controls, from each test's exported fields → `comparisons/bands_<Def>_<tests>.png / .csv` |
| `peak_j_figure.py` | Small multiples of local areal stretch J / interior mean at each protocol's peak (Def11–19, from the full-domain tracks; no exports needed). The lowest 3 × 2 mm region at Def11 is outlined, and its J/mean is printed per protocol. The color scale is fixed (0.8–1.2) so specimens compare directly. The strain window, interior edge and low-spot box come from `analysis_settings.json`; points tracked in fewer than half the frames are gray and left out of the mean. The interior used is drawn dashed, stated in the subtitle and recorded in the .csv (`interior_edge_z_mm`, `interior_edge_theta_mm`). `--window-mm L` uses another Gaussian width, `--box-px 75` the earlier 75 px box, `--edge-mm Z TH` another interior → `<test>/DIC_results/peak_J_over_mean_Def11-19[_w<L>mm \| _box<N>px][_edge<Z>x<TH>mm].png / .csv` (no suffix: the settings window and interior) |
| `compare_peak_j.py` | Cross-specimen comparison (`--samples N1739 N1690`: one to four samples, a map row each; window from `analysis_settings.json` or `--window-mm`): J/mean maps at Def11's peak for every test of the samples on one mm scale, with each test's lowest 3 × 2 mm region, plus per-test summaries (spread, interior share more than 10% below the mean, low-spot value; Def11 and the range over Def11–19). Uses the same computation as `peak_j_figure.py`, from the tracks; `--edge-mm Z TH` for another interior (recorded in the .csv) → `comparisons/peak_J_<samples>_w<L>mm[_edge<Z>x<TH>mm].png / .csv` |
| `sample_info.py` | The label and notes of a test, in a window filled in with what is saved → `sample` ("…, labeled healthy") and `notes` in `test_config.json` (the `info` step) |
| `sort_images.py` | Sorts a test's images into `DefNN_*` folders from its PhaseTiming log (`--dry-run` first; safe to re-run), then checks every folder for missing images, Data.csv rows and Coords1.csv rows (gaps in their times) and images without a row at their time → `Images/sort_summary.csv`; `--review` shows the table and asks whether to proceed (no: exits with an error) |
| `review_branches.py` | Loading branch of each deformation (`--defs`, default Def11–Def19), from the Stretch start to the rake peak (`dic2d.rake_peak`), shaded on the rake forces → `DIC_results/loading_branches.png`; `--review` asks whether they are right, and the wrong ones are clicked (start and end) → `loading_branches.json` (same file as the MATLAB version), read by `dic2d.branch_by_hand`: `--loading-only` tracking and the fits use it |
| `view_tracking.py` | Tracking report: the quality of each tracked deformation's loading branch (valid points and median ZNCC per frame, worst frame; Good ≥ 93% and 0.95, Problem < 80% or 0.90, Check between, as in the MATLAB version) → `DIC_results/tracking_report.csv` and a verdict line; then a window with the deformations side by side (tracked points and grid over the images, stretches, valid points and ZNCC of the frame shown), one slider and Play for all of them, opened at the worst frame of the worst one; "Open large" shows one in a big viewer. `--no-window` writes the report only; `--ask` (runDIC.sh's report step, before the figures and fits) asks after the window, when a deformation is marked Check or Problem, whether to proceed or stop the test (exits with an error, so its later steps are skipped) |
| `review.py` | Check windows for `sort_images.py`, `make_config.py` and `review_branches.py`: on with `--review` or `DIC_REVIEW=1` (`runDIC.sh`, `review = yes`) where there is a display; native backend on macOS, Tk on Linux |
| `status.py` | Status of every test from the files on disk: sorted, config, scale, widths, thickness, tracked deformations, J/mean figure, exports (`--notes` adds each test's notes) |

## Usage

```bash
.venv/bin/python strain_tracking/track.py --test N1739/test1 --def Def01 --tag test
```

```bash
.venv/bin/python strain_tracking/compare_cellscale.py N1739_LFPA/test1_08252026_PLC25_38-20/DIC_results/Def11_Equibiaxial.npz
```

```bash
.venv/bin/python strain_tracking/export_fields.py N1739_LFPA/test1_08252026_PLC25_38-20/DIC_results/Def11_Equibiaxial.npz
```

```bash
.venv/bin/python strain_tracking/animate_fields.py N1739_LFPA/test1_08252026_PLC25_38-20/DIC_results/Def11_Equibiaxial.npz
```

- **`--def`** accepts `Def13`, `Def1` (= Def01) or any unique folder prefix.
- **Subset size:** to change it, use `--radius 30 --tag r30`.
- **Time window:** to track an arbitrary window across all of `Images/`, use `--t0 … --t1 …`.
  The output goes to the `DIC_results/` of the deformation containing the first frame, or
  to wherever `--out` points.

## DIC settings

**Defaults, in physical units.** `analysis_settings.json` at the project root sets every size
in mm, the same for all samples and tests:

| Setting | Value | Used by |
|---|---|---|
| `target_px_per_mm` | 50 | images are resampled to this scale before tracking (area averaging, downsampling only) |
| `subset_mm` | 1.0 | subset side; radius = round((subset_mm × px/mm − 1) / 2) at the tracking scale (51 px) |
| `step_mm` | 0.2 | grid step, rounded to whole px (10 px) |
| `window_mm` | 1.5 | Gaussian strain window (equivalent box width) for F, J, E in every export and figure |
| `edge_mm` | z 1.5, θ 2.0 | interior for J/mean, spread and low spot: at least 1.5 mm from the left/right (axial) tine-tip lines and 2.0 mm from the top/bottom (circumferential) ones. User-set: edit it here, or pass `--edge-mm Z TH` to `peak_j_figure.py` / `compare_peak_j.py` for one run. `edge_effects.py` reports where the tine-pitch modulation of the normal stretch decays to the core level against it; on the tests so far that was mostly within 1.6 mm of the left/right lines and 2.0 mm of the top/bottom lines (a few tests up to 2.2–2.4 mm, by one or two 0.2 mm bins). Both values clear the strain-window floor (about 1.2 mm). Where a rake's outermost tine stops well short of a corner, that corner is under-stretched and can reach into the interior |
| `low_spot_box_mm` | [3, 2] | low-spot box (z, θ) |

**Resampling.** Each test's images are resampled to `target_px_per_mm` with its own
`scale_calibration.json` (factor = 50 / native px/mm, e.g. 0.80 at 62.8 px/mm), so every test
is matched with the same px as well as the same mm: 51 px subsets on a 10 px grid. Px-based parts
of the tracker (median-test tolerance, convergence, pre-blur, interpolation) and the displacement
noise in mm are then the same for every specimen, and tracking time no longer grows with zoom.
Saturated glints are detected on the native image; a resampled pixel is usable only if all of it
was clear. Tine tips, grid and excluded rectangles are converted to the resampled image, and
tracks are saved back in native image px, so configs, photos and CellScale points line up. Every
track records `px_per_mm` (native), `resample`, `track_px_per_mm`, `subset_mm` and `step_mm`;
`radius` is in the resampled px. On N1739 test1 Def11, tracking at 50 px/mm instead of native
62.8 px/mm changed peak displacements by a median of 0.20 px (3 µm) and J/mean by an RMS of
0.0017 (r = 0.999), with 99.9% instead of 99.6% of points valid per frame.

The effective strain resolution is about √(window_mm² + subset_mm²) ≈ 1.8 mm. `--radius` /
`--step` (resampled px) and `--win` (a square window in grid points) still override the settings.
`track.py` (single deformations, diagnostics) tracks at native resolution.

Tracks made before these settings used 61 px subsets on a 15 px grid with a 75 px box window.
That came from a resolution study on N1739 test1's Def01: it closed nearly all tracking holes
(99.9% valid per frame) and kept about 95% of the band contrast of a finer 41 px / 10 px
setting. Re-tracking N1739 test1 at the mm settings changed displacements by a median of
0.13–0.22 px and J/mean maps by an RMS of 0.0015–0.0025 (r ≥ 0.998)
(`comparisons/retrack_N1739_test1_mm_vs_px.png / .csv`). The Stage 2 local fits use a 45 px window (`--win 1`),
which the synthetic tests favored.

**Point validity, per frame.** A point is valid when it passes three checks:
- IC-GN converges (up to 100 iterations).
- ZNCC ≥ 0.75.
- The frame-to-frame displacement passes a normalized median test against the 8 neighbors
  (Westerweel & Scarano 2005; threshold 3, ε = 0.1 px).

A point that fails is handled as follows:
- It is flagged invalid, and its position for that frame is interpolated from its neighbors.
- It keeps its last good reference subset, and rejoins as a real measurement once it matches again.
- At a keyframe update, only the currently valid points are re-referenced.

**Glint masking.** A saturated glint (a bubble or a wet highlight) does not move
with the tissue and changes with load, so a subset that contains one can fail to match for a
whole protocol. With `mask_saturated` on (`--mask-saturated`, or `dic.mask_saturated: true` in
`test_config.json`), pixels of the reference (keyframe) image at raw gray ≥ 240, grown by 2 px,
get zero weight in each subset's ZNSSD. The mean, norm and Gauss–Newton update then use only the
usable pixels, and a subset needs at least half its pixels usable. Saturated pixels in the
current image are not masked by default (`TrackParams.sat_current`). Masking them as well made
each subset's pixel set change from frame to frame. Subsets next to the tines, whose highlights
flicker, then failed far more often: on N1739 test2 the edge ring fell from 100% to 85% valid.
Every output `.npz` records `mask_saturated`, `sat_level` and `sat_dilate`. With the option off,
the tracker is unchanged (bit-identical on test1).

In N1739 test2, the glint near (1016, 1251) px left about 24 points invalid in every protocol. On the
first 32 s of Def11 with the full-domain grid, the reference mask brings that area from 9% to
100% valid and the worst frame from 97.4% to 99.9%. Elsewhere, positions change by a median of
0.06 px, and the displacement field is slightly smoother.

**Excluded regions (optional).** `dic.exclude_px` in `test_config.json` lists rectangles in the
reference frame's pixels whose pixels get zero weight in every subset. They move with the
tissue: each protocol's rectangles are shifted by its link displacement, and they're shifted
again at every re-reference. Only subsets centered on a rectangle lose their match; their
neighbors track on the remaining pixels, and the strain field fills the gap. In N1739 test3, a
speckle particle near (1064, 849) px slides up to 0.47 mm over the tissue at peak load. It
produced a red/blue artifact in J (0.83–1.17); masking its strip removes the artifact (0.96–1.04
around it) with full strain coverage, and points away from it are unchanged. A saturated glare
patch whose edges shift with load is the other use: it pulls the subsets that overlap it along
with its edges (a red/blue pair in J), which a rectangle over the patch stops.
`artifact_check.py` finds such regions and writes the rectangles; `peak_j_figure.py --exclude` instead
blanks a region in the figure only.

**Protocol labels.** Figure titles name each deformation by its loading protocol
(`dic2d.protocol_label`, e.g. "Axial reduced, 0.75 : 1 (Def12)") and its programmed rake stretch,
read from the protocol table of the test's `.tst` file (`dic2d.programmed_stretch`; X = axial,
Y = circumferential, in %). Times in figures are measured from the start of the protocol.

**Image sampling.** Subsets are sampled bicubically (`cv2.remap`). Samples whose 4 × 4 stencil
leaves the image are set to NaN explicitly: a NaN border value inside `cv2.remap` leaks NaN into
some interior samples at certain subset sizes (63 and 65 px, for example), which would void
every subset.

**Strain and F.** A point needs ≥ 50% of its window valid. A point that is itself invalid
still gets a value from its neighbors if they surround it (support centroid within one grid
step), so values are interpolated, never extrapolated. Pass `min_frac=0.6, fill_holes=False`
to `deformation_gradient_field` for the stricter rule.

**Comparing specimens.** Spread measures (SD, fraction below a threshold, profile contrast)
shrink as the smoothing grows, and normalizing by the mean doesn't remove that. Both the subset
and the strain window smooth, so both are set in mm (`analysis_settings.json`): a window fixed in
px (the earlier 75 px box was 0.86 mm on N1690 test1 and 1.19–1.41 mm on N1739) or a subset
fixed in px leaves the smoothing zoom-dependent. The Gaussian window
(`deformation_gradient_field_mm`) uses σ = L/√12. It is exact on affine fields, gives the same
mean J as the box, and gives about 5% lower SD than a box of the same nominal width.

**Synthetic check.** A real frame was warped affinely up to E_zz = 0.45, saved as JPEG q92,
with σ = 2 gray-level noise added. The results:
- Mean strain recovered to 1e-4.
- Point-to-point standard deviation ≤ 0.0015.
- Position RMS < 0.09 px.
