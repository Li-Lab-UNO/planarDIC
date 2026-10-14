# Handoff: local strain-based constitutive parameter estimation

Sample N1739 LFPA, test1 (08/25/2026), CellScale planar biaxial, images at 5 Hz.
Attach this file at the start of a new thread. It gives the state of the project, the frozen
bulk baselines that local fits are compared against, and the plan for the local estimation.

Status: **Stage 0 (local data set and diagnostics) and Stage 1 (uniform-stress patch map, with and without a smoothness penalty) are complete, and the domain has been segmented by deformation-gradient jump. No equilibrium-based (Stage 2) or finite element fit has been run.**

---

## 1. Standing instructions (from the user)

- **Fitting data:** Def11–Def19 only, which are all of the stretch-ratio protocols (PhaseTiming Sets 2–10).
  - Never use the preconditioning protocols Def01–Def10 for estimation.
  - Don't suggest Def13 or Def19 as single-protocol test cases.
- **No fitting until asked:** don't run constitutive fits until the user asks. Ask before long batch runs.
- **Writing style:**
  - Write subscripts with underscores (C1_el, λ_z, E_θθ, t_zz), never HTML `<sub>`.
  - Use American spelling in replies, plot labels, comments and docs.
- **Axes:**
  - Image x = tester X = axial **z**, increasing to the right.
  - Circumferential **θ** runs along image y and increases **upward** (θ = −image y).
  - Tracking works in image coordinates. `dic2d.to_specimen(F)` converts F to (z, θ), and `dic2d.THETA_SIGN` holds the sign.
  - Only the signs of F_zθ, F_θz, E_zθ, u_θ and φ_1 depend on this; stretches and J do not.
  - Every export records its convention in `axis_convention`.
  - Export names use ASCII `theta` (`E_thetatheta`, `lambda_theta`, `u_theta`).
  - Spatial plots put θ on the vertical axis, measured from the bottom tine tips; z is measured from the left tine tips.
- **Zones:** a grid of patches over the region between the tine tips, set by the user (default 10 × 10). Horizontal strips are not needed and are not defined.
- **Figures:** for maps compared across protocols, use common color limits per column and keep a colorbar beside each panel.
- **Solution family:** the elastin-off solution (C1_el ≈ 0, γ ≈ 34°) is excluded as non-physiological. Keep fits in the elastin-active family (C1_el > 0.1 kPa).
- **Tine tips:** dropped as a stretch input.
- **Folder layout:**
  - Images live in `Images/DefNN_*` folders.
  - Code lives at the project root `DIC/`: `DIC/strain_tracking/` (DIC) and `DIC/constitutive_fit/` (fits, including this file), with one shared `DIC/.venv`. Sample folders (`N1739_LFPA/`, `N1690_RFPA/`) sit next to them.
  - Every script needs `--test SAMPLE/testN`, e.g. `--test N1739/test1` (or `DIC_TEST`); test-specific settings are in `<test>/test_config.json`.
  - Tracks and per-protocol results go in `<test>/DIC_results/`, named by deformation.
  - Constitutive results go in `<test>/constitutive_fit/`.

## 2. Paths

Test root: `N1739_LFPA/test1_08252026_PLC25_38-20/` (relative to the project folder)

| What | Where |
|---|---|
| Python | `DIC/.venv/bin/python` (use `python -m pip`); scripts run with `--test N1739/test1` |
| Test settings | `test_config.json` (ROIs, tine tips, widths, H, reference protocol) and `scale_calibration.json` |
| Raw forces and rake size | `ID_N1739_LFPA_PLC25 38-20_08.25.2026_test1Data.csv` (`Time_S`, `X_Force_mN`, `Y_Force_mN`, `X_Size_um`, `Y_Size_um`; latin1 encoding) |
| CellScale 4-point tracking | `Coords1.csv` (skip 1 row; X1..Y4; TL, TR, BR, BL = points 1, 3, 4, 2) |
| Phase timing | `Logs/*_PhaseTiming.csv` (Stretch and Recover start times; the loading branch ends at the rake peak, see §5) |
| DIC tracks, full domain | `DIC_results/Def1N_*.npz` |
| DIC tracks, original ROI | `DIC_results/Def1N_*_standard.npz` (sliced from the full-domain track when not tracked separately) |
| Bulk baseline, full domain (primary) | `constitutive_fit/baseline/` |
| Bulk baseline, original ROI (reference) | `constitutive_fit/standard_roi/baseline/` |
| Local data set and Stage 0 figures | `constitutive_fit/local/` |
| Stage 1 results | `constitutive_fit/local/stage1*/` (with and without the outer ring; with a smoothness penalty) |
| Segmentation | `constitutive_fit/local/segmentation/` |
| Thickness measurements | `N1739_LFPA/Thickness, 1739/Thickess measurements.xlsx` (sic); the values used are in `test_config.json` |

## 3. DIC tracks

Two sets of tracks exist for Def11–Def19. Both use the same settings and the same common reference: the first frame of Def11, at t = 760.0 s. Each protocol's rest frame is matched to that reference by IC-GN ("bridging").

| | Full domain (`_tine_region`) | Original ROI |
|---|---|---|
| Script | `track_protocols.py` (default `--grid tine_region`) | `track_protocols.py --grid standard`, or sliced from the full domain by `dic2d.load_track` |
| ROI setting (x0, y0, x1, y1, px) | [760, 755, 1285, 1365] | [790, 785, 1255, 1335] |
| Grid (cols × rows) | 36 × 41 = 1476 points | 32 × 37 = 1184 points |
| Extent of grid points (px) | x 760–1285, y 755–1355 | x 790–1255, y 785–1325 |
| Distance from the tine tips | 5 px left and right, 10 px top, 15 px bottom | 35–45 px |

- The original grid is rows 2:-2, cols 2:-2 of the full-domain grid. Its tracks agree with the full-domain ones to a median of 0.03 px.
- **Rest-state tine tips** (from `calibrate_scale.py`, px): axial rakes at x ≈ 755 (left) and 1290 (right); circumferential rakes at y ≈ 745 (top) and 1370 (bottom). The region between the tips is 8.52 mm (z) × 9.95 mm (θ).
- **Arrays:**
  - `X0` is the reference grid; `shape` is [rows, cols].
  - `pos` has shape (n_frames, n_points, 2), in image px.
  - `valid`, `zncc` and `fail` are per frame.
  - `X_start` and `bridge_zncc` describe the bridge.
  - `times`, `paths` and `keyframes` describe the frames.
- **Coverage (full domain):** valid points are ≥ 98.8% in every frame.
- **DIC settings:**
  - Subset 61 px (radius 30).
  - ZNCC ≥ 0.75.
  - Per-point reference updating.
  - Normalized median test for outliers.
  - Grid step 15 px (0.239 mm).
  - Strain window 75 px = 5 × 5 points = 1.19 mm (`dic2d.window_points(step)`).
- **Scale:** 62.80 px/mm ± 2% (`<test root>/scale_calibration.json`), or 15.92 µm/px.
- **Per-protocol outputs from the full-domain tracks** (comparison with CellScale, fields `.npz`/`.mat`, animations) exist for Def11 only, tagged `_tine_region`.

## 4. Code you can reuse (`DIC/strain_tracking/` and `DIC/constitutive_fit/`)

| Module | Relevant API |
|---|---|
| `dic2d.py` | `deformation_gradient_field(X0, x, valid, shape, win, min_frac=0.5, fill_holes=True)` → (F[rows, cols, 2, 2] in image frame, ok mask); `to_specimen(F)`; `THETA_SIGN`; `affine_fit`, `bilinear_F_center`, `green_lagrange`, `window_points`, `load_scale`; name constants `STRAIN_KEYS`, `STRETCH_KEYS`, `DISP_KEYS`, `F_KEYS`, `FORCE_KEYS`, `AXIS_CONVENTION` |
| `constitutive.py` | `PARAM_NAMES = (Cgr, C1_el, C2_el, C1_smc, C2_smc, C1_col, C2_col, gamma)` (gamma in degrees); `cauchy(F, P)` for a **general in-plane F** (shear allowed, incompressible, plane stress); `cauchy_biaxial(lz, lt, P)` closed form; `first_pk_membrane(F, P, H)` → N = H·t·F^-T (kPa·mm = N/m); Div N = 0 on the reference grid |
| `compare_stretch_inputs.py` | `build_table(sets, tracks="standard" \| "tine_region")` (protocol → folder mapping, loading-branch frames, tare); `fit(lz, lt, tzz, ttt, starts, rng, seeds=())`, the multi-start elastin-active fit; bounds `LO`/`HI`; `EL_MIN = 0.1`; `WIDTHS` from the test's `widths_mm` (test1: tine span 9.57 / 9.02 mm, specimen 15.4 / 13.4 mm). CLI: `--tracks`, `--table-only` (writes the stretch table, no fitting) |
| `baseline.py` | `load(tracks="standard" \| "tine_region")` (checks the checksum) → (baseline dict, data DataFrame); `evaluate(params, input="dic_affine", tracks=...)` → R2, RMS and per-protocol RMS on the frozen data, with an "outside 95%" flag per parameter. The default `tracks` is `"standard"` (original ROI), so local work must pass `tracks="tine_region"` |
| `track_tine_region.py` | full-domain tracking; `ROI`, `TINE_TIPS` (from the test's `test_config.json`), `SUFFIX` |
| `testcfg.py` | `get()` → the selected test (`root`, `images`, `fit_dir`, `local_dir`, `scale_file`, `setting(...)`, `folder("Def11")`); `select("N1739/test1")` in Python |
| `build_local_data.py` | Stage 0: builds `local_data.npz` and the Stage 0 figures. `--patches N_Z N_THETA` (default 10 10), `--baseline tine_region \| standard` (default `tine_region`) |
| `stage1_patches.py` | Stage 1: one stiffness factor per patch from the patch-mean F and the global stress; `load_inputs()` returns the patch-mean F, the baseline model stress at it and the global stresses. `--exclude-ring N` drops the outer N rings of patches from scales, statistics and diagnostics (the fits are independent per patch, so interior α do not change) |
| `stage1_smooth.py` | Stage 1 with a Tikhonov penalty on neighbor differences of log α, λ by leave-one-protocol-out cross-validation; `--exclude-ring N` leaves the outer rings out of the fit and the penalty |
| `segment_fjump.py` | deformation-gradient-jump segmentation (Witzenburg et al. 2016) on bilinear elements of 2 × 2 grid cells; F components scaled by default (`--unscaled` for the paper's Eq. 3); `--uniform-jump` null control |
| `export_fields.py`, `compare_cellscale.py`, `animate_fields.py` | per-protocol field export, CellScale comparison and MP4s for any track `.npz`; specimen frame with θ up |
| `check_preprocessing.py` | template for "refit variant X and compare it with the baseline and its 95% interval"; `--tracks tine_region` (default) or `standard` |

## 5. Bulk baselines (the comparison target)

The bulk baseline has been fitted two ways, which differ only in the DIC tracks behind the global stretches:

- **Full domain** (`baseline/`): tracks out to the tine tips. This is the **primary comparison target**, because its domain matches the local data set.
- **Original ROI** (`standard_roi/baseline/`): tracks over the inner ~87% of the gap between the tine tips. It is kept as a reference.

Everything else is identical. Both use the same 1492 loading frames, the same tared forces and the same fit settings (below). A parameter set must be scored on the baseline built from the same domain: the original-ROI parameters scored on the full-domain data give RMS 0.569 kPa rather than 0.431 kPa.

Each baseline folder holds:
- `baseline.json`: settings, parameters for every stretch input, bootstrap 95% intervals, per-protocol RMS, SHA-256 checksums.
- `baseline_data.csv`: the 1492 loading frames.
- `baseline_params.csv`.
- The figures.

**Superseded inputs.** Both frozen baselines, and every local result built on them (Stages 0–2, segmentation), used H = 2.533 ± 0.228 mm and tine-span widths w_θ = 9.83 mm, w_z = 9.06 mm. test1's `test_config.json` now has the corrected thickness H = 2.756 ± 0.249 mm and opposing-rake mean widths w_θ = 9.57 mm, w_z = 9.02 mm (`measure_tines.py`). All stresses scale with 1/(w·H): t_zz falls by 5.6% and t_θθ by 7.7% relative to the frozen values, so Cgr and the C1 constants move by about as much, while exponents and γ change little. The numbers below describe the frozen baselines until they are re-fitted.

`fit_4pt.py` fits the same model to the CellScale 4-point data on the same 1492 frames (outputs `constitutive_fit/fit_4pt_*`). Its best fit matches the baselines' CellScale 4-point row exactly. Its intervals come from a frame-level bootstrap (2000 resamples), so they are much narrower than the baselines' protocol-level intervals and should not be compared with them.

### Model and fit settings

- **Model:** Jadidi et al. 2019 bulk four-fiber model, Eq. 2. It has a neo-Hookean term (Cgr) plus exponential fiber families:
  - elastin along e_z;
  - SMC along e_θ;
  - two collagen families at ±γ from z, sharing C1/C2.
  - Fibers are tension-only (Macaulay bracket on I4 − 1). The material is incompressible, in plane stress.
- **Thickness (frozen baselines):** H = 2.533 ± 0.228 mm (corrected: 2.756 ± 0.249 mm).
- **Stress:** t_zz = F_z·λ_z/(w_θ·H) and t_θθ = F_θ·λ_θ/(w_z·H), with tine-span widths w_θ = 9.83 mm and w_z = 9.06 mm in the frozen baselines (now the opposing-rake means, 9.57 and 9.02 mm).
- **Loading branch:** from the logged Stretch start to the frame where the rakes peak (`compare_stretch_inputs.loading_end`, searched up to 2 s past the logged Recover start). The rake peak equals the logged Recover start in every protocol except Def12, where the rakes and forces keep rising for 0.5 s past the logged 864.0 s and peak at 864.5 s. Def12's timestamps also shift by 0.1 s at that moment (frames move from even to odd tenths of a second), which points to a timing glitch in the tester log.
- **Forces:** tared to the first frame of each protocol. At rest the load cells read −11 to −18 mN.
- **Objective:** unweighted sum of squared Cauchy residuals over both directions.
- **Optimizer:** multi-start trust-region-reflective (150 starts per input, then re-seeding from every input's best solution).
- **Uncertainty:** protocol-level bootstrap (200 resamples). In both fits, no dic_affine resample fell into the elastin-off solution.

### Primary input: DIC affine F (`dic_affine`)

| Param | Full domain | Bootstrap 95% | Original ROI | Bootstrap 95% |
|---|---|---|---|---|
| Cgr (kPa) | 6.075 | 4.24 – 6.17 | 6.252 | 4.43 – 6.45 |
| C1_el (kPa) | 4.106 | 3.76 – 8.07 | 3.583 | 3.10 – 7.46 |
| C2_el | 1.795 | 0.001 – 1.93 | 1.891 | 0.001 – 2.10 |
| C1_smc (kPa) | 10.466 | 9.51 – 14.93 | 10.047 | 8.99 – 14.82 |
| C2_smc | 10.457 | 9.01 – 10.99 | 10.747 | 9.13 – 11.49 |
| C1_col (kPa) | 1.103 | 0.65 – 1.92 | 1.318 | 0.74 – 2.10 |
| C2_col | 7.950 | 5.37 – 8.73 | 7.723 | 5.50 – 8.58 |
| gamma (deg) | 46.31 | 31.3 – 48.6 | 47.66 | 32.0 – 49.4 |
| R2 / RMS (kPa) | 0.9978 / 0.431 | | 0.9977 / 0.437 | |

Every full-domain parameter lies inside the original-ROI 95% interval. The largest relative shifts are in the two least constrained parameters: C1_el (+15%) and C1_col (−16%).

### Other inputs, full-domain baseline (same settings)

| Input | Cgr | C1_el | C2_el | C1_smc | C2_smc | C1_col | C2_col | γ | RMS (kPa) |
|---|---|---|---|---|---|---|---|---|---|
| CellScale 4-point | 5.91 | 3.55 | 1.81 | 9.27 | 10.03 | 1.83 | 7.15 | 49.4 | 0.492 |
| DIC 4-point | 6.06 | 3.54 | 1.71 | 9.99 | 9.97 | 1.30 | 7.22 | 46.9 | 0.438 |
| DIC region-mean λ | 5.43 | 4.72 | 1.55 | 11.62 | 10.89 | 1.13 | 7.97 | 43.8 | 0.424 |
| DIC region edges | 5.70 | 4.30 | 1.82 | 10.97 | 11.15 | 1.27 | 8.18 | 45.8 | 0.428 |
| rake | 4.17 | 4.19 | 1.18 | 9.83 | 6.70 | 1.35 | 6.12 | 46.1 | 0.425 |

CellScale 4-point and rake don't depend on the DIC domain, so they are identical in the original-ROI baseline. DIC 4-point differs only in the fourth digit.

### Findings from the bulk fits

- **Stretch input:** the tissue-based inputs (CellScale 4-point and the DIC measures) agree within the elastin-active family. Rake stretch gives softer C2 exponents.
- **Domain:** including the material near the tines lowers the global axial stretch by about 0.005 at the peak and raises the circumferential one by about 0.002. dic_affine peaks go from 1.357 / 1.173 (original ROI) to 1.352 / 1.175 (full domain).
- **Preprocessing:**
  - From `check_preprocessing.py` on both baselines; results are in `preprocessing_sensitivity.csv` and `.png` in each baseline folder. The two baselines give the same conclusions.
  - Taring matters. Without it, dic_affine gives Cgr −44% (full domain) or −43% (original ROI), C1_col +69% or +62%, C1_el +33% or +38%, and C2_el −22%, and RMS rises from 0.431 to 0.558 kPa (full domain) or from 0.437 to 0.564 kPa (original ROI). Cgr falls outside the 95% interval in both; C1_col falls outside only for the original ROI, whose interval is narrower. For CellScale 4-point, Cgr (−47%) and C1_col (+52%) both fall outside.
  - Taring to the mean of the first 5 frames instead of the first frame changes the parameters by at most 13% (C1_col), all inside the 95% interval.
  - Savitzky–Golay smoothing (11 or 21 frames) changes the dic_affine parameters by less than 0.5% and the CellScale 4-point parameters by at most 2%.
- **Identifiability:**
  - A single Def01 path is not identifiable.
  - Def11–19 together identify the parameters well, apart from the two near-equal minima (elastin on or off).
- **Student's results:** the student's S1 parameters fall inside these intervals once the family labels are aligned (her c1/c2 = collagen, c1^3/c2^3 = SMC, c1^4/c2^4 = elastin).

### Loading branch per protocol (full-domain DIC affine peaks, tared forces)

| Set | Folder | Frames | λ_z max | λ_θ max | F_z max (mN) | F_θ max (mN) |
|---|---|---|---|---|---|---|
| 2 | Def11_Equibiaxial | 191 | 1.351 | 1.174 | 954 | 887 |
| 3 | Def12_1-0.75_ratio | 144 | 1.262 | 1.172 | 376 | 658 |
| 4 | Def13_1-0.5_ratio | 101 | 1.171 | 1.171 | 217 | 561 |
| 5 | Def14_1-0.25_ratio | 101 | 1.082 | 1.168 | 108 | 501 |
| 6 | Def15_Equibiaxial | 191 | 1.352 | 1.175 | 954 | 883 |
| 7 | Def16_0.75-1_ratio | 191 | 1.346 | 1.131 | 756 | 436 |
| 8 | Def17_0.5-1_ratio | 191 | 1.342 | 1.086 | 643 | 251 |
| 9 | Def18_0.25-1_ratio | 191 | 1.339 | 1.041 | 576 | 152 |
| 10 | Def19_Equibiaxial_Save | 191 | 1.351 | 1.174 | 945 | 884 |

## 6. What is known about heterogeneity (Stage 0)

- **The strain field is real, not a code artifact.** It was validated with a null test, imposed synthetic bands, a brute-force strain check and independent template matching, all on Def01 (`validate_bands.py`).
- **Circumferential strain varies mainly with θ.** The z-averaged ratio E_θθ / global E_θθ has the same shape in every protocol. There are stiffer horizontal regions at θ ≈ 4–6 mm and along the bottom tine tips (ratio ≈ 0.4–0.8), and softer ones at θ ≈ 1.3–4.2 and 5.9–8.7 mm (≈ 1.1–1.3).
- **Axial strain forms oblique streaks.** λ_z is highest along streaks running from lower left to upper right, which repeat across protocols, and near the right edge and top-right corner. Its z-averaged profile is nearly flat in θ, but it varies by about ±3% along z.
- **The tines affect a margin of about 1–1.5 mm.** The shear angle changes sign from one tine to the next along every edge, and the corners show local peaks. Patch-level spread is largest in the edge and right-column patches.
- **Local stretch exceeds the machine (rake) stretch.** At every protocol's peak, 19–29% of the area exceeds the machine λ_z and 9–36% exceeds the machine λ_θ, in 5–9 separate regions. The local maximum is 0.08–0.16 above the machine stretch, far above DIC noise (CellScale agreement RMS 0.002–0.004 in strain). The exceedance is present in the interior as well as near the tines.
- **Local and global are consistent on the full domain.** At each peak, the area-weighted mean of the local stretch is within 0.007 of the full-domain dic_affine stretch, and 45–53% of the area lies above it.
- **Thickness:** there is no thickness map, so local thickness cannot be measured.
- **Shear:** the local F includes in-plane shear (F_zθ, F_θz). Local fits must use `constitutive.cauchy(F, P)`, not the diagonal closed form.

## 7. Constraints that shape the method

1. **Only the boundary resultants are measured.** F_z and F_θ (two numbers per frame) are the only force data. Local stress is never measured, so local parameters must come from an inverse method that couples the local F field with equilibrium and those two resultants.
2. **Thickness and stiffness are confounded.** Equilibrium involves N = H·t. Without a thickness map, only **membrane stiffness** (H_k × moduli for patch k) is identifiable. Exponents (C2_*) and γ are not confounded by H. Report patch results either as H·C products or as moduli at the nominal H (test1: 2.756 mm), stated as such.
3. **The measured domain and the force widths differ.**
   - The tracked region is 8.52 mm (z) × 9.95 mm (θ), between the tine tips. The bulk fits spread the forces over the tine-span widths, the means of opposing rakes (test1: w_θ = 9.57 mm, w_z = 9.02 mm; the rakes differ by 4–8%), and the specimen is 15.4 × 13.4 mm.
   - To stay comparable with the baseline, use **the same width assumption**. A section cut through the domain carries the fraction (cut length / tine-span width) of the measured force. This is the baseline's uniform-flux assumption, applied cut by cut.
   - Equations built from interior virtual fields that vanish on the domain boundary need **no** width assumption. They constrain only the *relative* patch stiffness. Keep them as a separate block so their effect can be reported on its own.
4. **No simple stress assumption is exact.**
   - The θ variation is close to horizontal, so for circumferential load the material acts mostly in series, and uniform t_θθ is a reasonable approximation.
   - The axial streaks are oblique, so neither uniform stress nor uniform stretch holds in z.
   - Any shortcut that pairs local stretch with the global stress is biased, most strongly in z.
5. **Patch values are not independent at the default resolution.**
   - Default patches are 0.852 × 0.995 mm (10 × 10, aspect 0.86), smaller than the 1.19 mm strain window and the 0.97 mm subset. Neighboring patch means share displacement data.
   - With 100 patches, even one factor per patch means 100 unknowns on top of the shared parameters. In Stage 1 each patch is fitted on its own and cross-validation rejects any smoothing (§8), but in Stage 2 the patches are coupled through equilibrium and only two resultants are measured, so a penalty (or fewer unknowns) is expected to be necessary there.
   - Near-square alternatives: 6 × 7 gives 1.42 mm squares (larger than the strain window); 10 × 12 gives 0.85 × 0.83 mm; 12 × 14 gives 0.71 mm squares.
6. **Edges and coverage.**
   - Grid cells extend to the tine tips, so the patches tile the domain with no gaps. Cells that straddle a patch edge are split by area.
   - The outer 2 rows and columns of points have truncated strain windows (`full_window = False`), and some corner points have no F. Corner patches have 75–80% F coverage at 10 × 10; all others are ≥ 98%.
7. **Data volume:** 1492 loading frames × 1476 points. Subsample frames (for example every 2nd–4th frame), then check that the result doesn't change.

## 8. Plan (staged; each stage gets a script and a comparison against `baseline.evaluate(..., tracks="tine_region")`)

### Stage 0: local data set and diagnostics (done)

`build_local_data.py` writes `constitutive_fit/local/`:

- **`local_data.npz`:**
  - F in the specimen frame for each of the 1492 baseline frames, on the 41 × 36 grid. Layout [row, col] = [θ, z], with row 0 at the bottom tine tips and col 0 at the left tine tips.
  - Masks `ok`, `valid_track`, `full_window`.
  - Cell geometry: `z_mm`, `theta_mm`, `cell_edges_*_mm`, `area_mm2`.
  - Patches: `patch_weights_mm2[i, j, row, col]` = area of cell (row, col) inside patch (i, j), with i the θ row from the bottom and j the z column from the left; `patch_edges_*_mm`; `patch_coverage`; patch-mean strain concentrations.
  - Tared forces, t, set and folder, and the dic_affine λ, all from the full-domain baseline.
  - Strain-concentration maps `conc_theta`, `conc_z`.
  - Checksums of the source tracks, the baseline data and the scripts.
- **`patches.json`:** patch counts, sizes, edges, coverage and strain concentrations.
- **Figures:**
  - `stage0_peak_maps.png`: λ_z, λ_θ, shear angle and J_2 at the peak of Def11 and each ratio protocol.
  - `stage0_band_profiles.png`: means along z vs θ, and along θ vs z.
  - `stage0_patches.png`: patch-mean strain concentrations, coverage, within-patch spread and patch-mean stretch histories.
  - `stage0_patch_layout.png`: the patch grid on the reference image, and a zoom on the DIC grid cells.

### Stage 1: patch map from local stretch and the global stress (done)

`stage1_patches.py` fits one factor α_k per patch (ladder level (a)): the patch-mean F goes through the full-domain baseline (dic_affine parameters), α_k scales Cgr, C1_el, C1_smc and C1_col, and α_k is fitted against the **global** t_zz and t_θθ of all 1492 frames (uniform-stress assumption; closed form, since the model is linear in those four parameters). α_z and α_θ (fits to one stress component each) are diagnostics. Uncertainty: 200-resample protocol bootstrap. Outputs in `local/stage1/` (all patches) and `local/stage1_interior1/` (outer ring of patches left out of scales and statistics): `stage1_maps.png`, `stage1_fits.png`, `stage1_diagnostics.png`, `stage1_patches.csv`, `stage1_results.npz`, `stage1_summary.json`.

| | All 100 patches | Interior 64 (outer ring out) |
|---|---|---|
| α range (area-weighted mean) | 0.06 – 1.97 (0.94) | 0.41 – 1.45 (0.94) |
| α_θ / α_z range | 0.06–2.83 / 0.06–2.21 | 0.44–1.82 / 0.37–1.85 |
| Patch RMS, median (α = 1 for comparison) | 2.67 kPa (4.04) | 2.11 kPa (3.21) |
| α 95% interval excludes 1 | 82 / 100 | 52 / 64 |
| r(α, E_θθ concentration) / r(α, E_zz concentration) | −0.55 / −0.62 | −0.46 / −0.67 |
| r(α_z, α_θ) | 0.33 | 0.23 |

Findings:
- **The extremes are tine effects.** α of 0.06–0.2 (right column and top-left corner) and 1.5–2.8 (bottom row) all sit in the outer ring next to the tine tips. Under uniform stress, any strain concentration reads as soft material, and the exponential model turns a moderate excess strain into a very small α.
- **Interior pattern (α_θ, the more defensible map):** stiffer rows at θ ≈ 4–6 mm across the full width (1.16–1.56), softer rows at θ ≈ 7–9 mm (0.54–0.87), and a stiff-to-soft gradient along the row at θ ≈ 1–2 mm (1.82 at the left to 0.67 at the right). This matches the Stage 0 E_θθ bands.
- **Interior pattern (α_z):** stiffer upper left (1.4–1.85) and right side below θ ≈ 6 mm (1.2–1.5), softer upper center (0.42–0.7). It tracks the oblique λ_z streaks, not the θ bands.
- **The two directions disagree in a structured way.** α_z / α_θ ranges from 0.46 to 2.30, while most interior intervals are only 2–20% wide. One factor per patch cannot fit both directions under uniform stress.
- **The fits are poor.** Patch RMS is 0.7–7.6 kPa against 0.43 kPa for the bulk fit. At a single patch the protocols follow different stress–stretch curves, which is the signature of the uniform-stress assumption failing. Stage 1 gives a ranking of stiffer and softer regions, not stiffness values.

**Stage 1 with a smoothness penalty** (`stage1_smooth.py`; outputs in `local/stage1_smooth/` and `local/stage1_smooth_interior1/`): minimize the stress misfit plus λ × the mean squared neighbor difference of log α, with λ chosen by leave-one-protocol-out cross-validation.

| | All 100 patches | Interior 64 (outer ring left out of the fit) |
|---|---|---|
| Cross-validated λ | 0 | 0 |
| Held-out RMS at λ = 0 | 3.41 kPa | 2.56 kPa |
| Held-out RMS at half the roughness (λ ≈ 0.3–0.4) | 4.44 kPa | 2.79 kPa |
| Held-out RMS at 10% of the roughness (λ = 10) | 8.48 kPa | 3.63 kPa |
| α range at λ = 0 → 50% → 10% roughness | 0.06–1.97 → 0.12–1.42 → 0.27–0.53 | 0.41–1.45 → 0.47–1.23 → 0.66–0.92 |

- Every amount of smoothing predicts held-out protocols worse, so the patch differences are consistent across protocols rather than noise. Within Stage 1, a penalty would be a prior, not something the data support.
- The smoothest surviving component is a gentle gradient, stiffer toward the lower left (about 0.9) and softer toward the upper right (about 0.66).
- With all patches in the fit, heavy smoothing pulls every α down to 0.27–0.53, because the tine-adjacent patches have very large model stress and dominate the stress-weighted fit. Penalized fits should leave them out or treat them separately.

### Segmentation by deformation-gradient jump (done)

`segment_fjump.py` follows Witzenburg, Dhume, Lake & Barocas (2016): bilinear elements on the tracked points (2 × 2 grid cells, 0.48 mm; 18 × 20 = 360 elements), F at each element-edge midpoint, the jump across each shared edge averaged over the last 5 loading frames of each protocol, each protocol normalized by its 99th-percentile jump and the nine summed. Links are removed in order of edge betweenness × jump (Girvan–Newman), and each partition and its modularity are recorded. Pieces under 2% of the elements (always a few in the tine-adjacent corners) are recorded but not counted as subdomains. Outputs in `local/segmentation/`: `fjump_segmentation.png`, `segmentation.npz` (summed and per-protocol jumps, labels for every partition), `segmentation.json`.

- **Component scaling (default).** With the paper's ||F_1 − F_2||^2, the larger axial stretches make F_zz dominate and the cuts follow only the z structure. By default each F component's squared jump is divided by its mean square over the interior links, per protocol, so every component counts equally (`--unscaled` gives the paper's form). Typical jumps before scaling: F_zz 0.016–0.033, F_θθ 0.013–0.017.
- **Checks.** Interior jumps at the peak are 230–650 times their rest-frame value, so noise does not drive the result. With every jump set equal (`--uniform-jump`), the first cut is a horizontal bisection; the real data cut differently, so the jumps, not the grid geometry, place the boundaries. The largest jumps sit along the top and right boundaries at the tine positions.
- **Partitions (scaled).** Modularity is 0.34, 0.57, 0.63 and 0.69 for 2–5 subdomains and levels off at 0.69–0.71 for 5–10. The 5-subdomain partition (sizes in elements; mean E_θθ / E_zz concentration):
  - lower right, z ≈ 6.7, θ ≈ 3.2 mm (88; 1.00 / 0.91);
  - lower left, z ≈ 2.5, θ ≈ 2.2 mm (85; 0.89 / 1.05);
  - upper right, z ≈ 6.0, θ ≈ 7.9 mm (81; 1.11 / 1.16);
  - mid-left band, z ≈ 2.7, θ ≈ 6.0 mm (58; 0.83 / 0.97), part of the stiff θ band;
  - upper-left corner, z ≈ 1.1, θ ≈ 7.9 mm (37; 1.19 / 0.77).
- **The tissue looks graded, not piecewise uniform.** Share of the element-level strain-concentration variance explained:

  | Partition | Regions | E_θθ | E_zz |
  |---|---|---|---|
  | Jump segmentation, scaled | 5 / 6 | 21% / 24% | 40% / 45% |
  | Jump segmentation, unscaled | 5 / 6 | 18% / 19% | 29% / 36% |
  | Horizontal strips (Stage 0 profile) | 4 / 5 | 38% / 40% | 3% / 3% |
  | Checkerboard | 4 / 9 | 6% / 25% | 29% / 53% |

  The segmentation does no better than a plain grid with a similar number of regions. The method assumes uniform subdomains separated by sharp interfaces; here the θ bands have soft edges spread over about 1 mm and the λ_z streaks are oblique. A smooth spatial description (penalized patches, or a few smooth basis functions in z and θ) may suit this tissue better than discrete subdomains.

### Stage 2: equilibrium-based patch identification (main method)

Use the Virtual Fields Method on the membrane form, in the reference configuration.

**Unknowns:**
- Shared global parameters: by default C2_el, C2_smc, C2_col and γ, which are not confounded by H.
- Per-patch scalar multipliers, with a smoothness penalty between neighbors whose weight is chosen by cross-validation or an L-curve.

**Parameterization ladder** (move up a level only if it is supported):

- **(a)** One membrane-stiffness factor per patch, α_k = (H_k/H)·(modulus scale), applied to Cgr, C1_el, C1_smc and C1_col together. Because the scaling is uniform, this is the same as an effective local thickness.
- **(b)** Two factors per patch: matrix/elastin (Cgr, C1_el) and fibers (C1_smc, C1_col).
- **(c)** Separate per-patch C1_smc and C1_col, with the exponents and γ shared.

**Equation blocks** (one equation per frame; the resultant blocks also run over cuts):

1. **Section-cut resultants (z direction).**
   - For vertical cuts at constant z across the domain: ∫ N_zz dθ over the cut (N already includes H) = F_z·(L_cut/w_θ).
   - Use first PK membrane tension `first_pk_membrane(F, P_k, H)`, where P_k is the patch's parameter set. A cell split between patches contributes its F to each patch, weighted by its area share.
2. **Section-cut resultants (θ direction).** The same with cuts at constant θ: ∫ N_θθ dz = F_θ·(L_cut/w_z).
3. **Interior equilibrium.**
   - Use virtual fields u* that vanish on the domain boundary, for example a sine or bubble basis per direction: Σ_points N : ∇u* · dA = 0.
   - This block needs no width assumption, and it is the one most sensitive to *relative* patch stiffness.

**Fit and weighting:**
- Least squares over all blocks, frames and protocols.
- Weight each block by its noise level, estimated from DIC noise in the rest frames passed through the model.
- Seed from the full-domain baseline. Keep C1_el above EL_MIN.

**Uncertainty:** protocol-level bootstrap, as for the baselines. Also run sensitivity checks on:
- frame subsampling;
- strain window (61 / 75 / 91 px);
- patch grid (for example 6 × 7 vs 10 × 10) and the smoothness weight;
- width assumption (tine span vs specimen);
- edge handling (with and without the patches next to the tines).

### Alternative to Stage 2: finite element model of the tracked domain (Barocas-group approach)

The Barocas group (GAIM: Raghupathy & Barocas 2010; Raghupathy et al. 2011; Nagel et al. 2014; Witzenburg & Barocas 2016; Pearce & Witzenburg 2024) fits a finite element model to the full-field displacements and the grip forces together; Nagel et al. showed that neither data source alone identifies a unique parameter set. Adapted to this test:

- Model only the region between the tine tips and impose the **measured boundary displacements** (the tracked points reach 5–15 px from the tips). Whatever the tines do to the tissue then enters as data, so tine-driven strain concentrations inside the domain are reproduced by the model instead of being read as soft material. Use boundary displacements, not F, because F is noisy where the strain window is truncated.
- Match the **reaction forces** on the boundary lines facing the rakes to the measured forces scaled by the domain's share of the tine span: 9.95 / 9.57 ≈ 1.04 of F_z and 8.52 / 9.02 ≈ 0.94 of F_θ (test1). These two fractions can also be estimated, so the interior displacements fix the relative stiffness and the forces fix only the scale, or bracketed with the specimen-width sensitivity check.
- Their samples are cruciforms gripped by the arms, which makes each arm carry the full grip force; a square sample on rakes has no such cut, which is why the two fractions are needed. (This description of their boundary handling is from the segmentation paper's text and the other papers' abstracts; the full GAIM papers were not read.)
- Spatial variation: subdomains (from the segmentation or a coarse grid) or penalized patches, with exponents and γ shared and λ chosen by protocol cross-validation.
- Needs a finite element solver (for example FEniCSx), which is not installed in this environment.

A full-specimen model with the tines is not needed for this; the measured boundary displacements replace it.

### Related work

- Bersi, Bellini, Di Achille, Humphrey, Genovese & Avril (2016), J Biomech Eng 138:071005: regional parameters of murine aortas from panoramic DIC with an inverse method based on the principle of virtual power (extension–inflation tests). Only the abstract was available; how variation between regions is constrained was not confirmed. In inflation tests the pressure constrains the local wall stress, which planar biaxial tests on rakes do not.
- Witzenburg, Dhume, Lake & Barocas (2016), IEEE TMI 35:29–41: deformation-gradient-jump segmentation (read in full; implemented in `segment_fjump.py`).
- Pearce & Witzenburg (2024), J Biomech Eng 146:121006: GAIM with an orthotropic constraint and full-field thickness from laser micrometry, which separates thickness from stiffness (constraint 2).

## 9. Verification built into each stage

- **Consistency with the baseline:**
  - With every patch factor tied (α_k = 1), the Stage 2 machinery must reproduce the full-domain dic_affine baseline within its bootstrap 95% interval.
  - Use `baseline.evaluate(..., tracks="tine_region")`, and also compare against the local version of the stress.
- **Synthetic recovery:**
  - Generate F fields for a known heterogeneous material (for example a patch pattern resembling the Stage 0 concentration maps) under the measured boundary stretches. Add noise matched to the DIC rest-frame noise.
  - Confirm that Stage 2 recovers the patch factors at the chosen grid and smoothness weight. Confirm also that Stage 1 shows the expected bias in z.
- **Scoring every local fit on the bulk data:** the patch-averaged response, evaluated with `baseline.evaluate(..., tracks="tine_region")` (for example from area-weighted homogenized parameters, or by refitting a single bulk model to the patch model's predicted resultants), should give an RMS no worse than 0.431 kPa.

## 10. Decisions for the user in the new thread

1. **Spatial description for the equilibrium-based fit:** penalized 10 × 10 patches, the 5-subdomain jump segmentation, a coarser grid, or a few smooth basis functions in z and θ.
2. **Method:** Stage 2 (virtual fields, no finite element solves) or a finite element model of the tracked domain with measured boundary displacements (needs FEniCSx or similar).
3. **Tine-adjacent margin:** own subdomains, kept in the model but out of the displacement objective, or excluded by moving the model boundary about 1–1.5 mm inward.
4. **Force fractions:** fixed from the tine-span widths (test1: 1.04 of F_z, 0.94 of F_θ) or estimated.
5. **Free per-patch parameters:** start at ladder level (a), (b) or (c).
6. **Reporting thickness:** H·C products, or moduli at the nominal H.

**Suggested first message in the new thread:**

> "Using @local_estimation_handoff.md, set up Stage 2 (virtual fields) on the interior 10 × 10 patches with a smoothness penalty chosen by protocol cross-validation, and test it on synthetic data first. Don't fit the measured data yet."
