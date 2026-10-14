# Constitutive fits

The bulk four-fiber family model (Jadidi et al. 2019, Eq. 2) is fitted to the stretch-ratio
protocols Def11–Def19 (PhaseTiming Sets 2–10), and local material properties are estimated
from the DIC strain fields. `local_estimation_handoff.md` documents the local method, its
status and its findings.

Every script needs a test (`--test SAMPLE/testN`, e.g. `--test N1739/test1`, or `DIC_TEST`; see the top-level `README.md`).
The thickness and widths default to that test's `test_config.json`. The scripts import
`dic2d` and `testcfg` from `../strain_tracking/`, and they read the tracks in `<test>/DIC_results/`.

## Outputs, per test (`<test>/constitutive_fit/`)

```
constitutive_fit/
├── params_by_stretch_input.csv / .png, stretch_inputs_data.csv, stretch_inputs_histories.png
│                            bulk fits for each global stretch input, full-domain tracks
├── fit_stress_stretch_<input>.png           data and model, every protocol, one plot per input
├── baseline/                frozen bulk baseline, full-domain tracks (primary)
├── standard_roi/            the same on the standard ROI (--tracks standard), with its baseline/
├── fit_4pt_*                bulk fit to the CellScale 4-point tracking
└── local/                   local estimation: local_data.npz, Stage 0 figures,
                             stage1*/, segmentation/, stage2_*/
```

## Scripts

### Bulk fits

| File | Purpose |
|---|---|
| `constitutive.py` | Four-fiber model: strain energy, plane-stress Cauchy stress for a general in-plane F, closed-form biaxial stress, and membrane first PK tension. `python constitutive.py` runs the self-checks |
| `compare_stretch_inputs.py` | Fits the model for each global stretch input: CellScale 4-point, DIC 4-point, DIC region-mean λ, DIC affine F, DIC region edges, rake. Fits stay in the elastin-active family (C1_el > 0.1 kPa), with a protocol-level bootstrap. `fit_stress_stretch_<input>.png` (one per input) plots the measured Cauchy stress (symbols: circles axial, triangles circumferential) and the model at the same stretches (solid lines) for every protocol on one plot; the plotted curves are lightly smoothed (5-frame running median, then an 11-frame Savitzky–Golay filter), the fit uses the raw data. All the plots of a sample share their axis limits. Options: `--tracks tine_region` or `standard`, `--width tine_span` or `specimen`, `--table-only`, `--plot-only` (redraws the plots from the saved results), `--plot-input` (default: one plot per input) |
| `baseline.py` | `freeze` stores a frozen baseline: parameters, bootstrap 95% intervals, per-protocol errors, the exact data set, settings and checksums. `evaluate` scores any parameter set against it. Pass `tracks="tine_region"` for the primary baseline |
| `check_preprocessing.py` | Sensitivity of a baseline to force taring and per-protocol Savitzky–Golay smoothing |
| `fit_4pt.py` | Bulk fit to the CellScale 4-point tracking, with a frame-level bootstrap and a γ profile → `fit_4pt_*` |

### Local estimation

| File | Purpose |
|---|---|
| `build_local_data.py` | Stage 0: builds `local_data.npz`. It holds F per frame on the full-domain grid in the specimen frame, tared forces, patches and checksums, plus the band diagnostics. No fitting |
| `stage1_patches.py` | Stage 1: one stiffness factor per patch from the patch-mean F and the global stress (uniform stress assumed). `--exclude-ring N` |
| `stage1_smooth.py` | Stage 1 with a neighbor smoothness penalty on log α. The penalty weight is chosen by leave-one-protocol-out cross-validation |
| `segment_fjump.py` | Splits the domain into subdomains wherever F jumps between neighboring elements (Witzenburg et al. 2016) |
| `membrane_fe.py` | Finite element solver for the incompressible plane-stress membrane on the DIC grid (Q1 elements, prescribed boundary displacements) |
| `stage2_setup.py` | Stage 2 equations and conditioning checks. The equations are section cuts plus interior virtual fields. No fitting |
| `stage2_synthetic.py`, `stage2_synthetic_b.py`, `stage2_synthetic_c.py` | Synthetic recovery tests on finite element data. They cover one factor per patch, separate matrix and fiber factors, and exponents that vary by patch. Finite element solutions are cached as `fe_*.npz` |
| `stage2_fit.py` | Stage 2 fit of patch stiffness factors on the measured data, with a bootstrap, a bulk score and a finite element forward check |
| `stage2_residuals.py` | Stage 2 misfit by protocol, load level and virtual-field mode, and cross-prediction between loading-ratio groups |
| `stage2_exponent_fit.py` | Fits a stiffness factor and an exponent factor per patch (`--synthetic NAME` runs the recovery check) |

## Usage

```bash
../.venv/bin/python compare_stretch_inputs.py --test N1739/test1 --tracks tine_region --table-only
```

```bash
../.venv/bin/python baseline.py --test N1739/test1 evaluate --params "6.24,3.61,1.88,10.1,10.7,1.31,7.75,47.7"
```

```python
import sys; sys.path[:0] = ["../strain_tracking", "."]
import testcfg; testcfg.select("N1739/test1")
import baseline
baseline.evaluate([...], tracks="tine_region")
```
