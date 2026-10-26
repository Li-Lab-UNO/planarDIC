# Constitutive fits

The bulk four-fiber family model (Jadidi et al. 2019, Eq. 2) is fitted to the stretch-ratio
protocols Def11–Def19 (PhaseTiming Sets 2–10), and local material properties are estimated
from the DIC strain fields in stages, each tested on synthetic data before the measured data (below).

Every script needs a test (`--test SAMPLE/testN`, e.g. `--test N1739/test1`, or `DIC_TEST`; see the top-level `README.md`).
The thickness and widths default to that test's `test_config.json`. The scripts import
`dic2d` and `testcfg` from `strain_tracking/`, and they read the tracks in `<test>/DIC_results/`. Commands are
run from the project folder.

## Outputs, per test (`<test>/constitutive_fit/`)

```
constitutive_fit/
├── params_by_stretch_input.csv / .png, stretch_inputs_data.csv, stretch_inputs_histories.png
│                            bulk fits for each global stretch input, full-domain tracks
├── fit_stress_stretch_<input>.png           data and model, every protocol, one plot per input
├── baseline/                frozen bulk baseline, full-domain tracks (primary)
├── standard_roi/            the same on the standard ROI (--tracks standard), with its baseline/
├── fit_4pt_*                bulk fit to the CellScale 4-point tracking
└── local/                   local estimation: stage0_*, stage1_*, stage1b_*, stage2_* (.csv / .json / .npz / .png)
```

## Scripts

### Bulk fits

| File | Purpose |
|---|---|
| `constitutive.py` | Four-fiber model: strain energy, plane-stress Cauchy stress for a general in-plane F, closed-form biaxial stress, and membrane first PK tension. `python constitutive.py` runs the self-checks |
| `compare_stretch_inputs.py` | Fits the model for each global stretch input: CellScale 4-point, DIC 4-point, DIC region-mean λ, DIC affine F, DIC region edges, rake. Fits stay in the elastin-active family (C1_el > 0.1 kPa), with a protocol-level bootstrap. `params_by_stretch_input.png` shows every parameter with its bootstrap interval, and each input's R² and RMSE in its key; `fit_stress_stretch_<input>.png` (one per input) plots the measured Cauchy stress (symbols: circles axial, triangles circumferential) and the model at the same stretches (solid lines) for every protocol on one plot; the plotted curves are lightly smoothed (5-frame running median, then an 11-frame Savitzky–Golay filter), the fit uses the raw data. All the plots of a sample share their axis limits. Options: `--tracks tine_region` or `standard`, `--width tine_span` or `specimen`, `--table-only`, `--plot-only` (redraws the plots from the saved results), `--plot-input` (default: one plot per input) |
| `baseline.py` | `freeze` stores a frozen baseline: parameters, bootstrap 95% intervals, per-protocol errors, the exact data set, settings and checksums. `evaluate` scores any parameter set against it. Pass `tracks="tine_region"` for the primary baseline |
| `check_preprocessing.py` | Sensitivity of a baseline to force taring and per-protocol Savitzky–Golay smoothing |
| `fit_4pt.py` | Bulk fit to the CellScale 4-point tracking, with a frame-level bootstrap and a γ profile → `fit_4pt_*` |

### Local estimation

Stages from the bulk fit toward local properties. Every stage uses the fixed eight-parameter model and the test's
bulk fit (DIC affine F); the finite element stages model the tracked domain with the measured boundary positions
prescribed, and compare strains in the interior (`edge_mm` of `analysis_settings.json`).

| File | Purpose |
|---|---|
| `stage0_identifiability.py` | Stage 0: profile of the bulk misfit for each parameter (the other seven refitted): the range within +1% and +5% RMS, whether it reaches a bound, and the parameter it trades off with |
| `stage1_uniform_fe.py` | Stage 1: a uniform material (the bulk fit) in the finite element model. Share of the interior and edge-band strain variance it explains, rest-frame noise bound, force resultants across mid-domain cuts; peak maps of every protocol |
| `stage1b_best_uniform.py` | Stage 1 check: the uniform material whose response shape (seven parameters, C1_smc fixed) best fits the interior strains of Def11, Def14 and Def18, scored on the other six protocols and on the bulk forces |
| `stage2_alpha_field.py` | Stage 2: a smooth stiffness factor α(x) on the whole response (Cgr and every C1; in a membrane the same as local thickness), log α a Legendre polynomial basis (total degree, or separate degrees in z and θ such as `2x6`). `--synthetic NAME` for recovery checks (`uniform`, `smooth`, `bands`, `bands_theta` positive control, `bands_shape` model-error control); `--groups` fits each loading-ratio group separately and cross-predicts (direction test), `--groups sequence` the equibiaxial protocols in run order |
| `membrane_fe.py` | Finite element solver for the incompressible plane-stress membrane on the DIC grid (Q1 elements, prescribed boundary displacements, optional stiffness factor per Gauss point) |
| `solid_fe.py` | Finite element solver for the same model as a 3D slab on the DIC grid: hexahedral elements, several layers through a thickness that can vary over the grid, near incompressibility by a mean-dilatation penalty, side faces prescribed through the full thickness, top and bottom faces free; Newton with a sparse direct solve, started from a membrane solution or a coarser layering. `python constitutive_fit/solid_fe.py` runs the self-checks (stress and forces against the energy, homogeneous stretch against the membrane model) |

`legacy/` holds the earlier local scripts (uniform-stress patch fits, deformation-gradient-jump segmentation,
virtual-field equations) and their notes. They are not maintained and write to the same `local/` folder.

## Usage

`./runDIC.sh` fits every test it runs (its `fit` step). Fit one test (every stretch input, bootstrap, plots → `<test>/constitutive_fit/`):

```bash
.venv/bin/python constitutive_fit/compare_stretch_inputs.py --test N1739/test1
```

Redraw the parameter and fit plots from the saved results, without fitting:

```bash
.venv/bin/python constitutive_fit/compare_stretch_inputs.py --test N1739/test1 --plot-only
```

Freeze the result as the test's baseline, then score another parameter set against it:

```bash
.venv/bin/python constitutive_fit/baseline.py --test N1739/test1 freeze
```

```bash
.venv/bin/python constitutive_fit/baseline.py --test N1739/test1 evaluate --params "5.71,4.07,1.69,9.88,10.5,0.991,7.93,45.6"
```

```python
import sys; sys.path[:0] = ["strain_tracking", "constitutive_fit"]
import testcfg; testcfg.select("N1739/test1")
import baseline
baseline.evaluate([...])
```
