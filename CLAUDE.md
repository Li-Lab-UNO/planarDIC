# Project conventions

## Writing
- American English everywhere: chat replies, code, comments, printed output, figure text, READMEs
  (center, centerline, color, fiber, gray, normalize, analyze, modeling, labeled).
- Write subscripts with underscores (lambda_z, E_ztheta, lambda_1), not Unicode subscripts.
- READMEs describe the current state; no changelog voice ("now", "previously", "was changed to").
- READMEs don't enumerate samples or keep per-sample status (samples keep being added): per-test caveats go in
  `notes` of test_config.json, and `strain_tracking/status.py` reports status from the files on disk.

## Analysis
- Analyze and fit Def11–Def19 only; Def01–Def10 are preconditioning.
- No constitutive fitting unless asked.
- Physical settings (subset, step, strain window, interior edge, low-spot box, resampling target) live in
  `analysis_settings.json`; the interior edge `edge_mm` is the user's choice, reported, never changed silently.

## Working
- Ask before long batch runs (tracking, multi-test exports); report progress as deformations finish.
- Never script PowerPoint or Keynote (AppleScript open/export/close) on the user's documents. Ask whether a
  deck is open before editing it with python-pptx.
