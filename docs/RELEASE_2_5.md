# ComPy 2.5.0 — 30 October 2026

ComPy now includes BRUIT-FM 1.3.0 as `compy_pressure`, with the original
campaign scripts, methods and tests under `pressure/`. The dated BRUIT-FM
reports remain linked to their archived source.

## Migration and complete runs

`Pressure_calibration.py` has been removed. Use the
[calibration migration table](../README.md#pressure-calibration). The tag
`v2.4.0` preserves commit `7ad79c2`, including the old calibration module
needed by the historical YV/XE scripts. The release stays in the 2.x series
at the owner's request despite this documented API removal.

Install from the source archive or Git checkout with
`python -m pip install '.[pressure,maps]'`, set `COMPY_RUN_DIR` and
`COMPY_PIPELINE_CONFIG`, then run:

```bash
python pipeline/run_pipeline.py CC11 --workers 4 --processes 4 --fetch
```

This executes preprocessing, cleaning, Rayleigh calibration, compliance,
inversion, figures and report tables in order. See
[pipeline/README.md](../pipeline/README.md) for input paths, CRUST1.0,
resuming and explicit alternate/nominal gain fallbacks.

The wheel installs the numerical modules, `compy_pressure`, its reference
data and eight console commands. The source archive additionally includes
the complete pipeline, examples and `pressure/` project folder.

## Release validation

- Python 3.10 and 3.13: **636 tests pass** on each interpreter.
- Ruff lint, formatting and `git diff --check` pass.
- The source archive and universal wheel build successfully.
- The wheel was installed in separate environments outside the checkout
  on Python 3.10 and 3.13. Both verify version 2.5.0, installed import paths,
  packaged reference data and all eight console commands, and run the
  offline synthetic inversion and pressure-QC workflow.
- Regression coverage includes calibration station/channel/source changes,
  child-process failure propagation, station scope and excluding gains from
  failed calibrations.

The PR's original real-data stage-D comparison remains recorded in
[PRESSURE.md](PRESSURE.md): three XE stations and YV.RR52 give gains
1.4–1.9% above the previous calibration. These are prior validation runs;
this release does not regenerate the archived full-network inversions.
Scientific selection policies and open supervisor questions remain recorded
with the campaign results.
