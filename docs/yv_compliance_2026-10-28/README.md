# YV broadband seafloor compliance, 28 October 2026

> **Note.** `scripts/d_calibrate.py` imports `Pressure_calibration`, which
> ComPy removed after 2.4.0. Rerun it with ComPy 2.4.0 (commit `7ad79c2`). The
> gains here are those of that method; the pipeline's stage D now uses
> `compy_pressure`.

This folder holds the compliance functions of the eight RHUM-RUM (YV) broadband OBS stations
(RR28, RR29, RR34, RR36, RR38, RR40, RR50, RR52), computed with ComPy 2.2 on this branch.
The inversion was not run.

- `report.html`: the full report (method, audit, gains, results, open issues). Open it locally.
- `results/<STA>/`:
  - `<STA>_compliance_{bruitfm,compy}.csv`: frequency, normalized compliance [1/Pa],
    Bendat & Piersol error, coherence, phase, and an in-band flag;
  - `.npz`: the same, plus stacked spectra and accepted window start times;
  - `<STA>_summary.json`: gates, rejected-window counts and monthly stacks;
  - `epoch/`: stable-epoch results for RR36, RR38, RR40 and RR50 (recommended for those
    stations);
  - `<STA>_tiskit_check.json`: independent TiSKitPy estimate on the same windows.
- `results/tables/`: the summary table, the report table and the sensitivity-step scan.
- `scripts/`: the full chain. Paths assume the IPGP Mac (data on `/Volumes/DATA_bruitfm`).
  1. `a_preprocess.py`: decimation and responses.
  2. `b_clean.py`: events, transients, tilt and coherent noise.
  3. `d_calibrate.py`: ComPy gain.
  4. `c_compliance.py [--end]`: compliance; `--end` keeps only data before a date (stable epoch).
  5. Checks and figures: `e_figures.py`, `g_tiskit_check.py`, `h_stability.py`,
     `f_legacy_compare.py`.
- `figures/`: all figures.

The intermediate waveforms (about 11 GB) stay on the Mac in
`~/Desktop/BRUIT-FM/compliance_yv_2026-10-28/work/`.
