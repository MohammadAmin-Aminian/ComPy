# ComPy end-to-end pipeline

From a broadband OBS network in an SDS archive to compliance, shear-velocity
models and serpentinization, unattended. The method is that of Aminian et al.
(2025, GJI, Fig. 2), with the corrections of ComPy 2.2–2.4 (`docs/AUDIT_2026-10.md`,
`docs/INVERSION.md`).

All network-specific values are in [`config.py`](config.py). Set
`COMPY_RUN_DIR` to the run folder, and set `COMPY_PIPELINE_CONFIG` to a file
that overrides any of its names. The defaults reproduce the XE (Pacific Array
Young ORCA) run in [`docs/xe_young_orca_2026-10-29`](../docs/xe_young_orca_2026-10-29/README.md).

## Stages

| Stage | Script | What it does | Output (in the run folder) |
|---|---|---|---|
| 0 | `fetch_inventory.py`, `fetch_catalog.py` | StationXML with Z, 1, 2 and pressure responses; USGS catalogue M ≥ 5.85 | `metadata/<NET>.xml`, `work/catalog_*.qml` |
| A | `a_preprocess.py` (`a_fast.py` for slow disks) | Decimation to 1.25 Hz with FIR filters. Clipped (≥ 99 % of 2²³ counts) or dead (σ < 5 counts) 10-min blocks are cut from all channels (±5 min). Responses are removed with no water level and an explicit pre-filter. | `work/<STA>/{counts,phys,qc}/` |
| B | `b_clean.py` | Earthquake spans (1.5 d per magnitude unit above 5.85), optional hourly-glitch removal, hourly tilt smoothed over 24 h and applied per sample, STA/LTA local events, horizontal-coherent noise removed in ≤ 24 h windows | `work/<STA>/clean/` |
| D | `d_calibrate.py` | DPG gain from Rayleigh-wave P/a with `compy_pressure`: native-rate SDS data, Mw ≥ 6 events by magnitude level, bins with squared coherence ≥ 0.95, median of ≥ 3 event gains and their spread | `work/<STA>/calibration/` |
| C | `c_compliance.py` | Non-overlapping 35-min windows. Gates: median coherence ≥ 0.9, or 0.8 if fewer than 500 windows pass; pressure 10–50 dB; acceleration −170 to −150 dB. Coherent band, stacked gravity-corrected compliance, B&P error, monthly stacks | `results/<STA>/<STA>_compliance_*.{npz,csv}`, `_summary.json` |
| E | `e_invert.py` | Bayesian Vs inversion with DRAM and tempering: free crustal Vp/Vs, gain sampled within ±5 %, CRUST1.0 Moho. Adds a prior-only run, the Vp/Vs 1.85 variant (`main`) and the H2 bracket | `results/<STA>/inversion/` |

`run_bdc.py` starts B, D and C for each station as soon as its stage A is
done, and `run_inversions.py` runs E as compliance results arrive. Several
inversion queues can run at once (`--reverse` takes stations from the end;
lock files give each station to one queue). `f_figures.py` draws the figures
of Aminian et al. (2025) for the network (map, PSD stages, tilt, calibration,
coherogram, selected windows, compliance, Vs heat maps, profiles and
serpentinization), and `g_report.py` collects the numbers of a report.
`a_download.py` replaces `a_fast.py` when the raw data must come from an FDSN
dataselect service.
Serpentinization is computed from the posterior draws with
`compy_inversion.serpentinization` (0 % at 4.2 km/s, 100 % at 2.2 km/s).

## Running

Install `.[pressure,maps]` (or `.[dev,pressure,maps]` for development).
For a complete run, set the archive, station, response and CRUST1.0 paths in
your configuration, then use the sequential runner:

```console
export COMPY_RUN_DIR=~/runs/xe
export COMPY_PIPELINE_CONFIG=/absolute/path/to/my_config.py
python pipeline/run_pipeline.py CC11 CC05 --workers 8 --processes 8 --fetch
```

This runs A → B → D → C → E → F → G and writes individual logs plus a JSON
run manifest under `logs/`. Omit `--fetch` when the inventory and catalogue
already exist. `--start-at D` resumes at calibration; `--stop-after C` ends
after compliance. Existing stage checkpoints are reused. Every selected
station must finish each requested stage. Missing calibration is recorded as
an explicit fallback to the alternative or nominal gain, with the existing
wider inversion prior; failed preprocessing or empty compliance stops the
run with a nonzero exit code. Figures and report tables use the same selected
stations. A completed workflow does not establish posterior convergence;
inspect the saved R-hat, ESS and fit diagnostics.

The original concurrent schedulers remain available for large campaigns:

```console
export COMPY_RUN_DIR=~/runs/xe   # holds metadata/, work/, results/, logs/
cd ComPy/pipeline
python fetch_inventory.py && python fetch_catalog.py
python a_preprocess.py CC11 CC05 --workers 8     # or: python a_fast.py ... (USB disks)
python run_bdc.py --jobs 6 CC11 CC05             # B, D, C per station
python run_inversions.py --processes 8           # E
```

CRUST1.0 (`crust1.bnds`, `.vp`, `.vs`, `.rho`) must be in `metadata/crust1`.
Stage E spawns processes. Run it from these scripts, not from stdin.

## Choices to revisit for a new network

- `FACTORS`: their product must equal the raw rate divided by 1.25 Hz (62.5 Hz: 5·5·2; 50 Hz: 5·4·2).
- `DEGLITCH`: on for instruments with a periodic transient (e.g. the hourly
  mass-position glitch of the INSU-IPGP BBOBS). The fold strength in
  `clean/deglitch.json` shows whether one is present.
- **Gain source.** Stage E inverts the ComPy Rayleigh gain where it exists. Otherwise it uses the BRUIT-FM gain with a wider prior (±15 %).
- **Stage D settings** (`CALIBRATION_BAND`, `CALIBRATION_TILT`, `CALIBRATION_DISPERSION`):
  - `CALIBRATION_TILT = "tiskitpy"` removes tilt noise from the vertical but
    needs horizontals with valid StationXML orientations;
  - `"none"` uses pressure and vertical only, which XE needs (its BH1 and BH2
    have a dip of −90°);
  - the default dispersion is RR52's, BRUIT-FM's default for every network.

  After ComPy 2.4.0, stage D is no longer `Pressure_calibration`'s method,
  so gains of new runs can differ from those of 2.4 runs.
- **Clipping thresholds** (`a_preprocess.CLIP_COUNTS`, `DEAD_STD_COUNTS`) assume 24-bit digitizers.
- **Stations without any calibration** are computed with the nominal response (`nominal`), and stage E samples the gain over 0.6–1.4. Report them separately.
- **The 500-window minimum** of the paper is applied when the results are collected (`f_figures.MIN_WINDOWS`). Stage C itself always returns its best attempt.
