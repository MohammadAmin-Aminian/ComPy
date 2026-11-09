# Data of the whole-archive audit

Small tables copied from the run folder `~/Desktop/ComPy/13_audit` (not in the repository:
per-window spectra cache, about 6 GB, and the restacked stage C results). Each file names
the script that made it. Stations are `NET.STA`; times are UTC; phases in degrees.

| File | Script | Content |
|---|---|---|
| `stations.csv` | `s05_anomaly_scan.py` | One row per inverted station: tier and quality, band, windows, gate bias (exclusive vs inclusive), \|Z/P\| changepoints with the sensor that stepped and the stage D event-gain ratio (`steps`, JSON), epoch shift, residual phase, pressure-corner factor, tidal ratio and implied compliance bias |
| `anomalies.csv`, `anomalies.md` | `s05_anomaly_scan.py` | Ranked anomaly records (`compy.anomaly` plus the cross-checks of s05) |
| `network_info.json` | `s05_anomaly_scan.py` | Per network: common mode of monthly deviations, synchronous steps, phase-sign test |
| `phase_metadata.csv` | `s02_phase_metadata.py` | Residual Z/P phase, model fits, instrument family and metadata corners and in-band phase of Z and P |
| `dpg_tide.csv` | `s03_dpg_tide.py` | Corner factor k fitted to the phase (two models), predicted and measured M2 ratio, K1 ratio, SNR, EOT20 amplitude |
| `level_steps.csv`, `dated_steps.csv` | `s04_level_steps.py` | Changepoints of each station's Z and P levels minus its network median (microseism and compliance band), and steps at given dates |
| `tilt_tides.csv` | `s06_tilt_tides.py` | Tidal lines (M2, M4, K1, S1) in hourly tilt angle and tilt-correction variance reduction |

Inputs: the whole-archive run (`00_inventory`, `01_preprocessed/*/counts`, `02_cleaned`,
`03_tilt`, `04_ppsd`, `06_compliance*`, `12_inversion_all_stations`) and the EOT20 tide
model (Hart-Davis et al. 2021, doi:10.17882/79489) from the pressure-atlas cache.
