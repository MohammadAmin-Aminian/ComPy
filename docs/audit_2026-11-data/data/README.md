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

## Research part (RESEARCH.md)

| File | Script | Content |
|---|---|---|
| `monthly_tide.csv` | `r01_monthly_tide.py` | M2, S2, K1 amplitudes in the pressure counts, 30-day windows every 15 days |
| `tide_at_steps.csv` | `r01_monthly_tide.py` | M2 change across every \|Z/P\| step, with a verdict (P / not P / unclear) |
| `tide_phase_monthly.csv` | RESEARCH.md §1 | Monthly M2 vs compliance-band phase and \|Z/P\| per station |
| `coherence_stations.csv`, `coherence_episodes.csv`, `coherence_hole.csv` | `r02_coherence_search.py` | Coherence per band per station; network-synchronous episodes; minimum coherence just above f_c |
| `earthquake_epochs.csv` | `r03_earthquake_epochs.py` | Tilt and \|Z/P\| changes around earthquakes, with z-scores against random times |
| `tide_constituents.csv` | `r04_tide_constituents.py` | O1, K1, M2, S2 ratios after deconvolution and the slope n of R(f) |
| `response_bracket.csv` | `r05_response_bracket.py` | Leak-corner factor from the tide vs from the phase, predicted phase, compliance bias |
| `dpg_leak_factors.csv` | RESEARCH.md §6 (D2) | Recommended leak-corner factor per DPG (tide where clear, else phase) |
| `yv_corrected_inversion.csv` | `r06_corrected_inversion.py` | YV stage E with the response correction vs the published run |
| `poor_fits.csv` | `r08_poor_fits.py` | 7D.J35A and J44A: band capped at 0.75 f_c, with and without the correction |
