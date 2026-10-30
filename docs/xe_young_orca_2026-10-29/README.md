# XE (Pacific Array Young ORCA) compliance, Vs and serpentinization, 29–30 October 2026

> **Note.** `scripts/d_calibrate.py` imports `Pressure_calibration`, which
> ComPy removed after 2.4.0. Rerun it with ComPy 2.4.0 (commit `7ad79c2`). The
> gains here are those of that method; the pipeline's stage D now uses
> `compy_pressure`.

This is the second network processed with ComPy, after YV (RHUM-RUM; Aminian et al. 2025). It is the validation run of
the end-to-end `pipeline/` of ComPy 2.4.

- Network: 30 Scripps broadband OBS (Trillium 240 + DPG, 50 Hz), 4386–4824 m water depth, April 2018 – May 2019, on
  ~40 Ma Pacific seafloor (Eilon et al. 2022).
- Data: the BRUIT-FM SDS archive (`/Volumes/DATA_bruitfm/SDS`), with responses from EarthScope (`metadata/XE.xml`).
- Report: `report/XE_young_orca_compliance_report.pdf`, a draft in the layout of Aminian et al. (2025), with its
  LaTeX source.

## Stations

| Outcome | Stations |
|---|---|
| Not processed | CC02 (DPG not working), CC10 (50 days), CC12 (17 days) |
| No 13-h span with all four channels usable (clipped or dead seismometer) | see `results/tables/report_numbers.json` (`status`) |
| Compliance from fewer than 500 windows (rejected) | EC04, WC01 |
| No Rayleigh-wave calibration; nominal response, gain sampled over 0.6–1.4 | CC09 (inverted, not pooled) |
| Inverted | listed in Table 1 of the report |

## Choices specific to this network

- **Clipping.** Many seismometers sat at full scale (2²³ counts) for weeks. Stage A cuts every 10-min block with a
  clipped or flat sample from all channels (±5 min) before decimation. This removed about a quarter of the record.
- **No glitch removal.** No hourly glitch: the fold strength was at most 1.76, below 2.0, in 584 blocks at 6 stations
  (`config.DEGLITCH = False`).
- **Gains.**
  - The ComPy Rayleigh-wave gain is the one used. It agrees with the BRUIT-FM microseism gain where both exist (median ratio ≈ 0.99).
  - BRUIT-FM's own Rayleigh calibration failed for XE (`unresolved_metadata`).
- **Variants.** Every inverted station has the final model (Vp/Vs free, gain sampled ±5 %, tempering), a prior-only run, and the Vp/Vs = 1.85 variant (`main`). CC04 and CC05 also have the H2 bracket.
- **Staging.** Raw files were copied to local disk station by station (`a_fast.py`) because parallel reads from the USB
  archive collapsed. A short attempt to download from EarthScope was throttled after a few hundred files and was abandoned.

## Main results

See the report. In short:

- Compliance:
  - 4.9–14.6 mHz;
  - coherence 0.93–0.97;
  - thousands of windows per station;
  - month-to-month changes of 1–3 %, except CC07, which drifts by about ±10–25 %.
- Vs: 80–170 m of sediment, then about 2.4 km/s at 0.3 km, 2.9 km/s at 1 km and 3.2 km/s at 2 km.
- Vp/Vs: lower-crustal Vs depends on it (≈ 3.5 km/s free, ≈ 4.05 km/s with 1.85, against CRUST1.0 4.05 km/s).
- Serpentinization: read through the Fig. 10c relation of Aminian et al. (2025), these normal-crust velocities would mean about 60 % serpentinization at 1 km.
- At 2 km, the corrected SWIR/CIR Vs of YV falls inside the Young ORCA range. Compliance-derived Vs alone therefore does not identify serpentinite.

## Reproduce

```console
cd scripts
python a_preprocess.py <STA> ... --workers 8   # or a_fast.py for the USB archive
python run_bdc.py --jobs 6 <STA> ...           # B, D, C
python run_inversions.py --variants final      # E (and --variants main)
python f_figures.py CC11 && python g_report.py && python h_latex.py && python i_html.py <out>
```

`scripts/config.py` holds every network-specific value. Intermediate waveforms (`work/`, about 60 GB) and posterior draws
(`results/*/inversion/*.npz`) stay on the Mac and are not in git.
