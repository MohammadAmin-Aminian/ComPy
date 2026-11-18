# Pacific Array Young ORCA (XE): latest results, November 2026

30 Scripps broadband OBS (Trillium 240 and Cox-Webb DPG at 50 Hz), 4386-4824 m of water on about 40 Ma Pacific
seafloor, April 2018 to May 2019. These results replace the run of 29 October 2026
([xe_young_orca_2026-10-29](../xe_young_orca_2026-10-29/README.md), kept as the record of ComPy 2.4-2.8).

| Step | Run | Code |
|---|---|---|
| Stages A-D (cleaning, tilt, Rayleigh-wave gain, compliance) | 13 November 2026, from raw data | ComPy 3.0 |
| Stage E (Vs inversion), final model `geo` | 16 November 2026 | ComPy 3.1 (`auto`, ordered Vs, Moho rule on the whole crust) |

The final inversion adds two priors to the default: the Moho from the complete Bouguer anomaly (crust 6.5 ± 1 km
where the anomaly is above 320 mGal, 7.5 ± 1 km below), and the sediment thickness from DSDP Sites 73-75 (σ 50 %).
No low-velocity zone is allowed at any station. CC07 uses only the windows before its seismometer dropped from the
release arm (22 November 2018; Eilon et al. 2021).

## Numbers

- 27 stations reached stage C and 26 have a compliance estimate; 15 have a Rayleigh-wave DPG gain. All 26 were
  inverted; 20 are interpreted.
- Median Vs of the 20 interpreted stations: **3.75 km/s at 2-5 km below the seafloor**, 3.10 km/s at
  0.7-1.5 km below the basement; sediment 87 m; Moho 7.4 km below the seafloor (12.0 km below sea level).
- Table of every station: [`tables/XE_final_model.csv`](tables/XE_final_model.csv) (Vs at fixed depths with
  90 % intervals, sediment, Moho, mantle Vs, fit, quality label).

## Figures

| File | Content |
|---|---|
| `figures/fig01_map.png` | Network and compliance status |
| `figures/fig03_psd_stages.png` | CC11 vertical and pressure spectra before and after cleaning |
| `figures/fig04_tilt.png` | CC11 hourly tilt |
| `figures/fig05_calibration.png` | CC11 Rayleigh-wave DPG gain |
| `figures/fig06_coherogram.png` | CC11 spectrograms and Z-P coherogram |
| `figures/fig07_windows.png` | CC11 selected windows |
| `figures/fig08_compliance.png` | Compliance of every station by line |
| `figures/CC11_inversion.png` | CC11 final inversion: fit, Vs, what the data constrain |
| `figures/profiles.png` | Vs profiles of every station |
| `figures/section.png` | Median Vs of the interpreted stations, west to east |
| `figures/map_results.png` | Maps of Vs, sediment and Moho |

The full reports (compliance of 13 Nov, inversion and gravity of 16 Nov) were sent with the runs; their figures and
tables are reproduced here.
