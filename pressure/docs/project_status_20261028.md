# Project status, 28 October 2026

Start here next session. The code is at **v1.2.0**, everything is on `main`,
and nothing is running.

## Done

| Deliverable | Result | Where |
|---|---|---|
| Pressure PPSDs | 846 unique pressure channel epochs, 30 networks | `~/Desktop/BRUIT-FM/fresh_run_2026-10-23/` |
| Seismometer PPSDs | 1,345 channel epochs at 378 pressure stations | `seismic_run_2026-10-26/` |
| QC tests 1–3 | all 846 channels; 456 good / 171 questionable / 189 bad | `qc_all_2026-10-26/`, `docs/reports/2026-10-26_all_station_qc/` |
| Calibrations | Rayleigh 54, microseism 167 (scored), local P 54 channels; methods agree (median +0.1 dB, Rayleigh vs microseism) | `docs/reports/2026-10-28_seismic_calibration/` |
| Final class (6 scores) | 439 good / 188 questionable / 189 bad | `seismic_run_2026-10-26/table/pressure_qc_calibration_table.xlsx` |
| Provisional steps 1–3 | policy frozen; metadata flags 2A, YD, Z4, X9, 4N; 396 accepted at 311 stations | `docs/reports/2026-10-28_provisional_decisions/` |
| Noise model (no gain) | 311 stations; bands, heat maps, station lines by basin, depth, sensor, season | `docs/reports/2026-10-28_provisional_decisions/noise_model/` |
| **Gain-corrected campaign** | 265 channels × calibration g; 464 good; X9 cleared; 421 accepted at 323 stations; noise model and all figures redone | `docs/reports/2026-10-28_gain_corrected/README.md`, `gain_corrected_2026-10-28/` |
| Package | installable, GPL-3.0, CLI commands, 503 tests, releases v1.0.0–v1.2.0 | `README.md`, `CHANGELOG.md` |
| Results report | methods, figures, worked example (one figure per QC test and calibration) | `docs/REPORT.md` |

Skipped by the owner's decision: network-wide tilt correction (pilot only:
RR52 −14 dB at 2–10 mHz) and the raw-data Godin pilot.

## Decisions waiting for the supervisor (Wayne Crawford)

1. **Thresholds** in `policy_v1_provisional.json`: QC bounds, class rules,
   calibration score (1.5/6 dB), the 7-day minimum and the 15 dB metadata flag.
2. **Gain correction**: whether to use the gain-corrected results; the method
   order (Rayleigh → microseism → local P); sharing within a gauge; the 0.1–10
   limit, which leaves 2A, YD and Z4 uncorrected although their estimates
   (×220–900) agree with their level offsets.
3. **Metadata networks** 2A, YD, Z4 (and 4N, 1T, X3): checking against the
   operators' calibration sheets means contacting outside groups, which needs
   approval first.
4. **The 198–228 channels marked *review*** and the accepted list.
5. **Open science points**: local P is about 0.44 above 4 Hz at the YV DPGs;
   APGs need their own model; Godin excess of +10–20 dB in the Pacific and
   Indian oceans.

## Rerunning after the decisions

Edit the rules in `gain_corrected_2026-10-28/decisions/` (or in
`provisional_decisions_2026-10-28/` for the uncorrected route), set
`--label approved` in `run.sh`, and run `bash run.sh`. A full rerun takes about
10 minutes.

## Housekeeping

- Remote: only `main`. Tags: `v0.1.0`, `v1.0.0`, `v1.1.0`, `v1.2.0`,
  `archive/local-p-check` (its useful change is already in main).
- The repository is private. Making it public is the owner's decision.
- BRUIT-FM will later be merged into ComPy (both GPL-3.0).
- The gh token lacks `workflow` scope, so `.github/workflows` cannot be pushed.
- The data drive is `/Volumes/DATA_bruitfm`. Connect it directly with its own
  power supply; it dropped twice on a hub.
