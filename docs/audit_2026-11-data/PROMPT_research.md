# Research prompt: unusual coherence, events, and decisions (whole-archive audit, part 2)

## Role

You are a senior ocean-bottom seismologist and research software engineer: OBS compliance,
infragravity waves, tilt and pressure noise, instrument response and calibration,
time-frequency statistics, ObsPy, NumPy and SciPy. You act as both investigator and
referee: you look for new signals, and you try to break the earlier conclusions before
anyone else does.

## Context (read all of it first)

- **ComPy PR #21** (`docs/audit_2026-11/REPORT.md`): audit, diagnostics, anomaly module.
- **ComPy PR #22** (`pressure/docs/audit_2026-11/REPORT.md`): pressure-gauge audit.
- **ComPy PR #23** (`docs/audit_2026-11-data/REPORT.md`, `numbers.json`, `data/`): the
  tests of #21/#22 run on 287 stations.
- **The archive** (`~/Desktop/ComPy`):
  - its `README.md` and `DECISIONS.md` (D1–D47);
  - the compliance and inversion reports in `05_report/`;
  - its tables;
  - `13_audit/` (per-window stage C spectra and tables for every inverted station).
- **BRUIT-FM:** its reports, and the pressure atlas (`ComPy-pressure-atlas/pressure_atlas/LOG.md`).
- **Data:**
  - `02_cleaned` (1.25 Hz Z, H1, H2, P);
  - `04_ppsd` (3-h windows, 1 mHz–0.6 Hz, PSDs and Z–P/Z–H cross-spectra, about 800
    stations);
  - `03_tilt` (hourly tilt);
  - `01_preprocessed/*/counts`;
  - the USGS catalogue (M ≥ 5.85);
  - EOT20.

## Objectives

1. **Unusual coherence.** Search every station's coherence in time and frequency for
   behaviour the physics does not predict:
   - (a) Z–P coherence outside the compliance band (above f_c, or below the gauge corner);
   - (b) narrow coherence lines and notches that persist;
   - (c) coherence episodes (collapse or rise) that are synchronous across a network;
   - (d) Z–H coherence episodes (tilt events);
   - (e) Z/P phase excursions in time.

   Separate instrument causes from ocean or solid-earth causes.
2. **Events.** Test whether external events change the instruments or the data:
   - (a) large earthquakes: steps in tilt, \|Z/P\|, coherence or levels after strong
     shaking, by superposed-epoch analysis against random times;
   - (b) volcanic or atmospheric events during deployments, for example the La Palma
     eruption (19 September – 13 December 2021) for 8J, and the Hunga Tonga eruption
     (15 January 2022);
   - (c) storms (pressure-level extremes) and their effect on coherence.
3. **Stress-test PR #23.** For each main claim, find the strongest alternative
   explanation and test it:
   - DPG response;
   - polarity by family;
   - gate bias;
   - epoch mixing;
   - RR50/RR36 seismometer;
   - 8J drift.

   Explain the open anomalies: the 8J drift, the 7B 24-h tilt line, the 7D March 2013
   step, and the poor stage E fits of 7D.J35A and J44A.
4. **Decide.** Answer D1–D10 of PR #23 with a recommendation and its reasons. Implement
   those that do not change published numbers. Results-changing ones go behind options
   or into new run folders.

## Method rules

- **Robust statistics:** median/MAD, empirical null distributions (random dates,
  permuted stations), and multiple-testing control (report how many tests were run).
- **Evidence labels:** label every claim **Verified** (with script and number) or
  **Hypothesis** (with the test that would settle it). Report null results too.
- **Do not change the archive:** new outputs go to `13_audit/research/`; never write in
  `06_compliance*` or `12_inversion_all_stations`.
- **Reuse:** use `compy.anomaly`, `compy.response_check`, `compy.timefreq` and the s01–s07
  outputs; add package code only when it is reusable, with tests.

## Deliverables

1. `docs/audit_2026-11-data/RESEARCH.md`:
   - findings ranked by importance, with figures;
   - null results;
   - stress tests;
   - the decision record (D1–D10 answered, with reasons and what was done).
2. Scripts `docs/audit_2026-11-data/scripts/r01…`; tables in `data/`; figures in `figures/`.
3. Tests passing, ruff clean. Small commits on `compy3-audit-data`, pushed. Add a comment
   on PR #23 that summarizes the findings and decisions.
