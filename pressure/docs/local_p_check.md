# Local-P pressure–vertical check: data, pilot and commands

**Current reviewed implementation: v2.** The development/pilot sections below
are historical. The documented noise-window safeguard is now enforced; see
[the closeout report](local_p_phase_report.md) for corrected LT02 numbers,
remaining scientific limits and safe commands. Legacy results are preserved,
not silently relabelled as v2. Helper outputs must use fresh paths below
`results/LOCAL_P/`; the exporter never updates a production master.

Calibration check 3 (`calibration_check_pressure_vertical_local_p_waves`).
Method, conventions and limits: [local_p_method.md](local_p_method.md).
No gain is applied to any data, no master is written and no score is assigned.

## 1. Available inputs

**Waveforms** (`~/Desktop/Data_Bruit_FM/Data`, read in place). Local P needs
pressure and vertical well above 1 Hz.

| Network | High-rate pairs on disk | Period | Local-P status |
|---|---|---|---|
| YV (RHUM-RUM) | 43 × BDH/BHZ, 62.5 Hz | 2013-01-01 → ≤ 2013-11-25 | eligible (43) |
| XO (AACSE) | LA33, LA34, LT02: HDH APG 125 Hz (separate APG logger) and EDH hydrophone 100 Hz, with HHZ 100 Hz; WD66, WD70: HDH DPG/HHZ 100 Hz | 2018-05-20 → 2018-12-31 | eligible (8 pairs) |
| 7D (Cascadia yr 1) | 15 × BDH/BHZ (BDH 40 or 50 Hz, BHZ 50 Hz) | 2011-10-19 → 2011-12-31 | eligible (15) |
| YL (Lau), XO.WD48 | 1 Hz LDH/LHZ only | – | `not_applicable_sample_rate` (usable band ≤ 0.4 Hz) |

- 7D BXH has archive files but no full response in any available StationXML
  (`no_full_pressure_response`). `INVENTORIES/7D.xml` high-rate epochs have no
  response stages, so 7D uses the full-response file fetched for the 7D
  Rayleigh run (`results/NETWORKS/7D_ALL_NO_TILT/_inputs/responses/…xml`).
- XO and 7D filenames carry location `00` while headers are blank. Files are
  indexed by name but identified by header; only that discrepancy is accepted.
- XO day files exist in both SDS and `waveforms/` layouts for some days; the
  SDS copy is used and the days are recorded.
- **XO.LT02 EDH is non-functional** during the whole 2018 record (std 2–15
  counts, ~20 distinct values per day); LA33/LA34 EDH are live.

**Catalogues.** None existed locally (only the M ≥ 5.5 teleseismic files).
USGS ComCat CSVs, with location/depth/magnitude uncertainties, were fetched
by `scripts/fetch_local_catalog.py` into `results/LOCAL_P/catalogs/` (each has a
sidecar JSON with the exact URL, time and SHA-256):

| File | Box (station extent + 320 km) | Rows |
|---|---|---|
| `usgs_local_XO_20180520_20190101_m2p5.csv` | M ≥ 2.5 | 1,835 |
| `usgs_local_7D_20111019_20120101_m2p0.csv` | M ≥ 2.0 | 33 |
| `usgs_local_YV_20121231_20131127_m2p5.csv` | M ≥ 2.5 | 103 |

No phase picks exist for these OBSs. Onsets are measured automatically (AIC)
and labelled as measurements. Slowness comes from TauP `iasp91`.

## 2. Pilot (XO.LT02 and XO.WD66)

**Selection, fixed before the pilot:** the two XO stations with published
references from Deng et al. (2022). LT02 is an absolute gauge (Paroscientific
APG) in 152 m of water at their Eq. 5 validation site, the steady-state regime.
WD66 is their DPG calibration example at 4941 m, the pre-multiple regime. All
eligible pairs of both stations ran with default settings (80 events maximum,
priority by magnitude and distance), at commit `fadfb06`.

Development runs preceded the pilot (`results/LOCAL_P/dev/`): LT02 EDH and
HDH with 15–25 events, and WD66 with 25 events. They led to the forward-model
revision described in the method note.

| Pair | Outcome | Accepted / processed | Main rejection reasons |
|---|---|---|---|
| LT02 EDH / HHZ | `unresolved` | 0 / 80 | pressure onset not measured 61 (dead channel) |
| LT02 HDH / HHZ (lag search ±150 s) | `conditional_gain_candidate` | 7 / 80 | lag refinement 63 (35 coarse lag > 3 s, 28 low onset SNR) |
| WD66 HDH / HHZ | `diagnostic_only` (pre-multiple) | 10 / 80 | onset not measured 59, onset mismatch 6, S too close 5 |

Pressure multipliers g (`g·p_nominal = p_model`), median over events:

| Pair | Band (Hz) | n | g | MAD | noise bracket | model interval | phase | polarity |
|---|---|---|---|---|---|---|---|---|
| LT02 HDH | 0.5–1 | 5 | 0.972 | 0.041 | 0.941–0.995 | 0.895–1.070 | 33° | +1 (3), indeterminate (2) |
| LT02 HDH | 1–2 | 7 | 1.021 | 0.040 | 0.990–1.042 | 0.910–1.173 | 31° | +1 (5), indeterminate (2) |
| WD66 HDH | 1–2 | 5 | 0.781 | 0.014 | 0.769–0.792 | 0.736–0.844 | 3° | +1 (5) |
| WD66 HDH | 2–4 | 8 | 0.776 | 0.037 | 0.727–0.804 | 0.726–0.832 | −1° | +1 (8) |
| WD66 HDH | 4–8 | 8 | 0.861 | 0.049 | 0.805–0.918 | 0.807–0.925 | −2° | +1 (8) |
| WD66 HDH | 8–16 | 10 | 0.949 | 0.137 | 0.846–1.014 | 0.890–1.020 | −3° | +1 (10) |

Interpretation, within the stated assumptions:

- **LT02 APG.** Both bands agree within their spread (scalar 0.998). An
  absolute quartz gauge is expected near 1, so the chain (responses, water
  model, forward estimator) gives a physically sensible result where a reference
  exists. This does **not** test the vertical absolute response, which is
  shared. With 5–7 events and a ±10–17% model interval (depth and c₀ near the
  resonance), it is a conditional candidate only. The ~30° residual phase is
  timing after onset alignment, not polarity.
- **WD66 DPG (diagnostic only).** Phase ≈ 0° and polarity +1 relative to the
  nominal response, which has negative sensitivity. The multiplier is ≈0.78 at
  1–4 Hz and rises above 4 Hz, where the nominal DPG response rolls off
  (|R|/S = 0.85 at 5 Hz, 0.58 at 10 Hz). This is frequency dependent.
  Context only: Deng et al. (Fig. 4) report WD66 C_p = 0.91, 0.91, 0.90, 0.85,
  0.79 for 0.02–0.07 Hz, and the saved Rayleigh multiplier is 0.909 on the
  1 Hz LDH channel (0.03–0.07 Hz). Different bands, channels and relations: no
  agreement or disagreement is claimed.

**APG clock drift (XO.LT02 HDH).** Lags refined per event (17 events with an
accepted refinement) follow a straight line:
`lag = −0.77 s − 0.4825 s/day × (t − 2018-05-24)` (−5.6 ppm), residual RMS 0.28 s.
At the dates of the M7 teleseisms that match the regions in Deng et al. Fig. 3
(our identification), the line predicts 92.8 s (S. Alaska, 30 Nov), 94.9 s
(New Caledonia, 5 Dec), 97.8 s (South Sandwich, 11 Dec) and 102.4 s
(Komandorski, 20 Dec). Deng et al. measured 93.04, 95.16, 98.04 and 102.69 s.
The local-P timing diagnostic thus independently reproduces their offsets
within 0.3 s. It confirms that the LT02 APG cannot be used without per-event
timing correction. (The LA33/LA34 APGs turned out not to drift; see §2.1.)

Outputs: `results/LOCAL_P/pilot/XO_LT02_WD66/` (campaign summary, per-pair
`summary.json`, `events.csv`, `spectra.npz`, figures, logs) and
`results/LOCAL_P/pilot/XO_LT02_WD66_export/`. All input hashes were unchanged.

### 2.1 APG timing follow-up (after the full campaign `all_v1`)

In `all_v1` the logger-name rule enabled a ±150 s lag search on all three XO
APGs. LA33 and LA34 then lost 80/80 and 76/80 events at lag refinement.

- **LA33/LA34 are synchronized.** Wide windows (−450…+150 s) around five
  M5.0–6.0 events show HDH and HHZ arrivals coincident; 2–10 Hz
  cross-correlation lags at LA33 are −0.02…0.00 s, and 0.5–2 Hz lags at LA34
  are 0.00, −0.04, −0.06 s for the three best events. They do not drift.
- **The real limit is the APG noise floor above ~2 Hz.** For those events the
  2–10 Hz SNR is 200–1500 on HHZ but 3–25 on HDH; at 0.5–2 Hz HDH reaches
  17–87. The pressure onset (AIC, 2–10 Hz, SNR ≥ 3) is therefore rarely
  measurable on these deep APGs.
- In LA34 the 1.2 s pulse recurs at 7.1 s, one water round trip
  (2h/c₀ = 6.15 s) later, as the first sea-surface multiple.
- Timing is now set per channel from recorded evidence in
  [`local_p_timing.json`](local_p_timing.json): LT02 uses a fitted drift model
  (`scripts/fit_local_p_clock_drift.py`: −0.4834 s/day, 16 events,
  RMS 0.09 s, valid 1 Jun–30 Dec 2018); LA33/LA34 use zero lag. Unlisted
  separate-logger channels are flagged `separate_logger_timing_unverified`.

Rerun of the three APG pairs (`results/LOCAL_P/campaigns/apg_timing_v2/`,
commit `ff45286`):

| Pair | Timing | Outcome | Accepted | Bands (g, MAD, n) |
|---|---|---|---|---|
| LT02 HDH | drift model | conditional gain candidate | 9 (was 7) | 0.5–1 Hz 0.972 ± 0.041 (5); 1–2 Hz 1.021 ± 0.044 (9); scalar 0.999 |
| LA33 HDH | none | unresolved | 0 | pressure onset not measured 67, mismatch 3, vertical onset 10 |
| LA34 HDH | none | unresolved | 0 | pressure onset not measured 54, mismatch 5, vertical onset 18 |

The single LA34 event accepted in `all_v1` (polarity −1) is gone; it was an
artefact of the unnecessary lag search, so that polarity is withdrawn. At LT02,
51 events still fail refinement: 32 for onset SNR < 3, and 19 where the two
AIC picks lock onto different arrivals despite a lag predicted to ~0.1 s.

## 3. Commands

Use the `bruitfm` environment from the checkout that holds this code. The
catalogues are local (git-ignored); recreate them in a fresh checkout with the
fetch commands in §3.4.

### 3.1 Plan only (no data read)

```bash
conda run --no-capture-output -n bruitfm python scripts/run_local_p_campaign.py \
  --campaign-dir results/LOCAL_P/campaigns/all_v1 --dry-run
```

Prints every pressure channel with its status: currently 66 eligible,
39 `not_applicable_sample_rate`, 15 `no_full_pressure_response` (7D BXH).

### 3.2 Full available dataset (resumable)

```bash
conda run --no-capture-output -n bruitfm python scripts/run_local_p_campaign.py \
  --campaign-dir results/LOCAL_P/campaigns/all_v1
```

Rerunning the same command resumes: completed pairs with an unchanged
fingerprint are skipped; changed or failed pairs are moved to `attempts/` and
rerun. Restrict with `--networks YV` or `--stations RR52 RR40`. Expected run
time is about 10 s per pair (≈ 15 min in total). Per-pair logs are in `logs/`;
the final table is `campaign_summary.csv` / `.json`.

### 3.3 Export and comparison (read-only)

```bash
conda run --no-capture-output -n bruitfm python scripts/export_local_p_results.py \
  --campaign-dir results/LOCAL_P/campaigns/all_v1 \
  --output-dir results/LOCAL_P/campaigns/all_v1_export
```

### 3.4 One pair, explicit inputs

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.local_p_check \
  --pressure XO.LT02..HDH --vertical XO.LT02..HHZ \
  --inventory ~/Desktop/Data_Bruit_FM/INVENTORIES/XO.xml \
  --waveform-dir ~/Desktop/Data_Bruit_FM/Data/XO/LT02 \
  --catalog results/LOCAL_P/catalogs/usgs_local_XO_20180520_20190101_m2p5.csv \
  --results-dir results/LOCAL_P/single/LT02_HDH --options-json lt02_timing.json
```

where `lt02_timing.json` holds `{"timing_drift": <the "drift" object of
XO.LT02..HDH in docs/local_p_timing.json>}`. The launcher reads that file
automatically; `--timing-lag-search S` remains for exploratory use.

Catalogue recreation (new files only; the script refuses to overwrite):

```bash
M=~/Projects/bruit-fm/results/NETWORKS
conda run --no-capture-output -n bruitfm python scripts/fetch_local_catalog.py \
  --inventory "$M/XO_ALL_GAUGES_NO_TILT/_inputs/responses/XO.LA33,LA34,LT02,WD48,WD66,WD70.xml" \
  --stations LA33 LA34 LT02 WD66 WD70 --start 2018-05-20 --end 2019-01-01 \
  --min-magnitude 2.5 --output results/LOCAL_P/catalogs/usgs_local_XO_20180520_20190101_m2p5.csv
conda run --no-capture-output -n bruitfm python scripts/fetch_local_catalog.py \
  --inventory "$M/7D_ALL_NO_TILT/_inputs/responses/7D.J28A,J29A,J30A,J35A,J36A,J38A,J43A,J44A,J45A,J52A,J53A,J54A,J55A,J63A,J67A.xml" \
  --start 2011-10-19 --end 2012-01-01 --min-magnitude 2.0 \
  --output results/LOCAL_P/catalogs/usgs_local_7D_20111019_20120101_m2p0.csv
conda run --no-capture-output -n bruitfm python scripts/fetch_local_catalog.py \
  --inventory ~/Desktop/Data_Bruit_FM/INVENTORIES/YV.xml \
  --start 2012-12-31 --end 2013-11-27 --min-magnitude 2.5 \
  --output results/LOCAL_P/catalogs/usgs_local_YV_20121231_20131127_m2p5.csv
```

ComCat is revised over time; a refetch can differ. Compare the SHA-256 in the
sidecar with the pilot values (XO `101995d9…`, 7D `b6237a9b…`, YV `a05a5376…`).

## 4. Output fields

Per pair, `summary.json` holds inputs with hashes, git commit, options,
`execution_status` (`completed`/`failed`/`skipped`), scientific `outcome`,
per-band results, `frequency_dependence`, `polarity_majority`,
`event_reason_counts`, the pressure epochs used, `score: null`.
`events.csv` lists every candidate with its status, reason, onsets, lags,
regime, windows and per-band values. `spectra.npz` keeps, per accepted event,
the observed and predicted cross-spectra, noise spectra, coherence, gains,
model-sensitivity ratios, phase and every QC mask.

The export writes `local_p_*` columns keyed by `Seed_id`,
`channel_start_time`, `channel_end_time` and `inventory_file`, with
`master_epoch_match` (`exact`, `channel_not_run`, `no_master_epoch`).
