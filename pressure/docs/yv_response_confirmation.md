# YV response confirmation: RR17, RR18, RR22 and RR56

21 October 2026. Branch `claude/yv-response-confirmation`, from main `bd34fe8`.
Part 1 is a metadata review that fits no gain. Part 2 is the owner-approved
pilot, which computes development gains from saved spectra into an isolated
folder. Neither part opens a waveform, repairs XML or changes a production
default, master, response snapshot or saved result.

## Bottom line

1. **The official record still contains the same values.** Current RESIF
   StationXML, retrieved 2026-10-21 11:35 UTC, is identical to the saved
   response snapshots for all eight P/Z channels. This covers epoch, sample
   rate, every stage gain, A0, pole, zero, decimation value and FIR
   coefficient. The only difference is the letter case of the overall
   sensitivity units (`PA`/`M/S`/`COUNTS` became `Pa`/`m/s`/`counts`).
2. **There is no second archive.** EIDA routing and the EarthScope federator
   both name only EPOS-France/RESIF. The official response therefore agrees with
   local `YV.xml` by descent, not as an independent measurement. No corrected
   official response exists to replace H2 or H4.
3. **The official metadata cannot identify hardware.** Every channel has a
   generic sensor description ("Hydrophone"/"Seismometer") with no
   manufacturer, model or serial number. Instrument-specific records
   (DEPAS/AWI R, C_hyd and S_d; Geomar and Güralp calibration sheets) remain
   the only way to confirm a physical response. Nobody was contacted.
4. **New internal-consistency evidence for RR56.** Its recorded vertical A0
   (608362) agrees with the **conjugate-paired** roots to +0.097%. That is the
   same offset as the correctly paired RR17/RR18 verticals. The duplicated
   roots as stored give −3.41%. So the normalization was computed from a paired
   root set, which supports H2 as a sign lost at export. RR22 points the same
   way, but more weakly: +0.29% paired against +1.51% stored.
5. **RR17/RR18 pressure poles are self-consistent.** Their recorded A0 matches
   the stored 0.747/0.766 rad/s poles to within 3×10⁻⁶. A pole ten times
   smaller would give A0 ≈ 1.00007. So any error lies in the inputs used to
   compute the pole (R, C, or the formula), not in a later edit of one field.
   The metadata cannot distinguish a wrong input from a true 1.3 s time
   constant. **H4 stays a post hoc hypothesis.**
6. **Pilot, development only.** The trial replay and all controls passed,
   including RR52 = 0.6787278140712512 exactly. Development gains on inspected
   events:
   * RR56 H2: 1.166 (4 events).
   * RR18 H4: 1.221 (4 events).
   * RR17 H4: 0.989 (10 events).
   * RR22 H2: 1.419 (5 events; phase only partly explained).
   None is a calibration. Validation on untouched events waits for the new
   computer.

## Shared facts (all four stations)

| Item | Evidence | Status |
|---|---|---|
| Saved inputs | All 8 snapshot hashes equal the values recorded in each station's `comparison.json`. `YV.xml` equals the recorded `58afc2d6…` | Confirmed |
| Snapshot vs source | Each snapshot channel is identical to the same epoch in `YV.xml` | Confirmed |
| Local source | `YV.xml` is a RESIF fdsnws-station 1.1.1 export, `Created` 2024-02-27 | Confirmed (file header) |
| Official epoch | One epoch per channel, with the same start/end as the saved epoch | Confirmed |
| Data centre | EIDA routing: `seisdata.epos-france.fr`; federator: `DATACENTER=EPOSFR` | Confirmed; single source |
| Units, rate | Hydrophone PA→V→COUNTS; vertical M/S→V→COUNTS; 50 Hz; analogue stage LAPLACE rad/s, A0 at 1 Hz | Confirmed, identical |
| Digitizer chain | Hydrophone 10066300 counts/V, vertical 373380 counts/V; nine FIR stages from 512 kHz to 50 Hz (the paper's count for 50 Hz), same coefficients | Confirmed identical; ADC values not independently verified (see below) |
| Unit-case change | BRUIT-FM upper-cases units before every check (for example [ppsd.py](../../compy_pressure/ppsd.py)) | Harmless to unit checks. Fresh XML would change snapshot hashes, which the provenance guards refuse |
| Hardware identity | No model, serial or calibration date in the official or local metadata | **Unknown** |

Primary source: [Stähler et al. (2016)](https://doi.org/10.5194/adgeo-41-43-2016),
local PDF SHA-256 `7102c298…`. Table 3, Appendix A2 and PDF page 14 (Table A1
(a)/(c), and the Eq. A3 paragraph) were checked against a rendered page image,
not only the text layer.

## Station-by-station evidence

"Confirmed" means a directly checked fact about the records. "Hypothesis"
means an interpretation the records do not establish.

### RR17 (DEPAS LOBSTER, 2205 m)

| Item | Saved = official (identical) | Stähler 2016 | Status |
|---|---|---|---|
| Epoch | 2012-10-02T15:50:20 to 2013-10-23T15:41:08 | Deployed 3 Oct 2012 (on-deck start allowed); recovered and end of record 23 Oct 2013; seismometer and hydrophone "good" | Confirmed consistent |
| Pressure stage | 191.0 µV/Pa; pole **+0.747005** rad/s; zero 0; A0 1.00704 | DEPAS: about 185 ± 8 µV/Pa; p₁ = −1/(R·C_total); station R, C not printed | Positive sign contradicts Eq. A6: **confirmed**. Magnitude: **unverified** |
| A0 vs pole | Matches the stored pole (+2.5×10⁻⁶); −p/10 would need A0 1.00007 | — | Confirmed self-consistent |
| Implied R·C_total | 1.339 s, so 6.4 nF at 210 MΩ or 2.7 nF at 500 MΩ | C_hyd 56.3 ± 3.5 nF; R 210 or 500 MΩ; C_amp not printed | Anomalous against the printed parameters; cause **unknown** |
| Overall sensitivity | 1918.1 counts/Pa; stage product 1922.66 (−0.24%); evalresp at 1 Hz +0.17% | Not printed | Confirmed small inconsistency (<0.25%) |
| Vertical roots | Table A1(a) values as stored (max relative difference 7.5×10⁻⁶); pairs valid | A1(a) | Confirmed match |
| Vertical gain | 2012 V/(m/s) × 373380 = 7.51241×10⁸ counts/(m/s) | Typical 1980 ± 15 V/(m/s) | Within about 2σ |
| H4 (pole −p/10) | Not in any record | Not printed | **Post hoc hypothesis** |

### RR18 (DEPAS LOBSTER, 4743 m)

| Item | Saved = official (identical) | Stähler 2016 | Status |
|---|---|---|---|
| Epoch | 2012-10-05T06:07:29 to 2013-10-29T16:25:04 | Deployed 6 Oct 2012; recovered and end of record 29 Oct 2013; both "good" | Confirmed consistent |
| Pressure stage | **214.0** µV/Pa; pole **+0.765683** rad/s; A0 1.0074 | DEPAS about 185 ± 8 µV/Pa | Sign contradicts Eq. A6: confirmed. Sensitivity is about 3.6 sd above the pool mean but is per-instrument: unverified, not a defect |
| A0 vs pole | Self-consistent (−2.1×10⁻⁶) | — | Confirmed |
| Implied R·C_total | 1.306 s, so 6.2 nF or 2.6 nF | As RR17 | Anomalous; cause **unknown** |
| Overall sensitivity | 2152.14 counts/Pa; stage product −0.095% | — | Confirmed small |
| Vertical | A1(a) as stored; 1946 V/(m/s); 7.26597×10⁸ counts/(m/s) | A1(a) | Confirmed match |
| H4 | Not in any record | — | **Post hoc hypothesis** |

### RR22 (DEPAS LOBSTER, 4920 m)

| Item | Saved = official (identical) | Stähler 2016 | Status |
|---|---|---|---|
| Epoch | 2012-10-08T15:49:42 to 2013-11-01T13:56:08 | Deployed 9 Oct 2012; recovered and end of record 1 Nov 2013; both "good" | Confirmed consistent |
| Vertical roots | −0.074016 ± 0.074016i, −471.239, **(−395.097 − 850.693i) twice**, −2199.11; zeros 0, 0 | Table A1(c): −395.1 **±** 850.69i | Unpaired in the **official** record: confirmed. Matches A1(c) after pairing (4.5×10⁻⁶): confirmed |
| A0 vs roots | Recorded 9.09162×10¹¹; stored roots +1.51%; paired +0.29% | — | Closer to paired: **supporting, weaker than RR56** |
| Vertical gain | 1944 V/(m/s); 7.25851×10⁸; evalresp at 1 Hz −1.56% from the stated value | 1980 ± 15 | Gain plausible; the 1 Hz mismatch follows from the unpaired roots |
| Pressure stage | 200.0 µV/Pa; pole +0.099011 (R·C 10.10 s, so 48.1 nF at 210 MΩ); A0 consistent | DEPAS values not printed | Sign contradicts Eq. A6; magnitude typical. Whether 200.0 is individual or nominal: **unknown** |
| H2 (pair restored) | Not in any record | Printed pairs | **Supported hypothesis**; station hardware not confirmed |

### RR56 (Geomar OBS, 4230 m)

| Item | Saved = official (identical) | Stähler 2016 | Status |
|---|---|---|---|
| Epoch | 2012-10-20T18:27:45 to 2013-06-29T20:55:33 | Deployed 21 Oct 2012; recovered 29 Nov 2013; **end of record 29 Jun 2013**; both "good"; Geomar | Confirmed consistent |
| Vertical roots | **(−0.074016 + 0.07347i) twice, (−502.655 + 596.903i) twice**; zeros 0, 0 | A1(a), which the caption applies to all 60 s stations except RR13/RR22 | Four unpaired roots in the **official** record: confirmed. Matches A1(a) after pairing: confirmed. Geomar sensor model not printed separately |
| A0 vs roots | Recorded 608362; stored −3.41%; **paired +0.097%** (the RR17/RR18 offset) | — | **Confirmed: normalization was computed from paired roots** |
| Vertical at 1 Hz | evalresp is +3.47% above the stated 7.32945×10⁸ | — | Consequence of the stored roots |
| Pressure stage | 200.0 µV/Pa; pole +0.102613 (R·C 9.75 s) | Geomar HTI-01-PCA: 199.5 µV/Pa and p₁ = 0.10774 rad/s | XML differs from the paper by +0.25% (S) and −4.8% (pole): **confirmed**; which is right: **unknown** |
| H2 (pairs restored) | Not in any record | Printed pairs | **Supported hypothesis**; strongest of the four stations |

## What the fixed hypotheses change in band

These numbers come from deterministic metadata arithmetic, using the existing
`hypothesis_terms`; no spectrum was read. "Factor" multiplies the saved |P|/|A|
amplitude ratio, so a corrected ratio is the saved ratio divided by the factor.
"Δphase" is the shift relative to H1, the completed trial's reflection.

| Station | Hypothesis | Factor at 0.03/0.05/0.07 Hz | Δphase vs H1 |
|---|---|---|---|
| RR22 | H2 | 1.000 / 0.999 / 0.999 | −130.2° (constant) |
| RR56 | H2 | 1.991 / 1.568 / 1.389 | +116.8 / +106.1 / +103.0° |
| RR17 | H4 | 3.800 / 2.510 / 1.943 | −54.2 / −53.8 / −49.8° |
| RR18 | H4 | 3.876 / 2.559 / 1.978 | −54.0 / −54.0 / −50.3° |

So RR22's defect is phase-only in band. RR56's and H4's corrections also change
the amplitude in a frequency-dependent way. Any pilot must therefore apply the
amplitude part, not only the all-pass phase factor that the completed trial
used.

## Other observations (not acted on)

* **Primary-source inconsistency.** Appendix A prints S_ADC = 3.62×10⁵
  counts/V and an overall LOBSTER sensitivity of "roughly 7.4×10⁵ counts
  (m s⁻¹)⁻¹" (checked on the page image). These disagree with each other:
  1980 × 3.62×10⁵ = 7.17×10⁸. The XML's 373380 counts/V reproduces the
  mantissa 7.4 (1980 × 373380 = 7.39×10⁸). A 3.1% vertical ADC difference would
  pass directly into g, so the ADC value needs a DEPAS record. The XML value is
  not established as wrong.
* The hydrophone ADC gain (10066300 counts/V) is not printed in the paper.
* Hydrophone overall sensitivities differ from their stage-gain products by
  0.10–0.24%. Calibration evaluates the full stage chain
  ([paired_spectra.py](../../compy_pressure/paired_spectra.py) `checked_response`),
  so the stated overall value is not what is applied.

## Saved-spectrum pilot (owner-approved 21 October 2026, run once)

The owner approved the proposal below on 21 October 2026. It ran once, on saved
spectra only, into
`results/EXPERIMENTS/RESPONSE_PILOT/20261021T115554Z_a8cde25f/`. Curated copies
are in [pilot_20261021T115554Z_a8cde25f/](yv_response_confirmation/pilot_20261021T115554Z_a8cde25f/)
(`pilot.json` `46160ac6…`, `stations.csv` `cad470cb…`). **All gains below are
development values from inspected events, not calibrations.** No XML, master,
snapshot, default, score or rating changed, and no gain was applied anywhere.

**Frozen plan.**

| Station | Variant | Basis | Label |
|---|---|---|---|
| RR22 | H2: restore Table A1(c) pair; pressure as H1 | A1(c); A0 closer to paired | Supported metadata hypothesis |
| RR56 | H2: restore Table A1(a) pairs; pressure as H1 | A1(a); A0 matches paired | Supported metadata hypothesis |
| RR17, RR18 | H4: pressure pole −p/10; vertical as recorded | No official or manufacturer value | **Post hoc**, development only |

Each station was run with H0, H1 and its single variant. There was no divisor
search and no other variant.

**Unchanged rules.** The inputs were the completed pole-only trial's saved event
spectra (`experiment.json` `446b8e36…`), with no reselection; the seven stations
have no recovered boundary events. Everything else is as before:

* 0.03–0.07 Hz band.
* coh² ≥ 0.95 for both median and mean, with at least 3 coherent bins.
* 20° phase screen with at least 3 valid bins.
* Arrival-overlap guard.
* At least 3 events, and the station value is the median.
* `absolute` objective.

**Code.** `calibration_experiment.apply_response_hypothesis` is new and
experiment-only. It multiplies the saved cross-spectrum by rp·conj(rz) and the
|P|/|A| amplitude ratio by |rp|/|rz|, where r = R_nominal/R_alternative. The
existing `evaluate` QC and fit then run unchanged. Coherence is invariant under
these deterministic per-bin factors. `scripts/run_response_pilot.py` binds and
hashes every input before and after the run (174 files).

**Replay before hypothesis.** At all seven stations, H0 and H1 reproduced every
per-event status and multiplier of the completed trial's `nominal_absolute` and
`pole_absolute` variants (atol 10⁻¹⁰), and its station gains. The run stops
otherwise. Phase RMS also reproduces the diagnosis:

* RR17: 78/56/4°.
* RR18: 83/53/4°.
* RR22: 73/110/21°.
* RR56: 153/115/9°.

### Results

| Station | Variant | Candidate / phase fail / overlap / coherence fail | Development g | Top bin removed | Event range | Phase RMS, H1 → variant |
|---|---|---|---:|---:|---|---|
| RR17 | H4 (post hoc) | 10 / 0 / 1 / 12 | **0.989** | 0.985 | 0.961–1.006 | 56.3° → 4.0° |
| RR18 | H4 (post hoc) | 4 / 0 / 0 / 20 | **1.221** | 1.200 | 1.154–1.267 | 52.8° → 4.0° |
| RR22 | H2 | 5 / 2 / 0 / 16 | **1.419** | 1.380 | 1.345–1.482 | 110.1° → 21.1° |
| RR56 | H2 | 4 / 0 / 1 / 16 | **1.166** | 1.151 | 1.127–1.169 | 115.3° → 8.6° |

| Control | Check | Result |
|---|---|---|
| RR52, identity under H2 and H4 | g = 0.6787278140712512 exactly; otherwise stop | **Passed exactly**; ran first; RR52 not refit |
| RR16, identity under H2 | trial g 1.078308697242057 | **Passed** (identical) |
| RR10, negative control under H4 | phase RMS must worsen | **Passed**: 2.7° → 17.3°. But all six events still pass the 20° screen and g moves 1.084 → 1.146 |

**Context.** The completed trial's 23 DEPAS/Geomar YV gains (INSU stations from
Table A1(d)/(e) excluded) have median 1.140 and range 1.015–1.592:

* RR56 (1.166) and RR18 (1.221) fall inside that range.
* RR17 (0.989) lies just below the lowest value.
* RR22 (1.419) is exceeded only by RR55 (1.592).

**A0 sensitivity.** The variants keep the recorded A0. If A0 were recomputed from
each variant's roots, g would change by these factors:

* RR56 H2: ×0.9990. This is the same 0.097% vertical offset that every
  correctly paired A1(a) station already carries.
* RR22 H2: ×0.9971.
* RR17 and RR18 H4: ×0.9921 and ×0.9918, about −0.7% beyond that offset.

All of these are smaller than the event scatter.

### Interpretation (development evidence only)

* **RR56 (strongest).** H2 removes the phase problem (115° → 9°; all coherent,
  non-overlapping events pass). The A0 evidence independently supports the
  paired root set, and g sits inside the DEPAS/Geomar range.
* **RR22.** H2 explains most of the phase but not all. RMS is 21° and 2 of 7
  coherent events still fail the screen, with the residual pattern noted in the
  diagnosis. H2 has no in-band amplitude effect here, so the high g (1.42)
  comes from the nominal amplitudes. If all of it were hydrophone sensitivity,
  the hydrophone would read about 141 µV/Pa rather than the recorded 200.0. That
  is far outside the paper's DEPAS 185 ± 8 µV/Pa, so other contributors
  (structure, vertical gain, water model, event windows) are not excluded.
  Unresolved.
* **RR17/RR18.** H4 removes the phase residual (about 55° → 4°), and RR17's
  event gains are tight (0.961–1.006). But H4 was chosen post hoc, and RR10
  shows that a wrong pole of this kind can still pass the phase screen while
  moving g by about 6%. Passing the screen is therefore not proof. Promotion
  needs the DEPAS R/C_hyd records.
* **Status.** No station is calibrated by this pilot. These events were
  inspected during the diagnosis.

## What remains, and when

1. **Independent validation (needs the new computer and raw data; ask first).**
   Freeze these four development gains and hypotheses. Evaluate them on M ≥ 6
   events never attempted at these stations, with RR52 as the control and
   unchanged windows and screens.
   * Do not spend those events on development.
   * RR56's epoch ends in June 2013, which limits its pool.
   * RR18 has only 4 development events.
2. **Instrument records.** The drafts below are not sent; the owner decides
   whether and to whom to send them.
3. **Outside this task.**
   * RR47 (H3) and RR53 remain as in the diagnosis.
   * No master/XML update follows from this pilot.

### Draft request to DEPAS (AWI Bremerhaven) — not sent

> Subject: RHUM-RUM (YV) hydrophone and logger calibration values for RR17, RR18 and RR22
>
> Dear DEPAS team,
>
> We are re-analysing RHUM-RUM (network YV, 2012–2013) pressure data and
> calibrating hydrophones against seismometers. Could you share the
> per-instrument values used to build the YV StationXML for:
>
> 1. **RR17 and RR18 hydrophones:** serial numbers, sensitivity S_d, C_hyd, and
>    the data-logger input resistance R and amplifier C_amp. Their published
>    poles (0.747 and 0.766 rad/s) are about ten times those of other DEPAS
>    stations.
> 2. **RR22:** the seismometer serial (6-pole CMG-OBS40T response) and whether
>    its hydrophone sensitivity of 200.0 µV/Pa is individual or nominal.
> 3. **The GEOLON-MCS vertical ADC sensitivity.** The StationXML uses 373380
>    counts/V; Stähler et al. (2016) print 3.62×10⁵ counts/V.
>
> Thank you very much.

### Draft note to RESIF / EPOS-France data centre — not sent

> Subject: YV (RHUM-RUM) StationXML: unpaired complex poles on DEPAS/Geomar seismometer channels
>
> Dear RESIF data centre,
>
> In fdsnws-station output for network YV (retrieved 2026-10-21), the analogue
> stage of YV.RR22 BHZ lists −395.097−850.693i twice. YV.RR56 BHZ lists
> −0.074016+0.07347i and −502.655+596.903i twice each. A real-valued response
> needs conjugate pairs, and Stähler et al. (2016, Table A1) prints them as ±
> pairs; RR56's recorded A0 matches the paired roots, not the stored ones. Our 2024
> export shows the same pattern at RR23, RR24, RR33, RR39, RR42 and RR53, which
> we did not re-query. Could the station metadata be reviewed?
>
> Thank you.

## Provenance

Official downloads are in
[official_20261021T113504Z/](yv_response_confirmation/official_20261021T113504Z/)
(`manifest.json` SHA-256 `9e395b09…`), separate from `YV.xml`:

| Query | Request UTC | HTTP | Bytes | SHA-256 |
|---|---|---|---:|---|
| EIDA routing (ORFEUS) | 2026-10-21T11:35:04 | 200 | 1178 | `1d9a5e9d…` |
| EarthScope fedcatalog | 2026-10-21T11:35:04 | 200 | 1234 | `6534ae7c…` |
| RESIF station, RR17 BDH,BHZ, level=response | 2026-10-21T11:35:05 | 200 | 48550 | `000fe5ea…` |
| RESIF RR18 | 2026-10-21T11:35:05 | 200 | 48550 | `ae47590a…` |
| RESIF RR22 | 2026-10-21T11:35:05 | 200 | 48831 | `8475a409…` |
| RESIF RR56 | 2026-10-21T11:35:05 | 200 | 48458 | `a9c5372e…` |

The exact URLs and full hashes are in the manifest. The query window
2012-01-01 to 2014-01-01 is the one recorded in the local export's ModuleURI.
There were no failures. The fetch ran with an earlier revision of the helper
(code SHA-256 `08751cba…`). Later edits added the band terms, LF-only CSV
line endings and stopped a 204 "no data" reply from being saved as an empty
file. Neither change affects these
HTTP 200 retrievals, so they were not repeated.

Comparison output:
[comparison_20261021/](yv_response_confirmation/comparison_20261021/)
(`evidence.json` `600635e2…`, `channels.csv` `7ab226df…`), written with helper
code `08e0ee59…`. It binds the frozen no-tilt workflow (`cbe5af4a…`), four
`comparison.json` files, eight snapshots, `YV.xml` (`58afc2d6…`) and the six
official files. All 20 consumed inputs were re-hashed before writing. `YV.xml`,
the snapshots and the owner's `gain_histogram.png` are unchanged.

## Reproduce

```bash
cd /Users/mohammadamin/Projects/bruit-fm-response-confirmation
PY=/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python
# Only if a fresh retrieval is wanted; always into a NEW folder:
PYTHONPATH=src:. $PY scripts/confirm_yv_responses.py fetch \
  --out docs/yv_response_confirmation/official_<UTC> --stations RR17 RR18 RR22 RR56
PYTHONPATH=src:. $PY scripts/confirm_yv_responses.py compare \
  --workflow /Users/mohammadamin/Projects/bruit-fm/results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json \
  --workflow-sha256 cbe5af4ae7782903d96ea5f6eeff3c07ca3c57e85b44f4f3d55c68f9f85486d2 \
  --official docs/yv_response_confirmation/official_20261021T113504Z \
  --out <new folder> --stations RR17 RR18 RR22 RR56
```

Pilot (saved spectra only; seconds; always into a NEW folder):

```bash
M=/Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. $PY scripts/run_response_pilot.py \
  --workflow $M/results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json \
  --workflow-sha256 cbe5af4ae7782903d96ea5f6eeff3c07ca3c57e85b44f4f3d55c68f9f85486d2 \
  --trial-root $M/results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343 \
  --trial-sha256 446b8e36b94c9952736b9a2c9b2353f0dfb0c95e05e8cc3719b9cea72fcc20cf \
  --out $M/results/EXPERIMENTS/RESPONSE_PILOT/<UTC>_<id>
```

All three commands refuse an existing output folder. `compare` and the pilot
also refuse changed inputs. The pilot stops if RR52 does not reproduce exactly
or if H0/H1 do not replay the trial.

An earlier run of the same pilot (`20261021T115425Z_6fbc1046`) was discarded
and deleted. It had identical gains and controls, but its A0-sensitivity field
assumed the recorded A0 matched the stored roots, which is false for the
RR22/RR56 verticals. The field was corrected and covered by a test before the
reported run.

## Verification performed

* `tests/test_confirm_yv_responses.py` has 12 tests. They cover root
  matching, identifying the root set A0 was computed from, field-level
  comparison, the Table A1 check, implied R·C, the query scope, and fetch
  provenance (200, 204, HTTP and network failures). They also cover exclusive
  outputs, refusal of a changed download, end-to-end hash binding, refusal of a
  changed snapshot or workflow, and detection of a corrected official response.
* `tests/test_calibration_experiment.py` gained 2 tests for the transform:
  exact identity with no mutation of inputs, phase/amplitude/coherence
  behaviour and equivalence to the trial's pole-reflection path.
* `tests/test_response_pilot.py` has 4 tests:
  * trial replay, then recovery of a known gain under H2, plus the A0 field;
  * a stop on a non-replaying trial outcome;
  * refusal of changed snapshots and trial manifests;
  * identity and negative-control logic, including RR52 exactness.
* Final full suite: **465 passed**, 126 warnings. Repo-wide Ruff check and
  format passed, and `git diff --check` is clean.
