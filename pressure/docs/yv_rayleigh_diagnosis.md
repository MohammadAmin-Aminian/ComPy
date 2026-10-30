# YV Rayleigh calibration: why stations have no estimate

Date: 19 October 2026. Branch `claude/yv-rayleigh-diagnosis`, based on main
`2fb1fc0`. Requested by [yv_rayleigh_claude_handoff.md](yv_rayleigh_claude_handoff.md).

**Status: diagnostic only.**

**Codex integration review, 21 October 2026:** the saved numerical findings below
were reproduced without refitting gains. Interpretations have been qualified:
H2/H3/H4 are response/polarity hypotheses, not approved instrument corrections;
failure of the current pilot does not prove data are irrecoverable. The original
Claude `evidence.json` and `stations.csv` remain byte-for-byte unchanged. The
hardened helper writes a new review with expanded provenance. See
[the integration review](reviews/yv_rayleigh_diagnosis_integration.md).

* No gain was fitted, applied or stored.
* No waveform was read.
* No response, StationXML, master, score, saved result or production default
  was changed.
* No calibration batch was run.

Every number below comes from saved products, or from in-memory evaluation of
fixed response hypotheses on saved spectra.

Evidence files:

* [yv_rayleigh_diagnosis/stations.csv](yv_rayleigh_diagnosis/stations.csv):
  one row per pressure epoch.
* [yv_rayleigh_diagnosis/evidence.json](yv_rayleigh_diagnosis/evidence.json):
  full metrics plus SHA-256 of 210 inputs.
* Helper: `scripts/diagnose_yv_rayleigh.py`. Tests:
  `tests/test_yv_rayleigh_diagnosis.py`.

## Bottom line

The 43 YV pressure epochs split as follows:

| Workflow | Estimates | Unresolved | Band exclusions |
|---|---:|---:|---:|
| Frozen no-tilt nominal baseline | 16 | 25 | 2 |
| Completed pole-only trial | 31 | 10 | 2 |

Neither set is accepted absolute calibration. Twelve stations have no estimate
in either workflow. Their bottlenecks fall into four groups.

| Group | Stations | Binding limitation | Status |
|---|---|---|---|
| Coherent signal, phase reduced by explicit response/polarity hypotheses | RR17, RR18, RR22, RR47, RR56 | Current phase screen; physical interpretation still conditional | Testable with saved spectra; needs metadata/polarity confirmation |
| Coherent signal, phase unexplained | RR53 | Unknown response or hardware difference | Needs new information |
| Current products do not support a P/Z estimate | RR12, RR31, RR42, RR57 | Flat pressure; DPG artefacts; reported failed verticals | No estimate under the current protocol; recoverability not established |
| Too few events | RR45 | Noisy vertical, 61 of 141 epoch days held locally | Needs the missing 2012 data; may still fail |
| Reference band | RR41 (and RR57) | 0.03–0.07 Hz band reaches the water-column resonance guard | Separate band protocol; RR57 also has no usable vertical |

Five findings are new relative to the existing record.

1. **Vertical StationXML is invalid at RR22, RR42, RR53 and RR56.** Their BHZ
   (and BH1/BH2) analogue stages contain complex poles without conjugate
   partners. A real-valued instrument cannot have such a response.
   * The defect is in the source `INVENTORIES/YV.xml` (SHA-256 `58afc2d6…`),
     not only in the snapshots. The same defect appears at non-local RR23,
     RR24, RR33 and RR39.
   * Stähler et al. (2016), Table A1, lists the same values as conjugate "±"
     pairs.
   * Restoring the published pairs (H2) explains RR56's phase. On the saved
     earthquake spectra the RMS drops from 115° to 9°. On independent
     microseism bins it drops from 102° to 3°.
   * H2 also explains most of RR22's phase: 110° → 21° on earthquakes and
     120° → 11° on microseism bins.
   * Calibration processing divides by this response unchecked.
2. **RR17 and RR18 pressure poles are about ten times larger than the network.**
   They are 0.747 and 0.766 rad/s, against 0.042–0.114 rad/s elsewhere.
   * Under Stähler's Eq. (A6), p₁ = 1/(R·C_total). RR17/RR18 imply
     R·C_total ≈ 1.3 s, against 8.8–23.9 s at every other DEPAS/Geomar
     hydrophone. Those other values correspond to about 47 nF effective
     capacitance with either published resistance (210 or 500 MΩ).
   * A fixed alternative, the pole reflected and divided by ten (H4), removes
     the phase residual: 56°/53° → 4°/4° RMS on all 15 coherent events.
   * H4 flattens the amplitude slope across the band: −0.76/−0.71 → −0.09/−0.03.
   * H4 explains independent microseism bins: about 28° → 3°.
   * H4 makes every normal DEPAS station worse, so it is not a generic
     improvement. The exact true pole is not identified by these data.
3. **RR47 is consistent with an opposite relative P/Z polarity under the assumed
   forward relation.** H3 gives 7° on all 7 coherent
   events (3° on microseism bins, as the earlier microseism check found). The
   data cannot tell which channel is inverted.
4. **Primary-source instrument records explain four stations without estimates**
   (Stähler et al., 2016, Table 3 and §3.1):
   * RR42 and RR57 carried the new 120 s seismometers, which "failed to level,
     yielding no useful seismometer data".
   * RR45's seismometer "recorded, but at high noise levels".
   * RR31's DPG "had high artefacts roughly every 9000 s".
   * The same paper lists RR12's hydrophone as "good", while the local archive
     is flat in 99.8% of its saved PPSD windows. That points to an archive or
     export question, not necessarily hardware.
5. **The examined local day-file names start on 2013-01-01 at every YV station.** The Rayleigh
   catalogue starts on 2012-12-29. October–December 2012 has therefore never
   been examined: 16–57% of each epoch, RR42 aside.

The immediate recorded blockers are input, coherence, phase or reference-band
screens, not uncertainty in the final gain's absolute scale. The gain-fit water
reference does not directly enter the frozen-bin coherence/phase masks. This
does **not** rule out structure, mode mixtures, vertical response or water
properties as physical contributors to observed disagreement. The dispersion
model also sets packet-window placement, and depth/sound speed enter the
RR41/RR57 resonance preflight. Changing preprocessing or windows is not the same
as applying a deterministic response factor to unchanged spectra.

No arithmetic defect was found in the calibration code. My independent response
evaluation reproduces every saved trial reflected phase to 6×10⁻¹⁴ degrees. The
one software gap is that responses are never checked for physical validity
(conjugate symmetry, positive real poles); see the response section below.

## Sources and how they were used

Read-only inputs. Hashes were checked before use; all 210 are recorded in
`evidence.json`, and a second run found the same values.

| Input | Identity |
|---|---|
| Frozen no-tilt workflow `results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json` | SHA-256 in evidence |
| Its 43 explicit `comparison.json` files and 583 saved event spectra (`*_spectra.npz`) | per-station digest in evidence |
| Pole-only trial `results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343/experiment.json` | `446b8e36…` (matches handoff) |
| 43 `YV_RRxx/trial.json`, `summaries.pole_absolute` only | SHA-256 in evidence |
| 82 per-channel response snapshots | each matches the hash recorded in its comparison |
| Source `Desktop/Data_Bruit_FM/INVENTORIES/YV.xml` | `58afc2d6…` (matches comparisons) |
| Six-check review `results/REVIEWS/SIX_CHECKS/final_20261019/review.json` (infragravity, expected-pressure metrics) | `241e9ab1…` (matches handoff) |
| Microseism YV run `bruit-fm-claude/results/EXPERIMENTS/MICROSEISM_YV/run_20261014T175553Z_30002f23` | `summary.json` `68fede25…` (matches [yv_microseism_check.json](yv_microseism_check.json)) |
| Saved RR12 pressure PPSD `results/NETWORKS/TEST/YV/RR12/…BDH…/ppsd.npz` | counted only |
| Raw folder `Desktop/Data_Bruit_FM/Data/YV/*/{BDH,BHZ}.D` | file names only, no content |
| Stähler et al. (2016), publisher PDF, DOI [10.5194/adgeo-41-43-2016](https://doi.org/10.5194/adgeo-41-43-2016) | SHA-256 `7102c298…`; Table 3 (pp. 48–49), §3.1 (pp. 49–50), Appendix A (pp. 56–57) |

Completed studies were read and not repeated:
[yv_heldout_results.md](yv_heldout_results.md),
[yv_heldout_arrival_review.md](yv_heldout_arrival_review.md),
[socrates_all_station_comparison.md](socrates_all_station_comparison.md),
[yv_rayleigh_saved_results_review.md](yv_rayleigh_saved_results_review.md),
[yv_microseism_check.md](yv_microseism_check.md) and
[paired_spectra.md](paired_spectra.md).

## Method

### Event outcomes

Event outcomes are copied from each saved comparison and kept in separate
categories:

| Category | Meaning |
|---|---|
| `selection_excluded_24h` | Within 24 h of an attempted mainshock; a selection rule, not a failure |
| `input_refused` | Missing, gapped, flat or misaligned window |
| `coherence_failed` | Event-level coherence screen failed |
| `phase_review` | Fewer than three bins with \|phase\| ≤ 20° |
| `arrival_overlap_review` | Another predicted arrival overlaps the window |
| `fit_bound_unresolved` | Gain grid minimum on its boundary |
| `candidate` | Event accepted |

Counts are per station. They are never added up into "failed stations".

### Coherent events

An event counts as coherent under the production rule: median and mean coh² of
at least 0.95, with at least three such bins. Coherence does not depend on
response metadata, so every hypothesis below is evaluated on the same coherent
events.

### Response hypotheses

These are five explicit alternatives evaluated without a numerical optimizer.
H4's divisor ten was chosen after inspecting the metadata/data pattern; it is
not an independently specified or validated calibration parameter.

| ID | Pressure response | Vertical response | Polarity |
|---|---|---|---|
| H0 | as recorded | as recorded | as assumed (sign −1) |
| H1 | positive real pole reflected to −p (the completed trial) | as recorded | as assumed |
| H2 | as H1 | unpaired complex roots restored to conjugate pairs | as assumed |
| H3 | as H1 | as H2 | reversed |
| H4 | positive pole moved to −p/10 | as H2 | as assumed |

For each hypothesis the helper computes the phase that nominal-response
processing would show if the hypothesis were true and the physical P/Z relation
real. Responses are evaluated as ratios of analogue-stage poles and zeros, in
the same s = +iω convention as ObsPy evalresp. It then reports:

* RMS residual against the saved phase over coherent bins;
* edge-bin residuals;
* the number of coherent events with at least three bins within 20°. This is
  the phase screen only; arrival overlap and fit bounds are not re-evaluated,
  so these are **not candidate events**;
* the level-free slope ln(r(0.07 Hz)/r(0.03 Hz)) of the observed/expected
  amplitude ratio after the hypothesis's amplitude factor.

The same hypotheses are scored on the saved microseism bins: 0.05–0.30 Hz,
coh² ≥ 0.5, convention D·conj(P), which is the negative of the Rayleigh phase.

### Controls

* **Reproduction.** H1 reproduces every saved `reflected_phase_deg` in the
  completed trial (maximum difference 5.7×10⁻¹⁴°).
* **Identity.** Where the metadata has neither a positive pole nor unpaired
  roots, H1, H2 and H4 equal H0. This holds at all eight INSU estimate
  stations (RR28–RR52: 2–8° under H0, H1, H2 and H4).
* **Negative controls.** H4 worsens every DEPAS station with a normal pole, for
  example RR10 3° → 17° and RR05 15° → 32°. H3 gives at least 160° at every
  station except RR47 (7°) and RR17, RR18 and RR53 (101–127°), none of which
  it explains.

### Order of discovery

The analysis was not blind, so the order matters.

1. H2 came from a metadata scan for physically invalid roots.
2. H4 came from the cross-station pole distribution together with RR17/RR18's
   saved amplitude multipliers. I had seen one RR18 event's phases before
   computing H4.
3. H3 for RR47 was already reported by the earlier microseism check.
4. The microseism bins were inspected in that earlier review but not used to
   build any hypothesis.

So the microseism agreement corroborates the hypotheses but is not a blind
test. All earthquake-spectrum agreement is development evidence, not
validation.

## Station-by-station: the 12 stations without any estimate

Coverage is the number of epoch days that have both BDH and BHZ day files
locally. "Attempted" excludes the 24 h selection skips.

### Coherent, phase reduced by metadata hypotheses (testable with saved data)

**RR17**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 296/387 | 24: 12 coherence, 11 phase, 1 input | 11 | H0 78°, H1 56°, **H4 4°** |

* The pressure pole is 0.747 rad/s (≈8.4 s corner). Stähler gives a nominal
  100 s DEPAS hydrophone corner.
* Under H1 the residual is flat at −58°/−52°. Under H4 it is −3°/−3°, and all
  11 events pass the phase screen.
* The amplitude slope goes from −0.76 to −0.09. The nominal ratio level of 2.53
  divided by H4's factor of 2.51 at 0.05 Hz gives about 1.0. That lies inside
  the range of other DEPAS stations, 0.76–1.01 (excluding RR22, RR53, RR55 and
  RR56).
* Microseism RMS: H1 29°, H4 3° (668 bins).
* Stähler Table 3 lists the seismometer and hydrophone as good.
* **Next action:** obtain the DEPAS R and C_hyd (or the calibration sheet) for
  this hydrophone.

**RR18**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 302/390 | 24: 20 coherence, 4 phase | 4 | H0 83°, H1 53°, **H4 4°** |

* The pressure pole is 0.766 rad/s.
* H4 gives residuals of 3°/0°, and all 4 coherent events pass the phase screen.
* The amplitude slope goes from −0.71 to −0.03. The level is 2.13 divided by
  2.56, about 0.83.
* Microseism RMS: H1 27°, H4 3° (678 bins).
* Coherence is also weak: 20 of 24 events fail, many between 0.85 and 0.95.
* **Next action:** same as RR17. Even with the phase resolved, only 4 events
  are coherent.

**RR22**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 305/390 | 24: 16 coherence, 7 phase, 1 input (0.0001 s grid) | 7 | H0 73°, H1 110°, **H2 21°** |

* The vertical metadata repeats −395.1 − 850.69i instead of Table A1(c)'s
  ± pair. That puts a nearly constant −130° on Z.
* Under H2 the residuals are +27°/+15° and 5 of 7 events pass the phase screen.
  This is the same residual pattern as other DEPAS stations with poles near
  0.1 rad/s under H1 (RR19 16°, RR21 14°, RR26 19°, RR48 21°).
* Microseism RMS: H1 120°, H2 11°.
* The 0.0001 s grid refusal is a single event and irrelevant to the phase
  (0.002° at 0.05 Hz).
* **Next action:** confirm the response at the data centre, then run the
  saved-spectra pilot.

**RR47**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 172/250 | 21: 14 coherence, 7 phase | 7 | H0 140°, H1 174°, **H3 7°** |

* Relative P/Z polarity is reversed against the network convention. All 7
  events pass the phase screen under H3.
* The data cannot say whether the pressure or the vertical channel is
  inverted. The vertical dip is −90° as recorded.
* Its epoch ends 2013-06-22; about 19 M≥6 events were never attempted.
* **Next action:** independent polarity evidence for one channel (P-wave first
  motion on Z, or deployment notes). Then run the saved-spectra pilot.

**RR56**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 180/253 | 21: 16 coherence, 5 phase | 5 | H0 153°, H1 115°, **H2 9°** |

* Geomar OBS. Both vertical pole pairs are recorded with the same sign. That
  adds +117° → +103° of phase and a 1.99–1.39 amplitude factor across the band.
* Under H2 all 5 events pass the phase screen and the slope goes from −0.33 to
  +0.03.
* Microseism RMS: H1 102°, H2 3° (669 bins).
* **Next action:** same as RR22.

### Coherent, phase unexplained (needs new information)

**RR53**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 303/378 | 24: 14 coherence, 10 phase | 10 | H0 69°, H1 29°, H2 80°, H3 101°, H4 63° |

* Geomar OBS with exactly RR56's defective vertical metadata, but repairing it
  makes RR53 worse, not better.
* Ten coherent events are highly consistent: amplitude multipliers 2.06–2.31
  and a phase spread of about 3°. This supports a repeatable measured transfer;
  it does not rule out coherent noise or shared processing/response errors.
* The ratio slope (+0.59) runs opposite to its peers.
* The saved infragravity pilot sits 32 dB above the model, the highest of 39
  YV stations. All others except RR31 and RR57 lie between +6 and +20 dB.
* No microseism-coherent bins exist.
* The actual RR53 seismometer and/or pressure response probably differs from
  RR56's metadata, but I cannot identify how.
* **Next action:** hardware identity and calibration from Geomar or RESIF. No
  data experiment yet.

### Current products do not support a P/Z estimate

**RR12**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 299/392 | 24: 24 input refused | 0 | — |

* Pressure is flat: 4,769 of 4,778 saved PPSD windows sit at the zero-power
  floor, and only 9 scattered windows do not. This agrees with the earlier
  finding of constant −8,268,016 counts.
* Stähler Table 3 lists RR12's hydrophone as good for 390 days.
* **Next action:** compare the local BDH files with the data-centre archive.
  This is new information and needs approval to fetch.

**RR31**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 167/199 | 21: 21 coherence | 0 | — |

* Coherence is at most 0.38 in every event, including M7.5–8.0 events.
* The infragravity pilot is −33 dB from the model (YV range +6…+20 dB), and
  the expected-pressure fraction is 0.05.
* Stähler: "high artefacts roughly every 9000 s… may be difficult to use";
  Table 3 lists the hydrophone as noisy.
* The response metadata is normal (INSU DPG, same as RR28).
* **Next action:** retain the unresolved pilot result. A dedicated artefact
  assessment/removal study would be a separate protocol; the paper notes that
  seismic signals occur between artefacts, so irrecoverability is not established.

**RR42**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — | 1/300 | 22: 22 input refused | 0 | — |

* Only 2013-129 is held locally.
* Stähler: a 120 s seismometer that failed to level, so no useful vertical;
  the hydrophone recorded 298 days.
* The vertical metadata also has RR22's unpaired-pole defect.
* **Next action:** none for Rayleigh P/Z. Fetching more pressure data would not
  supply a vertical.

**RR57**

| Nominal g | Trial g | Coverage | Status |
|---|---|---|---|
| — | — | 304/376 | Conservative band preflight exclusion; pointwise model ratio 0.898 at 5200 m |

The pointwise ratio and the preflight's maximum-phase-speed bound are different
calculations. The pointwise ratio is below 0.9; the conservative preflight still
excludes the requested band. Preserve both rather than treating this as a failed
pointwise threshold.

* Stähler: a 120 s seismometer that failed to level.
* The infragravity pilot is −18 dB (second-lowest) and there are no microseism
  coherent bins.
* **Next action:** keep excluded. A lower band would not supply a usable
  vertical.

### Too few events (needs new data)

**RR45**

| Nominal g | Trial g | Coverage | Attempted | Coherent | Phase RMS |
|---|---|---|---|---|---|
| — | — (2 trial candidates) | 61/141 | 13: 11 coherence, 2 phase | 2 | H0 50°, H1 8° |

* All eligible events were attempted at every magnitude level.
* Both coherent events pass the phase screen under H1.
* Stähler: the seismometer is noisy.
* 80 of 141 epoch days are not held locally; 78 of them are 2012-10-15 to
  2012-12-31.
* **Next action:** the owner decides whether to fetch the 2012 data. The noisy
  vertical may still limit coherence.

### Reference band

**RR41**

| Nominal g | Trial g | Coverage | Status |
|---|---|---|---|
| — | — | 168/246 | Band exclusion (ratio 0.938 at 5430 m) |

* Stähler Table 3 lists both sensors as good.
* **Next action:** lowest priority. A lower-band protocol would be separate
  development, not a fix.

### Notes on the 31 trial estimates

* The 15 trial additions (RR01, RR05, RR06, RR10, RR16, RR19, RR20, RR21,
  RR26, RR30, RR43, RR44, RR46, RR48 and RR55) all failed the nominal phase
  screen: H0 29–58° RMS, H1 3–21°.
* Their microseism bins agree: H0 10–23°, H1 2–25° (14 stations; RR01 has no
  coherent microseism bins).
* Stähler Eq. (A6) defines the hydrophone pole as p₁ = −1/(R·C_total), with a
  minus sign. This is primary-source support for the completed trial's
  reflection. It is not a new fix and does not certify those gains.
* **RR55** (trial estimate 1.592) has an extra metadata inconsistency.
  * Stähler states that RR45 and RR55 carried Geomar HTI-01-PCA hydrophones
    with p₁ = 0.10774 rad/s. RR55's XML pole is 0.043098; RR45 and RR56 carry
    0.102613.
  * Only RR53 matches the stated value.
  * This adds to the handoff's RR55 hardware caution; I did not interpret it
    further.

The full 43-row evidence is in
[stations.csv](yv_rayleigh_diagnosis/stations.csv).

## Classification

Categories follow the handoff.

| Category | Finding | Stations |
|---|---|---|
| (a) Observed screening or coverage limitation | 2012 part of every epoch absent locally (catalogue from 2012-12-29) | all 43 |
| | Bounded selection: 8 attempts per magnitude level (M>7, ≥6.5, ≥6.0), leaving roughly 50–60 M≥6 events unattempted at full-year stations | all with `attempt_limit_reached` |
| | 24 h mainshock rule: 18 skips per full station (selection, not failure) | all |
| | Coherence screen dominates event loss | RR18, RR31, RR45, and partly every station |
| (b) Demonstrated software defect | None in arithmetic: H1 reproduction exact; prior audits stand | — |
| | **Gap:** no physical-validity check on poles and zeros, so invalid vertical metadata and positive hydrophone poles pass silently | affects RR22, RR42, RR53, RR56 and all DEPAS/Geomar pressure |
| (c) Supported response-convention concern | Unpaired vertical poles, demonstrated invalid and contradicting Table A1 | RR22, RR42, RR53, RR56 |
| | Positive hydrophone pole (Eq. A6 negative; completed trial) | all 32 DEPAS/Geomar stations with saved spectra |
| | Pressure pole about 10× too large (H4 supported; exact value unknown) | RR17, RR18 |
| | Reversed relative polarity | RR47 |
| | Hydrophone pole inconsistent with the stated Geomar value | RR55 (and RR45/RR56 by 5%) |
| (d) Unresolved physical or hardware hypothesis | Unknown transfer-function difference | RR53 |
| | Flat archive despite "good" hydrophone | RR12 |
| | DPG artefacts | RR31 |
| | Failed verticals | RR42, RR57 |
| | Noisy vertical | RR45 |

The defective vertical poles also affect BH1/BH2. Any historical tilted-workflow
product at RR22, RR53 or RR56 inherits them.

## What can be addressed, with what

**With existing saved data, after owner approval:** RR17, RR18, RR22, RR47 and
RR56. One fixed response factor per station could be tested on the already-saved
event spectra, without reselection.

* This would turn diagnostic phase agreement into numerical, experimental gains.
* Those gains would still be development results, because the same spectra
  suggested the hypotheses.

**Needs new information first:**

* **Metadata:** the current data-centre StationXML for RR17, RR18, RR22, RR47,
  RR53 and RR56; DEPAS R and C_hyd values for RR17/RR18; Geomar seismometer
  and hydrophone identity for RR53; and polarity evidence for RR47.
* **Data:** the RR12 archive comparison, and the 2012 portion of all epochs.
  RR45 is the station most limited by the missing 2012 portion.

**Not resolved by the current protocol:** RR31 (DPG artefacts), RR42 and RR57
(reported failed verticals) and RR53 (hardware/response information missing).
This is not a proof that all archived records are irrecoverable. RR41 only under a
separate band protocol. Concluding that a station cannot be calibrated with
available information is a valid outcome here.

## Proposed next steps (owner approval needed for each)

Ordered by information gained per cost.

### Step 0: metadata confirmation (no waveform; needs approval for an outside query)

* Fetch the current RESIF/FDSN StationXML for YV RR17, RR18, RR22, RR42, RR47,
  RR53 and RR56, and a data-availability listing for RR12, RR42 and RR45.
* Compare them root-by-root against local `YV.xml` (`58afc2d6…`).
* Separately, ask DEPAS/AWI for the RR17/RR18 hydrophone R and C_hyd and
  Geomar for the RR53/RR56 seismometer poles and zeros.
* If a corrected official response exists, it replaces H2/H4 as the hypothesis.
  Its agreement with the saved phase then becomes a genuine, external check.

### Step 1: saved-spectra pilot (seconds; no raw data; needs approval)

**Factor.** One fixed factor per station, frozen before any amplitude is seen:

| Station | Hypothesis |
|---|---|
| RR22, RR56 | H2 (published Table A1 pairs) |
| RR47 | H3 |
| RR17, RR18 | The Step 0 pole if available; otherwise H4, explicitly labelled post hoc |

**Unchanged.** Saved events, 0.03–0.07 Hz band, coh² 0.95 median and mean,
20° phase screen, at least 3 events, 24 h rule, original (not weighted-log)
objective. No event reselection.

**Controls.**

* RR52 must reproduce g = 0.6787278140712512 exactly: no hypothesis changes it.
* RR16 must reproduce its trial g of 1.078309 under H2 (identity).
* RR10 under H4 is a negative control: phase must worsen, and it is not
  adopted.
* Input hashes are recorded before and after.

**Output.** An isolated experimental folder with gains labelled by hypothesis.
No master, XML or default change.

### Step 2: independent validation (raw data; after the new computer; needs approval)

* Freeze the Step 1 gains and evaluate them on M≥6 mainshocks never attempted
  at these stations, using the existing held-out runner plus the same factor.
* Use RR52 as the control; it had 14 isolated usable events and a +0.664%
  median in the completed study.
* Keep windows and screens unchanged.
* Expect coherence attrition. RR01/RR06 lost 39/37 of 50 events to coherence,
  and RR47/RR56 epochs end in June 2013 (about 19–20 unused events each).
* **Do not spend that unused event pool on development.** It is the only
  disjoint validation set.

### Optional software guard (production change; needs approval)

* Add a preflight that flags or refuses analogue stages with unpaired complex
  roots, or with positive real poles, unless an explicit hypothesis is
  supplied.
* `unpaired_roots` and its synthetic tests in this branch could serve as the
  regression.

**Not recommended now:**

* any RR53 data experiment before hardware information;
* a lower band for RR41;
* window relocation for coherence: this would be development on inspected
  cases, see [yv_heldout_arrival_review.md](yv_heldout_arrival_review.md);
* retrying RR12, RR31, RR42 or RR57.

## Limitations

* **Not blind.** H4's ÷10 is one fixed alternative suggested by the metadata
  pattern. The saved phases constrain RR17/RR18's true pole only roughly, to
  about 0.05–0.1 rad/s; H4 does not identify it.
* **Phase screen only.** The pass counts above omit arrival overlap, cluster
  and fit-bound checks, so they are not candidate counts. No gain was
  computed.
* **Coverage from file names only.** Day-file names do not reveal intra-day
  gaps. No waveform or raw byte was read or hashed.
* **Microseism comparison.** It reuses one week (2013-06-10) from a run in the
  `bruit-fm-claude` worktree; its summary hash matches the tracked record.
* **Primary source.** Stähler Table 3 gives RR45's end of record as 4 Jun 2013
  but its record length as 138 days. The XML epoch end of 2013-03-04 matches
  the 138 days.
* **Absolute calibration.** Absolute vertical calibration, Earth model, depth
  and published-agreement questions remain as in the existing record. This
  review does not address them.

## Reproduce (read-only; cannot start calibration or touch masters)

```bash
M=/Users/mohammadamin/Projects/bruit-fm
PYTHONPATH=src:. /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python \
  -m scripts.diagnose_yv_rayleigh \
  --workflow $M/results/NETWORKS/YV_RAYLEIGH_NO_TILT/_workflow/rayleigh.json \
  --trial-root $M/results/EXPERIMENTS/SOCRATES/all_20261013T195356Z_00cf3343 \
  --expect-trial-sha 446b8e36b94c9952736b9a2c9b2353f0dfb0c95e05e8cc3719b9cea72fcc20cf \
  --microseism-run /Users/mohammadamin/Projects/bruit-fm-claude/results/EXPERIMENTS/MICROSEISM_YV/run_20261014T175553Z_30002f23 \
  --data-root /Users/mohammadamin/Desktop/Data_Bruit_FM/Data/YV \
  --output-dir <new folder>
```

* The command refuses an existing output folder.
* It reads JSON, NPZ, CSV and StationXML snapshots, plus raw directory file
  names.
* It runs in about 1 s.

## Return handoff

* **Branch:** `claude/yv-rayleigh-diagnosis`, worktree
  `~/Projects/bruit-fm-yv-rayleigh`, from `2fb1fc0`. Not merged.
* **Files:** this report; `docs/yv_rayleigh_diagnosis/{evidence.json,stations.csv}`;
  `scripts/diagnose_yv_rayleigh.py`; `tests/test_yv_rayleigh_diagnosis.py`.
* **Tests run:**
  * `pytest tests/test_yv_rayleigh_diagnosis.py`: 7 passed.
  * Ruff check and format on both new files: clean.
  * The full suite was not rerun, because no existing code changed.
* **Preservation:**
  * Master CSV/XLSX and `YV.xml` hashes were the same before and after.
  * All 210 diagnosis inputs had identical hashes across two runs.
  * Codex's checkout was not switched or written; only its `results/` were
    read.
  * The Stähler PDF was downloaded into a session scratch folder, not the
    repository.
* **Approvals needed:**
  * Step 0: outside metadata query and data-centre listing.
  * Step 1: saved-spectra pilot.
  * Step 2: raw held-out validation.
  * The optional preflight guard.
  * Whether to fetch the 2012 data.
* **Open questions for the owner/Socrates:**
  * Is the published Table A1 sufficient authority for H2, or must the data
    centre confirm it?
  * Is a post hoc H4 acceptable even as labelled development if no
    calibration sheet exists?
  * How should RR47's polarity be assigned to a channel?
