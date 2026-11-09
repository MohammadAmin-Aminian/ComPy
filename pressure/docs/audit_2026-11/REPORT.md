# Pressure-gauge audit: QC tests, calibrations, anomalies and roadmap (November 2026)

This audit covers the pressure-gauge studies that moved into ComPy from BRUIT-FM
(`compy_pressure/` and `pressure/`).

- **Base:** branch `pressure-audit`, stacked on `compy3-audit` (PR #21), because it reuses
  `compy.anomaly` and the figure style.
- **Inputs:**
  - **history:** the BRUIT-FM repository at `8595452`: README, `integration/README.md`,
    `docs/REPORT.md` and the dated reports;
  - **code and docs:** ComPy's `compy_pressure` (40 modules, about 20 k lines) and
    `pressure/docs`;
  - **data:** the v1.2.0 master table of 846 pressure-channel epochs from 30 networks.
- **Defaults:** no gain, response or score changes. Every proposed change is listed under
  [Decisions](#9-decisions).

**Evidence labels**, as in the ComPy 3 audit:

- **Verified:** computed here from code or the archived table, with the script named.
- **Hypothesis:** consistent with the evidence, with a test that would settle it.

## 0. Where we started and what was checked

BRUIT-FM is a workflow for ocean-bottom pressure channels:

- it builds an inventory, then PPSDs;
- it runs three QC tests: expected levels (Webb 1998 / Brown et al. 2014), nearby
  stations, and infragravity (Godin et al. 2013);
- it runs three calibrations: teleseismic Rayleigh waves, microseism admittance and local
  P waves;
- it scores, classes and accepts channels, then builds a station-balanced noise model.

The v1.2.0 campaign scored 846 channels:

- 439 good, 188 questionable, 189 bad;
- 215 channels have at least one calibration;
- five networks are flagged for metadata;
- 396 channels are accepted into the noise model.

**Constraints:**

- No ObsPy, PPSD files or waveforms here, so the audit runs on the master table and on
  code reading.
- The new code is ObsPy-free and tested here and in CI.
- Analyses that need PPSDs or waveforms come with a command to run on the IPGP machine.

## 1. Executive summary

1. **The decisions are not in the package (verified).** Calibration scoring, the choice of
   pressure sign, the final class, the metadata-offset flags and station acceptance are
   computed by scripts in BRUIT-FM's dated report folders:
   - `docs/reports/2026-10-28_seismic_calibration/final_table.py`;
   - `.../provisional_decisions/{metadata_offsets,acceptance}.py`, duplicated under
     `gain_corrected/decisions/`.

   They are not in `compy_pressure` and have no tests. One visible consequence: 5 dead
   channels carry a "measured" offset of −3088 dB, which the offset step turns into
   multipliers of 2.4 × 10¹⁵⁴.
2. **"The methods agree" is not independent validation (verified).** Rayleigh and
   microseism agree (median +0.1 dB, 34 of 43 within 3 dB), but they share most of their
   modelling:
   - **the same water-layer theory**, written twice (`rayleigh_calibration.water_ratio`
     and `microseism_check.admittance_height`);
   - **the same homogeneous water column** (1028 kg/m³, 1500 m/s);
   - **one crust for every network:** the Rayleigh run uses RR52's dispersion.

   Local P is the only independent physics, and it agrees less well: the Rayleigh/local-P
   median is +2.9 dB (3 of 6 within 3 dB), and microseism/local-P +1.3 dB (8 of 16).
3. **The water-column assumptions bias deep calibrations (verified, analytic).** For a
   perfect gauge, a 30–70 mHz Rayleigh multiplier is biased as follows:

   | Error in the assumed model | ≤ 3 km | 4 km | 5 km | 5.4 km |
   |---|---|---|---|---|
   | Sound speed 1 % higher | ≤ 0.4 % | 0.9 % | 3.1 % | 4.8 % |
   | Phase velocity ±5 % | ≤ 0.3 % | 0.6–0.7 % | 2.2–2.5 % | 2.3–3.8 % |
   | In-situ density (column mean) | ≥ −0.6 % | −0.9 % | −1.1 % | −1.2 % |

   - **Why it grows with depth:** the least-squares fit gives the bins nearest the
     water-column resonance 27–53 % of the weight at 4–5.4 km.
   - **Mitigation:** a log-ratio fit halves the sound-speed bias.
   - **Impact:** irrelevant for the factor-of-2 QC scoring; relevant for compliance-grade
     (1 %) calibration.
4. **Gains depend on the network, not on depth (verified).** Microseism gains fall with
   depth across the table (Spearman −0.28, p = 2 × 10⁻⁴). Within network × sensor groups
   the trend vanishes (ρ = −0.006, p = 0.94; −0.3 dB/km), so gains are properties of
   instrument batches or metadata. A depth-dependent sensitivity of the gauges, or a
   depth-dependent model bias larger than about 1 dB/km, is ruled out.
5. **Test 1 (expected levels) mostly measures the band above 0.2 Hz (verified).** Webb
   enters as 4 discrete points at or below 0.2 Hz, while every bin above 0.2 Hz is
   compared with Brown's in-water hydrophone curve. At J38A, at most 4 of the 40 bins come
   from Webb. As a result, 18 channels pass Test 1 while sitting 6–32 dB from their
   neighbours; for example, all three 7D.J52A channels are 14–18 dB low.
6. **The infragravity summary overstates the excess (verified).** The report says the
   Pacific and Indian oceans sit +10 to +20 dB above Godin. The table's own Test 3
   residuals give:
   - good Pacific DPGs: +5.8 dB (IQR 2.2–11.1);
   - Indian DPGs: +5.7 dB (n = 5);
   - Indian hydrophones: −1.5 dB;
   - Atlantic: about −1 dB.

   Separately, the 36 YR hydrophones sit **−22 dB** below Godin although their metadata
   declare a usable band down to 4 mHz. This most likely means the response metadata lack
   a high-pass, and it is being scored as data quality.
7. **The local-P sign splits within 7D (verified; cause is a hypothesis).**
   - Rayleigh signs are uniform within each network: YV and 8J −1; 7D, 7K and YO +1.
   - Local-P polarity splits within 7D: 20 positive, 12 negative.
   - The split is the same in event count (median 7) and gain, and occurs in the same
     year and station families (G-stations 2012 have both).
   - This points to mixed wiring across Cascadia operators, which the Rayleigh sign
     choice ("the sign with more events") hides.
8. **New tools:**
   - `compy_pressure.water_layer`: sensitivities, bias and alternative fit;
   - `compy_pressure.table_audit`: 11 detectors, 206 ranked anomalies;
   - `compy_pressure.audit_plots`: 5 figures, plus the ranking;
   - a script, the data copy and 10 tests;
   - the ComPy 3 waveform diagnostics (coherograms including P–H) apply unchanged to
     pressure stations.

## 2. Audit findings (severity-ranked)

| # | Sev. | Where | Finding | Evidence | Status |
|---|---|---|---|---|---|
| P1 | High | BRUIT-FM `docs/reports/2026-10-28_seismic_calibration/final_table.py`; `…/provisional_decisions/metadata_offsets.py:66–114`, `acceptance.py` (and their copies under `gain_corrected/decisions/`) | **Decision logic outside the package:** scoring (e = \|20 log₁₀ g\|, 1.5/6 dB), sign choice, final class, network flags and acceptance; untested and duplicated | Verified (code location); 5 offsets of −3088 dB on nonworking channels (7D.FS01B, 7D.J06D ×2, XE.CC02, YV.RR12) stored as `measured`, multiplier 2.4e154 | Open: D1 |
| P2 | High | `rayleigh_calibration.py:54` (`water_ratio`), `microseism_check.py:70` (`admittance_height`), `local_p_check.py:479` | **Two "independent" calibrations share their physics:** the same water-layer equation, ρ, α, and (Rayleigh) RR52's dispersion for all networks. Their agreement validates the code, not the model | Verified (code); local-P agreement weaker (+2.9 dB, n = 6; +1.3 dB, n = 16) | Open: D2, D3 |
| P3 | Med–High | `rayleigh_calibration.py:68` (guard x < 0.9π/2), `:80–96` (`fit_multiplier`, least squares) | **Depth-dependent model bias:** the sensitivity of P/a to sound speed reaches 6–9 at 5–5.4 km inside the guard; near-resonance bins dominate the fit | Verified (`water_layer`, tests against finite differences) | Option: `water_layer.fit_gain(weights="log")`; D3, D4 |
| P4 | Med | `pressure_test.py:378–394` (4 Webb points vs every Brown bin) | **Test 1 weighting:** the score is dominated by ≥ 0.2 Hz against in-water hydrophone references; the infragravity and notch bands barely weigh | Verified (code); 18 `nearby_hidden` channels | Open: D5 |
| P5 | Med | BRUIT-FM `docs/REPORT.md` §1 | **IG-excess claim (+10–20 dB) not reproduced** by the per-channel Test 3 residuals (Pacific good DPG +5.8 dB) | Verified (`numbers.json`, `infragravity`) | Open: D6 (doc) |
| P6 | Med | Test 3 scoring (`infragravity_test.py:85`) | **Metadata scored as quality:** YR hydrophones are −22 dB with a declared band to 4 mHz (22 bad, 14 questionable) | Verified (table); cause is a hypothesis (§5.B) | Open: test |
| P7 | Med | sign selection (`final_table.py`), local-P polarity | **Polarity:** no cross-method or within-network consistency check. 7D local-P splits 20/12 while its Rayleigh sign is uniformly +1 | Verified | **Detector added**; D7 |
| P8 | Low–Med | `rayleigh_calibration.py:1731–1733` | **Gain uncertainty** is "event scatter, not a confidence interval" and omits the model terms of P3 | Code | Open: roadmap |
| P9 | Low | `microseism_check.py:70`, `rayleigh_calibration.py:54`, now `water_layer.py` | **Three copies of the water-layer formula** (the new one is ObsPy-free and tested against the Rayleigh one in CI) | Code | Open: D8 |
| P10 | Low | Hydrophone gains | Median g = 1.26 (Rayleigh, n = 23) and 1.29 (microseism, n = 44): hydrophone sensitivity overstated by about 2 dB. Not explained by the roll-off (Spearman with band low edge +0.37, p = 0.09; +0.21, p = 0.21) | Verified | Info |
| P11 | Low | APG microseism | 7D APGs give a median microseism g of 0.48 (n = 8). Quartz APGs are absolutely calibrated, so the method (shallow swell in 0.03–0.1 Hz) is the suspect | Verified | Info |
| P12 | Low | `group_noise_model.py:55, 94–128` | Checked: the local-season mapping (DJF winter north, JJA south), station balancing and the 30 dB dead-window rule are correct as described | Verified | — |

The ComPy 3 audit items M1 (homogeneous water column in the calibration) and M3
(pressure polarity) are the compliance-side view of P3 and P7.

## 3. New diagnostic plots

The plots live in `compy_pressure/audit_plots.py`; the examples are in `figures/`.

| Figure | Purpose | Reveals | Example |
|---|---|---|---|
| `plot_gain_consistency` | Same-channel multipliers by two methods, with the ±1.5 and ±6 dB score bands | Method failures; pairs that agree only because they share a model | `gain_consistency.png` |
| `plot_gain_depth` | Gains vs depth across the table and within network × sensor groups | Confounding vs physics (§1, item 4) | `gain_vs_depth.png` |
| `plot_water_layer_sensitivity` | Sensitivity of P/a to sound speed over depth × frequency, and multiplier bias vs depth for each assumption | P3. The saw-tooth beyond 5 km is bins leaving the resonance guard one by one | `water_layer_sensitivity.png` |
| `plot_infragravity_residuals` | Godin residual vs depth by basin and sensor | P5, P6 | `infragravity_residuals.png` |
| `plot_network_scores` | Network × six-check mean scores with counts | Whole-network failures (metadata) vs scattered ones | `network_scores.png` |
| `plot_anomaly_ranking` (ComPy) | Ranked anomalies | Triage | `anomaly_ranking.png` |

**Waveform figures for pressure stations.** These need PPSD or waveform data:

- `examples/compy3_waveform_diagnostics.py` on any P + Z (+ H1, H2) station gives
  spectrograms, P–Z and P–H coherograms and tilt.
- A sustained P–H coherence is the signature of a pressure channel picking up instrument
  motion. This is the test for the "5 Hz floor 20 dB high" seen at 1T.

## 4. Anomaly and correlation detection (`compy_pressure.table_audit`)

**Detectors:**

- gain outliers within network × sensor (robust z ≥ 4 and ≥ 3 dB);
- disagreement between methods (≥ 6 dB);
- polarity, per channel across methods and within networks;
- depth structure, across vs within groups;
- hydrophone roll-off;
- Godin basin excess and outliers within basin × depth class;
- channels hidden from Test 1 (pass Test 1, ≥ 6 dB from neighbours);
- offset outliers within a network (≥ 10 dB);
- invalid offsets (beyond ±120 dB).

Records use `compy.anomaly.Anomaly` and the same score. The ranked table is `anomalies.md`
(and `.csv`), with 206 records at σ ≥ 2.

| Detector | Count | Most affected |
|---|---:|---|
| method_disagreement | 61 | 7D 16, 7K 15, 8J 9, Z3 7 |
| infragravity_outlier | 54 | YR 17, 7D 9 |
| offset_outlier | 51 | 7D 18, YR 7 |
| nearby_hidden | 18 | 7D 12 (J52A, G20B, G11B, J61A on all their channels) |
| gain_outlier | 14 | 7D 7, 7B 2, XE 2 |
| offset_invalid | 5 | 7D 3, XE 1, YV 1 |
| polarity / network structure / basin | 3 | 7D; all; Mediterranean |

The detector finds the known whole-network metadata cases (2A, YD, Z4, X9, 4N) only
through method disagreement and offsets, because they are flagged already. The new items
are listed in §1.

## 5. Physical interpretation and discriminating tests

### A. Gains that fall with depth (§1, item 4)

| Hypothesis | Prediction | Test | Result |
|---|---|---|---|
| Hydrostatic pressure changes DPG or hydrophone sensitivity | Trend within each network with a spread of depths | Within-group slope | **Rejected:** −0.3 dB/km, p = 0.94 (159 channels in groups of at least 4) |
| Water-layer model bias with depth (P3) | ≤ 0.4 dB at 5 km for plausible errors | Analytic | Too small to cause the across-table trend |
| Network confounding (deep networks carry low-gain batches, e.g. YV DPG 0.62–0.77) | Trend across groups only | Same | **Supported** |

### B. YR hydrophones 22 dB below Godin (P6)

| Hypothesis | Prediction | Test | Result |
|---|---|---|---|
| The metadata lack the hydrophone's high-pass | Deficit growing toward low frequency, roughly f² per pole; the microseism band is fine | Frequency-resolved residual from the PPSDs | **Open:** `compy-pressure-infragravity` on YR, then plot residual vs f |
| Godin's model does not apply on the margin | Shallower water has *more* infragravity energy, not less | Physics | Rejected as the cause of a deficit |
| Real low infragravity energy at the site | Neighbouring YR DPGs would also be low | Compare YR DPGs at the same sites | Open: no YR DPG in the scored table |

### C. Mixed local-P polarity in 7D (P7)

| Hypothesis | Prediction | Test | Result |
|---|---|---|---|
| Low-SNR polarity estimates | Minority sign at low event counts | Events per sign | **Rejected:** median 7 events for both signs; gains 0.90 and 0.95 |
| Different wiring between Cascadia operators (LDEO, SIO, WHOI) | Sign follows the operator or sensor model | Join with StationXML equipment fields | Open: G-stations of 2012 have both signs, so it is not year-based |
| Response sign errors per station | Rayleigh, if run per channel without "majority" sign, flips too | Rerun Rayleigh with sign fixed per station | Open |

### D. 7K: three methods, three answers

Rayleigh gives 0.99 (n = 6), microseism 1.67 (n = 10) and the peer offset 0.56 (n = 12)
at about 2.0 km.

- **Rayleigh is the only one independent of peers and of the 0.03–0.1 Hz band.** It
  suggests the gains are right.
- **The peer offset (0.1–2 Hz vs the basin's good gauges)** is then a site effect
  (sheltered or shallow-slope microseisms).
- **The microseism bracket** is probably mis-modelled at 7K's depth, for example the
  mode-0 velocity of RR52's crust.
- **Test:** compute the microseism admittance with a 7K-specific layered model
  (`microseism_check.modal_phase_velocity` with local layers).

### E. Invalid offsets (P1)

Dead or constant data give a PSD near 0, whose log becomes −3088 dB. The fix belongs in
the ported offset step (D1): return NaN and a status `not_measurable` when the channel is
nonworking or has less than about 1 % live windows.

## 6. Implemented changes

| Commit | Change | Numbers changed? |
|---|---|---|
| `Add an ObsPy-free water-layer module…` | `compy_pressure/water_layer.py` | No (new module) |
| `Scan the pressure QC/calibration master table…` | `compy_pressure/table_audit.py`, `audit_plots.py`, `tests/test_pressure_audit.py` | No |
| `Add the pressure audit script, its data and outputs` | `pressure/scripts/pressure_audit_2026_11.py`, `data/` (master table without local paths, offsets; provenance README), figures, anomalies, `numbers.json` | No |
| `Add the pressure audit report` | This report, links, changelog | No |

## 7. Validation

- **No existing module changed.** The ComPy 3 default estimator remains bit-identical
  (PR #21), and every pressure result is untouched.
- **`water_layer.ratio`:**
  - equals the Zha & Webb expression on both branches (rtol 10⁻¹²);
  - the analytic sensitivities match central finite differences (rtol 10⁻⁵);
  - the guard reproduces `rayleigh_calibration`'s rule;
  - `fit_gain("ls")` equals `fit_multiplier` and `ratio` equals `water_ratio` (rtol 10⁻¹²,
    in CI where ObsPy is installed).
- **`table_audit`:** each detector is tested on a synthetic table with planted
  anomalies; the archived table reproduces 54 Rayleigh gains and the 5 invalid offsets.
- **Tests:** 9 passed and 1 skipped here (the ObsPy check); ruff check and format are
  clean.

## 8. Roadmap (pressure part of ComPy 3)

| Rank | Item | Impact | Effort |
|---|---|---|---|
| 1 | Port scoring, sign choice, metadata offsets and acceptance into `compy_pressure.decisions`, with tests, one copy and a NaN guard for dead channels (P1) | 5 | 2 |
| 2 | Independent calibration: an absolute-gauge comparison, or tides for APGs; tidal DPG time constant (Doran et al. 2019); local P promoted where it has ≥ 5 events | 5 | 3 |
| 3 | Station-specific water column (WOA or Mackenzie mean sound speed, in-situ density) and crust (CRUST1.0 or a local layered model) for Rayleigh and microseism; model uncertainty propagated into g | 4 | 2 |
| 4 | Robust multiplier fit (log-ratio, or inverse model-variance weights) and a per-event confidence interval | 3 | 1 |
| 5 | Test 1 by band (infragravity, notch, microseism, > 1 Hz), each scored separately, with Brown used only for hydrophones | 3 | 2 |
| 6 | Response-shape QC: frequency-resolved Godin and Webb residuals to detect missing poles (YR, 1T) before scoring quality | 3 | 2 |
| 7 | Polarity registry per station (sign, method, evidence), checked across methods, with metadata corrections proposed, never applied silently | 3 | 1 |
| 8 | One water-layer implementation (`water_layer`) used by Rayleigh, microseism and local P | 2 | 1 |
| 9 | Noise model per instrument batch and with gain-corrected inputs only after decision D1 | 2 | 2 |

## 9. Decisions

Each line is a yes/no question.

- **D1.** Port `final_table.py`, `metadata_offsets.py` and `acceptance.py` from the
  BRUIT-FM report folders into a tested `compy_pressure.decisions` module, with NaN for
  dead channels? (Y/N) Impact: reproducible decisions; the 5 invalid offsets become "not
  measurable"; no class changes expected.
- **D2.** Stop presenting Rayleigh–microseism agreement as independent confirmation in
  reports, and require local P or an independent reference for a "validated" gain? (Y/N)
  Impact: wording and acceptance rules for calibration.
- **D3.** Use station-specific sound speed, density and crust in the Rayleigh and
  microseism models? (Y/N) Impact: gains change by ≤ 0.4 % above 3 km and by about 1–5 %
  at 4–5.4 km; QC classes are unchanged.
- **D4.** Switch `fit_multiplier` to the log-ratio fit (or model-variance weights) for
  compliance-grade calibrations? (Y/N) Impact: halves the sound-speed bias at 5 km; gains
  change at the percent level for deep stations.
- **D5.** Score Test 1 per band, so that the infragravity and notch bands weigh as much
  as the microseism band? (Y/N) Impact: the Test 1 scores and some classes change.
- **D6.** Correct the "+10 to +20 dB above Godin" sentence to the per-channel values
  (Pacific good DPG +5.8 dB, IQR 2–11)? (Y/N) Impact: documentation only.
- **D7.** Record and check pressure polarity per station across methods, and investigate
  the 7D split against operator metadata? (Y/N) Impact: no gain change; affects
  phase-aware uses.
- **D8.** Make `rayleigh_calibration` and `microseism_check` call `water_layer.ratio`?
  (Y/N) Impact: none numerically (tested equal); one formula to maintain. The code hash
  of `rayleigh_calibration` recorded in provenance changes.
- **D9.** Move YR hydrophones from "bad" to "metadata error" pending the
  frequency-resolved test (§5.B)? (Y/N) Impact: 36 channels re-classed; noise model
  unchanged (not accepted either way).

## References

- Brown, S. et al. (2014). Hydrophone noise references (Figure 8), as digitized in
  `compy_pressure/reference_data`.
- Doran, A. K. et al. (2019). Calibration of differential pressure gauges through in situ
  testing. *Earth Space Sci.* 6. <https://doi.org/10.1029/2019EA000783>
- Godin, O. A. et al. (2013). *Geophys. Res. Lett.*, Equation 6.
  <https://doi.org/10.1002/grl.50418>
- Webb, S. C. (1998). Broadband seismology and noise under the ocean. *Rev. Geophys.* 36.
- Zha, Y. & Webb, S. C. (2016). Crustal shear velocity structure in the southern Lau
  basin constrained by seafloor compliance. *JGR* 121.
- BRUIT-FM v1.2.0 results report (`bruit-fm/docs/REPORT.md`) and the ComPy 3 audit
  (`docs/audit_2026-11/REPORT.md`).
