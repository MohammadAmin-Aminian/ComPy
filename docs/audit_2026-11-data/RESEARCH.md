# Research on the whole-archive results: unusual coherence, events, stress tests and decisions

Second part of the whole-archive audit ([REPORT.md](REPORT.md), PR #23). Brief:
[PROMPT_research.md](PROMPT_research.md). Scripts `scripts/r01–r07`; outputs in
`~/Desktop/ComPy/13_audit/research`; numbers in [`numbers.json`](numbers.json) →
`research`.

**Evidence labels:**

- **Verified:** computed here, with the script named.
- **Hypothesis:** consistent with the evidence, with the test that would settle it.
- **Null results are reported.**

## 1. Summary

1. **The tide is a pressure-only reference, and it settles which sensor changed**
   (r01, verified). The M2 amplitude in the raw pressure counts, fitted month by month,
   follows the pressure channel's low-frequency sensitivity and nothing else.
   - **YV.RR52, September 2013:** M2 is flat across the +7 % \|Z/P\| step (−1.0 ± 1.4 %),
     so the step is not the DPG.
   - **YV.RR50:** M2 falls gradually by about 10 % from May to July 2013, inside the
     "stable epoch" used for stage E, before the vertical steps on 30 July.
   - **YV.RR40:** M2 rises by up to 58 % in March 2013, with nothing visible in \|Z/P\|.
2. **8J hydrophones: the low-frequency response changes during the deployment** (r01,
   verified).
   - **Size:** M2 in the counts rises by a factor of 2–30 at most deep stations between
     September 2021 and April 2022.
   - **Same cause in the compliance band:** month by month, M2 tracks the compliance-band
     phase (Spearman +0.61) and moves against \|Z/P\| (−0.57); n = 591 station-months,
     28 stations.
   - **Not an external event:** the timing differs from station to station (0–240 days
     after deployment), which rules out La Palma or Hunga Tonga.
   - **Hypothesis:** the HTI-01/04-PCA/ULF hydrophone settles after deployment.

   One gain per record cannot correct 8J.
3. **The DPG response is now pinned by the tide** (r04, r05, verified).
   - **The ratio is flat.** The tide ratio R(f) is flat across O1, K1, M2 and S2 (slope
     −0.11, n = 99), as a scaled leak corner predicts.
   - **The tide predicts the phase.** The corner from the M2 tide alone predicts the
     measured compliance-band phase (Spearman 0.82, n = 119, median offset −1.8°).
   - **Response bias:** −3.0 % median at the band bottom (p10 −6.4 %, p90 −0.6 %); −2.0 %
     from the phase.
   - **A bug in PR #23's first version:** its "M2 and K1 differ by 5–12 dB" claim was an
     `s03` bug, now fixed.
4. **Unusual coherence: none instrumental across 678 stations** (r02, null result).
   - **The hole above f_c exists at 595 of 597 stations** (median minimum γ² 0.013). The
     two exceptions are compliance coherence extending just past f_c at about 1.1 km.
   - **Narrow coherence "lines"** are the microseism peak on the 1/8-octave grid, not tones.
     Water-column resonance is not supported (p = 0.085).
   - **All 17 network-synchronous episodes are natural:** storms (9H, January 2008, +21 dB),
     swell (6A, March 2011), surface waves (7B, M7.6 on 19 October 2020) and an array-wide
     bottom-current event (7K, 25–26 May 2019).
5. **Earthquakes do not damage the data or move the instruments for good** (r03, null
   result).
   - **No dose-response:** across 1,574 station–earthquake pairs, the size of tilt or
     \|Z/P\| changes does not grow with shaking (ρ = 0.02 and 0.08).
   - **The cases that stood out were transients:**
     - 7D, 11 April 2012: a 2-day dip after the M8.6 Sumatra sequence, surface waves the
       event masks missed;
     - XO.LA32, 30 November 2018: a gap and a 1-day transient.
6. **The poor stage E fits of 7D.J35A and J44A are partly the band top** (r08, verified).
   Their residuals are largest at 0.8–0.9 f_c.
   - **Capping the band at 0.75 f_c:** J44A goes from χ² 179/11 to 38/9, and J35A from
     246/13 to 167/9.
   - **Adding the response correction:** J35A improves further to 112/9.
   - **Something else is left at J35A.**
   - **Recommendation for the next run:** cap the band at 0.75 f_c, or add a model-error
     term near f_c.
7. **Decisions taken** (§6): D1, D2, D4, D5, D7, D8, D9 and D10 yes; D3 and D6 no. Code for
   D2, D5 and D9 is on this branch with tests. The D2 correction was tested by rerunning
   stage E on the eight YV stations (§5).

## 2. Stress tests of PR #23

| Claim in REPORT.md | Strongest alternative | Test | Outcome |
|---|---|---|---|
| The residual phase is the DPG response | EOT20 is wrong at the seafloor; nonlinear DPG | APGs match EOT20 within 1 dB at 47 of 54 (atlas, now fixed); R flat across four constituents (r04); tide-only corner predicts the phase (r05) | **Holds, stronger** |
| Phase and tide disagree by 3.4 dB; M2 and K1 differ | A second leak path | Recomputed with the response at each constituent | **Retracted:** K1 bug; the offset comes from the phase-only corner |
| RR50 and RR36 jumps are the seismometer | Pressure | Monthly M2 (r01) | RR36: holds. **RR50: refined:** vertical step on 30 July *and* a DPG change from May to July inside the stable epoch |
| YV September 2013 step: sensor unknown | — | Monthly M2 | **RR52: not the DPG** (M2 flat, ±1.4 %). RR28 and RR29: M2 too noisy (SNR 10). RR40: M2 falls, but a DPG excursion is ongoing |
| 8J drift is a sensitivity drift | Seismometers drifting together | Monthly M2 in the hydrophone counts | **The hydrophone**, and frequency-dependent: a response change, not a gain change |
| Polarity by DPG family | Low-SNR local-P | — | Unchanged (95 % agreement) |
| Gate bias +0.24 % | — | — | Unchanged |

## 3. Unusual coherence (r02)

**Inputs:** the per-window cross-spectra of `04_ppsd` (3-h windows, 1 mHz–0.6 Hz) for 678
stations. Per window: Z–P coherence in the infragravity band, above f_c, below the gauge
corner and in the microseism band; Z–H coherence; Z/P phase.

- **Synchronous episodes.** Network episodes were tested against 200 circular shifts per
  network (99.9 % level). This gave 17 bins in 9 networks, each traced:

| Network | When | Signature | Cause |
|---|---|---|---|
| 9H (Gulf of Cadiz) | 4 January 2008 | Infragravity coherence up at 87 % of stations | Storm: infragravity pressure +21 dB |
| 6A (off Namibia) | 13 March 2011 | Coherence above f_c up | Long swell: 30–60 mHz pressure +12 dB (not Tohoku: no fitting aftershock; tsunami periods far lower) |
| 7B (Pacific ORCA) | 19 October 2020 | Z/P phase excursion at 61 % | M7.6 Simeonof earthquake one hour earlier (surface waves in the window) |
| 7K (New England slope) | 25–26 May 2019 | Z–H coherence up and Z–P down at all 6 stations | Bottom current: horizontals +10–20 dB for about 12 h, verticals +7–10 dB, no permanent tilt |
| 2D, 1L, X3, YO, ZZ | single bins | Small (5–6 stations) | Not traced; at the null's edge |

- **Coherence above f_c.** The per-window coherence there is dominated by the estimator
  bias (few averages per window); stacked, it is about 0.01. Three stations with
  coherent energy just above f_c turned out to be the microseism band (from about 50 mHz)
  or energetic infragravity waves at about 1.1 km. No cross-talk.
- **Lines.** 138 narrow peaks at 117–152 mHz and on bin centres are the microseism peak.

## 4. Events (r01, r03)

- **Earthquakes.** Every station was paired with every M ≥ 5.85 event above a shaking
  proxy; tilt and \|Z/P\| changes (48 h and 4 days) were compared with random times at
  the same station.
  - **Excess:** 4.4 % of tilt pairs and 6.5 % of \|Z/P\| pairs are at p < 0.01, against
    1 % nominal.
  - **No dose-response:** |z| does not grow with shaking (ρ = 0.02 and 0.08), so the
    excess comes from non-stationary series and shared dates.
  - **Individual cases were transients:**
    - 7D, 11–13 April 2012: −10 to −17 % for 2 days after the M8.6 and M8.2 Sumatra
      earthquakes and a local M6.0;
    - XO.LA32, 30 November 2018: a 44-h gap starting an hour *before* the M7.1 Anchorage
      earthquake, then a 1-day transient.
  - **Finding:** the stage B earthquake spans may be too short for great earthquakes
    (multi-orbit surface waves). The effect on year-long stacks is below 0.1 %.
- **La Palma (19 September – 13 December 2021) and 8J:** the 8J response changes are not
  synchronous, so not caused by the eruption (§1.2).
- **Hunga Tonga (15 January 2022):** not tested here. The Lamb-wave pressure pulse lies
  below the 8J hydrophone corner.

## 5. Response correction tested on YV (r06)

The compliance stage E inverted was corrected for each station's leak-corner factor:

- from the tide where clear, otherwise from the phase;
- with `compy.response_check.correct_compliance` (compliance, errors, monthly stacks);
- stage E (variant `final`, ComPy 2.8.0, fixed seeds) rerun in
  `13_audit/research/inversion_corrected`.

| Station | k (source) | χ²/N before → after | Sediment [m] | ΔVs 0.5 km | ΔVs 1 km | ΔVs 3 km | ΔVs 5 km |
|---|---|---|---|---|---|---|---|
| YV.RR28 | 0.20 (tide) | 3.7/6 → 2.5/6 | 292 → 312 | +9.1 % (+0.3σ) | +2.2 % (+0.2σ) | -3.0 % (-0.4σ) | -7.8 % (-1.0σ) |
| YV.RR29 | 0.13 (tide) | 3.8/5 → 2.8/5 | 79 → 82 | +3.5 % (+0.2σ) | +1.9 % (+0.1σ) | -0.2 % (-0.0σ) | -2.8 % (-0.4σ) |
| YV.RR34 | 0.56 (phase) | 3.0/6 → 2.7/6 | 215 → 234 | +5.0 % (+0.2σ) | +2.4 % (+0.2σ) | -0.7 % (-0.1σ) | -4.7 % (-0.5σ) |
| YV.RR36 | 0.75 (phase) | 3.3/8 → 3.2/8 | 77 → 80 | +3.5 % (+0.2σ) | +2.9 % (+0.3σ) | +1.0 % (+0.1σ) | -1.2 % (-0.2σ) |
| YV.RR38 | 0.43 (tide) | 10.5/6 → 6.7/6 | 232 → 231 | -2.0 % (-0.2σ) | -1.4 % (-0.3σ) | -1.1 % (-0.4σ) | -0.9 % (-0.5σ) |
| YV.RR40 | 1.62 (tide) | 2.3/5 → 2.1/5 | 184 → 182 | -2.0 % (-0.1σ) | +1.0 % (+0.1σ) | +3.7 % (+0.5σ) | +6.8 % (+1.1σ) |
| YV.RR50 | 0.31 (tide) | 2.5/6 → 3.1/6 | 85 → 87 | +3.5 % (+0.2σ) | +1.4 % (+0.1σ) | -1.1 % (-0.2σ) | -4.1 % (-0.9σ) |
| YV.RR52 | 0.23 (tide) | 3.2/9 → 4.3/9 | 53 → 67 | +9.3 % (+0.7σ) | +4.5 % (+0.5σ) | +0.5 % (+0.1σ) | -3.8 % (-0.6σ) |

σ is half the 68 % width of the published posterior.

- **Where the true leak corner is lower than the metadata (k < 1, six of eight stations):**
  - compliance rises at the band bottom, so deep Vs falls by 1–8 % (up to 1σ);
  - the top 500 m becomes 3–9 % faster.
- **RR40 (k = 1.6) goes the other way:** +6.8 % at 5 km.
- **Fit:** χ² improves at RR28, RR29 and RR38 and worsens slightly at RR50 and RR52.
- **Meaning:** the correction is not negligible next to the current errors, which are
  probably too wide (χ²/N ≈ 0.4). It moves the deep crust most, the depths the paper
  interprets.
- **Status:** these runs stay in the research folder until the full
  corrected rerun of all DPG networks (D2) is made in a new run folder.

## 6. Decision record (PR #23, D1–D10)

The user asked me to decide. **Implemented** = done on this branch or in the run folder.

| # | Decision | Answer | Why | Done |
|---|---|---|---|---|
| D1 | The residual phase is the DPG response | **Yes** | Gauge type decides (APG vs DPG); the tide predicts the phase (ρ 0.82) | Recorded |
| D2 | Correct the DPG response per station | **Yes**, with k from the tide (k from the phase as a check) | Two methods agree; the bias is −3 % median, −6 % at YV | `correct_compliance`, `scale_leak_poles` (tested); correction table `data/dpg_leak_factors.csv`; YV stage E rerun (§5). The full rerun of all DPG networks is the next run, in a new folder |
| D3 | Rerun stable epochs for 21 stations | **No** (flag instead) | Most are 8J (D4). XE.CC07 alternates ±13 % with no pressure change (regimes, not epochs). 7D.G19B fails the phase test (exclude it). 7D.J68A is a vertical step: use the pre-step epoch when it is next rerun | Flags in `data/stations.csv`; roadmap item 2 |
| D4 | 8J provisional | **Yes** | Time-varying hydrophone response (§1.2) | Archive D48 |
| D5 | Stage D sign from the phase | **Yes**, opt-in `PRESSURE_SIGN = "phase"` | 107/108 agree; saves the second stage D pass | `pipeline/d_calibrate.py` + test; default unchanged |
| D6 | Exclusive gate by default | **No** | Median 0.24 %, smaller than the response and epoch systematics; the inclusive gate keeps results comparable with the paper. Add it as a systematic term instead | — |
| D7 | Reword D40 | **Yes** | RR36: vertical. RR50: vertical step plus a DPG change inside the stable epoch | Archive D48 |
| D8 | Fix the atlas tide check | **Yes** | 7 APGs lost to a NaN | Atlas branch commit `94ba4e2`: 47 of 54 within 1 dB |
| D9 | Save per-window levels in stage C | **Yes**, on by default (new file only) | The anomaly scan needs them | `pipeline/c_compliance.py` + test; bit-identical to the s01 levels on RR52 |
| D10 | Merge #21–#23 as 2.9.0 | **Yes** | Defaults unchanged, CI green, 3 stacked PRs | Merged and tagged after this push |

## 7. Still open

- **Rerun with the correction:** all DPG networks with the response correction (D2), in a
  new run folder.
- **8J:** a time-dependent response from monthly tides and phase, or excluding 8J from
  interpretation.
- **Event spans:** lengthen them after great earthquakes (stage B).
- **7B:** the S1 (24 h) line in the 7B tilt angle.
- **J35A:** it still does not fit after the band cap and the correction (χ² 112/9). It had
  \|Z/P\| steps in March–April 2012 and a transient in April 2012; the next check is a
  stack without March–April 2012.
- **Band top for the next run:** at most 0.75 f_c, or a model-error term near f_c
  (results change).
