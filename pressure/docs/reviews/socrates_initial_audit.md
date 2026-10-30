# Socrates initial scientific audit of BRUIT-FM

| Item | Value |
|---|---|
| Auditor | Socrates (independent scientific review). Codex develops; Sentinel reviews software. |
| Date | 2026-10-12 |
| Revision | 2, same day. Adds four primary sources supplied by the user: Zha & Webb (2016), Deng et al. (2022), Ruan et al. (2014) and Aminian et al. (2025). Changes are marked "Rev. 2"; §9 summarises them. |
| Branch / baseline | `review/sentinel` at `007348e` (Sentinel handoff); Sentinel's base `faa39ff`; `main` at `4387fbe` |
| Checkout | separate clone `~/Projects/bruit-fm-socrates`; Codex's checkout and the shared env were not modified |
| Software checks | `pytest`: 141 passed (HOME redirected); `ruff check .` and `ruff format --check .` passed |
| Reproduction | `PYTHONPATH=src python docs/reviews/socrates_validation.py` (read-only; prints V1–V10) |
| Data access | No raw waveform was read. Read-only: committed docs/JSON/CSV evidence, the `Profondeur` column of the local master CSV, four GMRT point-depth queries |
| Code changes | None. Corrections below are proposals for Codex and decisions for the user |

Labels used below:

- **Confirmed**: demonstrated by derivation, numerical test or primary source.
- **Concern**: plausible, with evidence, but not proven to change a result.

Severity is the impact on the project's scientific deliverables if left unaddressed:
High, Medium or Low.

---

## 1. Summary

### Overall scientific readiness

**The spectral and physical core is correct. The calibration products are not yet
ready to be applied or scored.** I verified the following independently:

- the water-column reference;
- spectral normalisation;
- FFT/evalresp phase conventions;
- the pressure PPSD path;
- the infragravity dispersion solver;
- the least-squares gain.

The provisional gains are reproducible. They are not yet calibrations, for four reasons:

1. **YV metadata error.** A metadata error, not a software bug, decides most YV
   outcomes (SOC-001).
2. **Documented artefacts not handled.** A documented instrument artefact affects
   all nine YV DPG stations and is not handled (SOC-002).
3. **Scalar-gain model fails at some stations.** A single scalar gain demonstrably
   fails at least one station (B01W).
4. **Uncertainty underestimated.** The real uncertainty is several times the event
   scatter now reported (§6).

For the WP2 noise-model deliverable, the PPSD engine is sound. The statistics built
on it still need:

- instrument usable-band censoring;
- glitch handling;
- exclusion of dead-channel floor values.

### Five highest-priority scientific issues

| Rank | ID | Issue | Severity |
|---|---|---|---|
| 1 | SOC-001 | DEPAS/Geomar hydrophone metadata place the RC high-pass pole in the right half-plane, contrary to the published hardware transfer function. **None of the 23 YV stations with pole ≥ 0.06 rad/s produced a candidate**: 22 unresolved and 1 band-excluded. The cause is a predicted phase error of ≥ 23.6° against a ±20° mask, not instrument failure. | High |
| 2 | SOC-002 | The hourly (3620 s) INSU levelling glitch, 15 dB above vertical noise at 20–500 s (Stähler et al., 2016), is not handled anywhere. It contaminates INSU vertical PPSDs, tilt training and Rayleigh windows. | High |
| 3 | SENT-004 / SOC-006 | The scalar-gain model can fail. B01W is quantitatively a single-pole high-pass: fc ≈ 0.08 Hz (τ ≈ 2 s), in-band gain ≈ 0.28. Averaged over Zha's band, this model gives 0.17–0.19, reproducing Zha's 0.19 (Rev. 2). Zha & Webb themselves suspected "a leaky valve in the DPG", and Deng et al. (2022) find frequency-dependent DPG coefficients at many stations. A two-parameter (gain + corner) model with an amplitude–phase consistency test is needed. | High |
| 4 | §6 + SOC-004/005 | Event scatter (often < 2%) is reported as the uncertainty. The combined systematic budget is about ±3–6% for a well-behaved deep station and larger near the band top. The unweighted fit also puts most weight on the most model-sensitive bins. | Medium |
| 5 | SOC-007 | PPSD-based noise statistics (the WP2 low/high models) have no instrument self-noise or usable-band censoring, down to 0.001 Hz. DEPAS long-period self-noise, C07W's reduced-gain days and RR12 floors would enter the low-noise model. | Medium (High for WP2 final products) |

### What is verified as correct

| Item | How verified |
|---|---|
| P/a_up = ρ·tan(rH)/r, r = ω√(α⁻² − c⁻²) (`rayleigh_calibration.py:37`) | **V1**: independent first-order Lagrangian acoustic-gravity ODE; the g→0 limit matches to 1e-13. Identical to Zha & Webb (2016) Eq. A5 and Deng et al. (2022) Eq. 5, read in the primary texts (Rev. 2). |
| Gain direction: g multiplies pressure; Zha's coefficient = 1/g | Zha & Webb (2016) Fig. A2: the published coefficient multiplies the theoretical ratio (Rev. 2) |
| One-sided density scale 2/(f_s Σw²), cross-power averaged before ratios | Welch (1967); Bendat & Piersol (2010); Sentinel's Parseval tests |
| NumPy e^{−iωt} FFT ↔ evalresp s = 2πif; dip > 0 → negate | **V7** (causal high-pass phase matches to 1e-3°); Sentinel's analytic-derivative test (reviewed, sound) |
| ObsPy PPSD pressure path: `output="VEL"` on a Pa-input response returns counts/Pa unchanged; no ω² | **V7**; ObsPy 1.5.1 `spectral_estimation.py` source read |
| IG dispersion x·tanh x = ω²H/g (`infragravity_test.py:41`) | Matches `brentq` to 5 decimals at 2–5 km, 0.002–0.02 Hz |
| Depth inputs (StationXML elevation) | Equal to Zha Table 1 at all 25 YL stations; agree with GMRT within 0.8% at RR03/RR11/RR52/RR54 (SOC-010) |

---

## 2. Literature

### Read for this audit (full text or source code)

| Source | Used for |
|---|---|
| Stähler, S. C., et al. (2016). Performance report of the RHUM-RUM ocean bottom seismometer network around La Réunion, western Indian Ocean. *Adv. Geosci.* 41, 43–63. doi:10.5194/adgeo-41-43-2016. Full PDF: §2–3, Table 3, Appendices A1–A3 | Instrument pools; hydrophone and DPG transfer functions (Eqs. A5–A6); INSU glitch; DPG calibration practice; depths; clock skews |
| ObsPy 1.5.1 `obspy/signal/spectral_estimation.py` | PPSD internals: nfft, taper, dB averaging, hydrophone handling, `dtiny` floor |
| TiSKitPy 2.3.1 package source (`Pressure_noise_model`; module list including `PeriodicTransient`) | Independent Brown et al. (2014) curves; available glitch-removal tool |
| All BRUIT-FM modules, method documents and Sentinel's audit, handoff and evidence on this branch | Code review |
| Zha, Y., & Webb, S. C. (2016), Crustal shear velocity structure in the Southern Lau Basin constrained by seafloor compliance, *JGR Solid Earth* 121, 3220–3237, doi:10.1002/2015JB012688. Full text, including Table 1 and Appendix A (Rev. 2) | Eq. A5 reference; coefficient direction; band 0.04–0.10 Hz; M ≥ 7 events with P/a coherence > 0.9; Zha et al. (2014a) Lau velocity model; nominal DPG leak time constant 600–1000 s; YL sensors Güralp CMG3T + DPG; B01W "leaky valve" |
| Deng, H., An, C., Cai, C., & Ren, H. (2022), *GRL* 49, e2021GL096952, doi:10.1029/2021GL096952. Full text (Rev. 2) | Eq. 5 is identical to Eq. A5 (q = 1/c), validated with an APG up to c₀/4h; DPG coefficients vary with frequency; body-wave calibration with q < 10 s/deg is insensitive to the Earth model; ρc₀v holds only for the first crest |
| Ruan, Y., Forsyth, D. W., & Bell, S. W. (2014), *JGR Solid Earth* 119, 6357–6371, doi:10.1002/2014JB011162. Full text (Rev. 2) | H1-type D/P estimator; coherence > 0.95; DPG scale factor constrained below 0.1 Hz; microseism-band (0.1–0.2 Hz) D/P strongly sediment-dependent; ρ = 1030 kg/m³, α = 1.5 km/s; nonzero transfer-function phase flags response problems |
| Aminian, M. A., et al. (2025), Shallow crustal structures of the Indian ocean derived from compliance function analysis, *GJI* 242, 1–16, doi:10.1093/gji/ggaf253. Full text (Rev. 2) | Published RR28–RR52 DPG gains (Table 1); hourly glitch removal (Deen et al., 2017); tilt and event exclusion; ρw ≈ 1040 kg/m³ in the gravity term |

### Abstracts or metadata only (full text not read)

| Source | Status |
|---|---|
| Godin, O. A., et al. (2013), *GRL* 40, 2159–2165, doi:10.1002/grl.50418 | Full text blocked (HTTP 403 at Wiley, Ifremer and the WHOI archive). **Eq. 6 and its 2000 Pa²/Hz normalisation remain unverified against the primary text.** |
| Doran, A. K., et al. (2019), *Earth Space Sci.* 6, 2663–2670, doi:10.1029/2019EA000783 | Abstract: in-situ DPG sensitivity **and time constant** can differ from nominal by a factor of two or more |
| Deen, M., et al. (2017), First observation of the Earth's permanent free oscillations on ocean bottom seismometers, *GRL* 44(21), 10988–10996 | Abstract: removal of RHUM-RUM INSU glitches by subtracting averaged glitch signals |
| Janiszewski, H. A., et al. (2023), Broad-band ocean bottom seismometer noise properties, *GJI* 233, 297–315 | Abstract: OBS noise grouped by sensor type and depth; tilt depends on sensor, compliance on depth |

### Relied on from standard knowledge (not re-read here)

- Welch (1967), *IEEE Trans. Audio Electroacoust.* 15, 70–73.
- Bendat & Piersol (2010), *Random Data*, 4th ed. (coherence, transfer-function
  estimators and their errors).
- Carter, Knapp & Nuttall (1973), *IEEE Trans. Audio Electroacoust.* 21, 337–344.
- Harris (1978), *Proc. IEEE* 66, 51–83.
- Peterson (1993), USGS OFR 93-322.
- McNamara & Buland (2004), *BSSA* 94, 1517–1527.
- Cox, Deaton & Webb (1984), *J. Atmos. Oceanic Technol.* 1, 237–246.
- Crawford, Webb & Hildebrand (1991), *JGR* 96, 16151–16160.
- Webb (1998), *Rev. Geophys.* 36, 105–142.
- Crawford & Webb (2000), *BSSA* 90, 952–963.
- Ringler & Hutt (2010), *SRL* 81, 972–983.
- Ryan et al. (2009), GMRT, *G-cubed* 10, Q03014.
- Lamb, *Hydrodynamics*: the linearised acoustic-gravity equations behind V1.

### Not read in this session

- Mohammad-Amin Aminian's PhD thesis (not available to me here).
- Zha's thesis.
- Webb, Zhang & Crawford (1991) full text.
- Brown et al. (2014) full text.
- Rebeyrol et al. (2024), cited in Aminian et al. (2025).

I make no claim about their contents. The thesis's calibration chapter (Figure 3.13)
should be reviewed against SOC-001, SOC-002 and SOC-006 in the next pass.

### Missing from `docs/references.bib`

- **Stähler et al. (2016).** This is the most important omission: it documents the
  YV hardware, responses and artefacts that SOC-001 and SOC-002 depend on.
- Cox et al. (1984).
- Webb (1998).
- Crawford & Webb (2000).
- Crawford et al. (1991).
- Deen et al. (2017).
- Janiszewski et al. (2023).
- Welch (1967).
- Bendat & Piersol (2010).
- Ringler & Hutt (2010).
- Brown et al. (2014) and Webb et al. (1991): currently cited only inside
  `pressure_references.json`.
- Ryan et al. (2009), for any bathymetry check.
- **Deng et al. (2022)** (Rev. 2): frequency-dependent DPG calibration and
  Earth-model-independent body-wave calibration.
- **Rebeyrol, S., Ker, S., Duval, L., & Crawford, W. C. (2024)**, "Revisiting the OBS
  seafloor compliance signal removal with a stationarity and stacking-based approach:
  the BRUIT-FM toolbox", *GJI* 239(1), 386–401 (Rev. 2). It carries the project's own
  name, as listed in Aminian et al. (2025); I have not read it.
- Harmon et al. (2022), *Seismica* 1(1): the hourly tilt-rotation method used in
  Aminian et al. (2025).

---

## 3. Sentinel handoff: verdicts

The three status classes are: corrections already made by Sentinel, Sentinel claims
that needed scientific validation, and new findings (§4).

### Corrections already made by Sentinel

| ID | Change | Socrates verdict |
|---|---|---|
| SENT-002 (`2d4e8d6`) | Spectra and true QC status kept when the gain fit hits the grid boundary | **Scientifically valid.** The phase guard is decided before an amplitude-only fit, which is the correct precedence. No acceptance rule changes. |
| SENT-003 (`41a86f5`) | Tilt horizontals take the paired vertical's band code | **Valid.** CleanRotator rotates Z only in the sensor basis, so N/E azimuths are not needed. |
| SENT-011 (`3f663c5`) | Independent convention tests | **Valid.** Bin-centred tones on alternate bins isolate sign, scale and phase. By construction they cannot detect SENT-010; Sentinel states this. |

### Sentinel claims that needed scientific validation

**SENT-001 (C07W): supported and strengthened.**

- The state ratio is 1.909 against C08W (95% CI 1.886–1.960) and 1.970 against C09W
  (1.935–1.987). The event-gain ratio at 0.03–0.07 Hz is 1.963 (V4).
- The same factor of about 2 in the ambient 0.12–0.3 Hz band and in coherent Rayleigh
  waves at 0.03–0.07 Hz means a frequency-independent change. That points to
  electronic gain, not mechanics (mass centring or feedback faults change response
  shape).
- Stähler et al. (2016, App. A1) record that single-ended wiring halves Trillium 240
  sensitivity. Intermittent loss of one leg of a differential output is therefore a
  documented-type mechanism for an exact factor of 2. **This is a hypothesis**; the
  deployment logs are needed.
- Rev. 2: Zha & Webb (2016, §2.1) state that the YL compliance stations carried Güralp
  CMG3T sensors. Their balanced outputs make the same halving mechanism possible.
- Rev. 2: Zha's C07W coefficient (1.99 ± 0.01) rests on only 11 events, against 16–19 at
  most stations. Its value matches our normal-state events (1/g ≈ 1.90, the network
  offset), not the reduced state (≈ 3.7). Zha's accepted events therefore appear to come
  from the normal state; the paper does not say why.
- Additional impact: C07W LHZ PSDs for the 251 reduced days are about 6 dB low and
  must not enter noise statistics.

**SENT-004 (B01W): confirmed, with a physical model.** A first-order high-pass,
H = G·(if/fc)/(1 + if/fc), fits amplitude and phase **jointly** in all four events
(V3):

| Parameter | Value |
|---|---|
| fc | 0.073–0.091 Hz (τ = 1.8–2.2 s) |
| G | 0.260–0.299 |
| Amplitude rms | 3–9% |
| Phase rms | 5–8° |

- Amplitude and phase are linked by the minimum-phase relation |H|/G = cos φ, so this
  is two independent observables agreeing.
- The same fit on B04W and A02W, for the same events, finds no in-band corner. It
  recovers their station values: 1.35 against 1.354, and 1.43 against 1.436.
- The leak time constant is far shorter than nominal, which is consistent with Doran
  et al. (2019).
- Zha's 0.19 is a band average of this shape, not a calibration.
- Rev. 2, primary-source corroboration:
  - Zha & Webb (2016, §2.3) call B01W "anomalously low (~0.2) … which may have a leaky
    valve in the DPG".
  - They give the nominal DPG capillary-leak time constant as 600–1000 s (Cox et al.,
    1984). The fitted τ ≈ 2 s is about 300–500× shorter, which is what a leaking valve
    would produce.
  - Averaged over Zha's 0.04–0.10 Hz band, the fitted model gives 0.17–0.19 (V3),
    reproducing the published 0.19 ± 0.02. A physical model now explains the published
    number.
  - Ruan et al. (2014, §2) likewise treat a significantly nonzero P–Z transfer-function
    phase as a sign of instrument-response problems.

**SENT-004 (RR55): not established as the same mechanism.**

- Sentinel reports R/T 0.24–0.27 with a +52° to +76° phase lead. If R/T is flat, a
  single-pole high-pass cannot produce that lead, because it would raise R/T about
  2.5× across the band.
- Stähler et al. (App. A2) state that RR55 carried a **Geomar HTI-01-PCA**
  hydrophone, with nominal pole 0.10774 rad/s and 199.5 µV/Pa. The metadata pole is
  0.0431 rad/s, a DEPAS-type value, so the response may describe the wrong
  hydrophone.
- I could not reproduce RR55, because Sentinel's numbers are scratch output. It needs
  the V3 joint test.

**SENT-005 (timing tolerance): agree.**

- 360·f·Δt = 0.005° at 0.07 Hz for 200 µs, as stated.
- Pressure and Z share one digitizer clock, so OBS clock skew (−1.3 s to +3.4 s over
  about 13 months; Stähler Table 3) is common-mode and does not enter P/Z
  cross-phase.
- A tolerance derived from the phase budget is justified.

**SENT-006 (−5.4% against Zha): narrowed.** New evidence:

| Candidate cause | Finding |
|---|---|
| Depth | Our depths equal Zha's Table 1 at all 25 stations (V2) → **eliminated** |
| Incompressible-water reference | Would shift 1/g by +4.5% to +15% (our band) or +9.6% to +40% (Zha's band), strongly depth-dependent; the observed offset is nearly depth-independent (r = +0.21) → **eliminated** |
| In-situ density | Moves the offset by −0.6% |
| ρ = 1000 | Moves it by +2.8% |
| SENT-010 | About −0.6% |
| Estimator bias (V6) | ±1–3% |

A **depth-independent constant** is the most likely cause: a density convention,
nominal pressure or vertical sensitivity versions, or the estimator. The offset does
not by itself indicate a hardware change.

Rev. 2: Zha & Webb's Appendix A is now read.

**Confirmed by the paper:**
- The reference formula is identical (Eq. A5).
- The coefficient multiplies the theoretical ratio (Fig. A2), so comparing with 1/g
  is correct.
- The fit band is 0.04–0.10 Hz.
- Only M ≥ 7.0 events were used (21 in total), with "P/a coherence > 0.9".
- Theory used the Zha et al. (2014a) Lau average Rayleigh model, not the RR52 model
  used here.

**Not stated in the paper:**
- the water density and sound speed;
- whether coherence was magnitude or squared;
- the fit objective, beyond "minimizing the difference";
- the response or sensitivity version.

Implications:
- **Slower Lau phase velocities.** These raise Zha-style coefficients by only about
  +0.6% in their band (V2-style calculation).
- **Unstated choices are the most likely explanation.** The remaining −5% is most
  plausibly in the unstated choices. For example, ρ = 1000 (+2.8%) together with
  scalar nominal sensitivities (+1.8%, Sentinel) would account for about 4.6%. This
  is a hypothesis.
- **Event set and threshold also differ.** Our M ≥ 6 fallback and γ² ≥ 0.95 threshold
  make the comparison non-equivalent.

**SENT-007 (floors): agree.**

- −3076.5 dB = 10·log₁₀(2.225e-308), ObsPy's `dtiny`.
- RR12's flat value of −8,268,016 counts lies near the negative full scale of a 24-bit
  converter (−8,388,608), which is consistent with a railed channel. Such windows call
  for a **nonworking-channel** classification, not a level comparison.
- Floors must also be excluded from noise-model statistics (SOC-007).

**SENT-010 (taper-then-differentiate): mechanism confirmed (V5).**

- The bias of iω·FFT(w·v) against FFT(w·a) scales as (Δf/f)² times spectral slope.
- For velocity power ∝ f⁰, f⁻², f⁻⁴ and f⁺², the gain bias is −0.2%, **+0.6%**,
  +1.6% and −1.0%.
- The f⁻² case (white acceleration) reproduces Sentinel's real-data +0.65%.
- The bias is largest at 0.03 Hz, where R/T is 0.967–1.023. Because it is event- and
  station-dependent, it also inflates scatter.
- **Correction:** response-correct and differentiate over the full window, or a longer
  padded record, before segmenting. Validate with V5 on synthetic data and on
  re-processed saved events.

**SENT-012 (default sign −1): agree.**

- Physically, correct polarity gives P/a phase ≈ 0°, so a per-station sign diagnostic
  from high-coherence events is justified.
- On YV it is valid only after SOC-001: RR17/RR18 (pole ≈ 0.75 rad/s) carry a 119–152°
  phase error that would masquerade as a polarity inversion.

### Outside scientific scope

SENT-008, -009, -013, -014 and -015 are documentation and software items, outside
scientific scope. I found nothing that contradicts them.

---

## 4. New scientific findings

### SOC-001: DEPAS/Geomar hydrophone metadata have the wrong pole sign. This, not instrument failure, decides most YV outcomes.

- **Class / severity:** Confirmed (metadata error with demonstrated consequence) / High.
- **Location:** source StationXML for 34 YV pressure channels (listed in
  `docs/yv_rayleigh_independent_validation.json`). The consequence appears at the
  phase mask, `src/bruitfm/rayleigh_calibration.py:1238`. Code is not at fault.
- **Physics:**
  - Stähler et al. (2016, Eqs. A5–A6) give the hydrophone response as an RC high-pass:
    G = S·(RC·2πif)/(1 + RC·2πif). This has a zero at 0 and a **left-half-plane**
    pole, p₁ = −1/(R·C_total).
  - A right-half-plane pole is unstable and cannot describe a passive RC circuit.
  - The published text itself prints "p₁ = 0.10774 rad s⁻¹" without the minus sign,
    which is a plausible origin of the metadata error.
- **Evidence:**
  1. **Response amplitude at fixed bins is unchanged (V7).** A +p pole
     leaves |R| identical to 1e-6 and shifts phase by exactly −2·atan(p/ω).
     PPSD amplitudes are unchanged by this reflection alone. Recomputed scalar
     gains can change if phase QC admits different bins or events; this is not a
     promise that a reflected-pole run preserves station gains.
  2. **The deconvolved pressure therefore leads by +2·atan(p/ω).** That is 25–63° at
     0.03 Hz and 11–29° at 0.07 Hz for the metadata values 0.042–0.115 rad/s.
  3. **The pole magnitudes are physical RC values.** R = 500 MΩ or 210 MΩ (Stähler's
     two logger versions) with C ≈ 56 nF give about 0.036 and 0.085 rad/s; C_amp in
     series raises them further.
  4. **Outcome cross-tabulation (V8), 43 YV stations:**

     | Group | Candidates | Unresolved | Band-excluded | Notes |
     |---|---:|---:|---:|---|
     | Hydrophone, pole ≥ 0.06 rad/s (predicted error ≥ 23.6° at every in-band frequency) | 0 | 22 | 1 | — |
     | Hydrophone, pole 0.042–0.052 rad/s | 8 | 2 | 1 | Bins ≥ 0.039–0.047 Hz can pass ±20° |
     | INSU DPG (no positive pole) | 8 | 1 | 0 | The unresolved one is RR31, whose DPG has documented artefacts every ~9000 s (Stähler §3.1) |

  5. The 8 hydrophone candidates' gains are therefore **sub-band estimates**, from
     about 0.04–0.07 Hz only, selected by a metadata artefact.
- **Impact:**
  - YV's "25 unresolved" is mainly a metadata phase artefact.
  - Any per-station sign inference on DEPAS stations is unreliable until the poles are
    corrected.
  - Future compliance and coherence work is unaffected in amplitude but not in
    transfer-function phase.
- **Recommended correction (user decision):**
  1. Report the sign error to RESIF/DEPAS.
  2. As a separate, labelled re-analysis, never overwriting saved products, evaluate
     the phase QC with each positive-real pressure pole reflected into the left
     half-plane. This is justified by Eq. A5.
  3. Re-derive the sign assumption afterwards.
  4. Obtain per-instrument R and C_hyd values from DEPAS.
- **Validation:** after reflection, residual P/a phase should be ≈ 0° ± 10° across
  0.03–0.07 Hz at high coherence, as at the INSU stations. Pole magnitudes should match
  1/(R·C_total) from the DEPAS calibration sheets.

### SOC-002: the INSU hourly levelling glitch is not handled.

- **Class / severity:** Confirmed omission / High for noise models, Medium for gains.
- **Location:** global. A repository search finds no handling of `glitch`, `3620` or
  `PeriodicTransient`. Affected paths: `ppsd.py:98` (vertical PPSDs),
  `rayleigh_calibration.py:472`/`:531` (tilt training) and event windows.
- **Evidence (Stähler et al. 2016, §3.1):**
  - Stations RR28, RR29, RR31, RR34, RR36, RR38, RR40, RR50 and RR52 record a glitch
    every 3620 s.
  - Each lasts about 1200 s, is about 1.5 µm in the 20–500 s band, exceeds vertical
    noise by 15 dB, and is comparable to intermediate teleseismic body waves.
- **Impact:**
  - **PPSD:** every 3 h window contains about 3 glitches. INSU vertical PSDs at
    20–500 s describe the glitch, not the seafloor, which would bias a WP2 low-noise
    model high.
  - **TiSKitPy training:** every 2 h window (0.001–0.01 Hz) contains a glitch.
  - **Event windows:** about 83% of 1800 s windows overlap one. The glitch is
    incoherent with pressure, so it lowers coherence and, for weaker
    (M6–6.5-fallback) events, biases R low and therefore g high. YV used `tiskitpy`
    tilt correction, so both paths apply.
- **Recommended correction:**
  - Remove the glitches with TiSKitPy `PeriodicTransient`, already a dependency, or
    with template subtraction following Deen et al. (2017).
  - Alternatively, mask glitch intervals.
  - Store cleaned products separately.
- **Rev. 2, evidence from Aminian et al. (2025):**
  - The authors' own earlier pipeline (ComPy) removes this transient. They describe it
    as a serial interrogation of the mass-centre positions, inducing "up to 20 min of
    noise" and contaminating "a broad and extended frequency range". The removal uses
    Deen et al. (2017) stacking and subtraction.
  - In their Fig. 3e at RR52, tilt, event and glitch processing lower the median vertical
    PSD near 0.004 Hz by roughly 10–15 dB (read from the figure).
  - BRUIT-FM therefore lacks a step the project's own published method applies.
  - V11: for all eight INSU stations analysed in that paper, BRUIT-FM's candidate gains
    exceed the published ones by +0.3% to +6.1% (median +3.0%).
  - That direction is what uncleaned, pressure-incoherent vertical noise predicts. Other
    method differences could also contribute: tilt method, M ≥ 6 against M > 7 events,
    window picking, and SENT-010.
  - These two estimates share author and method lineage, so this is a reproducibility
    check, not independent validation.
- **Validation:** compare INSU vertical PSDs (20–500 s) before and after; compare
  gains from glitch-free sub-windows and from M ≥ 7 only, and recheck V11 after
  cleaning.

### SOC-003: constant ρ = 1028 kg/m³ and a homogeneous ocean bias the reference low by 0.25–1.3%.

- **Class / severity:** Confirmed model bias / Low.
- **Location:** `rayleigh_calibration.py:848` (default) and `:37`.
- **Evidence (V1).** I solved the linearised acoustic-gravity equations for the
  Lagrangian pressure a moving gauge records: free surface δp = 0, floor ξ_z = u,
  neutral compressional stratification ρ₀(z). Relative to `water_ratio`:

  | Effect | 0.9 km | 2 km | 2.8 km | 4 km | 5 km |
  |---|---:|---:|---:|---:|---:|
  | Gravity alone | +0.05% | +0.1% | +0.2% | +0.3–0.4% | +0.3–0.6% |
  | Gravity + in-situ density | +0.25% | +0.5% | +0.7% | +0.7–1.1% | +0.6–1.3% |

  The g → 0 limit reproduces `water_ratio` to 1e-13. The formula's physics is right;
  this is a small, depth-dependent density convention.
- **Recommended correction:** use a depth-mean in-situ density (TEOS-10 with
  climatological T/S), or carry it as a systematic term. Report the value used in
  every comparison with published coefficients.
- **Rev. 2:** the literature uses inconsistent conventions:
  - Ruan et al. (2014): 1030 kg/m³.
  - Aminian et al. (2025): ρw ≈ 1040 kg/m³ in the gravity-attraction term.
  - Zha & Webb (2016): density not stated.
  - BRUIT-FM: 1028 kg/m³.

  Adopt one documented convention across BRUIT-FM's P/a and compliance work.
- **Validation:** V1.

### SOC-004: the unweighted absolute-unit fit weights the most model-sensitive bins most heavily.

- **Class / severity:** Concern / Medium.
- **Location:** `fit_multiplier`, `rayleigh_calibration.py:77`.
- **Evidence:**
  - Minimising Σ(gR − T)² weights bins by T². T grows toward the band top and with
    depth.
  - These are the bins where the reference is least certain (V10):
    - d ln T/d ln H = 1.04–1.23 at 2 km, 1.16–2.63 at 4 km, and up to 6.5 at 5 km and
      0.07 Hz.
    - A steep arrival (apparent velocity 20 km/s, e.g. body-wave energy) changes T by
      +8% at 4 km and +22% at 5 km at 0.06 Hz.
  - Hann bins are not independent: adjacent-bin power correlation is 0.44, so "≥ 3
    valid bins" means fewer independent measurements.
  - Rev. 2: Deng et al. (2022, Fig. 4f) show the same physics. Near c₀/4h, the ratio
    F = p/(ρha) becomes strongly sensitive to the horizontal slowness q (= 1/c). They
    exploit that to invert phase velocity at 0.06–0.07 Hz in 4.9 km of water. In a
    calibration, that sensitivity is model error. For YV stations at 4–5 km
    (c₀/4h = 0.075–0.094 Hz), the top of the 0.03–0.07 Hz band is in this regime.
- **Recommended correction:**
  - Fit ln(T/R) with inverse-variance weights, using
    var(ln|H|) ≈ (1 − γ²)/(2Kγ²) (Bendat & Piersol, 2010).
  - Report the sensitivity of g to removing the band's top bin.
- **Validation:** refit saved spectra with both objectives and tabulate differences
  per station. This is read-only.

### SOC-005: estimator noise bias of ±1–3% is uncorrected and uncounted.

- **Class / severity:** Concern / Medium.
- **Location:** `event_spectra` and the fit, `rayleigh_calibration.py:594-639`.
- **Evidence (V6),** Monte Carlo of the production QC: 9 segments, per-bin γ² ≥ 0.95,
  event median and mean ≥ 0.95, ≥ 3 bins.

  | Noise case | Accepted events | Bias in g from R = √(Spp/Saa) |
  |---|---:|---|
  | Pressure noise only, true γ² = 0.95 | ~52% | −2.8% |
  | Pressure noise only, true γ² = 0.98 | — | −1.2% |
  | Acceleration noise only, true γ² = 0.95 | ~52% | +2.2% |
  | Acceleration noise only, true γ² = 0.98 | — | +0.9% |
  | Balanced noise | — | ≈ 0 |

  R is the geometric mean of |H1| and |H2|, so it is robust only when the noise is
  balanced. The sign of the bias is not known a priori.
- **Recommended correction:**
  - Estimate per-channel noise PSDs from the pre-event (quiet) window and form a
    noise-corrected ratio √((Spp − Npp)/(Saa − Naa)).
  - Report R, |H1| and |H2| together; their spread bounds this bias.
- **Validation:** V6 on synthetic data; on real events, compare corrected and
  uncorrected ratios.

### SOC-006: the scalar-gain model needs an explicit adequacy test.

- **Class / severity:** Concern, demonstrated at B01W / High.
- **Location:** `fit_multiplier`; QC flag `scalar_gain_frequency_independence`.
- **Evidence:**
  - B01W is a single-pole high-pass, not a scalar (SENT-004, V3).
  - Doran et al. (2019) find DPG sensitivity and time constant deviate from nominal
    by ≥ 2×.
  - Stähler (App. A3): seven of the nine INSU DPG responses copy RR28's land
    calibration, and the corner depends on oil viscosity (temperature, pressure).
  - DEPAS hydrophone corners depend on R (210 or 500 MΩ by version) and C_hyd (±6%).
  - Rev. 2, Deng et al. (2022, §4 and Figs S4–S7): narrow-band (0.01 Hz) DPG
    calibration at 29 Alaska stations finds coefficients that are often
    frequency-dependent. At WD66 the coefficient falls from 0.91 to 0.79 across
    0.02–0.07 Hz; at WD70 it ranges from 0.91 to 1.71. Some stations are nearly
    flat, e.g. WD48 ≈ 0.95.
  - Rev. 2: Zha & Webb (2016, App. A) note that a response shift below the Rayleigh band
    "is not detected with this method".
  - Scalar inadequacy is therefore common, not exceptional.
- **Recommended correction:**
  - Fit (G, fc) jointly to complex P/a. Accept a scalar only when fc is
    unconstrained below the band and the amplitude–phase (minimum-phase) residuals
    are small.
  - Never extrapolate a corrected response outside the fit band.
- **Validation:** V3 method. It recovers controls with no in-band corner and B01W's
  corner consistently across events.

### SOC-007: noise-model statistics need usable-band censoring and exclusions.

- **Class / severity:** Concern / Medium now, High for the final WP2 low/high models.
- **Location:** `ppsd.py:105–175` (to 0.001 Hz, no instrument floor);
  `pressure_test.py:312` and `infragravity_test.py:94` (percentiles of raw values).
- **Evidence:**
  - At 0.001 Hz the DEPAS CMG-OBS40T (60 s corner) and hydrophones (58–100 s corner)
    are about 25–50 dB below midband.
  - Stähler reports strong DEPAS long-period self-noise.
  - Response-corrected PSDs there describe instrument noise.
  - Floors (SENT-007), C07W reduced days (SENT-001) and INSU glitches (SOC-002) would
    also enter the statistics.
- **Recommended correction:** define each channel's usable band from its response and
  self-noise model (Ringler & Hutt, 2010), as in Janiszewski et al. (2023). Censor
  outside it, exclude floor and flagged intervals, and publish raw and tilt- or
  compliance-corrected vertical statistics separately.
- **Validation:** per-channel PSD percentiles against self-noise curves; stability of
  the low-noise percentile under these exclusions.

### SOC-008: no epicentral-distance bounds; packet isolation remains heuristic.

- **Class / severity:** Concern / Low–Medium.
- **Location:** `rayleigh_calibration.py:939–953` (candidate filter: M ≥ 6, depth
  < 150 km, Mw, no distance limit).
- **Evidence:** T(f) assumes a single plane fundamental-mode wave. Regional events
  (< ~20°) and near-antipodal ones (> ~160°) violate this. 20 of 48 YV candidate
  windows contain predicted P/S, and V10 quantifies their effect at depth.
- **Recommended correction:** add distance bounds, e.g. 20–160°, chosen and justified
  explicitly. Isolate the fundamental mode by group-velocity windowing or
  time-variable filtering.
- **Rev. 2, an independent calibration path.** Deng et al. (2022, §4) calibrate DPGs
  with teleseismic body waves of near-vertical incidence (ray parameter < 10 s/deg),
  in narrow frequency bands.
  - With q ≈ 0, F depends only on water depth and sound speed. The calibration is then
    **independent of the Earth model** and of Rayleigh-mode isolation.
  - Applied to the same stations, it would test the Rayleigh gains independently. It
    would also give the frequency dependence (SOC-006).
- **Validation:** compare gains with and without regional and steep-arrival-dominated
  windows, using saved spectra where possible.

### SOC-009: the RR52 Indian Ocean Earth model is reused for the Lau back-arc.

- **Class / severity:** Concern / Low.
- **Evidence:**
  - The table has c = 3.92–4.05 km/s with no water layer. Lau back-arc phase
    velocities are expected to be lower.
  - The effect on T is ≤ 1.3% (Sentinel). For c = 3.6 km/s it is 0.1–0.5% in
    0.03–0.07 Hz and 0.3–1.3% in 0.04–0.10 Hz (V2-style calculation), so it is minor.
    It still must not be described as a validated regional model.
- **Recommended correction:** use a regional model that includes the water layer.
  `rayleigh_dispersion.py` already accepts a fluid top layer.

### SOC-010: depth inputs. A positive check, with a caution about published tables.

- **Class / severity:** Concern / Low.
- **Evidence:**
  - StationXML depths match GMRT (Ryan et al., 2009; queried 2026-10-12) within 0.8%:

    | Station | StationXML | GMRT |
    |---|---:|---:|
    | RR03 | 4340 m | 4327 m |
    | RR11 | 3473 m | 3470 m |
    | RR52 | 2918 m | 2895 m |
    | RR54 | 3375 m | 3392 m |

  - Stähler Table 3 instead gives RR11 = 3941 m and RR54 = 2499 m. These appear to be
    typos and should not be used.
  - RR52's three sources span 2880–2918 m. Using 2880 m would lower RR52's g by 1.7%.
- **Recommended correction:** check all stations against GMRT/GEBCO and propagate a
  depth uncertainty using V10's d ln T/d ln H.

### SOC-011: the expected-level references are screening aids only.

- **Class / severity:** Concern / Low.
- **Evidence:**
  1. **Brown digitisation.** The repository digitisation differs from TiSKitPy's
     independent Brown et al. (2014) curves by up to 3.4 dB, with 8 of 34 values
     beyond the assumed 2 dB tolerance (V9). The tolerance in
     `pressure_references.json:21` is optimistic.
  2. **Brown is not a seafloor model.** Its curves describe IDC hydroacoustic
     stations (moored, SOFAR-channel hydrophones), not the seafloor.
  3. **Godin is unverified and regional.** Godin's normalisation is unverified here
     (§2) and was derived for New Zealand. Its level is plausible: 34.6 dB at 4 km and
     0.005 Hz, within Webb's 10–50 dB bounds.
- **Recommended correction:**
  - Keep these as screens with ≥ 3–4 dB tolerances.
  - Verify Godin Eq. 6 and Webb's bounds against the primary texts before using them
    in any score.

### SOC-012: physics for the two pending calibration checks.

- **Class / severity:** Concern / Low, for future work.
- **Microseisms (Calibration 2).**
  - The first water resonance α/(4H) is 0.075–0.19 Hz for 2–5 km, inside the
    microseism band. The scalar P/a reference has poles there.
  - A modal (acoustic-elastic) model is required, not the Rayleigh formula.
- **Microseisms, added in Rev. 2.** Ruan et al. (2014, Figs 4–6) show that the
  Rayleigh-wave D/P ratio at 0.1–0.2 Hz depends strongly on sediment thickness and
  shear velocity. In their inversion, the DPG scale factor is constrained mainly below
  0.1 Hz, where sediment sensitivity is small. A microseism-band check is therefore
  not independent of an assumed sediment model. It needs a modal model and, where
  possible, a scale factor fixed from below 0.1 Hz.
- **Local and teleseismic P (Calibration 3).** Rev. 2 corrects my earlier proposal.
  - I had proposed the impedance relation p = ρα·v_z/cos θ_w before the surface
    reflection returns. Deng et al. (2022, §2) note that p = ρc₀v follows only from the
    first transmission. It does not reduce to ρha at low frequency and matches
    simulations only for the first crest.
  - The appropriate reference is Deng et al. Eq. 5, which equals Zha Eq. A5 with q set
    to the ray parameter. They validate it for body waves with an APG up to slightly
    above c₀/4h.
  - With q < 10 s/deg it is nearly Earth-model-free (SOC-008). It requires
    pressure-sensor bandwidth at the chosen frequencies; the INSU DPG passband is
    0.002–30 Hz (Stähler).

### SOC-013: an independent network-level check of the Rayleigh gains exists.

- **Class / severity:** Concern, opportunity / Low.
- **Evidence (V8):**
  - Median candidate gains are 0.698 for INSU DPGs and 1.140 for DEPAS hydrophones.
  - These imply nominal INSU pressure levels 4.27 dB above DEPAS.
  - Stähler et al. (2016, §3.3.3) independently report DEPAS hydrophone levels
    "approximately 5 dB lower … than on INSU stations in the period range of 12–40 s",
    the same band as the calibration. They attributed this to instrument type.
- **Recommended use:**
  - Formalise this as a level test. Compare quiet-time 0.03–0.07 Hz pressure PSDs of
    neighbouring INSU/DEPAS pairs (nearby test) before and after candidate gains.
  - In the nearby test generally, separate a frequency-independent offset
    (20·log₁₀ g) from spectral shape, and depth-normalise the IG band by cosh(kH).

### Minor (no action needed)

- **ObsPy dB-domain averaging.** ObsPy averages octave bins in dB (geometric mean).
  log(χ²) gives about −0.3 dB bias for about 15 degrees of freedom. This is consistent
  with Peterson/McNamara conventions and is already handled when compared with
  linear-mean Godin pilot spectra.
- **No-overlap Welch.** Non-overlapping Hann segments in `event_spectra` discard
  edge-weighted data; 50% overlap would nearly double the effective degrees of
  freedom (Welch, 1967). This is an efficiency point, not a bias.
- **K = 9 significance.** For K = 9 the 95% zero-coherence level is 1 − 0.05^(1/8) =
  0.31. The 0.95 threshold is very conservative, and acceptance changes sharply around
  true γ² = 0.95 (V6: about 52% accepted).
- **Godin-pilot window.** The Godin-pilot window of 12.42 h (the M2 period) is good
  practice against tidal leakage.
- **Post-earthquake exclusion (Rev. 2).** The quiet-training coda rule
  `max(3600, (M−5.5)·43200)` s gives 18 h after an M7. Aminian et al. (2025),
  following Deen et al. (2017), exclude (Mw − 5.85)·36 h, which is 41.4 h after an
  M7. BRUIT-FM's rule is less conservative for the 2 h tilt-training window; justify
  the choice or align it.
- **Tilt stability (Rev. 2).** Aminian et al. (2025, Fig. 4) show the RR52 tilt angle
  drifting from about 1° to 2.3° over months, with a jump on 2013-05-02. A 2 h training
  window 6–48 h before an event is adequate except across such jumps. Check events
  near them.

---

## 5. Validation performed

All validation is in `docs/reviews/socrates_validation.py`. It is read-only and
deterministic; the seeds are fixed.

| ID | Test | Result |
|---|---|---|
| V1 | Water reference against the acoustic-gravity ODE | Exact in the g→0 limit; +0.05% to +1.3% physical corrections |
| V2 | Model variants against the −5.4% Zha offset | Depth eliminated; incompressible water eliminated; ρ conventions ±2.8% |
| V3 | B01W single-pole high-pass, joint fit | fc 0.073–0.091 Hz, G 0.26–0.30; controls have no corner. Rev. 2: the model's 0.04–0.10 Hz coefficient is 0.17–0.19, against Zha's 0.19 ± 0.02 |
| V4 | C07W state ratio | 1.91 / 1.97 / 1.96, consistent with a frequency-independent factor of about 2 |
| V5 | Taper-then-differentiate bias | −1.0% to +1.6% depending on slope; +0.63% for white acceleration |
| V6 | R, H1, H2 bias under production QC | ±1–3% depending on noise placement |
| V7 | evalresp / NumPy conventions; Pa through `VEL`; +p poles | All verified |
| V8 | YV outcome against the metadata pole | 0 of 23 candidates for p ≥ 0.06; INSU/DEPAS gains imply 4.27 dB |
| V9 | Brown digitisation against TiSKitPy | Up to 3.4 dB difference |
| V10 | Depth and apparent-velocity sensitivity of T | Tabulated |
| V11 (Rev. 2) | BRUIT-FM INSU gains against Aminian et al. (2025) Table 1 | All 8 higher, +0.3% to +6.1% (median +3.0%) |

Checks outside the script, with their sources:

- Stähler PDF extraction (Table 3, Appendices A).
- GMRT point depths (four queries).
- A repository-wide search for glitch handling.
- Review of ObsPy's PPSD source.

---

## 6. Indicative uncertainty budget (per station gain, scalar model valid)

| Term | Size | Source |
|---|---|---|
| Event scatter (currently reported) | often < 2% | saved results |
| Depth (±1%) | ±1.0–1.2% at 2 km; ±1.2–2.6% at 4 km | V10 |
| Density convention | +0.25% to +1.3% (one-sided) | V1 |
| Estimator noise bias | ±1–3% | V6 |
| Taper/differentiation | −1.0% to +1.6% | V5 / SENT-010 |
| Phase velocity / Earth model | ≤ 1.3% | Sentinel; SOC-009 |
| Steep-arrival contamination | 0 to +8% (4 km), larger at 5 km, at band top | V10 |
| Vertical absolute calibration | unknown; enters 1:1 | flag `vertical_absolute_response_unverified` |
| Response shape (B01W-type) | unbounded if not tested | SOC-006 |

Even for a well-behaved deep station, the combined systematic uncertainty is likely
±3–6%. That is two to three times the event scatter, before vertical calibration is
counted. Gains should be reported with such a budget, or with bootstrap intervals that
include these terms, before any acceptance or score is assigned.

---

## 7. Missing scientific validation, in order

1. **Pole-reflected phase QC** on YV hydrophones (SOC-001), as a separate labelled
   product, and RR55's metadata identity.
2. **Glitch removal** for the nine INSU stations, then a PSD and gain comparison
   (SOC-002).
3. **Two-parameter response fits** with an amplitude–phase consistency test at all
   stations (SOC-006); apply it first to B01W, RR55, RR22 and RR56.
4. **Independent vertical calibration** (relative gains, e.g. tilt-free teleseismic
   amplitudes against neighbours, as in SENT-001), because every gain inherits it.
5. **Weighted and noise-corrected refits** of saved spectra (SOC-004/005), and the
   uncertainty budget (§6).
6. **Distance bounds and fundamental-mode isolation** (SOC-008).
7. **Usable-band censoring** for noise statistics (SOC-007).
8. **Primary-source verification** of Godin Eq. 6 and Webb's bounds, plus a review
   of the thesis calibration chapter. Zha & Webb Appendix A was verified in Rev. 2;
   its unstated ρ, α and response version should be requested from the authors if the
   offset matters.
9. **Physics-based methods for Calibrations 2 and 3** (SOC-012), and the INSU/DEPAS
   level test (SOC-013).
10. **Near-vertical body-wave calibration** in narrow bands, following Deng et al.
    (2022). It is an Earth-model-independent check of the Rayleigh gains and of their
    frequency dependence (Rev. 2).
11. **Gain stability through time**, by splitting events into deployment halves.
    Aminian et al. (2025) report stable INSU gains; C07W shows why this must be
    tested per station (Rev. 2).

---

## 8. Decisions for the user, and proposals for Codex

These are not applied.

**Decisions for the user:**

- whether to commission the SOC-001 pole-reflected re-analysis, and the report to
  RESIF/DEPAS;
- the glitch treatment (SOC-002);
- adopting a two-parameter response model (SOC-006);
- the fit weighting and noise correction (SOC-004/005);
- the distance bounds (SOC-008);
- the density convention (SOC-003).

**Proposals for Codex, after approval:**

- a glitch-cleaning stage that writes new products: TiSKitPy `PeriodicTransient`, or
  porting the ComPy step used in Aminian et al. (2025);
- a near-vertical body-wave calibration diagnostic following Deng et al. (2022)
  (Rev. 2);
- a reflected-pole phase diagnostic, labelled and never overwriting;
- full-window response correction before segmenting (SENT-010);
- a weighted log fit and noise-corrected ratio reported alongside the current
  estimator;
- usable-band metadata per channel;
- the missing references.

**For Sentinel:**

- SENT-001 and SENT-004 B01W are supported and quantified here.
- SENT-004's RR55 equivalence is **not** established (§3).
- SENT-006 can record that depth and incompressible-water explanations are eliminated.
- Rev. 2: the handoff's "sources I could not read" list can drop Zha & Webb (2016).
  Their Appendix A and Table 1 now confirm the formula, the coefficient direction, the
  band and the event policy.
- Rev. 2: the B01W "leaky valve" reading is the authors' own (§2.3).

---

## 9. Revision 2: what the supplied papers changed

Codex integration-review qualifications (2026-10-12): V4's intervals are IID
day-bootstrap intervals. They omit serial dependence within multi-day amplitude
episodes and are descriptive, not calibration confidence intervals. The two
amplitude states are reproducible evidence; the specific hardware failure still
needs independent confirmation. Agreement with Aminian (2025) does not establish
glitch removal as the cause of the remaining gain offset.

| Item | Before | After reading the primary texts |
|---|---|---|
| Water-column reference | Verified by derivation and ODE | Also identical to Zha & Webb Eq. A5 and Deng et al. Eq. 5 |
| Coefficient direction (1/g against Zha) | Inferred from the repository | Confirmed by Zha & Webb Fig. A2 |
| SENT-006 offset | Primary method unknown | Band, events, coherence and velocity model known. ρ, α, objective and response version remain unstated, so the offset stays unattributed. |
| SENT-004 B01W | Physical model fitted (V3) | Corroborated by Zha & Webb's "leaky valve" note. The model reproduces their 0.19. Nominal τ is 600–1000 s, against τ ≈ 2 s fitted. |
| SENT-001 C07W | Mechanism hypothesis based on the Trillium | YL sensors are Güralp CMG3T. Zha's 11-event coefficient matches the normal state. |
| SOC-002 glitch | Stähler et al. only | The project's own 2025 pipeline removes it. BRUIT-FM INSU gains are +3.0% (median) above the published values, consistent in direction. |
| SOC-006 scalar gain | B01W, plus Doran et al. abstract | Frequency-dependent DPG coefficients are common (Deng et al., 29 stations) |
| SOC-012 local P | Proposed the ρα·v impedance relation | Corrected: use Deng Eq. 5 with the ray parameter. The impedance relation holds only for the first crest. |
| New independent check | — | Near-vertical body-wave calibration (Deng et al.), nearly Earth-model-free |
| Microseism calibration | Modal model needed | Also sediment-dependent (Ruan et al.); scale factor best fixed below 0.1 Hz |

No conclusion of the first revision was reversed, except the local-P recommendation
in SOC-012. The five highest-priority issues and the overall readiness assessment are
unchanged.
