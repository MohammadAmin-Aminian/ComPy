# Local-earthquake P-wave pressure–vertical check: method note

Calibration check 3 of the six-check master table. This note records physics,
conventions, applicability limits, fixed defaults and disclosed development
revisions. It is not a blind preregistration or a validated field calibration.
The reviewed v2 implementation and bounded validation are documented in
[the phase closeout report](local_p_phase_report.md).

## 1. Sources and what each supplies

| Source | Sections/pages read | Used for |
|---|---|---|
| Deng, An, Cai & Ren (2022), GRL 49, e2021GL096952 | §2 Eq. 1–6 (pp. 2–3); §3 LT02 test (pp. 3–4, Fig. 3); §4 DPG calibration (pp. 4–5, Fig. 4); §5 (pp. 5–7); §6 (p. 7) | Full plane-wave water-column relation (Eq. 5), its q-dependence, validation band, first-crest statement |
| An, Cai, Zheng, Meng & Liu (2017), GRL 44, 10,272–10,281 | §2 Eq. 1a–6 (pp. 10,273–274); §4.1 (pp. 10,276–277); §5 (pp. 10,278–279) | Explicit z-up sign convention; incompressible limit p = ρha; validity conditions; multiplier convention; SIO polarity |
| Zha & Webb (2016), JGR 121, 3220–3237 | §2.3 (p. 3224), Table 1; Appendix A Eq. A1–A6 (pp. 3234–3236) | Same relation with phase velocity; recorded-over-predicted gain convention; frequency-independence assumption |
| Aminian et al. (2025), GJI 242 | §4.1.5 (pp. 7–8), Table 1 | Multiplier convention of the repository Rayleigh gains (RR52 0.67) |
| Auger (2024), M1 internship report, 15 pp. (bruit-fm.org/Documents/Auger2024.pdf, SHA-256 `19859509…a7b342`) | Eq. 1–6 (pp. 3–8), data (p. 5), Figs 6–11 (pp. 6–10), annexes | Prior local-event attempt (Mayotte 1T, Lucky Strike 4G); its limitations |
| SOC-012 Rev. 2 (`docs/reviews/socrates_initial_audit.md`) | SOC-012 and correction table | Use Deng Eq. 5 with the ray parameter; p = ρc₀v only for the first crest |

The Deng et al. supporting information (Text S1 synthetics, Table S1, Figs
S4–S9) was **not available locally and was not read**. Statements here about it
come only from the main text. The prompt calls Auger (2024) an M2 report; the
document itself is an M1 report dated April 2024.

## 2. Conventions

- **Units.** Pressure in Pa, vertical velocity in m/s, acceleration in m/s²,
  depth in m, slowness in s/km internally (SI s/m in formulas), frequency in Hz.
- **Vertical.** Positive **upward**. Responses are evaluated to the declared
  input units (Pa; m/s) and divided out in the frequency domain exactly once.
- **Pressure.** Physical excess pressure, **compression positive**.
- **Fourier.** NumPy forward transform `X(f)=Σx·e^{-i2πft}`, so `d/dt ↔ +i2πf`
  and `A(f)=i2πf·V(f)`. Both model ratios used (p/a below resonance, p/v before
  the first multiple) are **real**, so their phase is convention-free.
- **Gain.** The repository convention (`fit_multiplier`; Aminian 2025;
  An 2017; Deng's C_p): **g is a pressure multiplier**, `g·p_nominal = p_model`,
  so `g = |H_model| / |H_obs|`. Zha & Webb (2016) Table 1 report the inverse
  (recorded/predicted); comparisons must invert it.
- **Polarity** is reported separately from the positive amplitude gain: the
  phase of `H_obs/H_model` near 0° → +1, near ±180° → −1, otherwise
  indeterminate. A nominal response with negative sensitivity (e.g. XO.WD66 HDH,
  −1008 counts/Pa) already carries a sign; polarity is relative to that response.

## 3. Physics

### 3.1 Plane wave in a uniform water layer

Water depth h, sound speed c₀, density ρ, horizontal slowness q (ray
parameter for body waves, 1/c for surface waves), free surface above. Write
`η = √(1 − q²c₀²) = cos i₂` (water incidence angle i₂, Deng p. 3) and
`θ = ωhη/c₀`. For a standing field with p = 0 at the surface (z up, seafloor
z = 0), `p(z) = P sin(k_z(h−z))`, `k_z = ωη/c₀`. Linear momentum
`ρ∂v_z/∂t = −∂p/∂z` gives at the seafloor

    p / a_up = ρ tan(k_z h) / k_z = ρ h · F,     F = tan θ / θ          (A)

This is **Deng Eq. 5–6** and **Zha Eq. A5** (`r = ω(1/α² − 1/c²)^{1/2}`,
`P/a = ρ sin(rH)/(r cos(rH))`). Limits and properties:

- `θ → 0`: `F → 1 + θ²/3`, so p = ρha (**An Eq. 6**, Filloux). Positive sign:
  upward seafloor acceleration raises bottom pressure (An Eq. 1b–1d, z up).
- Poles at `θ = (2n+1)π/2`, i.e. `f_n = (2n+1)c₀/(4hη)`; zeros at `θ = nπ`.
  At a pole the seafloor vertical motion vanishes for finite pressure. F < 0
  between the first pole and first zero. The model is used only below the first
  pole with a 10% margin (`θ ≤ 0.9·π/2`, same as `rayleigh_calibration.water_ratio`).
- **The fluid–solid reflection coefficient R cancels** (Deng Eq. 3: the factor
  `[1 + R e^{−2iθ}]^{-1}` multiplies p and a alike). For one plane wave with
  known q, Eq. A is independent of the sub-seafloor structure. Different q,
  scattered energy or lateral heterogeneity break this.
- Deng's printed potentials (Eq. 1–3) are not internally consistent in the
  direction of z (the free-surface term vanishes at z = −h while φ₁ propagates
  toward +z). Magnitudes are unaffected. The sign above is derived from first
  principles and checked by the independent simulation (§7). This is a likely
  source of the unexplained π shift in Auger (2024, p. 9).

### 3.2 Slowness and incidence for P waves

`q[s/km] = ray_param[s/rad] / 6371 km = ray_param[s/deg] / 111.195`.
For any P-type wave, `q ≤ 1/v_P,min`. With a conservative `v_P,min = 5 km/s`
and c₀ = 1.5 km/s: `q c₀ ≤ 0.30`, `η ∈ [0.954, 1]`, `i₂ ≤ 17.5°`. Converted
phases generated by the same incident wave keep the same q in a laterally
uniform medium. The q-uncertainty for a P window is therefore **bounded** and is
propagated as a model range, never fitted. Deng (§4) use teleseismic waves with
q < 10 s/deg (0.09 s/km) so that F is insensitive to q.

Auger (2024, pp. 10–11) obtains best-fitting water incidence angles of 45°
(local), 0° (low-frequency events) and 85° (teleseisms). `sin i₂ = qc₀` makes
45° (q = 0.47 s/km) and 85° impossible for P waves. The angle there absorbs
other mismatch, so **incidence is not a free parameter here**.

### 3.3 Two time regimes

Eq. A is the **steady state** after the water-column reverberations
(`φ₂, φ₃, …` in Deng Eq. 2) have built up. The round-trip time is
`τ = 2hη/c₀`. Before the first surface multiple returns (`t − t_P < τ`), only
the upgoing transmitted wave φ₁ exists at the seafloor, and exactly

    p = (ρ c₀ / η) · v_up           (first transmission, t − t_P < τ)    (B)

(Auger Eq. 4; Deng §2 p. 3; SOC-012). Eq. B is **not** a general full-window
relation: it ignores the surface reflections, does not reduce to ρha at low
frequency, and applies only to the pre-multiple segment ("first crest" in Deng
and Saito 2017; at local-P frequencies that segment contains several cycles).

| Regime | Condition (window T after onset) | Model | Use here |
|---|---|---|---|
| A: steady state | `T ≥ N_rt · τ_max` (N_rt = 10) and `f ≤ 0.9·f₀(η_min)` | Eq. A | **Conditional gain candidate** |
| B: pre-multiple | window ends before `t_P + τ_min` | Eq. B | **Diagnostic only**, never a gain |
| Mixed | neither | none valid | Unresolved, `window_regime_mixed` |

`τ_min` uses the shallower depth bound, `η_min` and the largest c₀; `τ_max` uses
the deeper bound, η = 1 and the smallest c₀. N_rt = 10 lets a reverberation with
|R| = 0.8 decay to 0.8¹⁰ ≈ 0.11 (|R| ≈ 0.3 for sediment, up to ~0.8 for bare basalt).

Consequences for the available data (c₀ = 1500 m/s):

| Depth | τ | f₀ = c₀/4h | Feasible regime for local P (S−P of a few to ~30 s) |
|---|---|---|---|
| 152 m (XO.LT02) | 0.20 s | 2.47 Hz | A below ~2.2 Hz with windows ≥ 2 s |
| 2.6–2.9 km (7D, YV.RR52) | 3.5–3.9 s | 0.13–0.14 Hz | A needs ≥ 39 s windows and resolves poles every 0.26 Hz; local P is above f₀ → only B |
| 4.5–4.9 km (XO.LA33/34/WD66, YV) | 6.0–6.6 s | 0.08 Hz | only B (window ≤ ~5 s) |

Deng observed real-data agreement only up to slightly below c₀/4h. At
2.5–3 Hz (above resonance, LT02) only the first half of each packet was
reproduced (§6). Their validation used near-vertical teleseismic body waves,
not local P.

### 3.4 Applicability of Deng's body-wave method to local P

| Issue | Teleseismic body waves (Deng) | Local/regional P (this check) |
|---|---|---|
| q | < 0.09 s/km, F insensitive | up to ~0.2 s/km; bounded by `η ≥ 0.954`, reported as range |
| Window vs reverberation | long coda windows, steady state | S−P limits windows; steady state only in shallow water |
| Frequency | below c₀/4h (0.02–0.07 Hz deep, 1.5–2 Hz LT02) | local P energy 1–10 Hz: above f₀ in deep water |
| Wavefront | plane | curvature; require hypocentral distance ≥ 10 km |
| Contamination | isolated phases | S, P coda with other q, overlapping events |
| Sensor band | long-period | needs pressure bandwidth to ≥ several Hz; OBS coupling resonances above ~10 Hz |
| Timing | Deng had to shift the LT02 APG by 93–103 s (separate logger) | same issue; onset/lag diagnostics required |

## 4. Estimator

1. **Events.** Catalogue events (provenance, location/depth errors kept) with
   epicentral distance ≤ 300 km, hypocentral distance ≥ 10 km, inside both
   channel epochs and waveform coverage. TauP `iasp91` predicts first P, first S
   and the P ray parameter (receiver at the seafloor). Candidates are ranked
   deterministically by magnitude and distance; the number processed is capped.
2. **Onsets.** An AIC onset (Maeda 1985) is measured independently on the
   vertical and on the pressure channel within a search window around the
   predicted P. These are measurements, labelled `auto_aic`, not catalogue picks.
   Their difference is a timing diagnostic: beyond 0.1 s the event is rejected
   unless a lag is explicitly estimated (`--timing-lag-search S`): a coarse lag
   maximizing the 0.5–2 Hz envelope correlation within ±S, refined by the
   onset difference (≤ 3 s) and rounded to a sample. A sign flip does not move
   an onset, so polarity remains testable; the phase slope (timing) does not.
   A third mode predicts the lag from a fitted logger clock-drift line and
   refines it the same way. The mode per channel is a recorded decision
   based on timing evidence only (`docs/local_p_timing.json`).
3. **Windows.** Signal window starts 0.2 s before the vertical onset. Its
   length is the regime limit, capped so that it ends at least 1 s before
   `t_P + 0.85·(t_S − t_P)`. Noise has the same length and ends **at least** 1 s
   before the signal and no later than predicted P minus its search uncertainty.
   V2 takes the earlier of these endpoints, rounded down on the working grid.
   If necessary it moves the noise window earlier; it never moves the signal
   window or relaxes a threshold. Insufficient padding is rejected. Event
   records save noise start/end, earliest predicted P and the guard result.
   Events with another catalogue P arrival in the expanded noise/signal span
   are rejected. V1 documented this guard but omitted its enforcement; v1
   outputs remain archived and explicitly marked legacy/unverified on export.
4. **Forward model.** The whole event segment of vertical counts is converted to
   upward velocity with the nominal complex response, multiplied by the regime
   model (`p/v = i2πf·ρhF` in steady state; `ρc₀/η` pre-multiple) and converted
   back to counts with the nominal pressure response: the **predicted pressure
   record**. Observed pressure, predicted pressure and vertical are band-limited
   by the same zero-phase cosine weight (0.1–0.2 Hz rise; steady-state top
   0.85–0.9 of the *lowest possible* first pole `c₀,min/(4h_max)`; pre-multiple
   top 0.9–1.0 of 0.8·Nyquist).
5. **Spectra.** DPSS multitaper (NW = 2.5, K = 4) on identical signal and noise
   windows of observed pressure, predicted pressure and vertical. The K tapers
   give independent averages; coherence from a single periodogram is never
   reported. The gain per bin is the ratio of predicted to observed
   pressure–vertical cross-spectra, `g₁ = |S_pred,v|/|S_obs,v|` (upper bound:
   vertical noise), `g₂ = |H₁,pred|/|H₂,obs|` (lower bound: pressure noise), and
   the noise-subtracted `g = √((S_pred−N_pred)/(S_obs−N_obs))` (primary).
6. **Valid bins.** squared coherence (observed pressure, vertical) ≥ 0.8; both
   signal/noise power ratios ≥ 10; `|R(f)|/|S| ≥ 0.1` on both responses;
   `f ≤ 0.8·Nyquist`; `f ≥ 2W` (`W = NW/T`); at least 2W inside the flat part of
   the passband weight. Polarity is classified from `arg(S_obs,v/S_pred,v)`.
   Model ranges from depth ± 5%, η ∈ [η_min, 1], c₀ ∈ [1450, 1540] m/s and
   ρ ∈ [1020, 1050] kg/m³ scale the gain into a sensitivity interval per bin.
7. **Bands.** Octave bands 0.5–1, 1–2, 2–4, 4–8 and 8–16 Hz. Per event and
   band: median gain over valid bins (≥ 3 bins). Per station and band: median
   and spread over events. No single scalar is forced: a scalar summary is
   given only when band medians agree within their event spread, otherwise the
   result is reported as frequency dependent.
8. **Station outcome** (scientific, never a score):
   `conditional_gain_candidate` (regime A, ≥ 3 events in a band);
   `diagnostic_only` (regime B results, or regime A with < 3 events);
   `unresolved` (events processed but no valid band);
   `no_suitable_events`; `not_applicable` (e.g. sample rate too low).
   Execution status (`completed`, `failed`, `skipped`) is separate.

The defaults are fixed for routine runs and not tuned to reproduce reference
gains. The development revisions below were made after preliminary observations,
so this work must not be described as fully preregistered blind validation.
V2 repairs the documented noise-window contract without changing QC thresholds.
Outputs keep every rejection reason and the protocol version.

**Revision after the first development runs (XO.LT02, before the pilot).**
The original design compared a windowed spectral ratio pointwise with F(f),
with a guard `|F(f±W)/F(f) − 1| ≤ 5%`. At LT02 the guard excluded every bin
with local-P signal (local P sits under the noise below ~1 Hz, while F varies
by >5% across ±W above ~1 Hz). The underlying problem is an estimator bias: a
windowed estimate of a frequency-varying transfer function is not the pointwise
model. Steps 4–5 (forward model through identical windows) replace the guard.
The pre-multiple model is flat, so its results are unaffected in principle.
For transparency: preliminary per-bin LT02 values (≈1.1–1.2 at 1.5–2 Hz) had
been printed during debugging before this change. No threshold was changed.
The lag-refinement bound (§4 step 2) was set to 3 s after the envelope lag
proved coarse at the few-second level; this affects timing only.

## 5. Shared assumptions and limits of independence

- **Vertical absolute response.** The gain is relative to the nominal vertical
  response: a vertical sensitivity error ε biases g by (1+ε). The Rayleigh and
  microseism checks share this, so **agreement with them is not independent
  proof** of absolute pressure calibration.
- **Shared water parameters.** ρ and c₀ enter both A and B; h enters A (and
  Rayleigh through ρH). Depth is the inventory elevation estimate, not verified
  bathymetry.
- **Different bands.** Rayleigh gains are 0.03–0.07 Hz; current conditional
  microseism gains use eligible bands ending at or below 0.10 Hz (the wider
  range is diagnostic). Local P is 0.5–16 Hz. Deng (§4) and An (§4.1) find frequency-
  dependent DPG coefficients, and responses differ above a few Hz (DPG roll-off,
  anti-alias FIR). Cross-band disagreement is not by itself a failure.
- **Coupling and tilt.** OBS–sediment coupling resonances can distort vertical
  motion above ~10 Hz. Tilt dominance is not established here; it must not be
  declared negligible for every instrument/band from coherence alone.
- **Plane-wave layered model.** Lateral heterogeneity, rough bathymetry over
  the Fresnel zone and scattered P coda with other q are not modelled.
- **Timing.** A separate pressure logger (XO APGs) can carry clock offsets.
  Lag estimation removes timing but also removes phase as an independent check.

## 6. Auger (2024) reading notes

- Events from Mayotte (1T, 4 km distance) and Lucky Strike (4G, 40 km) with
  hydrophones; teleseisms and low-frequency events at MOUM. These networks are
  not in the local archive.
- Time-domain peak-to-peak first-arrival amplitudes, picked manually in
  Snuffler, converted with **scalar sensitivity only** (no full response),
  acceleration = 2πf·velocity at one frequency.
- Compared with Eq. B (first transmission), with the full tan form and
  fitted incidence angles. An unexplained π phase difference remains.
- Takeaway: the report supports the expectation that first arrivals follow
  Eq. B, but it is not a calibration and it fits a physically excluded parameter.
  It is a lead, as `project_context.md` already says.

## 7. Validation design

- Closed-form limits: F → 1 + θ²/3; poles; q conversions; consistency with
  `rayleigh_calibration.water_ratio` (shared physics, not independent).
- **Independent reference:** a 1-D staggered-grid finite-difference acoustic
  simulation of a water layer over a half-space. A plane wave with slowness q
  reduces exactly to 1-D with per-layer vertical modulus `K/η²`. It never calls
  the closed-form code. Long windows must reproduce Eq. A away from poles,
  gated windows must reproduce Eq. B, and changing the half-space impedance
  (R) must not change either ratio.
- Gain recovery and polarity on synthetic data processed through the full
  estimator, including a known scale factor and sign.

## 8. Open scientific questions

1. Is the pre-multiple relation (Eq. B), strictly time-gated, acceptable as a
   *calibration* reference for deep stations? It is kept diagnostic here.
2. What is the DPG/hydrophone response accuracy above ~2 Hz, where the nominal
   responses roll off (e.g. RR52 BDH 0.84 at 5 Hz)?
3. How should frequency-dependent gains be combined with Rayleigh-band gains
   for the final score?
4. XO APG clock offsets: should lag estimation be accepted, given that it
   removes the phase check?
