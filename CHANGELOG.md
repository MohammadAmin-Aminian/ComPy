# Changelog

## 2.2.0 — 2026-10-28

Audit of the chain up to the compliance function. The full list is in
`docs/AUDIT_2026-10.md`.

### Numerical and scientific corrections

- The wave-attraction (gravity) term is no longer multiplied by coherence in
  `estimate`. Before this fix, compliance was biased low by (1-γ)·a_w/p, most
  strongly at the lowest frequencies.
- The theoretical (Bendat & Piersol) uncertainty now uses all selected
  windows. The number of averages counts distinct samples only, so overlapping
  windows are not double counted, and includes the Hann 50 % overlap factor.
  It previously used only the first window.
- `remove_coherent_noise` cleans the vertical alone, in 12–24 h windows,
  with transfer functions estimated outside event spans. The legacy hourly
  cleaning inside `Rotate` also altered the pressure channel.
  `Rotate(clean_horizontals=False)` now disables it.
- Example 0 removes responses with `water_level=None` and an explicit
  pre-filter. The default 60 dB water level gave errors of 22 dB in the
  compliance band when applied at the original sampling rate.
- Physical constants are now shared (`GRAVITY`, `GRAVITATIONAL_CONSTANT`,
  `FREE_AIR_GRADIENT`), and the gravity ratio uses the same g as the
  dispersion solver.

### New

- `stacked_compliance` stacks cross-spectra over windows and returns the
  compliance, its phase, the Bendat & Piersol error, and per-window
  statistics.
- `SelectionCriteria` with the `GJI2025` preset applies the published gates:
  coherence, pressure PSD 10–50 dB and acceleration PSD -170 to -150 dB.
- `coherent_band` gives the band where the mean coherence reaches the
  threshold.
- New module `compy_cleaning`:
  - unattended detection and subtraction of periodic transients;
  - smoothed per-sample tilt rotation using TiSKitPy angles;
  - STA/LTA spans for local events.
- `calculate_spectral_ratio(return_details=True)` returns per-event gains and
  their spread. The water density is a named parameter.

### Tests

- The always-passing gravity-ratio check (absolute tolerance 1e-8 on
  4e-11 values) is replaced by strict tests of the formula and its sign, and
  of the bias of the legacy and stacked estimators at low coherence.
- New tests check:
  - the stacked error against window scatter;
  - phase, gain scaling and gate reasons;
  - coherent-noise removal;
  - the rotation coefficients against TiSKitPy;
  - glitch recovery and STA/LTA detection.

### Migration

Legacy function signatures and default gates are unchanged. With the gravity
term no longer multiplied by coherence, and with the error computed over all
windows, the legacy outputs change slightly.

## 2.1.0

- Add optional ArviZ multiple-chain diagnostics: rank-normalized folded split-Rhat,
  bulk/tail ESS, retained-sample quantiles and explicit warnings.
- Save portable numeric run records with settings, seeds, dependency versions,
  source hashes and artifact checksums; refuse existing destinations.
- Add an offline four-chain noisy synthetic inversion example with a biased
  starting model and diagnostic plots. Existing sampler return tuples are unchanged.


## 2.0.1

Fixed short-record date labels and pressure PSD dB scaling; the download example now updates and attaches decimation response stages.

## 2.0.0 — 2026-10-07

### Numerical and scientific corrections

- Solve finite-depth gravity-wave dispersion without assuming the first bin is
  DC, and without unbounded secant iterations; share one validated solver.
- Stabilize the wave-attraction correction and convert acceleration/pressure
  into displacement/pressure before combining terms.
- Use consistent mean-Welch auto/cross spectra and amplitude coherence in the
  compliance estimators; retain explicit squared coherence for calibration gates.
- Unify amplitude-coherence uncertainty and handle zero coherence explicitly.
- Correct signs in the Mackenzie seawater sound-speed polynomial.
- Fit pressure gain analytically; preserve the sign and units of the gain.
- Honor the supplied calibration dispersion model and handle the removable
  zero vertical-wavenumber singularity.
- Rebuild event processing to include the final/single event, avoid hardcoded
  2501-sample requirements, and reject empty/invalid calibration sets.
- Initialize the inversion likelihood/prediction at iteration zero; honor custom
  models, require uncertainties and use log-domain Metropolis acceptance.
- Keep Vp/density fixed, evaluate roughness on a fixed metre grid and preserve
  all saved outputs after rejected proposals. Hold the half-space fixed and
  conserve total finite-layer thickness.
- Reject out-of-bound proposals rather than redraw or reflect them; use a local
  seedable RNG and count acceptance over actual proposals.
- Correct generic model construction and prevent hidden layer-count fallbacks.

### Reliability and maintenance

- Use copied sample-aligned windows with no duplicated boundaries, no skipped
  subsecond samples and no infinite overlap loop. Validate sampling grids/gaps.
- Raise contextual rotation failures; remove outlier edits that changed returned
  angles independently of the actual rotation. Report after/before variance.
- Update TiSKitPy imports, constructor arguments and time-bound tuples.
- Remove import-time network access and NumPy-incompatible warning suppression.
- Remove machine-specific automatic plot exports; make bathymetry inputs explicit.
- Restore the overwritten daily coherogram under a distinct name.
- Add installable packaging, version metadata, regression tests, formatting,
  CI, a standalone synthetic example, contributor guidance and scientific notes.

### Migration

Expect changed numerical results. V2 uses mean rather than mixed mean/median
spectral estimates, revised gravity terms and a corrected sampler. Reassess
selection thresholds, priors and regularization. Existing station results
have not been rerun. Public legacy imports and inversion return-tuple lengths
are preserved, but uncertainties are mandatory, `final_plot` requires an
image directory, and empty selection now raises instead of returning NaNs.
