# Changelog

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
