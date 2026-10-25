# Changelog

## Unreleased

- Sampler tuning (opt-in; defaults reproduce 2.1 chains draw for draw):
  per-layer `sigma_v`/`sigma_h`, Robbins-Monro step-size adaptation frozen after
  `adapt_until`, `initial_state` for over-dispersed chain starts that keep the
  prior of `starting_model`, and a `diagnostics` dict with tuned steps and
  per-move acceptance.
- Add a converged four-chain synthetic example and a Vs-depth posterior heat map.
- Explicit uniform Vs prior (`vs_bounds`) independent of the starting model, and
  `vary_thickness=False` for the fixed-geometry minimum-structure parameterization
  of Crawford et al. (1991).
- `invert_occam`: Occam's inversion (Constable et al. 1987) in ln Vs / ln compliance
  with a damped, monotone second phase; an independent cross-check of the sampler.
- `compliance_kernels` and `resample_model`: numerical Vs sensitivity kernels on an
  interface-exact slicing of the model.
- Tests: the forward model is checked against an independent Aki & Richards
  propagator (agreement ~1e-6), plus kernels, Occam and the new sampler options.
- Deep-crust synthetic example (14 layers to 13 km, 3.5-18.5 mHz, 4200 m water)
  with prior-versus-posterior resolution analysis.
- `alpha_bounds`: the smoothing weight is Gibbs-sampled from its Gamma conditional
  (hierarchical Bayes) instead of being fixed.
- `data_covariance`: full data covariance (e.g. a pressure-gain error common to all
  frequency bins) through a Cholesky-whitened likelihood.
- `interface_move="adjacent"`: move one interface at a time in free-interface runs.
- `compy_diagnostics.rank_rhat_ess`: rank-normalized split R-hat and bulk ESS on
  NumPy/SciPy alone; the examples no longer need ArviZ for diagnostics.

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
