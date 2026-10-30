# Changelog

## 2.4.0 — 2026-10-30

Integration release: every open pull request and branch is merged, the
inversion is carried through to serpentinization, and an end-to-end pipeline
is validated on a second network.

### Merged

- PR #6 (2.2, compliance-chain audit and YV compliance) and PR #7 (2.3,
  Bayesian Vs inversion), including the `inversion-improvements` branch.
- PR #5 (W. Crawford, `Pressure_calibration.py` and README). It was written
  against the pre-2.0 code; its conventions are ported onto the 2.x function:
  - explicit `zchan`/`pchan` channel codes;
  - `filt_freq1`/`filt_freq2` pre-filter arguments checked against the band;
  - Google-style docstrings;
  - no network/station constants at module level;
  - disba/Numba installation notes and a corrected README example.
- `refactor/professionalize-compy` (PR #4, closed on 2026-10-07 as superseded)
  was re-checked and is not merged. Its packaging metadata is older than main
  (version 0.2.0, Python ≥ 3.9), and its tests target the legacy
  `Roughness`/`liklihood` behaviour found defective in the 2.3 audit (V1, V2).

### Added

- `compy_inversion.serpentinization` and
  `Posterior.serpentinization_quantiles`. The degree of serpentinization of a
  peridotite crust is computed from Vs, linear between 4.2 km/s (0 %) and
  2.2 km/s (100 %) as in Aminian et al. (2025, Fig. 10c; Miller & Christensen
  1997), per posterior draw so that its credible interval follows.
- `pipeline/`: config-driven end-to-end scripts.
  - Stage 0: inventory and catalogue.
  - Stage A: decimation, responses, and removal of clipped or dead spans.
  - Stage B: cleaning, with optional glitch removal.
  - Stage D: Rayleigh-wave DPG calibration.
  - Stage C: compliance.
  - Stage E: inversion.
  - Schedulers that chain the stages per station. `a_fast.py` stages a slow USB archive through local disk.
- `calculate_spectral_ratio(..., zchan=, pchan=, filt_freq1=, filt_freq2=)`.
- From the XE validation run (`docs/xe_young_orca_2026-10-29`):
  - `pipeline/f_figures.py`: the figures of Aminian et al. (2025) for any network;
  - `pipeline/g_report.py`: report numbers;
  - `pipeline/a_download.py`: FDSN staging;
  - inversion queues with lock files;
  - a nominal-gain fallback for stations without calibration, with a 0.6–1.4 gain prior;
  - the paper's 500-window minimum applied when the results are collected.

### Fixed

- **F2.** The legacy serpentinization plots used (9000 − 2 Vs)/30, i.e. 0 % at 4.5 km/s and 100 % at 3.0 km/s. They now use the published caption relation: Vs = 3.5 km/s gives 35 %, not 67 %.
- **F1.** The legacy heat-map depth ticks assumed 167-m rows for profiles averaged into 100-m rows, which labelled depths 5/3 too deep. The ticks now come from `LEGACY_HEATMAP_ROW_M`, and the CRUST1.0 reference arrays keep their own row depth.
- CI was failing on the 2.3 branch:
  - the YV report script used Python 3.12 f-string syntax;
  - files added since 2.2 were not `ruff format`-ed.
- `compy.__version__` was 2.2.0 while the package said 2.3.0.
- Docstrings:
  - `remove_coherent_noise`: the sequential, conditioned removal does not depend on the channel order;
  - `legacy_wave_attraction_correction`: it gives the 2.0 estimator, not 2.2.

## 2.3.0 — 2026-10-29

Shear-velocity inversion, audited and rebuilt on top of the 2.2 compliance
chain and the `inversion-improvements` branch. Details, audit table and
validation: `docs/INVERSION.md`.

### Audit of the earlier inversion code (documented, legacy kept)

- v1 roughness `sqrt(sum(np.gradient(vs, 2)))` equals √((V_bottom − V_top)/2)
  and applied no smoothing to the published models.
- The v1 likelihood used the residual norm, not its square.
- v1 coupled Vp (Poisson 0.25) and density (Gardner) to Vs, while the paper's text says
  they were fixed.
- v1 proposals near the bounds and thickness moves were asymmetric, and the thickness
  moves did not conserve total depth.
- Pre-2.0 compliance also used a free-air constant of `3.07 * 10e-6` (3.07e-5).
- In v2, the default `alpha=0.25` is a 4 m/s Vs-jump prior, far stiffer than the data. This is now
  documented in the README.

### New module `compy_inversion`

- `ComplianceData`: a full data covariance assembled from the random error, the covariance of
  monthly (or other) sub-record estimates, a pressure-gain error common to all
  frequencies with the exact gain sensitivity of the gravity-corrected
  estimator, and a relative floor.
- `CrustalModel`: the Zha & Webb (2016) parameterization:
  - a sediment of free thickness with Hamilton (1976, 1979) Vs and Vp profiles and an optional Vs scale;
  - crustal layers doubling from 100 m to a CRUST1.0 Moho, with an optional Moho shift;
  - a mantle half-space;
  - Vp/Vs 1.85 with Christensen & Shaw density, Brocher (2005), or Vp and density fixed from a reference;
  - a hierarchical second-difference roughness prior.
- `sample` / `dram_chain`: a delayed-rejection adaptive Metropolis sampler.
  - The covariance is learned from the recent half of the chain, with a global scale tuned to 23 % acceptance; both are frozen after burn-in.
  - Delayed rejection uses the Tierney & Mira acceptance.
  - The smoothing weight is Gibbs-sampled.
  - Parallel tempering is optional; chains run in parallel processes.
- `Posterior`: rank-normalized R-hat and bulk ESS, χ², Vs-depth quantiles, and save to npz/json.
- `vs_kernels`, `sensitivity_depth` and `slice_model`: Vs sensitivity at constant
  Vp/Vs and the depth holding a fraction of it.
- `crust1_column` / `crust1_model`: read CRUST1.0 at a station.
- `Forward`: the `raydep` propagator evaluated at the seafloor only, compiled with
  Numba when installed (optional extra `fast`), about 140 times faster.

### Gain, epochs, figures (same release)

- `sample(gain_bounds=...)`: pressure-gain factor(s) sampled with a uniform prior.
- `ComplianceData.join`: several epochs of one station with their own gains.
- `ComplianceData.from_result`: build directly from a `ComplianceResult`.
- `CrustalModel(vp_vs_range=...)`: one free crustal Vp/Vs.
- `Posterior.trace_chi2`: the misfit at every iteration, including burn-in.
- Documented two errors in the publication figure code:
  - heat-map depths are 5/3 too deep relative to their labels;
  - the serpentinization formula does not match the Fig. 10 caption.

### Tests

New tests: the forward model against `raydep` and the static half-space,
the gain sensitivity, the covariance, DRAM on a truncated Gaussian, the
hierarchical marginal, tempering on a bimodal target, kernels, CRUST1.0
parsing, the fixed-Vp option and synthetic recovery.

## Inversion work (branch inversion-improvements, 2026-10-25/26)

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
- `brocher_elastic` and `elastic_from_vs`: optionally tie Vp and density to Vs
  (Brocher 2005), as in Zheng et al. (2024); Zha & Webb (2016) used a regional relation.
- `legacy_wave_attraction_correction`: converts compliance estimated with ComPy before
  2.0 to the 2.0 estimator. The legacy code added the wave-attraction term
  (acceleration/pressure) to displacement/pressure before the free-air factor, so the
  term was effectively dropped; published values are low by up to ~25-30 % near 7 mHz
  at 2.9 km water depth (RR52) and ~11 % at 4.1 km (RR50), decreasing with frequency.
- `invert_tempered`: parallel tempering (replica exchange) around the sampler, with
  independent ensembles for R-hat. For posteriors where single chains do not mix between
  equally good shallow layerings when the top kilometre is below resolution.

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
