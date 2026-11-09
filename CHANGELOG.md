# Changelog

## Unreleased (branch `compy3-audit`)

Audit for ComPy 3: [docs/audit_2026-11/REPORT.md](docs/audit_2026-11/REPORT.md).
**Published results are unchanged:** the default estimator is bit-identical to 2.8.0, and
all 61 archived compliance files are recomputed from their stacks to 2e-16.

### Added

- `compy.anomaly`: automatic detection and ranking of anomalies and correlations:
  - changepoints, tidal and solar periodicity of series and of accepted windows,
    pressure-level dependence and regimes;
  - network common mode, synchronous steps and PSD-level outliers;
  - epoch mixing and inversion gain pull;
  - instrument-response models fitted to the residual Z/P phase, with the compliance bias
    they imply;
  - a Monte Carlo of the coherence-gate selection bias.
- `compy.timefreq`: per-window spectra, coherence, admittance and tilt direction from
  plain arrays.
- `compy.results`: a loader for archived stage C and E products.
- `compy.noise_models`: Peterson NLNM/NHNM and Webb pressure bounds, without ObsPy.
- `compy.plotting.diagnostics`: 13 diagnostic figures (coherograms, error budget, monthly
  stability, phase diagnostic, tidal maps, common mode, anomaly ranking and more).
- `stacked_compliance(..., gate="exclusive")`: a per-bin coherence gate without the
  upward selection bias of the default gate. It is opt-in, and `pipeline/c_compliance.py`
  reads `COMPLIANCE_GATE` from the run configuration.
- `stacked_compliance_from_spectra`: the same estimator on precomputed spectra.
- Example scripts:
  - `examples/compy3_diagnostics.py`, for archived results;
  - `examples/compy3_waveform_diagnostics.py`, for cleaned MiniSEED.

### Pressure-gauge audit (branch `pressure-audit`)

[pressure/docs/audit_2026-11/REPORT.md](pressure/docs/audit_2026-11/REPORT.md). No gain,
response, score or class changes.

- `compy_pressure.water_layer` (ObsPy-free):
  - the water-layer P/a ratio and its analytic sensitivities to sound speed, phase
    velocity and depth;
  - the calibration bias when the water column or crust is wrong;
  - a log-ratio multiplier fit;
  - the depth-mean in-situ density.
- `compy_pressure.table_audit`: anomaly scan of the QC/calibration master table:
  - gain outliers and method disagreement;
  - polarity;
  - depth confounding;
  - Godin outliers;
  - channels hidden from Test 1;
  - invalid peer offsets.
- `compy_pressure.audit_plots`: five figures.
- `pressure/scripts/pressure_audit_2026_11.py`, with a copy of the BRUIT-FM v1.2.0 master
  table (local paths removed).

### Whole-archive audit (branch `compy3-audit-data`)

[docs/audit_2026-11-data/REPORT.md](docs/audit_2026-11-data/REPORT.md): the tests of the
two audits above, run on the whole-archive data (287 inverted stations, 17 networks).
**No default changes:** every archived stage C result (657 files) is reproduced
bit-for-bit through `stacked_compliance_from_spectra`.

- `compy.response_check` (NumPy only, tested):
  - fits the pressure-response corner factor to the residual Z/P phase;
  - predicts the tidal amplitude ratio and the compliance bias it implies;
  - fits tides in hourly counts, with gaps and a spline drift;
  - derives the pressure sign from the compliance phase.
- `docs/audit_2026-11-data/scripts/` (s01-s07): stage C restack (inclusive and exclusive
  gates), phase vs instrument metadata, DPG tide test, which sensor stepped, anomaly scan
  with stage D cross-checks, tilt vs tides, figures.
- `examples/compy3_waveform_diagnostics.py`: the docstring says to read tilt from stage A
  data (on cleaned data the azimuth splits into two bands).

### Fixed

- `stacked_compliance` no longer mixes windows with different frequency axes silently,
  and counts rejection reasons of custom criteria instead of raising `KeyError`.
- `compy.streams` imports ObsPy lazily, so `compy.compliance` works without ObsPy.

## 2.8.0 — 2026-11-05

**Stage E results change:** every variant of `pipeline/e_invert.py` now requires Vs to
increase across the Moho. Without this rule, the last crustal layer was faster than the
mantle in 7 to 81 % of the draws at eight YV stations. The library default is unchanged,
so published scripts repeat.

### Added

- `compy.inversion.literature_mcmc`: Monte Carlo inversion with the parameterizations of
  Zha & Webb (2016) (`ZhaWebb2016Model`) and Zheng et al. (2024) (`Zheng2024Model`), run by
  `sample_literature_model`.
- `compy.inversion.minimum_structure`: `invert_minimum_structure`, the minimum-structure
  inversion of Crawford (2004), which smooths Vs itself, with `crawford_thin_layers`.
- `compy.inversion.impedance`: a forward model by surface-impedance recursion, written
  independently of `propagator` and tested against the exact half-space and a
  matrix-exponential solution.
- `pipeline/e_alternative_inversions.py` runs these on one station's stage C result.
- Stage E variants `final_brocher` (Vp from Vs by Brocher 2005) and `final_crust1` (Vp and
  density fixed from CRUST1.0), to compare with the free Vp/Vs of `final`.
- `CrustalModel(require_moho_jump=True)` and `invert_occam(moho_jump=True)`, both off by
  default.
- `docs/ALTERNATIVE_INVERSIONS.md`, with these methods compared on eight YV stations.

### Fixed

Fixes from running every OBS network of the archive. Results of stations that
already worked are unchanged, sample for sample.

- Stage A: a piece ending before the next 0.8 s grid point no longer fails its
  chunk ("startime is larger than endtime"); invalid day files such as day 366
  of a non-leap year are skipped; sampling rates within 1e-5 of a nominal value
  are snapped to it so segments can merge; channels exceeding 2^23 counts in over
  half of their blocks (32-bit digitizers, absolute gauges) use a 32-bit clip
  level; data outside every StationXML epoch loses only its segment, listed under
  `no_response` in the QC file.
- Stage B: non-finite samples become gaps; gaps left between cleaning windows
  can be written to MiniSEED.
- Version-tagged names replaced, old names still work with a DeprecationWarning:
  `invert_compliance_beta` -> `invert_compliance_broad_bounds`,
  `plot_inversion_v2` -> `plot_inversion_summary`,
  `plot_inversion_beta` -> `plot_inversion_with_start`,
  `plot_inversion_serpentinization1` -> `plot_inversion_serpentinization_panels`.
- Unused variables and imports removed; mutable default argument removed.
- Every public function and class has a docstring. `compy_pressure` checkpoints that compare
  the package's code checksum will see new code once.

## 2.7.0 — 2026-11-02

A correctness pass before full-data reprocessing: a code review of the
package, a test of every real-data example, and a pipeline run from raw data.
Results of the existing workflows are unchanged.

### Fixed

- **`correct_tilt`, `remove_coherent_noise` and `estimate_tilt_angles`
  failed on float32 traces.** TiSKitPy needs complex128 cross-spectra. Float32
  traces are what stage A writes, and what FLOAT32 MiniSEED holds. Inputs are now
  copied to float64; float64 input gives bit-identical results.
- **`correct_tilt` per-hour cleaning (`clean_horizontals=True`, the default)
  failed at 1.25 Hz:** TiSKitPy's default spectral segment does not fit a 1-h
  window. When it does not fit, the segment falls back to one sixth of the
  window (with a warning), and `spectral_window_s=` sets it explicitly.
- **Pipeline stage F:**
  - it crashed on a reference station without stage B/C output (`fig04_tilt`);
  - it crashed on a short record without a monthly stack (`fig06`, IndexError) or without clean windows (NameError);
  - Fig. 3 now states the gain and its source, and says "uncalibrated" for the nominal response.
- **One rule for the gain source** (`e_invert.result_label`: Rayleigh,
  BRUIT-FM, nominal) is now shared by stages E, F and G, which each had a copy.
- `plot_coherogram_daily` took every channel's sampling rate from `st[3]`.
- `compy.compliance` no longer imports Matplotlib or the FDSN client at load
  time. They are imported only where used.
- Stale provenance string `compy.Rotate` in `compy_pressure/tilt.py` summaries.

### Added

- `correct_tilt(..., skip_failed=True)` keeps windows that TiSKitPy cannot
  rotate (e.g. tilt above its 20° limit, as on XE CC11 before it clipped). They
  stay uncorrected, with NaN angles and a warning, instead of aborting the
  record. The default still raises.
- Example 1 runs unattended with the method of stage B (periodic template with
  a fold-strength test, outside earthquake spans) on every channel.
  `--interactive` keeps TiSKitPy's PeriodicTransient, which prompts for clip
  levels.
- Example 2 follows the paper: hourly rotation without per-hour cleaning, then
  `remove_coherent_noise` in windows of up to 24 h.

### Verified

- 830 tests pass on Python 3.10 and 3.13; ruff 0.16.9 clean.
- Examples run with a fresh wheel install:
  - synthetic, reproducible and Vs inversions;
  - pressure QC;
  - FDSN download (RESIF YV.RR52);
  - transients (RR52);
  - rotation (XE CC11);
  - Vs inversion of CC11;
  - calibration (CC11; correctly unresolved with one usable event).
- `run_pipeline.py --fetch` was run from raw SDS (CC11) through stages 0–G.

## 2.6.0 — 2026-11-02

One `compy` package with consistent names, inversions started from CRUST1.0,
a reorganised repository, and an end-to-end check of the pipeline from raw
data. Every ComPy 2.x import keeps working, with a deprecation warning, and
numerical results are unchanged.

### Verified for this release

- The PR #16 restructure gives bit-identical results to 2.5 on XE CC11. Compared were:
  - stacked compliance and its errors;
  - the data covariance and forward model;
  - a seeded DRAM run;
  - the legacy solver;
  - serpentinization.
- 824 tests pass on Python 3.10 and 3.13.
- `pipeline/run_pipeline.py --fetch` was run from one month of raw CC11 SDS data
  through stages 0, A, B, D, C, E, F and G in a new run folder. The figure
  stage was also run on the full XE results.

### Fixed

- **Pipeline stage F** crashed when the reference station had no Rayleigh-wave
  calibration (`KeyError: 'compy'`), for example when its record contains no
  M ≥ 7 event. It now uses the same gain as stages C and E (`compy`, then
  `bruitfm`, then `nominal`; `f_figures.result_label`). It skips the window
  figures when the station has no compliance at all. Regression test:
  `tests/test_pipeline_figures.py`.

### Documentation

- The README no longer opens with the all-station Vs heat map.


### One `compy` package with consistent names

- The nine top-level modules are now one package, one module per processing
  step: `compy.streams`, `compy.cleaning`, `compy.preprocessing`,
  `compy.compliance`, `compy.inversion` (`dram`, `metropolis`, `occam`,
  `propagator`, `models`), `compy.diagnostics`, `compy.plotting` (`style`,
  `spectra`, `inversion`, `qc`) and `compy.numerics`. `compy_pressure` is
  unchanged.
- Misspelled, non-PEP 8 and misleading names are renamed, e.g.
  `Calculate_Compliance` → `compliance.calculate_compliance`,
  `Rotate` → `preprocessing.correct_tilt`, `calc_norm_compliance` →
  `inversion.propagator.normalized_compliance`, `liklihood` →
  `inversion.metropolis.likelihood`, `final_plot` →
  `plotting.inversion.plot_publication_figure`, `ffplot.psd` →
  `plotting.spectra.plot_psd`. Names that come from Crawford's MATLAB code
  (`raydep`, `argdtray`, `gravd`, `dtanh`) keep their names. Full table in
  `docs/USER_GUIDE.md`.
- Duplicates are merged: one `rms`, one `compliance_uncertainty`, one plot
  style (`set_plot_style`, `set_plot_style_large`, which set the same
  `rcParams` as the two former `plt_params`), and the `ffplot` wrappers of
  inversion functions now point to the originals.
- Commented-out code, empty `# %%` cells and loop-progress `print` calls are
  removed. No equation, default or order of operations changed: a seeded run
  of the Metropolis and DRAM samplers, the propagator and the kernels gives
  bit-identical arrays before and after.
- Compatibility: every ComPy 2.x import keeps working (`import inv_compy`,
  `from compy_processing import estimate`, `compy.Calculate_Compliance`, ...)
  and warns with the new location. `compy_numerics` and the other one-to-one
  modules are the new module objects, so private names and monkeypatching
  still work. The aliases will be removed in 3.0.
- `save_run` records the SHA-256 of every file in the package, keyed by path
  (`compy/diagnostics.py`, ...).

### Inversion starts from CRUST1.0 at the station

- `compy_inversion.crust1_start(model, column)` maps the CRUST1.0 cell of the
  station (`crust1_column`) onto a `CrustalModel`: sediment thickness,
  thickness-weighted Vs of each crustal layer, mantle Vs and, with free
  Vp/Vs, the crystalline-crust Vp/Vs, clipped to the prior.
- `sample(..., start=, start_spread=0.1)`: the first chain starts at that
  model and the others are dispersed around it by 10 % of each prior width
  (`dispersed_starts`), so R-hat still detects unmixed chains; sampled gains
  start at 1. The start is recorded in `settings["start"]`. Without `start`,
  chains still start from prior draws.
- Pipeline stage E and `examples/4_Vs_inversion.py --crust1` start every run,
  including the prior-only run, from CRUST1.0. The results in `docs/results/`
  were computed before this change, with prior-draw starts.
- Legacy sampler: `inv_compy.crust1_starting_model(column, sublayers=1)`
  builds the `[thickness, rho, Vp, Vs]` starting model from CRUST1.0, and the
  `invert_compliance*` functions accept a `crust1_column` as `starting_model`.
- `inv_compy.Model_V2` is deprecated. Its hard-coded templates differ from the
  CRUST1.0 cell of their station:

  | Station | Template versus the station's CRUST1.0 cell |
  |---|---|
  | RR28 | crust of RR38 (700/1540/4750 m) instead of 660/1460/4510 m; mantle Vs 4190 instead of 4510 m/s |
  | RR29 | sediment 10 m instead of 80 m; crust of RR38; mantle Vs 4190 instead of 4510 m/s |
  | RR34 | sediment 10 m instead of 410 m; crust of RR38; mantle Vs 4190 instead of 4500 m/s |
  | RR36 | crust of RR38 (700/1540/4750 m) instead of 690/1530/4710 m |
  | RR38 | half-space Vs 4370 m/s (RR36's mantle) instead of 4190 m/s |
  | RR40 | crust of RR38 instead of 680/1510/4640 m; mantle Vs 4190 instead of 4360 m/s |
  | RR50 | crust of RR38 (Vs 2700/3700/4050 m/s) instead of 1180/1180/4570 m with Vs 2540/3590/3940 m/s; half-space Vs 4370 instead of 4190 m/s |
  | RR52 | crust matches; half-space Vs 4370 instead of 4310 m/s |

  In the legacy sampler the starting model also sets the Vs prior, so
  inversions run with these templates used these priors.

### Repository organisation

No code or numerical changes. Update local paths as follows:

| Before | After |
|---|---|
| `_Example/` | `examples/` (new index: `examples/README.md`) |
| `_Images/` | `docs/images/` |
| `docs/yv_compliance_2026-10-28/` | `docs/results/yv_compliance_2026-10-28/` |
| `docs/yv_inversion_2026-10-29/` | `docs/results/yv_inversion_2026-10-29/` |
| `docs/xe_young_orca_2026-10-29/` | `docs/results/xe_young_orca_2026-10-29/` |

- README rewritten as a step-by-step workflow that follows station XE.CC11;
  the API reference moved to `docs/USER_GUIDE.md`; new `docs/README.md` index.
- `environment.yml` for conda; keywords, classifiers and project URLs in
  `pyproject.toml`; DOI, ORCID, abstract and keywords in `CITATION.cff`.
- GitHub issue and pull-request templates; Dependabot for workflow actions.
- Ruff pinned to 0.16.9: newer releases also format Python code blocks in
  Markdown, and the lint result should not change with each ruff release.
- The session handoff note `docs/HANDOFF_2026-10-29.md` is removed.

## 2.5.0 — 2026-10-30

ComPy becomes the complete package for seafloor pressure gauges: all of
BRUIT-FM moves in, and its calibration replaces `Pressure_calibration.py`.
This release stays in the 2.x series at the owner's request. It includes the
breaking calibration API change below; existing scripts must use the migration
guide or pin ComPy 2.4.0.

### Breaking changes

- **`Pressure_calibration.py` is removed.** Use the Rayleigh-wave calibration
  of `compy_pressure` instead (README, "Pressure calibration", has a
  migration table). Scripts that import the old module need ComPy 2.4.0
  (commit `7ad79c2`), the last version that has it:
  - `docs/yv_compliance_2026-10-28/scripts/d_calibrate.py`;
  - `docs/xe_young_orca_2026-10-29/scripts/d_calibrate.py`.
- **Pipeline stage D (`pipeline/d_calibrate.py`) now uses `compy_pressure`.**
  New runs can give gains that differ from 2.4 runs: on XE CC11, WC03 and
  EE04 and YV RR52 they are 1.4–1.9 % higher (`docs/PRESSURE.md`). Its outputs keep their
  names and keys (`compy_gain.json`, `<STA>_compy_calibration.npz` and
  `.png`), so stages C and G are unchanged. Stage F's calibration figure no
  longer assumes M ≥ 7 events or a 0 Hz first bin. New `config.py` settings:
  `CALIBRATION_BAND`, `CALIBRATION_TILT` and `CALIBRATION_DISPERSION`.

### New: `compy_pressure` and `pressure/`, all of BRUIT-FM

BRUIT-FM 1.3.0 (its last release; the repository is kept as a backup with
its dated reports), moved by its `integration/export_to_compy.py`:

- `compy_pressure/`, the package, holds all 40 modules and their reference
  data:
  - inventory and master table;
  - response-corrected PPSDs and plots;
  - QC tests against the Webb/Brown limits, nearby gauges and the Godin
    et al. (2013) infragravity model;
  - Rayleigh-wave, microseism and local-P calibrations;
  - tilt-corrected vertical PPSDs;
  - percentile noise models;
  - maps and network figures;
  - whole-drive coverage runs and quality review.
- `pressure/` is BRUIT-FM's project folder:
  - `tests/` (503 tests);
  - `scripts/` (59 run and review scripts);
  - `docs/` (methods, validation notes, reviews and policies; the dated
    reports are linked);
  - BRUIT-FM's README and changelog.

  It is laid out like BRUIT-FM's checkout, so its relative paths still hold.
- Commands: `compy-pressure` (stages `inventory` … `rayleigh`, `all`, `view`),
  `compy-pressure-inventory`, `-tilt`, `-local-p`, `-microseism`, `-paired`,
  `-noise-model` (group noise model) and `-approved-noise-model`.
- New extras: `pressure` (pandas, openpyxl) and `maps` (Cartopy). `dev` now
  includes pandas and openpyxl.
- `compy_pressure.rayleigh_network.calibration_master`: a calibration without
  PPSDs, with each channel's StationXML as its response snapshot. It is used
  by example 3 and stage D.
- The Rayleigh network lock and the master-table lock also work on Windows.
- `_Example/3_DPG_Calibration.py` is rewritten for `compy_pressure`.
- New `_Example/5_Pressure_QC.py`: inventory, PPSDs, plots and QC tests 1 and
  3 on a synthetic gauge (offline, seconds) or on real data.
- Guide: `docs/PRESSURE.md`.

### Fixed (in `compy_pressure`, since the move)

- Calibration masters validate the station, channel selection and complete
  StationXML source set before reusing a saved table.
- `pipeline/run_pipeline.py` runs A → B → D → C → E → figures/report in one
  command, saves a run manifest and stops on failed stages. Stage A now exits
  unsuccessfully if any preprocessing chunk fails.
- Source distributions now include the pipeline and BRUIT-FM's scripts,
  methods and tests, so the complete workflow is available from the release
  source archive as well as a Git checkout.

- The tilt-cleaned Rayleigh calibration read the channel responses from the
  PPSD folder (`<PPSD folder>/inventory.xml`) instead of the response
  snapshot. It therefore failed for a calibration without PPSDs (explicit
  snapshots, as with `calibration_master` and stage D) once tilt cleaning was
  on. This was inherited from BRUIT-FM 1.3.0; it is now fixed, with a test.

### Tests

- `pressure/tests/`: 503 tests on synthetic StationXML and waveforms. They
  pass on Python 3.10 and 3.13; pytest finds them through `testpaths` and
  `pythonpath` in `pyproject.toml`.
- Removed with the old module:
  - `tests/test_calibration.py`;
  - the Mackenzie sound-speed and plotting-helper tests.
- The gain-without-quantization test now checks `compy_pressure`'s fit.
- New `tests/test_pressure_calibration.py` covers:
  - `calibration_master` (snapshots, channel filter, reuse, changed
    StationXML);
  - the band check;
  - tilt training from explicit snapshots;
  - stage D's catalogue conversion and figure arrays.
- The wheel check also runs `compy_pressure` outside the checkout: packaged
  reference data, console scripts and `5_Pressure_QC.py`.
- Reproduction check against BRUIT-FM (90 gain-corrected epochs, the
  gain-corrected noise model): every output file is identical. Details are in
  BRUIT-FM's `integration/README.md`.

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
