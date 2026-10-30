# Calibration assumptions and checks to finish

The Rayleigh implementation is treated as a **working computational method**.
Its station gains are provisional calibration candidates. The network run does
not alter responses, waveforms, existing spectra or the six scientific scores.
Every result records `implementation_status: working_method_with_provisional_assumptions`
and machine-readable `review_flags`; the master links these under `rayleigh_*`.
An unresolved comparison is not an instrument failure.

The method follows the thesis/manuscript pressure–acceleration fit and the
published Zha/Webb water-column reference. BRUIT-FM recomputes from raw counts
using its own archive/spectral processing and actual TiSKitPy CleanRotator.
The historical Earth-model layers were reused from the earlier DPG code.
Agreement with 0.67 is not independent certification. Synthetic tests inject
known gains deliberately; production fitting has no hard-coded 0.67 target.

## Recorded flags and remaining work

| Flag | Check, test or improvement needed |
|---|---|
| `pressure_polarity_assumed` | Verify pressure/vertical hardware and response sign. The YV run uses one explicit −1 assumption, not per-event sign selection. |
| `vertical_absolute_response_unverified` | Independently verify vertical seismometer absolute sensitivity; its error propagates into pressure gain. |
| `depth_from_inventory` | Verify bathymetry/deployment depth. Inventory elevation is an estimate; compare depth sensitivity. |
| `earth_model_needs_station_validation` | Check the shared RR52 historical phase/group model against station/path structure or alternative models. The model is not independently validated for all YV stations. |
| `water_properties_assumed` | Test homogeneous density/sound speed against appropriate water-column information. |
| `sensor_axes_assumed_orthogonal` | Verify horizontal 1/2 geometry and metadata azimuths. Horizontal codes follow the paired vertical (BHZ → BH1/BH2, LHZ → LH1/LH2); general N/E pairing is not implemented. |
| `quiet_screen_and_coda_heuristic` | Inspect quiet training for uncatalogued transients and validate the ad hoc coda exclusion formula. |
| `tilt_stationarity_needs_testing` | Vary preceding training intervals; check that the rotation preserves Rayleigh motion and remains appropriate during events. |
| `rayleigh_window_mode_needs_review` | Review broad arrivals/overlaps and verify isolated fundamental-mode packets, rather than selecting the largest peak or maximizing coherence. |
| `qc_rules_need_field_validation` | Validate window/segment lengths, coherence/phase criteria, 24 h separation, depth limit and event-count rule. None is an instrument pass/fail standard. |
| `magnitude_fallback_needs_testing` | Review smaller-event signal quality; magnitude alone does not establish suitable Rayleigh motion. |
| `scalar_gain_frequency_independence` | Check frequency-dependent residuals. A scalar fit does not determine the DPG capillary time constant or justify extension outside the fit band. |
| `uncertainty_incomplete` | Include response, depth/model, water, event/window/estimator and selection uncertainties; event scatter is not total uncertainty. |
| `pressure_hardware_needs_confirmation` | Resolve generic or hydrophone labels before applying DPG-specific interpretation. A metadata label is not a hardware inspection. |
| `clock_alignment_needs_review` | Investigate reported sample-start offsets; justify any correction explicitly rather than rounding timestamps or interpolating gaps. |
| `local_tilt_method_needs_to_be_tested` | Only when explicitly selecting the retained experimental eigenvector method: comparative field validation is still pending. |

Flags remain present even when automatic event checks pass. They are review
items, not scores. The current software does not automatically mark them resolved.
When `--tilt-method none` is explicitly selected, tilt-training, stationarity and
horizontal-axis flags are omitted because no rotation is attempted. Pressure/Z
response, physical model, event-window and QC flags remain. See
[optional tilt correction](tilt_methods.md#optional-no-tilt-correction).

## Parameters: adjustable versus currently fixed

Network CLI defaults use the original figure's **0.03–0.07 Hz** band, pressure
sign −1 and **all eligible events**. Three remains the existing minimum for a
station estimate, not a sampling cap. The historical three-candidate target and
eight attempts per magnitude level remain available through explicit `bounded`
sampling and describe the completed saved YV/YL pilots. The standalone
runner retains its earlier 0.02–0.06 Hz default. Both use the same method.

Standalone adjustable inputs include band, density (1028 kg/m³), sound speed
(1500 m/s), depth, squared-coherence threshold (0.95), pressure sign, event
target/budget, tilt method, catalogue, dispersion table and reviewed UTC windows.
The network API also accepts these scientific options; its CLI exposes the most
common band/sign/count/catalogue/model controls.

Currently fixed pilot choices are source depth <150 km; a 24 h separation;
30-minute windows centered on the predicted 40 mHz group arrival; 200 s Hann
segments with at least eight; phase ±20°; two-hour quiet training searched
6/12/24/48 h before origin; and coda exclusion
`max(3600, (M-5.5)*43200)` seconds. These are documented heuristics, not hidden
calibration answers. Their values are now saved in each station result.
Gain-search bounds 0.1–3.0 with step 0.001 are checked against the analytical
least-squares solution; a boundary minimum stays unresolved.

The actual TiSKitPy constructor/training filters and response prefilter are
recorded in per-event tilt metadata; see [tilt details](tilt_methods.md).
The resonance exclusion margin and display axes are separate from instrument
acceptance. An unsuitable water-reference band is flagged without silently
changing the selected band to produce a gain.

## Checks already completed, with their limits

- [x] Compare the water formula with an independently solved acoustic ODE.
  This checks the equation numerically under its boundary/model assumptions.
- [x] Cross-check grid and analytical gains, reciprocal sensitivity and boundary
  handling; recompute real event gains directly from saved spectra.
- [x] Recover constructed gains, rotations and sign reversals from synthetic raw
  records. The synthetic generator shares the water formula; this is not a
  fully independent Earth/ocean-model test.
- [x] Compare BRUIT-FM's projection with actual TiSKitPy `apply()` and preserve
  its input traces. This verifies integration, not complete field-noise removal.
- [x] Rederive ratios, transfers, coherence, phase and masks from all 481 saved
  network event spectra; independently solve all 48 candidate gains and 16
  station/median-curve gains. Matches to roundoff, but reuses saved spectral
  inputs and model assumptions. [Initial review](yv_rayleigh_saved_results_review.md).
- [x] Validate epoch binding, gap/alignment guards, restricted master writes,
  extra-sheet preservation, network restart and source-change invalidation.
- [x] Visually inspect RR52 example plots against manuscript Figure 5 / thesis
  Figure 3.13. Network-wide physical packet validation remains pending.
- [x] Cross-check seven representative raw-count event estimates using separate
  count-domain SciPy CSD and direct reads. Matches saved spectra/gains to roundoff;
  the same response metadata/evalresp and saved projection are still shared.
  [Independent validation](yv_rayleigh_independent_validation.md).
- [ ] Independently verify absolute responses and physical forward conventions,
  including the positive-real Laplace pressure poles at 34 stations. RESIF repeats
  the questioned RR03/RR22 metadata; this does not establish physical correctness.
- [ ] Inspect each candidate station's waveforms, phases, quiet training and
  reference applicability before applying a response correction.
  All 48 broad arrival plots and 48 quiet intervals have now been screened;
  20 windows contain predicted P/S arrivals, 46 quiet intervals pass independent
  continuity, and two have fractional-grid boundaries. Pure-mode selection and
  noise stationarity are not yet established. Future Rayleigh reads reject those
  fractional boundaries before merging; original outputs are retained.
- [ ] Review the existing evalresp FIR-normalization warnings against the source
  StationXML and hardware response. Snapshot agreement alone does not validate
  the underlying metadata. Warnings are retained, not suppressed as a fix.
- [x] Remove the default first-three-success sampling cap; the completed all-YL
  run exhausted eligible events under existing selection/QC safeguards.
- [ ] Test held-out events; assess sensitivity to windows, bands, models and tilt
  training. More attempted events alone do not validate the selected windows.
- [ ] Produce an uncertainty budget and agree scientific acceptance/classification
  rules before assigning the six test scores or combining accepted noise spectra.

The earlier RR52 depth sensitivity example moved the current paper-band median
from 0.679 to about 0.636/0.724 when assumed depth changed by −/+5%. Those are
illustrative perturbations, not measured depth uncertainty. Small event scatter
cannot establish corresponding absolute accuracy.

[Network commands and outcomes](yv_rayleigh_calibration.md) ·
[Physical method and source review](rayleigh_calibration_review.md).

[Repository and saved-result correctness audit](project_audit.md): 130 tests pass;
primary gain bounds and scientific acceptance rules remain unchanged.
