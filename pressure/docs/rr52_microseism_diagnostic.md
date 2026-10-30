# RR52 microseism pressure–vertical diagnostic

## Outcome

The saved RR52 pressure/upward-displacement spectra support a reproducible
**admittance diagnostic**, not an independent pressure calibration. No pressure
gain is identifiable from these P/Z spectra alone. The observable combines the
pressure scale, vertical scale, and a frequency-dependent coupled
water–sediment–crust response. RR52's Rayleigh pressure multiplier
`gP = 0.6787278140712512` is therefore kept fixed and shown only as a comparison.
No gain was fitted, tuned, or applied to source data, and no score or master field
was changed.

The bounded result is recorded in the
[curated JSON](rr52_microseism_diagnostic.json) and
[figure](figures/rr52/microseism_pressure_vertical_diagnostic.png). The complete
new diagnostic is local at
`results/DIAGNOSTICS/RR52_MICROSEISM/run_28de4707b751/`. It contains only derived
CSV/JSON/PNG products; it does not replace the saved paired run.

## Primary-source theory

[Ruan, Forsyth, and Bell (2014)](https://doi.org/10.1002/2014JB011162) define
microseism Rayleigh-wave admittance as

```text
eta(f) = uz(f) / deltaP(f)
       = <uz(f) deltaP*(f)> / <deltaP(f) deltaP*(f)> .
```

This is exactly the saved BRUIT-FM transfer direction and cross-spectral
convention: `H_dp = S_dp / S_pp`, in m/Pa, with upward displacement. It is not
the conventional long-period compliance of a seabed deforming under external
infragravity-wave loading. In the 0.1–0.2 Hz microseism band, pressure and
displacement are properties of the same propagating Rayleigh mode.

The forward prediction is the ratio of the mode's vertical-displacement and
normal-stress eigenfunctions at the seafloor. Ruan et al. model a fluid layer
with local fixed depth, density 1030 kg/m³, compressional speed 1.5 km/s and zero
shear speed over elastic sediments, crust and mantle. Thus neither a hydrostatic
`rho*g` conversion nor coherence alone predicts D/P. Water depth, sediment
thickness and shear speed, and deeper structure affect its amplitude and shape.
The fundamental mode should have zero phase in an elastic isotropic model;
higher modes can change the sign and mix distinct admittances.

Ruan et al. used the 0.1–0.2 Hz band where fundamental-mode microseisms were
highly coherent and restricted their analysis to signals with reported
coherence above 0.95. They tilt-corrected the vertical first. They did not infer
a DPG scale from this band alone: at every sediment thickness/velocity grid
point they also adjusted a constant pressure scale, constrained mainly below
0.1 Hz where their model was less sediment-sensitive, while prescribing water
depth from acoustic ranging.

[Saikia et al. (2020)](https://doi.org/10.1007/s11001-020-09402-0) make the
identifiability condition explicit in a later application. Their fitted
prefactor accounts for pressure **and vertical** amplitude relative to the modal
prediction. It is a DPG gain only if the seismometer gain is assumed correct;
they also fixed water depth from multibeam bathymetry, removed tilt before the
P/Z transfer, screened days by average squared coherence above 0.8, and fitted
sediment structure separately in 0.1–0.2 Hz.

With correction multipliers defined by `P_true = gP P_nominal` and
`D_true = gZ D_nominal`, the saved transfer obeys

```text
H_nominal(f) = (gP / gZ) eta(f; water, sediment, crust).
```

Consequently, RR52 P/Z spectra cannot distinguish `gP` from `gZ` or from the
amplitude of an unconstrained structural prediction. A unique pressure
multiplier requires at least an independently verified vertical amplitude
response, verified water depth, and a station-specific modal structure fixed or
constrained by independent data. Fitting a free modal curve to force agreement
with the Rayleigh multiplier would be circular.

## Saved RR52 measurement

The analyzer read the master-linked run `run_bfbef17f2c13` for
`YV.RR52.00.BDH` / `YV.RR52.00.BHZ`, 10–17 June 2013. All five saved artifacts
exist. Their SHA-256 values, all recorded raw-file sizes/mtimes, the source
inventory hash, and both response-snapshot hashes verified before reuse. The
individual window spectra reproduce the saved weekly auto spectra, complex cross
spectrum, transfer and squared coherence. The older curated validation run's
`spectra.npz` is byte-identical, but it was not substituted for the current
master-linked path.

The week retains 167 hourly, nonoverlapping segments. In 0.1–0.2 Hz:

- 361 bins are valid; median squared coherence is `0.802142`.
- Median nominal `|D/P|` is `4.9857904e-8 m/Pa`, with median phase `+4.119°`.
- Only 29 bins (8.03%) exceed the Ruan `0.95` coherence reference; they occur
  between 0.10056 and 0.12056 Hz, not across the full band.
- Re-averaging the stored hourly spectra by UTC day gives six of seven days
  above the later Saikia mean-coherence reference of 0.8. Their complex
  daily-transfer stack has median nominal amplitude `5.0287852e-8 m/Pa`.
- Displaying the fixed Rayleigh multiplier changes the weekly median to
  `7.3457877e-8 m/Pa` and the six-day stack to `7.4091338e-8 m/Pa`, exactly by
  division by `0.6787278140712512`. These are reference transforms, not fitted
  microseism gains.

The two published coherence rules are diagnostic references only. This pilot
does not reproduce either study's duration, windowing, event rejection, tilt
correction or structural inversion, so their thresholds are not adopted as a
BRUIT-FM pass rule. The weekly coherence also falls with frequency through much
of 0.1–0.2 Hz; the result does not establish a clean fundamental-mode band over
the full interval.

## Sensor band, tilt and quality guidance

The saved StationXML describes a DPG and Trillium 240 OBS at 62.5 Hz. Relative
to each response's 1 Hz sensitivity, full nominal response amplitudes at
0.1/0.15/0.2 Hz are `1.0036/1.0054/1.0059` for BDH and
`0.9968/0.9969/0.9969` for BHZ. The band is therefore well represented by the
saved nominal responses and frequency grid. This is computational support, not
an independent usable-band or sensitivity certification. The inventory elevation
suggests 2918 m water depth, but bathymetry/acoustic depth and actual sensor
hardware remain unverified.

Tilt is not optional for a source-faithful calibration inference. Ruan et al.,
[Bell, Forsyth, and Ruan (2015)](https://doi.org/10.1785/0120140054), and Saikia
et al. first removed horizontal-correlated tilt from the vertical or explicitly
evaluated its ordering with pressure correction. The saved RR52 product contains
only P/Z spectra and cannot measure the tilt contribution. It is legitimate for
this model-free diagnostic, but not for estimating a pressure gain. A future
tilt-corrected comparison needs synchronized horizontal spectra or waveforms in
a new, traceable derived product.

The exact OBSIC material behind the brief is the official
[OBSIC quality-assessment method](https://obsic-metrics.whoi.edu/static/images/perform_quality_assessment-description.pdf).
It labels passive channel-hours from the PSD level at 0.2 Hz and whole-curve PSD
curvature, then aggregates good hours. Its own limitations include false labels
and omitted metrics such as anomalous clock drift. It is a performance overview,
not a P/Z transfer-function model, calibration coefficient, or acceptance rule
for this check.

## Reproduction and next constraint

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.microseism_diagnostic
```

Each invocation creates a fresh local result directory. The analyzer has no
waveform processing, master-update, gain-fitting or score-writing route.

The next defensible extension is to obtain an independent RR52 water-depth and
shallow-structure constraint, verify the vertical amplitude response, and then
build a station-specific coupled modal prediction. In parallel, a new bounded
P/H/Z product can quantify tilt before the transfer is interpreted. Until those
constraints exist, the saved P/Z curve is scientifically useful for response
phase/coherence/shape diagnostics but cannot certify the pressure multiplier.
