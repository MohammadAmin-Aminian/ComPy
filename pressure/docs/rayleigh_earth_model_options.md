# Optional Earth-layer model and YL gain comparison

The standalone and network calibrations accept either an existing dispersion CSV
or an explicit Earth-layer CSV, never both. The Earth-layer option computes
fundamental Rayleigh phase/group curves using disba 0.7.0. Numba accelerates the
solver; it is not an Earth model. Existing tabulated-input commands remain valid.

```bash
# Standalone option, replacing --dispersion-file:
--earth-model-file path/to/verified_layers.csv

# Network option, replacing --rayleigh-dispersion-file:
--rayleigh-earth-model-file path/to/verified_layers.csv
```

The layer CSV must have this exact header:

```csv
thickness_km,vp_km_s,vs_km_s,density_g_cm3
```

Every finite layer has positive thickness; the last layer has zero thickness and
represents the half-space. Velocities/densities are checked, and only the first
layer may have zero S velocity (fluid). The generator requests 160 increasing
periods, covering at least 0.02–0.1 Hz and the fit band. It requires complete
phase/group solutions and converts km/s to m/s. It does not infer missing layers,
change density or silently extrapolate a supplied dispersion table.

Generated curves are saved as `model_dispersion.csv` in each new station run.
The result records model/curve hashes, solver versions, algorithm, mode, units
and period limits. Network reuse depends on model contents and solver provenance.
One generated curve is shared across a network invocation, rather than solved
again per station. Scientific applicability still needs validation; a successful
solver is not an independently verified regional model.

Optional installation on another environment:

```bash
python -m pip install -e '.[dispersion]'
```

The active environment now has disba 0.7.0 and Numba 0.68.0, with the existing
NumPy/SciPy versions retained. A homogeneous-half-space test independently solves
the Rayleigh secular equation and checks phase/group velocities and units. A
synthetic calibration test verifies generated-curve provenance, known-gain
recovery, network wiring and master preservation. Regenerating the original
cached RR52 curves at their original periods gives zero maximum phase difference
and about 1.55e-11 m/s maximum group difference.

## Selecting a model

The bundled historical layer file is
`src/bruitfm/reference_data/rr52_historical_rayleigh_layers.csv`. It can be used
to test the generator or extend that same assumed model to 0.1 Hz. It is not a
verified Zha/Lau model.

For reproducing YL, obtain and verify numerical average Lau Basin layers or
dispersion from Zha et al. (2014), the model referenced in
[Zha & Webb (2016), Appendix A](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688).
The retrieved papers establish the reference, but no verified numerical model
has yet been supplied to this option. No preset is labelled Zha based on guessed
layers or a fit to the desired gains.

## Why the reported gains look different

BRUIT-FM fits `g * observed_P_over_a = theoretical_P_over_a`. Zha's Figure A2
describes the observed ratio being about 1.5 times the theoretical ratio. Thus
compare the published coefficients with **1/g**, subject to confirming the same
nominal instrument responses. The arithmetic comparison for the original YL pilot:

| Station | Our pressure multiplier g | Reciprocal 1/g | Zha coefficient | Relative difference |
|---|---:|---:|---:|---:|
| A02W | 0.704107 | 1.420239 | 1.50 | −5.32% |
| A08W | 0.741370 | 1.348854 | 1.43 | −5.67% |
| A10W | 0.585030 | 1.709313 | 1.82 | −6.08% |
| A12W | 0.639401 | 1.563965 | 1.67 | −6.35% |

A06W has only two accepted events, so no station coefficient is compared here.
These are method comparisons, not uncertainty bounds or calibration acceptance.
Source: [published Table 1 and Figure A2](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688).

Remaining differences have not been attributed to a single cause. Our saved
pilot used three events, 0.03–0.07 Hz, a borrowed model, automatic 30-minute
windows, additional phase/timing rules and the median of per-event gains. The
paper uses many more large events, a Lau Basin model, 0.04–0.1 Hz and an averaged
spectral comparison. Exact spectral processing/window choices and nominal
response versions still need reconciliation. The sampling cap has since been
removed in code, but the saved pilot has not been rerun.

A bounded saved-spectrum diagnostic held those original candidate events fixed
and used the same historical layers extended to 0.1 Hz. Applying the current
coherence/phase masks in 0.04–0.1 Hz gives pressure-multiplier medians of 0.723537,
0.750673, 0.589205 and 0.636017 respectively. This band change does not uniformly
close the discrepancy. These are diagnostics, not newly accepted station gains.
The calculation wrote separate validation output and preserved master/workflow
hashes. [Numerical comparison](yl_gain_convention_comparison.json).

No new real-data calibration batch, response correction or gain application was
performed for the model-option change.
