# ComPy user guide

Reference for the ComPy modules: input conventions, the compliance estimators,
the Bayesian and legacy inversions, pressure-gauge QC and the plotting tools.
The step-by-step workflow with figures is in the [README](../README.md).

- [Input conventions](#input-conventions)
- [Compliance estimators](#compliance-estimators)
- [Bayesian Vs inversion (`compy_inversion`)](#bayesian-vs-inversion-compy_inversion)
- [Legacy forward model and inversion (`inv_compy`)](#legacy-forward-model-and-inversion-inv_compy)
- [Deep-crust inversion tools](#deep-crust-inversion-tools)
- [Pressure-gauge QC and calibration (`compy_pressure`)](#pressure-gauge-qc-and-calibration-compy_pressure)
- [Plotting](#plotting)
- [Testing and reproducible runs](#testing-and-reproducible-runs)

## Input conventions

| Operation | Vertical channel | Pressure channel | Other inputs |
|---|---|---|---|
| `stacked_compliance`, `Calculate_Compliance*` | response-corrected **displacement, m** | nominal response-corrected pressure, **Pa** | positive water depth, m; multiplicative pressure gain |
| `compy_pressure` Rayleigh calibration | raw counts in an SDS archive; response removed internally to acceleration | raw counts; response removed internally to pressure | StationXML with full responses, ComCat GeoJSON catalogue |
| `Rotate`, `rotate_smoothly` | physically consistent seismic channels | passed through the cleaning workflow | `*1`, `*2`, `*Z` seismic channels |
| `calc_norm_compliance` / inversion | measured normalized compliance, **1/Pa** | not used directly | Hz, m; model columns below |
| `compy_pressure.rayleigh_dispersion` | not used | not used | layer CSV in **km, km/s, km/s, g/cm³**, or phase/group velocities in m/s |

`*Z` denotes the vertical component and `*H` the pressure channel (normally
`BDH`). ComPy expects exactly one of each for spectral processing. Process one
station and one response epoch at a time, and merge trace fragments and resolve
gaps before analysis. The core helpers reject masked samples, nonfinite data,
duplicate trace IDs, unequal sample rates and off-grid start times.

Model columns for the inversion and the elastic forward calculation:

```text
[finite-layer thickness (m), density (kg/m³), Vp (m/s), Vs (m/s)]
```

The last row is an elastic half-space; its thickness does not enter the
propagator (use zero in custom models). Do not pass a disba model directly to
the elastic compliance solver.

Pressure gain is defined as `calibrated_pressure = gain_factor * nominal_pressure`,
so the pressure PSD scales with `gain_factor**2`. Apply the gain once;
instrument-response removal is a separate step.

Calibration thresholds use **magnitude-squared coherence**; compliance outputs
use **amplitude coherence**, its square root.

## Compliance estimators

### Stacked compliance (recommended)

```python
import compy
from compy_processing import GJI2025, stacked_compliance

pair = cleaned.select(channel="BHZ") + cleaned.select(channel="BDH")
windows = compy.split_stream(pair, 2100)  # 35-min, non-overlapping
result = stacked_compliance(
    windows, depth, gain, band=(0.007, 0.019), criteria=GJI2025, nseg=1024
)
result.compliance, result.uncertainty, result.phase_deg, result.n_windows
```

`stacked_compliance` returns k·|S_pz|/S_pp from cross-spectra stacked over the
selected windows, with the gravity correction applied, together with the
Bendat & Piersol random error, the phase, and the per-window mean and standard
deviation. `SelectionCriteria` with the `GJI2025` preset applies the published
gates: median coherence ≥ 0.8 (or 0.9), pressure PSD 10–50 dB and
acceleration PSD −170 to −150 dB.

The gravity-wave wavenumber solves `ω² = g k tanh(kH)` with a bounded root
solver. The wave-attraction correction converts acceleration to displacement
units before it is combined with Z/P.

### Window-based estimators

```python
curves, coherence, windows, frequency, full_frequency, scatter = (
    compy.Calculate_Compliance_beta(
        displacement_pressure_stream,
        depth=4000,
        gain_factor=gain,
        time_window=2,
        f_min_com=0.007,
        f_max_com=0.02,
        nseg=4096,
        plot=False,
    )
)
```

| Function | Return tuple | Processing-window step |
|---|---|---|
| `Calculate_Compliance_beta` | curves, amplitude coherence, selected streams, selected frequencies, full frequencies, standard deviation across curves | 5 minutes |
| `Calculate_Compliance` | same first five items, peak-to-peak scatter, theoretical uncertainty estimate | 1 minute |

`curves` has shape `(accepted_windows, frequency_bins)`. Beta returns
frequencies from 0.001–0.1 Hz; the older estimator returns 0.005–0.025 Hz.
Providing `depth` avoids station-service calls; otherwise it is read from the
vertical channel's elevation. Spectra use a Hann window with 50 % overlap and
`nseg = 4096` by default; the frequency spacing is `sampling_rate / nseg`.
The theoretical uncertainty is `|η|·sqrt(1−γ²)/(γ·sqrt(2n))` with amplitude
coherence γ.

Figures of every processing step for one station are in the
[README workflow](../README.md#workflow-step-by-step).

## Bayesian Vs inversion (`compy_inversion`)

Parameterization and sampler follow Zha & Webb (2016); the CRUST1.0 Moho and the
±5 % gain range follow Aminian et al. (2025). Method and validation:
[INVERSION.md](INVERSION.md).

```python
import compy_inversion as ci

column = ci.crust1_column(-20.4723, 68.1094, "crust1.0/")
data = ci.ComplianceData.from_npz(
    "RR52_compliance_compy.npz",
    subsets=monthly_estimates,  # (n_months, n_freqs), optional
    gain_rel=0.0,               # gain is sampled below
    floor_rel=0.01,
)
model = ci.CrustalModel(moho_depth=column["moho_depth"], sediment_max=1000.0)
post = ci.sample(data, model, n_chains=8, n_iter=80000, n_burn=30000,
                 thin=10, gain_bounds=(0.95, 1.05))
prior = ci.sample(data, model, prior_only=True, n_chains=8, n_iter=80000, n_burn=30000)
post.summary()["max_rhat"], post.summary()["chi2"]
vs_q = post.vs_quantiles(depths)  # (5 quantiles, len(depths)) in m/s
post.save("RR52_final")           # .npz draws and .json summary
```

- **Data.** `ComplianceData` holds a full covariance: random error, monthly
  scatter, an optional gain error shared by all frequencies, and a floor.
  `ComplianceData.join` inverts epochs with different gains together.
- **Model.** `CrustalModel`: Hamilton sediment of free thickness over crustal
  layers that double in thickness to the Moho, above a mantle half-space.
  Vp/Vs is 1.85 by default, free within bounds, Brocher (2005), or Vp fixed
  from a reference.
- **Sampler.** Independent DRAM chains with a Gibbs-sampled smoothing weight;
  `temperatures=` adds parallel tempering for multimodal posteriors.
- **Results.** `Posterior` gives R-hat, ESS, χ² and Vs-depth quantiles;
  `vs_kernels` and `sensitivity_depth` give the depth sensitivity.
- **Serpentinization.** `ci.serpentinization` and
  `post.serpentinization_quantiles(z)` apply the Miller & Christensen (1997)
  relation (0 % at 4.2 km/s, 100 % at 2.2 km/s) to every posterior draw.

The forward model is the `raydep` propagator, compiled with Numba when it is
installed (about 140 times faster). Chains run in spawned processes: call
`sample` from a script guarded by `if __name__ == "__main__":` or from a
notebook, or pass `processes=1`.

## Legacy forward model and inversion (`inv_compy`)

```python
import numpy as np
import inv_compy as inv

model = np.array(
    [[100, 2200, 3000, 1500],
     [1000, 2800, 6000, 3500],
     [0, 3300, 8000, 4500]],
    dtype=float,
)
f = np.array([0.005, 0.007, 0.010, 0.015])
prediction = inv.calc_norm_compliance(4000, f, model)

chain, profiles, prior, misfit, predictions, likelihood, acceptance = (
    inv.invert_compliance_beta(
        measured_compliance, f, 4000,
        starting_model=model, s=measurement_uncertainty,
        iteration=10000, sigma_v=5, sigma_h=5, alpha=0.25,
        seed=0, return_profiles=False,
    )
)
```

`chain` has shape `(layers, 4, iterations)`; `predictions` has shape
`(iterations, frequencies)`; `misfit` and `likelihood` have shape
`(1, iterations)`. For long runs keep `return_profiles=False`: 100 000 profiles
on a 30-km metre grid need about 24 GB.

Density and Vp are fixed; Vs and finite-layer thickness are sampled with a
symmetric Gaussian random walk, and the half-space stays fixed. The data term
is `sum(((observed - predicted)/s)**2)`; the roughness term is the squared
second derivative on a 1-m grid. `s` must be finite and strictly positive.
`seed` controls a local RNG. `Model_V2` provides station templates with
3/6/9/12 subdivisions for the eight RHUM-RUM stations. The names
`invert_compliace*`, `liklihood`, `Roughness` and `Comliance_uncertainty`
remain available as aliases for existing notebooks.

## Deep-crust inversion tools

Opt-in additions to `inv_compy`; all defaults reproduce earlier chains draw
for draw.

| Addition | Purpose |
|---|---|
| `sigma_v`/`sigma_h` arrays, `adapt_until` | Per-layer steps tuned during burn-in to ~30 % acceptance, then frozen |
| `initial_state`, `vs_bounds` | Dispersed chain starts under an explicit prior |
| `vary_thickness=False` | Fixed-geometry minimum-structure grid (Crawford et al. 1991) |
| `alpha_bounds` | Smoothing weight sampled from the data (Malinverno & Briggs 2004) |
| `data_covariance` | Correlated errors, e.g. a pressure-gain error common to all frequencies |
| `interface_move="adjacent"` | One interface moved at a time in free-interface runs |
| `invert_occam` | Occam minimum-structure inversion (Constable et al. 1987) |
| `invert_tempered` | Parallel tempering (Geyer 1991; Sambridge 2014) |
| `compliance_kernels`, `resample_model` | Numerical Vs sensitivity kernels |
| `compy_diagnostics.rank_rhat_ess` | Rank-normalized split R̂ and bulk ESS without ArviZ |

The propagator agrees with an independent Aki & Richards solver to about 1e-6.
Compliance senses structure to roughly λmax/6 of the longest water wave
(Crawford et al. 1991): about 9 km at 4.2 km water depth and 3.5 mHz. Compare
posterior and prior-only widths to see which depths the data resolve.
Compliance measured with ComPy before 2.0 lacks the wave-attraction term; apply
`compy_numerics.legacy_wave_attraction_correction` before inverting it.

```bash
python _Example/deep_crust_inversion.py --output deep-crust   # ~25 min, 2 cores
python _Example/plot_deep_posterior.py deep-crust
```

## Pressure-gauge QC and calibration (`compy_pressure`)

`compy_pressure` runs QC and calibration on every pressure gauge in a waveform
archive. It reads an SDS archive and its StationXML and writes one master table
that scores every pressure channel:

- response-corrected PPSDs;
- spectral QC against the Webb (1998) and Brown et al. (2014) seafloor limits,
  against nearby gauges at a similar depth, and against the Godin et al. (2013)
  infragravity model;
- three independent gain estimates: Rayleigh waves, microseisms and local P
  waves;
- percentile noise models by ocean basin, depth, sensor type and season.

Runs are parallel, checkpointed and restartable. The package comes from
[BRUIT-FM](https://github.com/MohammadAmin-Aminian/bruit-fm), where it processed
846 pressure channels from 30 networks.

```bash
python -m pip install -e '.[pressure]'
python _Example/5_Pressure_QC.py --output qc-example    # synthetic gauge, offline
compy-pressure inventory --inventory-dir $XML --master-dir run/master
compy-pressure ppsd --data-dir $SDS --inventory-dir $XML --master-dir run/master \
    --results-dir run/results --network YV --pressure-only --workers max
compy-pressure pressure --data-dir $SDS --inventory-dir $XML --master-dir run/master \
    --results-dir run/results --network YV
```

Commands, outputs and the Python interface: [PRESSURE.md](PRESSURE.md).

### Migrating from `Pressure_calibration` (ComPy ≤ 2.4.0)

| `calculate_spectral_ratio` | `compy_pressure` |
|---|---|
| raw stream and inventory | SDS archive at the native rate, StationXML file or folder |
| `mag=7`, `event_spans` | ComCat GeoJSON catalogue; Mw ≥ 6 with a magnitude schedule |
| `f_min`, `f_max` | `band` |
| `coh_trsh`, `mean_trsh` (0.97) | squared coherence ≥ 0.95 per bin (`threshold`) |
| `zchan`, `pchan` | pairing from the StationXML; `calibration_master(channels=...)` |
| `filt_freq1`, `filt_freq2` | none: responses are removed in the frequency domain |
| `water_density` | `density` (1028 kg/m³); depth from the StationXML elevation |
| `return_details=True` | `comparison.json`: per-event gains, median and scatter |

Scripts that import `Pressure_calibration` run with ComPy 2.4.0 (commit `7ad79c2`).

## Plotting

`ffplot` contains PSD, coherence, spectrogram and publication plotting
routines. `ffplot.compl` is a quick compliance diagnostic without the gravity
terms; use `stacked_compliance` for quality selection and corrections.
Export figures explicitly with `matplotlib.pyplot.savefig(...)`.
`final_plot(..., image_dir=...)` expects `RR36.png`, `RR38.png`, `RR40.png`,
`RR50.png`, `RR52.png` and `Rift_valley.png` in `image_dir`.

## Testing and reproducible runs

```bash
MPLBACKEND=Agg python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m build
```

The regression suite covers dispersion-relation residuals, the homogeneous
elastic half-space limit, identical-layer splitting, known spectral transfer
functions, gravity units, calibration scaling, event handling, sampler state
consistency, RNG isolation, window boundaries and MiniSEED round trips. GitHub
Actions runs it on Python 3.10–3.13, builds the wheel and exercises it outside
the checkout. Details: [VALIDATION.md](VALIDATION.md).

A complete offline multi-chain experiment (noisy synthetic observations, a
biased starting model, four seeds) records models, observations, settings,
software versions and checksums, and reports rank-normalized split R̂ and
bulk/tail ESS:

```bash
MPLBACKEND=Agg python _Example/reproducible_inversion.py --output synthetic-results
```

See [REPRODUCIBLE_RUNS.md](REPRODUCIBLE_RUNS.md) for the run record format.
