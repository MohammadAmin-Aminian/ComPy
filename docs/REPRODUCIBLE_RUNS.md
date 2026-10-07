# Reproducible inversion and chain diagnostics

## Offline benchmark

After installing the development environment (`python -m pip install -e '.[dev]'`),
run from the repository root:

```bash
MPLBACKEND=Agg python _Example/reproducible_inversion.py --output synthetic-results
```

The example uses a three-layer known model, 12 frequencies, 3% independent
Gaussian measurement noise, a biased shallow starting velocity and four seeded
chains of 1,000 saved states each. It discards the first 250 states for summaries.
All chains share the same starting model because ComPy currently derives prior
bounds from that model. Changing initial models can change the posterior target,
which makes their joint convergence diagnostics inappropriate.

The short run demonstrates recording and diagnostics; it is not a converged
scientific inversion or an independent benchmark of the elastic solver. The
observations are generated with ComPy's own forward model. An improved data fit
need not imply recovery of a unique true velocity structure.

## Saved files

| File | Contents |
|---|---|
| `arrays.npz` | True/starting models, observations, uncertainties, frequency, all chains, predictions, weighted residual norms and acceptance rates |
| `metadata.json` | Settings, units, seeds, software/dependency versions, source hashes, UTC timestamp and chain diagnostics |
| `manifest.json` | SHA-256 checksums of the NPZ and metadata files |
| `diagnostics.png` / `diagnostics.svg` | Residual traces, observation/prediction comparison and shallow velocity traces |

Existing output paths are refused. Numeric run records are staged before publication,
so failed writes do not leave a partially written run. Load arrays without pickle:

```python
import numpy as np

with np.load("synthetic-results/arrays.npz", allow_pickle=False) as result:
    models = result["models"]  # chains × layers × 4 model columns × iterations
    predictions = result["predictions"]  # chains × iterations × frequencies
```

The sampler's misfit is a weighted residual **norm**, not chi-square. Square it
to obtain the summed squared normalized residual. Acceptance rates cover all
proposals, including burn-in. A checksum detects altered files; it does not prove
scientific correctness. Numerical reproducibility depends on the recorded
software/platform environment; UTC timestamps naturally differ between runs.

## Apply the helpers to your own inversion

```python
from compy_diagnostics import save_run, summarize_chains

# Draws from independent chains targeting the same posterior:
# shape = (chains, iterations, parameters).
report = summarize_chains(draws, parameter_names, burnin=500)
save_run(
    "my-run",
    {"models": models, "observations": observations, "uncertainty": uncertainty},
    settings={"seeds": seeds, "water_depth_m": depth, "alpha": alpha},
    diagnostics=report,
)
```

The helpers are separate from the legacy sampler; existing return tuples are
unchanged. Add complete sampler options, input identifiers and units to `settings`.
For a runtime-only installation, add the optional dependency with
`python -m pip install -e '.[diagnostics]'`.

## Interpretation

ArviZ computes rank-normalized folded split-Rhat and bulk/tail effective sample
sizes. Reports flag Rhat ≥ 1.01 and ESS below 100 per chain. Constant/stuck chains
produce undefined diagnostics (`null` in JSON), with an explicit warning. Duplicate
chains are also flagged. Five/50/95-percent quantiles describe retained samples;
they are not reliable posterior intervals until sampling has been assessed.

Inspect **all** sampled variables and relevant derived quantities. Fixed density,
Vp, half-space properties and conserved total finite depth are excluded from the
example's parameter diagnostics; individual finite-layer thicknesses, velocities
and the data misfit are included. These checks do not replace trace inspection,
Monte Carlo error assessment, sensitivity analysis, independent forward-solver
comparison or validation on actual station waveforms.

References:

- [Stan: convergence diagnostics](https://mc-stan.org/learn-stan/diagnostics-warnings.html).
- [Vehtari et al. (2021): rank-normalization, folding and localization](https://doi.org/10.1214/20-BA1221).
- [ArviZ diagnostics implementation](https://github.com/arviz-devs/arviz/blob/v0.22.0/arviz/stats/diagnostics.py).
