# ComPy documentation

## Guides

| Document | Content |
|---|---|
| [User guide](USER_GUIDE.md) | Input conventions, compliance estimators, inversion API, pressure QC, plotting |
| [Pipeline guide](../pipeline/README.md) | Running a whole network: stages A–G, configuration, restarts |
| [Examples](../examples/README.md) | Offline synthetic examples and per-step scripts for real data |

## Method and validation

| Document | Content |
|---|---|
| [Inversion](INVERSION.md) | Parameterization, priors, sampler, forward model and its validation |
| [Pressure tools](PRESSURE.md) | `compy_pressure`: QC tests, three calibration methods, noise models |
| [Validation](VALIDATION.md) | Test coverage and scientific references |
| [Reproducible runs](REPRODUCIBLE_RUNS.md) | Run records, array shapes and chain diagnostics |
| [Compliance audit (October 2026)](AUDIT_2026-10.md) | Review of the compliance chain against the literature and the fixes it led to |

## Results

| Folder | Content |
|---|---|
| [XE, Pacific Array Young ORCA](results/xe_young_orca_2026-10-29/README.md) | 30 stations: compliance, Vs, serpentinization, PDF report |
| [YV compliance, RHUM-RUM](results/yv_compliance_2026-10-28/README.md) | Compliance of eight Indian Ocean stations, gains, stability |
| [YV inversion, RHUM-RUM](results/yv_inversion_2026-10-29/README.md) | Vs posteriors, variants and paper-style figures |
| [Synthetic benchmark](benchmark/README.md) | Settings and diagnostics of the CI inversion benchmark |

Each results folder holds the scripts that produced it, its figures and its
tables. Large intermediate files (waveforms, posterior draws) are not
distributed.

## Releases

- [Changelog](../CHANGELOG.md)
- [Release notes 2.5](RELEASE_2_5.md)
