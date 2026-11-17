# PressureQC's project folder

Moved from PressureQC 1.3.0. This folder is laid out like PressureQC's checkout,
and the code it holds is ComPy's `compy_pressure` package (`../compy_pressure/`):

- [`docs/`](docs/): methods, validation notes, reviews, policies and the data
  behind them; the results report is [`docs/REPORT.md`](docs/REPORT.md).
- [`scripts/`](scripts/): the run and review scripts of the campaign.
- [`tests/`](tests/): the package's tests, run by ComPy's `pytest`.

The dated reports (PDFs, figures and scored tables, 136 MB) stay in
PressureQC, which is kept as a backup: [docs/reports](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports).
The ComPy user guide is [`docs/PRESSURE.md`](../docs/PRESSURE.md).

Run the scripts from this folder with both ComPy and it on the path:
`cd pressure && PYTHONPATH=..:. python scripts/plot_network_median_ppsd.py --help`.
The shell scripts expect the conda environment named in each file.

The documents are PressureQC's, unchanged apart from their links. In ComPy the
names are:

| PressureQC | ComPy |
|---|---|
| package `compy_pressure`, `compy_pressure/` | `compy_pressure`, `../compy_pressure/` |
| `compy_pressure` | `compy-pressure` |
| `compy_pressure-inventory` | `compy-pressure-inventory` |
| `compy_pressure-tilt` | `compy-pressure-tilt` |
| `compy_pressure-local-p` | `compy-pressure-local-p` |
| `compy_pressure-microseism` | `compy-pressure-microseism` |
| `compy_pressure-paired` | `compy-pressure-paired` |
| `compy_pressure-group-noise-model` | `compy-pressure-noise-model` |
| `compy_pressure-noise-model` | `compy-pressure-approved-noise-model` |

PressureQC's README and changelog are here as
[`HISTORY_README.md`](HISTORY_README.md) and
[`HISTORY_CHANGELOG.md`](HISTORY_CHANGELOG.md).
