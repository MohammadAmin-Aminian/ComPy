# BRUIT-FM's project folder

Moved from BRUIT-FM 1.3.0. This folder is laid out like BRUIT-FM's checkout,
and the code it holds is ComPy's `compy_pressure` package (`../compy_pressure/`):

- [`docs/`](docs/): methods, validation notes, reviews, policies and the data
  behind them; the results report is [`docs/REPORT.md`](docs/REPORT.md).
- [`scripts/`](scripts/): the run and review scripts of the campaign.
- [`tests/`](tests/): the package's tests, run by ComPy's `pytest`.

The dated reports (PDFs, figures and scored tables, 136 MB) stay in
BRUIT-FM, which is kept as a backup: [docs/reports](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports).
The ComPy user guide is [`docs/PRESSURE.md`](../docs/PRESSURE.md).

Run the scripts from this folder with both ComPy and it on the path:
`cd pressure && PYTHONPATH=..:. python scripts/plot_network_median_ppsd.py --help`.
The shell scripts expect the conda environment named in each file.

The documents are BRUIT-FM's, unchanged apart from their links. In ComPy the
names are:

| BRUIT-FM | ComPy |
|---|---|
| package `bruitfm`, `src/bruitfm/` | `compy_pressure`, `../compy_pressure/` |
| `bruitfm` | `compy-pressure` |
| `bruitfm-inventory` | `compy-pressure-inventory` |
| `bruitfm-tilt` | `compy-pressure-tilt` |
| `bruitfm-local-p` | `compy-pressure-local-p` |
| `bruitfm-microseism` | `compy-pressure-microseism` |
| `bruitfm-paired` | `compy-pressure-paired` |
| `bruitfm-group-noise-model` | `compy-pressure-noise-model` |
| `bruitfm-noise-model` | `compy-pressure-approved-noise-model` |

BRUIT-FM's README and changelog are here as
[`BRUITFM_README.md`](BRUITFM_README.md) and
[`BRUITFM_CHANGELOG.md`](BRUITFM_CHANGELOG.md).
