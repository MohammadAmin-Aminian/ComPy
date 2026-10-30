# BRUIT-FM

**Quality control and calibration of ocean-bottom pressure gauges, from a raw
SDS waveform archive and its StationXML to a scored table of every channel.**

[![Checks](https://github.com/MohammadAmin-Aminian/bruit-fm/actions/workflows/checks.yml/badge.svg)](https://github.com/MohammadAmin-Aminian/bruit-fm/actions/workflows/checks.yml)
![Version](https://img.shields.io/badge/version-1.3.0-2b6cb0)
![Python](https://img.shields.io/badge/Python-3.13%2B-3776AB?logo=python&logoColor=white)
![ObsPy](https://img.shields.io/badge/built%20on-ObsPy-276749)
[![License: GPL v3](https://img.shields.io/badge/license-GPL--3.0-blue)](LICENSE)

> **BRUIT-FM is now part of [ComPy](https://github.com/MohammadAmin-Aminian/ComPy)**
> The code is ComPy's `compy_pressure` package and its
> commands are `compy-pressure*`. The tests, scripts and documentation are in
> ComPy's `pressure/` folder. This repository is kept as a backup and holds the
> dated reports (`docs/reports/`). Version 1.3.0 is the last release; see
> [integration/README.md](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/integration/README.md).

BRUIT-FM finds every network, station and channel in a seismic data archive.
It computes response-corrected noise spectra (PPSDs) and runs three spectral
quality tests and three independent calibration checks on each pressure
channel. Every result is written back to one master table. The runs are
restartable, parallel and checkpointed, so a whole 18 TB drive can be
processed unattended on one workstation.

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_all_station_qc/report/img/class_map.png?raw=true" width="85%"
       alt="Map of 846 ocean-bottom pressure channels coloured by QC class">
  <br><em>846 pressure channels from 30 networks, coloured by QC class
  (v1.0.0 run, provisional thresholds). See the <a href="docs/REPORT.md">results report</a>.</em>
</p>

## What it does

| Step | What you get | Module |
|---|---|---|
| **Inventory** | A master table (CSV + Excel) of every channel epoch in your StationXML: coordinates, depth, sensor, response, sample rate, data coverage | `bruitfm.inventory` |
| **PPSDs** | ObsPy PPSDs with the full instrument response removed: 3 h windows, 50 % overlap, down to 1 mHz. Saved as `.npz` with plots and provenance | `bruitfm.ppsd`, `bruitfm.plot_ppsd` |
| **Test 1: expected levels** | Fraction of frequency bins inside the Webb (1998) and Brown et al. (2014) seafloor pressure-noise limits | `bruitfm.pressure_test` |
| **Test 2: nearby stations** | Median daily difference from neighbouring gauges at a similar depth | `bruitfm.nearby_test` |
| **Test 3: infragravity** | Residual from the depth-dependent Godin et al. (2013) infragravity-wave model | `bruitfm.infragravity_test`, `bruitfm.godin_pilot` |
| **Calibration 1: Rayleigh waves** | Pressure gain from teleseismic Rayleigh-wave pressure/vertical ratios, 30–70 mHz, with optional tilt cleaning | `bruitfm.rayleigh_calibration`, `bruitfm.rayleigh_network` |
| **Calibration 2: microseisms** | Gain bracket from the water-layer admittance of coherent microseism energy | `bruitfm.paired_spectra`, `bruitfm.microseism_check` |
| **Calibration 3: local P waves** | Gain from the pressure/vertical ratio of local P arrivals, by frequency band | `bruitfm.local_p_check` |
| **Tilt correction** | ComPy-cleaned vertical PSDs, before and after correction | `bruitfm.tilt` |
| **Noise model** | Station-balanced percentile noise curves (p5–p95) by ocean basin, sensor type, water depth and local season, drawn as percentile bands, heat maps and all station lines | `bruitfm.group_noise_model` |
| **Maps and network figures** | Station maps with coastlines; PPSD medians by network, depth, sensor type and ocean basin | `bruitfm.maps`, `bruitfm.plot_network_ppsd` |

Nothing is applied silently. The calibration checks estimate a gain and record
it, but they never change a response or the data. Every stage writes its inputs,
settings and code version beside its result.

## Installation

BRUIT-FM needs Python 3.13 or later. A conda environment is the easiest way to
get ObsPy and Cartopy:

```bash
conda create -n bruitfm -c conda-forge python=3.13 obspy cartopy
conda activate bruitfm
git clone https://github.com/MohammadAmin-Aminian/bruit-fm.git
cd bruit-fm
pip install -e '.[all]'          # core + Rayleigh (TiSKitPy), dispersion, maps
```

Optional extras:

| Extra | Adds | Needed for |
|---|---|---|
| `rayleigh` | `tiskitpy==2.3.1` | Rayleigh calibration with tilt cleaning |
| `dispersion` | `disba==0.7.0` | Rayleigh dispersion from a layered Earth model |
| `maps` | `cartopy>=0.24` | Station maps with coastlines |
| `tilt` | [ComPy](https://github.com/MohammadAmin-Aminian/ComPy) from GitHub | `bruitfm-tilt` |
| `dev` | `pytest`, `ruff` | Running the test suite |

Check the installation:

```bash
python -c "import bruitfm; print(bruitfm.__version__)"
python -m pytest            # about 500 tests on synthetic StationXML, no data needed
```

## Input data

- **Waveforms** in an [SDS archive](https://www.seiscomp.de/seiscomp3/doc/applications/slarchive/SDS.html):
  `<root>/<year>/<net>/<sta>/<chan>.D/<net>.<sta>.<loc>.<chan>.D.<year>.<doy>`.
- **StationXML** with full response stages, one or more files in a folder. A
  file without response stages is reported, not guessed.

## Quick start

Each stage resumes from its checkpoints, so the same command can be rerun after
an interruption. `--workers max` uses every CPU.

```bash
DATA=/path/to/SDS
XML=/path/to/stationxml
OUT=~/bruitfm_run

# 1. Master table of every channel epoch
bruitfm inventory --inventory-dir $XML --master-dir $OUT/master

# 2. PPSDs and plots for every pressure channel of one network
bruitfm ppsd  --data-dir $DATA --inventory-dir $XML --master-dir $OUT/master \
              --results-dir $OUT/results --network YV --pressure-only --workers max
bruitfm plots --data-dir $DATA --inventory-dir $XML --master-dir $OUT/master \
              --results-dir $OUT/results --network YV --pressure-only

# 3. The three QC tests
for test in pressure nearby infragravity; do
  bruitfm $test --data-dir $DATA --inventory-dir $XML --master-dir $OUT/master \
                --results-dir $OUT/results --network YV
done

# 4. Rayleigh calibration (USGS catalogue, RR52 dispersion, TiSKitPy tilt cleaning)
bruitfm rayleigh --data-dir $DATA --inventory-dir $XML --master-dir $OUT/master \
                 --results-dir $OUT/results --network YV \
                 --rayleigh-catalog events.json --rayleigh-tilt-method tiskitpy

# Open the saved plots of one station
bruitfm view --data-dir $DATA --inventory-dir $XML --master-dir $OUT/master \
             --results-dir $OUT/results --network YV --station RR52 --view-kind ppsd
```

`bruitfm all` runs every stage in order, from the inventory through the PPSDs, plots, the three QC tests, the Godin pilot and the paired spectra. Use
`--station`, `--seed-id` or `--channels 'BH?'` to narrow a run.

Other commands:

| Command | Purpose |
|---|---|
| `bruitfm-paired`, `bruitfm-microseism` | Paired pressure/vertical spectra for one week, then the microseism admittance check |
| `bruitfm-local-p` | Local P-wave calibration for one pressure/vertical pair and a local catalogue |
| `bruitfm-tilt` | Tilt-corrected vertical PPSDs for every station with a three-component seismometer |
| `bruitfm-group-noise-model` | Percentile noise models by basin, sensor type, depth and season from a scored table (e.g. `--select-column final_class --select-value good`) |
| `bruitfm-noise-model` | Approval-gated noise model from a formal review document |

All commands accept `--help`. They are also available as `python -m bruitfm.<module>`.

### From Python

```python
from pathlib import Path

from bruitfm.inventory import load_master_dataset
from obspy.signal import PPSD

master = load_master_dataset(Path.home() / "bruitfm_run/master")
pressure = master[master.ppsd_status.eq("complete") & master.Seed_id.str.endswith("DH")]
ppsd = PPSD.load_npz(pressure.ppsd_npz_file.iloc[0], allow_pickle=True)
ppsd.plot()
```

## Outputs

```text
master/master_inventory.{csv,xlsx}     one row per channel epoch, every stage's status and result
results/<net>/<sta>/<epoch>/
    ppsd.npz, inventory.xml, processing.json
    plots/                             PPSD, spectrogram, all PSD curves
    pressure_expected/  nearby/  infragravity/  rayleigh/  tilt/
```

## Results

The [**results report**](docs/REPORT.md) covers 846 pressure channel
epochs and about 1,300 seismometer channels from 30 networks. It includes
methods, figures and the provisional scores. The full dated reports, with the
scripts that built them, are in [`docs/reports/`](docs/reports/).

## Documentation

| Topic | Document |
|---|---|
| Results of the v1.0.0 campaign | [docs/REPORT.md](docs/REPORT.md) |
| Running every stage on a network | [docs/network_workflow.md](docs/network_workflow.md) |
| Processing settings | [docs/rr52_methods.md](docs/rr52_methods.md) |
| From-scratch pipeline over a whole drive | [docs/from_scratch_pipeline.md](docs/from_scratch_pipeline.md) |
| Rayleigh / microseism / local P methods | [docs/rayleigh_calibration_review.md](docs/rayleigh_calibration_review.md), [docs/microseism_pressure_vertical_check.md](docs/microseism_pressure_vertical_check.md), [docs/local_p_method.md](docs/local_p_method.md) |
| Tilt correction | [docs/tilt_methods.md](docs/tilt_methods.md) |
| Calibration assumptions | [docs/calibration_assumptions.md](docs/calibration_assumptions.md) |
| Moving the core into ComPy (handoff, questions, steps) | [integration/README.md](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/integration/README.md) |
| Full development history | [docs/development_log.md](docs/development_log.md), [CHANGELOG.md](BRUITFM_CHANGELOG.md) |
| References | [docs/references.bib](docs/references.bib) |

## Project layout

```text
src/bruitfm/     the package (inventory, archive, ppsd, tests, calibrations, network runner, maps)
scripts/         run scripts and one-off reviews used during development
tests/           pytest suite on small synthetic inventories
docs/            methods, validation notes, results report and dated reports
integration/     manifest, export and comparison tools for the move into ComPy
```

## Development

```bash
pip install -e '.[dev,all]'
python -m pytest
ruff check . && ruff format --check .
```

GitHub Actions runs the same checks, builds the wheel and smoke-tests it outside
the checkout. See [CONTRIBUTING.md](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/CONTRIBUTING.md).

## Citing

If you use BRUIT-FM, please cite it using [CITATION.cff](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/CITATION.cff) (GitHub's
"Cite this repository" button), along with the methods it implements: Webb (1998),
Brown et al. (2014), Godin et al. (2013) and the calibration papers listed in
[docs/references.bib](docs/references.bib).

## License and relation to ComPy

BRUIT-FM is released under the [GNU General Public License v3.0](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/LICENSE), the
same license as [ComPy](https://github.com/MohammadAmin-Aminian/ComPy). It
started as ComPy's pressure-calibration component and was developed here on
its own. Version 1.3.0 moved all of it back into ComPy, replacing ComPy's old
`Pressure_calibration.py`:

- the code became the `compy_pressure` package;
- the tests, scripts and documentation went to ComPy's `pressure/` folder.

The move was done by `integration/export_to_compy.py` and checked against
BRUIT-FM's outputs. Development continues in ComPy; this repository is kept
as a backup and as the home of the dated reports. The details are in
[integration/README.md](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/integration/README.md).

## Acknowledgements

Developed at the Institut de Physique du Globe de Paris (IPGP) by
Mohammad-Amin Aminian, with Wayne Crawford. Data from the FDSN networks listed
in the results report, accessed through EarthScope (IRIS) and the operators' data
centres.
