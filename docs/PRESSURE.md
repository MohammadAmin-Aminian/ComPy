# Pressure-gauge QC and calibration (`compy_pressure`)

`compy_pressure` checks and calibrates ocean-bottom pressure gauges across a
whole waveform archive. It starts from an SDS archive and its StationXML, and
ends with one master table that scores every pressure channel. It is the
core of [BRUIT-FM](https://github.com/MohammadAmin-Aminian/bruit-fm), moved
into ComPy unchanged apart from its name (see [Origin](#origin-and-updates)).

For one station and one teleseismic event at a time, ComPy's older
`Pressure_calibration.calculate_spectral_ratio` is still available; this
package does not modify it. `compy_pressure` is the archive-scale successor:
it adds three spectral QC tests, two more calibration checks, restartable
parallel runs and a provenance record for every result.

## Installation

```bash
python -m pip install -e '.[pressure]'        # pandas and openpyxl
python -m pip install -e '.[pressure,maps]'   # plus Cartopy, for station maps
```

TiSKitPy and disba, which the Rayleigh and microseism checks use, are already
ComPy dependencies. Python 3.10 to 3.13 is supported and tested.

## Inputs

- **Waveforms** in an [SDS archive](https://www.seiscomp.de/seiscomp3/doc/applications/slarchive/SDS.html):
  `<root>/<year>/<net>/<sta>/<chan>.D/<net>.<sta>.<loc>.<chan>.D.<year>.<doy>`.
- **StationXML** with full response stages, one or more files in a folder. A
  channel without response stages is reported, never guessed.
- For the Rayleigh calibration, an earthquake catalogue in USGS ComCat
  GeoJSON; for the local-P check, a local catalogue (ComCat CSV or GeoJSON).

## What it computes

| Step | Result | Command (`compy-pressure ...`) | Module |
|---|---|---|---|
| Inventory | Master table (CSV + Excel) of every channel epoch: coordinates, depth, sensor, response, sample rate, coverage | `inventory` | `inventory` |
| PPSDs | ObsPy PPSDs with the full response removed: 3 h windows, 50 % overlap, down to 1 mHz | `ppsd`, `plots` | `ppsd`, `plot_ppsd` |
| QC test 1 | Fraction of median PSD bins inside the Webb (1998) and Brown et al. (2014) seafloor limits | `pressure` | `pressure_test` |
| QC test 2 | Median daily difference from nearby gauges at a similar depth | `nearby` | `nearby_test` |
| QC test 3 | Residual from the depth-dependent Godin et al. (2013) infragravity model | `infragravity` | `infragravity_test` |
| Calibration 1 | Gain from teleseismic Rayleigh-wave pressure/vertical ratios, 30–70 mHz, with optional TiSKitPy tilt cleaning | `rayleigh` | `rayleigh_calibration`, `rayleigh_network` |
| Calibration 2 | Gain bracket from the water-layer admittance of coherent microseisms | `paired`, then `compy-pressure-microseism` | `paired_spectra`, `microseism_check` |
| Calibration 3 | Gain from the pressure/vertical ratio of local P arrivals | `compy-pressure-local-p` | `local_p_check` |
| Noise model | Station-balanced percentile noise curves (p5–p95) by ocean basin, sensor type, water depth and local season | `compy-pressure-noise-model` | `group_noise_model` |
| Figures | Network PPSD medians by depth, sensor and basin; station maps | — | `plot_network_ppsd`, `maps` |

Nothing is applied silently. The calibration checks estimate a gain and record
it, but they never change a response or the data. Each stage writes its
inputs, settings and code fingerprint next to its result, and reruns reuse a
finished result only when all of these still match.

The pressure gain follows ComPy's convention:
`calibrated_pressure = gain * nominal_pressure`.

## Quick start

```bash
DATA=/path/to/SDS
XML=/path/to/stationxml
OUT=~/pressure_run
COMMON="--data-dir $DATA --inventory-dir $XML --master-dir $OUT/master --results-dir $OUT/results"

compy-pressure inventory --inventory-dir $XML --master-dir $OUT/master
compy-pressure ppsd  $COMMON --network YV --pressure-only --workers max
compy-pressure plots $COMMON --network YV --pressure-only
for test in pressure nearby infragravity; do compy-pressure $test $COMMON --network YV; done
compy-pressure rayleigh $COMMON --network YV --rayleigh-catalog events.geojson \
    --rayleigh-dispersion-file compy_pressure/reference_data/rr52_historical_rayleigh_dispersion.csv
compy-pressure view $COMMON --network YV --station RR52 --view-kind ppsd
```

`compy-pressure all` runs the inventory, check, PPSD, plot, QC and
paired-spectra stages in order; the Rayleigh calibration runs only when asked.
Narrow a run with `--station`, `--seed-id` or `--channels 'BH?'`. Rerunning
the same command resumes from its checkpoints. Every command takes `--help`
and is also available as `python -m compy_pressure.<module>`.

When `--master-dir` and `--results-dir` are omitted, `exports/` and
`results/` in the current directory are used.

From Python, [`_Example/5_Pressure_QC.py`](../_Example/5_Pressure_QC.py)
runs the inventory, PPSD, plot and QC stages. Without arguments it does so on
a synthetic gauge in a few seconds; with `--data-dir`, `--inventory-dir` and
`--network` it runs on real data. To read the results:

```python
from obspy.signal import PPSD
from compy_pressure.inventory import load_master_dataset

master = load_master_dataset("pressure_run/master")
gauges = master[master.ppsd_status.eq("complete") & master.Seed_id.str.endswith("DH")]
PPSD.load_npz(gauges.ppsd_npz_file.iloc[0], allow_pickle=True).plot()
```

## Outputs

```text
master/master_inventory.{csv,xlsx}     one row per channel epoch, every stage's status and result
results/<net>/<sta>/<epoch>/
    ppsd.npz, inventory.xml, processing.json
    plots/                             PPSD, spectrogram, all PSD curves
    pressure_expected/  nearby/  infragravity/  rayleigh/
results/_workflow/<stage>.json         checkpoints and status counts of each stage
```

Updates to the master are locked between processes, staged, and backed up in
`master/backups/` before each write.

## Origin and updates

The package is exported from BRUIT-FM by its
`integration/export_to_compy.py`. `compy_pressure.BRUITFM_VERSION` and
`compy_pressure.BRUITFM_COMMIT` record the source. The export:

- renames `bruitfm` to `compy_pressure` in every import and message;
- leaves out the lines that BRUIT-FM marks `# bruitfm-only`:
  - the raw-data Godin pilot stage;
  - the checkout-relative default folders;
- reformats the copy with ComPy's ruff settings and checks it for Python 3.10;
- runs the 225 copied tests (`tests/pressure/`).

Until BRUIT-FM is frozen as the campaign record, fix bugs in BRUIT-FM and
export again. A change made only here is overwritten.

```bash
cd ~/Desktop/BRUIT-FM/bruit-fm
python integration/export_to_compy.py --force --run-tests --install-into ../ComPy   # update
python integration/export_to_compy.py --check ../ComPy                              # 0: in sync
```

`--check` also verifies the package entries of ComPy's `pyproject.toml`. Once
BRUIT-FM is frozen, this package becomes the maintained copy and the note in
`compy_pressure/__init__.py` can go.

What stays in BRUIT-FM:

- **The campaign record:** the dated reports, run folders and review scripts
  of the 846-channel campaign, and the raw-data Godin pilot.
- **Tilt correction:** ComPy has its own (`Rotate`).

The BRUIT-FM-only `bruitfm-tilt` command is not exported.

## Open scientific questions

These were put to the supervisor on 29 October 2026 and are not yet answered.
The code runs as it stands; an answer may change a default.

| Question | Code affected |
|---|---|
| Gain-corrected or uncorrected noise model as the result; gain limit 0.1–10 | campaign scripts in BRUIT-FM (no change here) |
| Shelf stations (< 200 m): depth-dependent reference for QC test 1, or a graded score; split the < 1 km depth class | `pressure_test`, `group_noise_model` |
| A separate model for absolute pressure gauges (APGs) | `group_noise_model` |
| Godin Eq. 6 amplitude B<sub>p</sub>: published 2.0·10³ Pa²/Hz or fitted (about 4·10³) | `infragravity_test.PRESSURE_NORMALIZATION` |

## Validation

- The tests in `tests/pressure/` run on synthetic StationXML and waveforms;
  no data are needed. They pass on Python 3.10 (pandas 2.3, NumPy 2.2) and
  3.13 (pandas 3.1, NumPy 2.5).
- Reproduction check (29 October 2026): BRUIT-FM and the installed ComPy wheel
  were run on the same inputs.
  - The QC tests ran on 90 gain-corrected pressure epochs (YV, ZM, 1L, X9).
  - The noise model was rebuilt from the provisional acceptance table.
  - Every output table, JSON file and array matches. Results and method are in
    BRUIT-FM's `integration/README.md`.
- Campaign results (846 pressure channels, 30 networks): BRUIT-FM's
  `docs/REPORT.md`.
