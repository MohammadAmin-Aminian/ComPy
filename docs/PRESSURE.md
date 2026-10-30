# Pressure gauges: QC, calibration and noise models (`compy_pressure`)

ComPy holds all of
[BRUIT-FM](https://github.com/MohammadAmin-Aminian/bruit-fm). BRUIT-FM checks
and calibrates ocean-bottom pressure gauges across a whole waveform archive.
It starts from an SDS archive and its StationXML and ends with one master
table that scores every pressure channel. BRUIT-FM 1.3.0 was its last
release; that repository is kept as a backup and holds the dated reports of
its 846-channel campaign.

| Where | What |
|---|---|
| `compy_pressure/` | the package (installed with ComPy) |
| `pressure/tests/` | its 503 tests, run by ComPy's `pytest` |
| `pressure/scripts/` | BRUIT-FM's run and review scripts |
| `pressure/docs/` | BRUIT-FM's methods, validation notes, reviews, policies and their data; the results report is [`pressure/docs/REPORT.md`](../pressure/docs/REPORT.md) |
| [`pressure/README.md`](../pressure/README.md) | the BRUIT-FM → ComPy name table |

`compy_pressure`'s Rayleigh-wave calibration replaces ComPy's old
`Pressure_calibration.py`, which was removed after ComPy 2.4.0. The README section "Pressure
calibration" has a migration table.

## Installation

```bash
python -m pip install -e '.[pressure]'        # pandas and openpyxl
python -m pip install -e '.[pressure,maps]'   # plus Cartopy, for station maps
```

TiSKitPy and disba, which the Rayleigh and microseism checks use, are ComPy
dependencies. Python 3.10 to 3.13 is supported and tested.

## Inputs

- **Waveforms** in an [SDS archive](https://www.seiscomp.de/seiscomp3/doc/applications/slarchive/SDS.html):
  `<root>/<year>/<net>/<sta>/<chan>.D/<net>.<sta>.<loc>.<chan>.D.<year>.<doy>`.
  A file name with location `00` also matches a StationXML location `""`.
- **StationXML** with full response stages, in one file or a folder. A
  channel without response stages is reported, never guessed.
- For the Rayleigh calibration, an earthquake catalogue in USGS ComCat
  GeoJSON; for the local-P check, a local catalogue (ComCat CSV or GeoJSON).

## What it computes

| Step | Result | Command | Module |
|---|---|---|---|
| Inventory | Master table (CSV + Excel) of every channel epoch: coordinates, depth, sensor, response, sample rate, coverage | `compy-pressure inventory` | `inventory` |
| PPSDs | ObsPy PPSDs with the full response removed: 3 h windows, 50 % overlap, down to 1 mHz | `compy-pressure ppsd`, `plots` | `ppsd`, `plot_ppsd` |
| QC test 1 | Fraction of median PSD bins inside the Webb (1998) and Brown et al. (2014) seafloor limits | `compy-pressure pressure` | `pressure_test` |
| QC test 2 | Median daily difference from nearby gauges at a similar depth | `compy-pressure nearby` | `nearby_test` |
| QC test 3 | Residual from the depth-dependent Godin et al. (2013) infragravity model | `compy-pressure infragravity` | `infragravity_test` |
| Calibration 1 | Gain from teleseismic Rayleigh-wave pressure/vertical ratios, with optional TiSKitPy tilt cleaning | `compy-pressure rayleigh` | `rayleigh_calibration`, `rayleigh_network` |
| Calibration 2 | Gain bracket from the water-layer admittance of coherent microseisms | `compy-pressure paired`, then `compy-pressure-microseism` | `paired_spectra`, `microseism_check` |
| Calibration 3 | Gain from the pressure/vertical ratio of local P arrivals | `compy-pressure-local-p` | `local_p_check` |
| Tilt | Tilt-corrected vertical PPSDs (ComPy's `Rotate`) for every three-component station | `compy-pressure-tilt` | `tilt` |
| Noise models | Station-balanced percentile curves (p5–p95) by ocean basin, sensor type, depth and local season | `compy-pressure-noise-model`; approval-gated: `compy-pressure-approved-noise-model` | `group_noise_model`, `noise_model` |
| Campaign tools | Whole-drive coverage runs, quality review, raw-data Godin pilot | `python -m compy_pressure.coverage`, `.quality_review`, `compy-pressure pilot` | `coverage*`, `quality_review`, `godin_pilot` |
| Figures | Network PPSD medians by depth, sensor and basin; station maps | — | `plot_network_ppsd`, `maps` |

Nothing is applied silently. The calibration checks estimate a gain and record
it; they never change a response or the data. Each stage writes its inputs,
settings and code fingerprint next to its result. Reruns reuse a finished
result only when all of these still match.

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

`compy-pressure all` runs these stages in order: inventory, check, PPSD,
plots, the three QC tests, the Godin pilot and paired spectra. The Rayleigh
calibration runs only when asked. Narrow a run with `--station`, `--seed-id`
or `--channels 'BH?'`. Rerunning the same command resumes from its
checkpoints. Every command takes `--help` and also runs as
`python -m compy_pressure.<module>`. Without `--master-dir` and
`--results-dir`, the outputs go to `exports/` and `results/` in the current
directory.

Two examples run from Python:

- [`_Example/5_Pressure_QC.py`](../_Example/5_Pressure_QC.py) runs the
  inventory, PPSD, plot and QC stages. With no arguments it uses a synthetic
  gauge and finishes in seconds.
- [`_Example/3_DPG_Calibration.py`](../_Example/3_DPG_Calibration.py) runs the
  Rayleigh calibration of one station. It needs no PPSDs:
  `rayleigh_network.calibration_master` uses each channel's StationXML as its
  response snapshot.

## The Rayleigh-wave calibration

The method has five steps:

- **Events.** Teleseisms at most 150 km deep from the catalogue, by magnitude
  level, largest first. Smaller levels, down to Mw 6, are used only when too
  few events are accepted.
- **Windows.** One window per event, centred on the Rayleigh-wave group
  arrival predicted by the dispersion file.
- **Spectra.** The full responses are removed from the native-rate pressure
  and vertical channels. With `tilt_method="tiskitpy"` (the default), tilt
  noise is removed from the vertical using the horizontals. This needs
  horizontals with valid orientations. `"none"` uses pressure and vertical
  only.
- **Fit.** In the band (default 0.02–0.06 Hz; 0.03–0.07 Hz for
  `compy-pressure rayleigh`), the bins whose magnitude-squared coherence is at
  least 0.95 are compared with the water-layer prediction for the station
  depth. The prediction uses the phase velocity, a homogeneous 1028 kg/m³
  water column and a sound speed of 1500 m/s. A gain is fitted analytically
  per event and checked against a 0.1–3 grid.
- **Station gain.** The median of at least three accepted events, with their
  scatter (`event_multiplier_std`, a scatter, not a confidence interval).
  `comparison.json` records every event's status, gain, spectra and figures,
  with review flags for the assumptions: pressure sign, depth from the
  StationXML, Earth model and water properties.

The dispersion file gives the Rayleigh phase and group velocities; RR52's,
from a layered model, is the BRUIT-FM default for every network.
`--rayleigh-earth-model-file` computes them from another layered model with
disba.

### Pipeline stage D

`pipeline/d_calibrate.py` (stage D of the end-to-end pipeline) runs this
calibration on the station's SDS waveforms (`config.SDS`) and StationXML
(`config.INVENTORY`). It first converts the saved USGS QuakeML catalogue to
GeoJSON. It writes the same outputs as in 2.4:

- `work/<STA>/calibration/compy_gain.json`: `gain`, `gain_std`,
  `event_gains`, `events_used`, `band_hz`, `depth_m`, or `error`;
- `<STA>_compy_calibration.npz` and `.png`.

Stages C and G are therefore unchanged; stage F's calibration figure accepts
the new frequency grid. Three settings in `config.py`
control it:

- `CALIBRATION_BAND`: default 0.02–0.06 Hz;
- `CALIBRATION_TILT`: `"none"` for XE, whose StationXML gives BH1 and BH2 a
  dip of −90°; `"tiskitpy"` for networks with valid horizontals;
- `CALIBRATION_DISPERSION`.

The full run is in `work/<STA>/calibration/compy_pressure/`.

Check against the 2.4 stage D (`Pressure_calibration`), 30 October 2026. The
new stage D ran on the original SDS data and StationXML of the XE run, in a
separate folder:

| Station | New stage D (`compy_pressure`) | 2.4 stage D | Difference |
|---|---|---|---|
| XE.CC11 | 0.711 ± 0.011 (14 events) | 0.698 ± 0.012 (12) | +1.9 % |
| XE.WC03 | 0.790 ± 0.010 (9) | 0.779 ± 0.008 (9) | +1.4 % |
| XE.EE04 | 1.165 ± 0.015 (12) | 1.146 ± 0.015 (12) | +1.7 % |

The new gains are consistently 1.4–1.9 % higher. That is within the ±5 % gain
prior of the inversion (stage E). The likely cause is the Rayleigh phase
velocity behind the water-layer prediction: RR52's dispersion here, against
the old module's own model. BRUIT-FM's microseism method gave 0.702, 0.778 and
1.159 for these gauges. These XE gauges are BRUIT-FM's first Rayleigh-wave
results: its campaign could not pair their channels because of the wrong
horizontal dips.

## Outputs

```text
master/master_inventory.{csv,xlsx}     one row per channel epoch, every stage's status and result
results/<net>/<sta>/<epoch>/
    ppsd.npz, inventory.xml, processing.json
    plots/                             PPSD, spectrogram, all PSD curves
    pressure_expected/  nearby/  infragravity/  rayleigh_calibration/
results/_workflow/<stage>.json         checkpoints and status counts of each stage
```

Updates to the master are locked between processes, staged, and backed up in
`master/backups/` before each write.

## Running BRUIT-FM's scripts and tests

`pressure/` is laid out like BRUIT-FM's checkout. Its scripts and tests
address `docs/`, `scripts/` and `tests/` relative to it, and the package
source as `../compy_pressure`. Run them from that folder, with ComPy and it
on the path:

```bash
cd pressure && PYTHONPATH=..:. python scripts/plot_network_median_ppsd.py --help
python -m pytest pressure/tests        # from the ComPy checkout
```

The shell scripts expect the conda environment named in each file. The
documents in `pressure/docs/` are BRUIT-FM's, unchanged apart from their
links, and still use BRUIT-FM's names; `pressure/README.md` maps them to
ComPy's.

## Origin

The move was made by BRUIT-FM's `integration/export_to_compy.py` from BRUIT-FM
1.3.0. `compy_pressure.BRUITFM_VERSION` and `compy_pressure.BRUITFM_COMMIT`
record the source. The export:

- renamed `bruitfm` to `compy_pressure`;
- adapted the paths to this layout;
- pointed the Markdown links at their new places, or at BRUIT-FM v1.3.0 for
  the dated reports;
- copied data files byte for byte;
- formatted and linted the copy with ComPy's ruff settings.

ComPy is the maintained copy from now on. BRUIT-FM's
`export_to_compy.py --check` would now report ComPy's later changes (for
example `calibration_master`) as differences.

Reproduction check (29–30 October 2026): BRUIT-FM and the installed ComPy
wheel ran on the same inputs. They ran QC tests 1–3 on 90 gain-corrected
pressure epochs (YV, ZM, 1L, X9) and rebuilt the gain-corrected noise model.
Every output table, JSON file and array matches, and the per-epoch QC tests 1
and 3 also match the archived campaign results. The method and numbers are
in BRUIT-FM's `integration/README.md`.

## Open scientific questions

These were put to the supervisor on 29 October 2026 and are not yet answered.
The code runs as it stands; an answer may change a default. Record the answers
here.

| Question | Code affected | Answer |
|---|---|---|
| Gain-corrected or uncorrected noise model as the result; gain limit 0.1–10 | campaign scripts and run folders | |
| Shelf stations (< 200 m): depth-dependent reference for QC test 1, or a graded score; split the < 1 km depth class | `pressure_test`, `group_noise_model` | |
| A separate model for absolute pressure gauges (APGs) | `group_noise_model` | |
| Godin Eq. 6 amplitude B<sub>p</sub>: published 2.0·10³ Pa²/Hz or fitted (about 4·10³) | `infragravity_test.PRESSURE_NORMALIZATION` | |
