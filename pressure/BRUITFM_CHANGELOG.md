# Changelog

All notable changes to BRUIT-FM. The project follows
[Semantic Versioning](https://semver.org/): the major number changes when the
command line, the master-table columns or saved-file formats change in an
incompatible way.

## [1.3.0] — 2026-10-30

The last BRUIT-FM release: everything moves into ComPy, which becomes the
complete package; this repository is kept as a backup. `integration/` exports
all of BRUIT-FM as ComPy's `compy_pressure` package and its `pressure/`
project folder. Only the dated reports and the export's own tests stay here.
Scientific results are unchanged: the reproduction check below found every
output identical.

### Added
- `export_to_compy.py`: a full move.
  - The 40 modules go into `compy_pressure/`; `tests/`, `scripts/` and
    `docs/` go into ComPy's `pressure/`, laid out like this checkout.
    Every relative path still holds there.
  - Lines and blocks marked `# bruitfm-only` are left out.
  - Every `bruitfm.` mention is renamed.
  - References to the package's source become `../compy_pressure`, and
    package modules reach `pressure/` from their new place.
  - Markdown links are pointed at the new places, or at this release for the
    dated reports.
  - Data files are copied byte for byte.
  - The copy is reformatted and linted with ComPy's ruff settings for Python
    3.10.
  - The source version and commit are recorded in
    `compy_pressure/__init__.py`.
  - `--python` runs the moved tests on another interpreter.
  - `--install-into COMPY` replaces ComPy's copy, only folders an earlier
    export made.
  - `--check COMPY` reports any drift, including missing `pyproject.toml`
    and pytest settings.
- `tests/test_core_boundary.py` tests the export: what moves and what stays,
  the markers, the layout rules, the link rewriting, and Python 3.10 syntax
  and standard-library names in every moved file.
- `bruitfm.common.PROJECT_DIRECTORY`: one definition for the default
  `exports/` and `results/` folders of the core stages. It is the current
  directory in ComPy; BRUIT-FM still uses the checkout.

### Changed
- All code, tests and scripts run on Python 3.10, ComPy's lowest supported
  version:
  - no `datetime.UTC`;
  - no `hashlib.file_digest`.
- The master-table lock works on Windows too (`msvcrt`). It is unchanged on
  macOS and Linux. `coverage_run` imports on Windows and refuses there to
  run a campaign (it needs `fcntl` locks).
- `network`: an unknown diagnostic stage now raises instead of running the
  Godin pilot.
- `local_p_check` provenance: `code.src_or_scripts_dirty` is now
  `code.dirty` with `code.dirty_paths` (the paths checked). In ComPy these
  are the package's own folder. No code reads this field.
- Manifest: everything moves except `[stays]`; the extras have version
  floors (pandas 2.2, openpyxl 3.1).

### Reproduction check
- 90 gain-corrected pressure epochs (YV, ZM, 1L, X9), QC tests 1–3 with the
  28 October parameters: BRUIT-FM and the installed ComPy wheel give the same
  output.
- The gain-corrected noise model rebuilt by both from
  `decisions/acceptance.csv` is identical.
- In the moved layout, all 503 moved tests pass on Python 3.10 and 3.13.
- Details: `integration/README.md`.

### Added before the export (integration preparation)
- `integration/`: preparation for the move into ComPy. `compy_manifest.toml`
  lists the 22 core modules, their tests and what stays;
  `export_to_compy.py` exports the core as a stand-alone package (default
  `compy_pressure`) and runs its tests; `compare_outputs.py` compares two
  result folders for the reproduction check; `README.md` is the handoff with
  the owner's decisions and the supervisor's open questions.
- `tests/test_core_boundary.py`: the core may not import or name campaign-only
  modules.
- `docs/reports/2026-10-29_noise_model_report/`: the noise-model report sent to
  the supervisor, with a Godin Eq. 6 comparison with and without tanh(kH)^(−3/2).
  Version 2 adds Figure 15: Eq. 6 with the published B_p = 2.0e3 Pa²/Hz and
  with B_p fitted per depth class (2.3–5.0e3 Pa²/Hz).

### Changed before the export
- `complete_samples` and `write_csv` moved from `godin_pilot` to the new
  `bruitfm.common`; `station_distribution` moved from `noise_model` to
  `plot_network_ppsd`. The old import paths still work.
- `bruitfm.network` imports the Godin pilot only when the "pilot" stage runs.
- The paired-spectra provenance hashes `common.py` instead of `godin_pilot.py`.
- Tests load the reference data through the package instead of the
  repository path, and the script-based microseism tests skip when the
  campaign scripts are absent.
- `docs/reports/2026-10-28_gain_corrected/`: the whole campaign redone after
  multiplying each pressure channel by its calibration gain (Rayleigh, then
  microseism, then local P; same-gauge gains shared; 0.1–10 limit). It holds
  the corrected QC tests, final classes, metadata flags, acceptance, noise
  model and network figures, with a full report in its README.

## [1.2.0] — 2026-10-28

### Added
- Noise-model figures: a heat map of the station-balanced PSD probability density
  with the percentile curves (`noise_model_<grouping>_heatmap.png`), and the
  median PSD of every station used as a line (`noise_model_<grouping>_lines.png`).
  The NPZ files now also hold the full distribution (`frequency_edges_hz`,
  `db_edges`, `station_balanced_probability`).
- Provisional steps 1–3 in `docs/reports/2026-10-28_provisional_decisions/`:
  a frozen provisional policy, the level offset of every pressure channel from
  comparable good gauges with network metadata flags (2A, YD, Z4, X9; 4N), the
  accepted/review/rejected list (396 accepted at 311 stations), and the noise
  model rebuilt from the accepted stations. `docs/REPORT.md` §8–§10 are updated.

## [1.1.0] — 2026-10-28

### Added
- `bruitfm.group_noise_model` (`bruitfm-group-noise-model`): station-balanced
  percentile pressure noise models (p5, p10, p50, p90, p95) from a scored
  channel table, grouped by ocean basin, sensor type, water depth and local
  season. It keeps one channel per station and sensor, uses only the usable
  response band and drops dead windows.
- Provisional noise model from the 352 stations classed *good*, in
  `docs/reports/2026-10-28_noise_model_provisional/` and §8 of
  `docs/REPORT.md`.
- GPL-3.0 license (the same as ComPy, into which BRUIT-FM will be merged).
- The historical formula/signal-processing review (YV checks, gain plots and
  scripts) and the local-P science review, merged from their archived branches.
- The calibration results and a worked example in `docs/REPORT.md`, with one
  figure and caption for each QC test and calibration method.

## [1.0.0] — 2026-10-28

First public release: the full pipeline from a raw SDS drive to a scored table
of every pressure channel, used for the all-network campaign described in
[docs/REPORT.md](docs/REPORT.md).

### Added
- **Whole-drive processing.** Networks, stations and channel epochs are
  discovered automatically. PPSD, plot and test stages run in parallel
  (`--workers max`), are checkpointed and resume after an interruption. A
  `--channels` filter selects channel codes with wildcards (PR #8, #9, #13).
- **QC references.** Webb (1998) limits at ≤ 0.2 Hz and digitised Brown et al.
  (2014) limits at ≥ 0.2 Hz for Test 1. Godin et al. (2013) Eq. 6 with and
  without the tanh(kH) factor for Test 3. The nearby-station test also reports
  the RMS of daily differences (PR #10).
- **Station maps with coastlines** through `bruitfm.maps.add_land`, with the
  optional `maps` extra (PR #12).
- **Tilt correction.** `bruitfm.tilt` computes ComPy (`compy.Rotate`)
  tilt-cleaned vertical PPSDs in 6-hour windows at 1 Hz and saves the
  uncorrected and corrected PSDs side by side (PR #13).
- **Calibration runners.** The YV all-methods runner (Rayleigh, microseism,
  local P) and the network Rayleigh stage, with bounded event sampling and a
  choice of pressure sign (PR #4).
- **Network figures.** Median PSD and PPSD-distribution plots by network
  (PR #3).
- **Packaging.** Console commands (`bruitfm`, `bruitfm-inventory`,
  `bruitfm-tilt`, `bruitfm-local-p`, `bruitfm-microseism`, `bruitfm-paired`,
  `bruitfm-noise-model`), `bruitfm.__version__`, the `all` and `tilt` extras,
  project URLs, `CITATION.cff`, `CONTRIBUTING.md` and this changelog.
- **Reports.** Dated campaign reports with the scripts that built them, in
  `docs/reports/`.

### Changed
- The README now covers installation and usage. The step-by-step development
  history it used to hold is in
  [docs/development_log.md](docs/development_log.md); handoff notes moved to
  `docs/internal/`.
- The Rayleigh code hash is recorded at import time, not at write time (PR #5).

### Fixed
- Verticals coded `3` (dip ±90°) are accepted by the tilt and calibration
  vertical selection.
- PPSDs of tilt-corrected streams use a complete unit response dictionary.

## [0.1.0] — 2026-10-21

Reviewed deployment checkpoint (PR #1): inventory, SDS reader, response-corrected
PPSDs and plots, the three pressure QC tests, the Godin pilot, paired
pressure/vertical spectra, Rayleigh, microseism and local-P calibration
diagnostics, coverage tools and approval-gated noise-model code. It was
validated on YV, YL, 7D and XO.

[1.3.0]: https://github.com/MohammadAmin-Aminian/bruit-fm/releases/tag/v1.3.0
[1.2.0]: https://github.com/MohammadAmin-Aminian/bruit-fm/releases/tag/v1.2.0
[1.1.0]: https://github.com/MohammadAmin-Aminian/bruit-fm/releases/tag/v1.1.0
[1.0.0]: https://github.com/MohammadAmin-Aminian/bruit-fm/releases/tag/v1.0.0
[0.1.0]: https://github.com/MohammadAmin-Aminian/bruit-fm/releases/tag/v0.1.0
