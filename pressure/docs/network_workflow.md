# Network workflow through the screened Godin stage

For the full YV command sequence, with each completed stage on its own line,
see [the README walkthrough](../BRUITFM_README.md#run-every-completed-stage-for-the-yv-network).
The commands below demonstrate automatic discovery across all networks;
the README versions explicitly select `--network YV` for calculations.

`src/bruitfm/network.py` runs the completed calculations for every channel epoch
with waveform files in the selected archive. It discovers network, station,
location and channel identifiers from SDS filenames and validates actual
MiniSEED headers during processing. No YV station list is hard-coded.

## Inputs and outputs

Defaults on this machine:

```text
~/Desktop/Data_Bruit_FM/
  INVENTORIES/                 Matching StationXML files
  Data/                       Waveform archive
    YV/<station>/<channel>.D/ YV example; other layouts are searched recursively
```

SDS filenames must have the usual `NET.STA.LOC.CHA.D.YEAR.DAY` form. Add another
network's waveform files under `Data` and its matching StationXML under
`INVENTORIES`. The same commands discover it automatically. The network code
comes from the filename and trace header, so directory names do not determine
the sensor identity. Missing StationXML is reported explicitly.

Master results remain in `exports/master_inventory.csv` and
`exports/master_inventory.xlsx`. Network products are stored separately:

```text
results/NETWORKS/TEST/
  <network>/<station>/
    <seed-id>__<response-epoch>__<inventory-hash>/
      ppsd.npz
      inventory.xml
      processing.json
      pressure_expected/
      nearby/
      infragravity/
      godin_pilot/<pilot-start>__<days>days/
    plots/
  _workflow/
    plan.json
    ppsd.json
    plots.json
    pressure.json
    nearby.json
    infragravity.json
    pilot.json
```

Existing RR52 products under `results/TEST` remain available. When network
processing completes, the matching master rows link to their new products.
Metadata, manually entered scores and extra Excel sheets are preserved.

## Commands in order

Run these in Terminal, **one command at a time**, waiting for each to finish.
Finish copying the source files before starting calculations; the processing
manifest requires a stable source set for reliable restart.

First enter the project and limit numerical-library threads per worker:

```bash
cd /Users/mohammadamin/Projects/bruit-fm
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
```

**1. Synchronize inventories without resetting existing results.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network inventory
```

This creates the master if needed. Otherwise, it appends newly discovered
source/channel epochs and missing metadata columns, backs up both files and
preserves all existing cells. An unchanged inventory leaves both files
untouched. Existing epoch metadata are not automatically refreshed; review
changes to already registered response epochs explicitly.

**2. Check available networks, stations, channels and matching inventory rows.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network check
```

This saves `plan.json` and prints the current counts and source size. Its
availability check uses filenames; actual header and response validation
happens during computation. Inventory channels without files are listed in
the plan and do not receive invented spectra or calibration failures.
Counts change as new data are added.

The completed-copy check found **43 YV stations, 172 channels/epochs and
45,593 files (370.08 GB)**, with no unmatched waveform identifiers. The
inventory also lists 55 channel identifiers without files in this archive.
The [availability snapshot](yv_archive_plan.json) records both lists.

**3. Calculate PPSDs for all available channels and their response epochs.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network ppsd --workers 4
```

Uses the full available recording period by default: three-hour windows,
50% overlap and the 0.001 Hz target. Sample rate and response handling are
selected per channel. Pressure spectra use Pa²/Hz; seismic spectra use
acceleration units. The full response is applied once. Four independent
channel/epoch processes share no Excel writes; only the parent updates the
master. Each worker receives its channel's file list rather than the entire
archive index. The existing PPSD checkpoint mechanism is retained.

**4. Save the three PPSD plot types for every processed channel.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network plots
```

Creates distributions, time spectrograms and all individual PSD curves.
Seismic distribution and curve plots include NLNM/NHNM. All three paths and
`ppsd_plot_status` are stored in the matching master rows. Figures are closed
after saving so a batch does not accumulate plot windows or large curve arrays.

**5. Compare pressure PPSDs with the Webb and Brown reference limits.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network pressure
```

Uses only pressure-channel epochs, saves the three color-coded comparison
views and appends `pressure_expected_*` results. The existing provisional
reference bounds and nominal response-band choices are unchanged.

**6. Compare pressure spectra with nearby stations at similar depth.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network nearby
```

Reuses saved pressure NPZs, compares common qualifying UTC days and appends
`nearby_*` diagnostic fields and file links to CSV/Excel. Saves spectral
comparisons, daily residual heatmaps, neighbor-selection plots and candidate
audits. Preferred screening uses 50 km/5% depth difference and matching declared
sensor labels; exploratory screening allows 100 km/10%. At least two distinct
reference stations and 30 jointly supported days per frequency are required.
The reference pool includes all registered ready pressure products, even if
targets are restricted to one network. These settings are configurable proposals;
unsupported stations remain unresolved and scores stay blank.
See [method and options](nearby_comparison.md).

Finish the PPSD batch first; this stage refuses master writes while its workflow
summary is unfinished. Diagnostics run sequentially using cached daily spectra.

**7. Calculate the depth-dependent Godin comparison from saved pressure PPSDs.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network infragravity
```

Uses each row's own depth estimate and response band. Saves model, PPSD-overlay
and residual plots and `infragravity_*` diagnostics. A channel with no eligible
full smoothing bins is recorded as `no_comparable_bins`, not as a failed
calibration. Missing depth or invalid responses produce an explicit processing
error while other epochs continue.

**8. Calculate the screened full-window Godin pilot for every pressure epoch.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network pilot
```

Selects a seven-day interval near each epoch's middle, avoiding deployment
startup data. Uses three full 12.42-hour Hann spectra per day, count-based
50–60 mHz screening, full-response correction and linear-power averaging.
Saves raw and corrected spectra, decisions, comparisons and three plots;
appends `infragravity_pilot_*` results. Multiple locations and response epochs
are handled independently. If the epoch cannot supply a week, the error is
recorded; `--pilot-days` can select a shorter interval deliberately.
When no matching PPSD frequency bins exist, a header-only comparison CSV is
saved and the plot states that limitation.

These commands process every network present in `Data`. Add `--network YV`
to commands 2–8 to restrict a run to YV. Add `--station RR01` or a complete
`--seed-id` to limit the scope further. Repeated network/station/SEED options
select several values. Use `--data-dir`, `--inventory-dir`, `--master-dir` and
`--results-dir` when moving the workflow to another disk or checkout.

## Open saved plots

Network batch commands save figures without opening hundreds of windows.
After the relevant stage, open one station's saved plots:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR01 --view-kind ppsd
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR01 --view-kind pressure
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind nearby
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR01 --view-kind infragravity
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR01 --view-kind pilot
```

Replace RR01 with the desired station. `view` opens existing PNGs and changes
neither spectra nor master results.

## Restart, errors and one-command operation

Rerun a stage's command to continue. A stage reuses an epoch only when its code,
settings, scientific inputs, source fingerprint and expected output files
agree. Otherwise it retries that epoch. Failure in one channel does not stop
other channels. Failed/partial results return a nonzero exit status, with
details in the stage JSON and matching master status/error fields.
Successful retries clear the previous processing error.

For PPSDs, Ctrl-C requests that workers stop after the current bounded time
chunk and save a checkpoint. Pending jobs are cancelled; completed master
updates are saved. For plots and diagnostics, finished epochs are checkpointed
and the interrupted epoch is retried on the next run.

`--force` retries saved stages while retaining PPSD checkpoint resume. If the
waveform source set or PPSD settings changed, an existing checkpoint refuses to
mix them. `--rebuild-ppsd` explicitly recalculates the selected network PPSD
products; existing original RR52 files in the separate directory are retained.
Keep alternative configurations in a separate `--results-dir` when comparing
them. Run only one batch writer against a master pair at a time.

After the input check, the complete sequence can also be run as one command:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network all --workers 4
```

This runs inventory → check → PPSDs → plots → pressure limits → nearby comparison → Godin comparison
→ Godin pilot → paired pressure–vertical spectra (seven raw days per pressure epoch;
see [paired spectra](paired_spectra.md)). The list below predates the paired stage:
→ Godin pilot. The same automatic discovery applies to subsequently added
networks. Individual station scripts remain available for detailed settings
and interactive work.

## Validation and scientific status

A workflow recheck passed the full 62-test suite and style/format checks.
All 14 commands in the README YV walkthrough were parsed against the actual
CLI without executing production calculations. A read-only check of 156 saved
YV pressure/seismic PPSD products found matching channel identifiers and sample
rates, finite spectra, unique timestamps and the required NPZ/response/provenance
files. This checks the products available at inspection time; it does not
establish complete archive coverage or validate calibration acceptance.
Production spectra and masters were left untouched by this review.

The suite passes **62 tests**. Network checks cover two discovered networks,
multiple pressure locations and response epochs, two PPSD workers, cached
restart, corrupted-channel isolation and retry, cooperative stopping,
inventory append preservation, missing prerequisites and saved-plot viewing.
Nearby checks additionally verify paired daily offsets, distinct-station counts,
joint coverage, exploratory outcomes, reference-dependent restart and CSV/Excel
preservation. [Isolated YV comparison](nearby_validation_summary.json) verifies
five completed pressure products without changing the production master.

A real-data verification used one interior week from RR01 pressure and vertical
motion, RR52 pressure and RR41 pressure, covering **50, 62.5 and 100 Hz**.
Each channel produced **110 PPSD windows**, all plotting stages completed,
and the three pressure pilots each retained **20 of 21** Hann windows.
The saved-PPSD Godin comparison correctly found no fully eligible smoothing
bins for RR41, while its finer raw pilot grid provided eligible frequencies.
This bounded verification used a separate master under
`results/network_validation`; the project master and original RR52 products
were preserved. The full archive has not been computed by this verification.

The [verification summary](network_validation_summary.json) records the
channel settings and stage outcomes. Example verified plots:
[RR01 seismic PPSD](figures/network_validation/RR01_BHZ_ppsd.png),
[RR41 empty-bin Godin diagnostic](figures/network_validation/RR41_BDH_godin_residuals.png)
and [RR41 full-window pilot](figures/network_validation/RR41_BDH_godin_pilot.png).

Some YV response evaluations emitted FIR-normalization warnings from evalresp.
The workflow retains ObsPy's existing response evaluation and records the
source inventories; these messages should be reviewed with the instrument
metadata when interpreting calibration.

The calculations are diagnostics. No new calibration score, correction factor
or modified StationXML response is assigned. All prior scientific limitations
of the [pressure comparison](pressure_comparison.md),
[nearby comparison](nearby_comparison.md),
[Godin comparison](infragravity_comparison.md) and [pilot](godin_pilot.md) apply
to every station, including depth and usable-frequency validation.

## Provisional Rayleigh calibration stage

The explicit `network rayleigh` stage now runs the working earthquake method
for available pressure epochs, using TiSKitPy and the paper's 0.03–0.07 Hz band.
It checkpoints stations and appends only `rayleigh_*` fields; scientific scores
and input responses remain unchanged. It is not included in `network all`.
[Commands, current partial execution and saved plots](yv_rayleigh_calibration.md)
and [assumptions/checklist](calibration_assumptions.md) distinguish a completed
calculation from independently validated field calibration.
