# BRUIT-FM handoff

**Fresh chat: start with [CODEX_NEXT_CHAT.md](CODEX_NEXT_CHAT.md).** README progress
runs through step 46. Reviewed engineering fixes and optional scientific trial
tools are merged into main; original calibration defaults stay unchanged.
The completed 89-station trial yields 61 provisional estimates, including C07W's
three-event normal-state subset, versus 46 preserved baseline gains.
173 tests pass; no gains, responses, masters, saved results or scores were changed.
[Comparison](../socrates_all_station_comparison.md) ·
[Contributions](../contributions.md). Earlier notes below are historical context.

Resume this project in the local `bruit-fm` checkout.
Python environment: conda `bruitfm` (Python 3.13, ObsPy 1.5.1).
See the [completed-step log](../research_log.md) and
[processing details](../rr52_methods.md) for the full settings and results.
The README lists forty-six completed steps in the original numbered format,
followed by exact YV stage commands and saved-plot viewing commands. It distinguishes
completed RR52 results from whole-network execution progress in `_workflow`.

## Data and current results

- Waveforms on the current machine: `~/Desktop/Data_Bruit_FM/RR52`.
- Matching StationXML: `Data_Bruit_FM/INVENTORIES/YV.xml`.
- Four channels: `YV.RR52.00.BDH` (DPG pressure), `BH1` and `BH2`
  (horizontal motion), and `BHZ` (vertical motion).
- Master metadata: `exports/master_inventory.csv` and `.xlsx`, containing
  16,784 channel-epoch rows. Processing updated the four RR52 rows and preserved
  original metadata and calibration results. Master-file backups are in
  `exports/backups/`.
- Full available RR52 data were processed: 5,249 PSD windows per channel.
- NPZ spectra, per-epoch response snapshots, and processing provenance are in
  `results/TEST/YV/RR52/`, with paths recorded in the master table.
- Original waveform files are read as digital counts. ObsPy PPSD applies the
  full response internally once. Pressure PSDs reference `Pa^2/Hz`; seismic
  PSDs reference acceleration `m^2/s^4/Hz`. Stored PSD amplitudes are in dB.

## Processing settings and timing

- PSD windows: 10,800 seconds (3 hours).
- Overlap: 50%, so a new PSD begins approximately every 5,400 seconds (90 minutes).
- Each time-plot column represents a 3-hour average, despite 90-minute spacing.
- Sampling rate: 62.5 Hz; Nyquist: 31.25 Hz.
- Frequency target: 0.001 Hz. Lowest saved period-bin center corresponds to
  approximately 0.000954 Hz.
- Internal FFT length: 131,072 samples, approximately 35 minutes, with ObsPy's
  Welch averaging inside the 3-hour windows.
- Histogram bins: -250 to 100 dB, in 1 dB steps; the wider lower bound preserves
  the very quiet vertical spectra.
- Final calculation used four worker processes, one per channel. Time chunks
  within each channel run sequentially; only the parent writes CSV/Excel.
- Plot generation currently runs sequentially.

## Code and commands

- `src/bruitfm/inventory.py`: build the StationXML master table and export it;
  update matching rows with processing results.
- `src/bruitfm/archive.py`: index daily SDS MiniSEED headers and read bounded
  time intervals, preserving gaps and handling overlapping daily files.
- `src/bruitfm/ppsd.py`: calculate response-corrected spectra, checkpoint and
  resume, save NPZ/provenance, and update the master files.
- `src/bruitfm/plot_ppsd.py`: plot saved distributions, time spectrograms, and
  all individual PSD curves. All 5,249 curves are included for each channel.
- Seismic distributions and curve plots include Peterson NLNM/NHNM reference
  models. Pressure plots use pressure units. Spectrograms show frequency versus
  UTC window start time, with PSD amplitude represented by color.
- Twelve PNGs (three per channel) are in `results/TEST/YV/RR52/plots/`.
- `src/bruitfm/pressure_test.py` now compares saved BDH spectra with supplied
  Webb points and digitized Brown hydroacoustic curves. Run
  `conda run -n bruitfm python -m bruitfm.pressure_test` to display/save three plots
  and update master diagnostics. Details are in [pressure comparison](../pressure_comparison.md).
- Pressure-comparison CSV/JSON/PNG files are beside the NPZ under
  `pressure_expected/`. Status is `comparison_only`; no calibration score is assigned.
- `src/bruitfm/infragravity_test.py` implements the Godin depth-dependent
  pressure prediction. Run `python -m bruitfm.infragravity_test` in the `bruitfm`
  environment to display/save three model comparison plots.
  CSV/JSON/PNG files are under `infragravity/`, with diagnostics in master CSV/Excel.
  RR52 has two supported bins and +16.13 dB median residual under default choices;
  status is `comparison_only`, with the third test score blank.

From the project directory, display/save all plots:

```bash
conda run -n bruitfm python -m bruitfm.plot_ppsd
```

Select only one type, if desired:

```bash
conda run -n bruitfm python -m bruitfm.plot_ppsd --plot-type spectrogram
```

Resume the existing calculation with its current settings:

```bash
conda run -n bruitfm python -m bruitfm.ppsd \
  --data-dir ~/Desktop/Data_Bruit_FM/RR52 --workers 4
```

Use the plotting script for viewing existing results. Running the inventory
builder again creates fresh exports, so preserve the current exports or choose
another output directory when rebuilding inventory metadata.

## Verification and implementation details

- Four saved NPZs reload correctly; spectra are finite, timestamps unique, and
  no values fall outside the widened histogram range.
- Original CSV values and Excel cells were compared with the original backups.
- Source waveform sizes/mtimes and the source StationXML hash remained unchanged.
- The expanded test suite passed 62 tests; Ruff checks and formatting passed.
- Plot verification checked every curve against saved PSD values, frequency/time
  alignment, and seismic-only reference models. PNGs were visually reviewed.
- Restart uses a correctly initialized PPSD followed by `add_npz()`, addressing
  an ObsPy 1.5.1 response-preloading problem with `load_npz(metadata=...)`.
- Fractional-sample phases at daily boundaries are accommodated by advancing
  read boundaries to avoid skipped windows. Source timestamps are unchanged on
  disk; in-memory merging follows ObsPy's sample-grid rules.
- Five commits record Steps 6–10; Git history and the research log preserve the
  implementation and validation details.

## Screened full-window Godin pilot

`src/bruitfm/godin_pilot.py` calculates one raw-pressure pilot near mid-recording
by default. Three full 12.42-hour Hann spectra per day, screening before
response correction in 0.05–0.06 Hz, and linear-power means follow the source
procedure. Segment-mean removal and periodic Hann are explicit choices.

RR52 supplied 21 complete windows; 20 retained, one rejected. The FFT has
2,794,500 samples and 0.0223654 mHz spacing. The primary retained mean has
a median +17.79 dB model residual in the nominal eligible band. A matched
diagnostic compares against 110 same-interval three-hour PPSD windows.
The two eligible matched bins are 1.79 and 2.32 dB lower than that PPSD mean.
Neither diagnostic establishes instrument calibration acceptance.

Pilot NPZ/CSV/JSON/PNG outputs are under the pressure epoch's separate
`godin_pilot/<interval-start>__<days>days/` folder. Master CSV/Excel row 13834
receives new `infragravity_pilot_*` fields. All old fields/scores remain intact.
Open saved plots with `python -m bruitfm.godin_pilot --plot-only` in `bruitfm`.
Calculate another interval with `--start` and `--days`. Same-interval reruns
replace that pilot folder; different intervals get separate folders.
See `docs/godin_pilot.md` for details.

## Whole-network workflow

`python -m bruitfm.network` now runs `inventory`, `check`, `ppsd`, `plots`,
`pressure`, `nearby`, `infragravity`, `pilot`, `view`, or `all`. Defaults discover networks
in `~/Desktop/Data_Bruit_FM/Data`, using `~/Desktop/Data_Bruit_FM/INVENTORIES` and
the existing `exports` master. No network/station/channel list is hard-coded.
Only actual files are processed; unmatched inventory identifiers are reported.

New network outputs are under `results/NETWORKS/TEST/<net>/<sta>/`, keeping
original RR52 files under `results/TEST` available. Batch processing links
matching master rows to its new products. Scores and extra sheets are preserved.
`inventory` appends missing source/channel epochs; it never resets existing
results. Four PPSD workers are the default; diagnostics/plots run sequentially.
Parent-only master writes and saved-stage reports permit restart and retry.

The Godin pilot remains seven days near each epoch's middle; PPSDs use the full
available period. A missing week/depth/response is recorded explicitly.
`view --station RR01 --view-kind pilot` opens saved plots without computation.
All exact commands are in `docs/network_workflow.md`.

Real-data verification in isolated `results/network_validation` used four
channels at 50/62.5/100 Hz, 110 PPSD windows each and three pressure pilots.
All three retained 20/21 Hann windows. All 62 tests pass. This isolated verification
preserved the project master and did not compute the full archive.
The subsequent production YV run completed all 172 PPSDs and base plot jobs,
all 43 expected-pressure/nearby/Godin comparisons, and attempted all 43 pilots.
Two pilot jobs failed: RR12 had no continuous windows in its selected interval;
RR42's interval was outside waveform coverage. These are processing issues,
not calibration failures. See [workflow snapshot](../yv_workflow_status.json).
The user finished copying the archive. The final preflight found 43 stations,
172 channels/epochs and 45,593 files (370.08 GB), all matched to inventory.
Another 55 inventory channel identifiers have no files in the provided folder.
Rerun `check` whenever adding more files or networks.

## Nearby-station comparison

Step 17 implements `nearby_test.py` and `network nearby`. It reads saved pressure
spectra once, builds covered UTC-day medians and compares supported bands to
distinct-station peers. Defaults: preferred 50 km/5% depth difference and matching
declared type/model; exploratory 100 km/10%; 30 common days, 12 hours/day and
two distinct simultaneous reference stations. These are proposed settings,
not published acceptance cutoffs. Generic sensor labels require review.

Three PNGs and numerical/audit CSV/JSON outputs are beside each NPZ in `nearby/`.
Both masters receive `nearby_*` diagnostics. All six scores are preserved;
the second status changes only if its score is blank. One pair is context only.
Restart fingerprints cover the ready reference pool, not just the target.
Wait for the current PPSD batch to finish before starting another master writer.

Isolated copies of RR43/RR44/RR47/RR48/RR52 pressure produced four exploratory
comparisons and one without suitable peers in that selected pool. RR44 used
three references over 153 days, median -1.99 dB. Original NPZs and scores were
preserved; the production master was not patched by verification.
See [method](../nearby_comparison.md) and [verification](../nearby_validation_summary.json).

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network nearby
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind nearby
```

## Possible next work — not yet approved or implemented

The first, second and third scientific comparisons are implemented as diagnostics.
Next: verify reference applicability, water depth and instrument bands before
assigning acceptance scores. Next code work is paired pressure–vertical
cross-spectra, coherence and transfer functions for the remaining checks.
A separate short-window product is an optional later task; preserve the current
3-hour results and their master references when adding configurations.

Example settings checked against RR52's FFT/smoothing support:

| Window | Spacing at 50% overlap | Example minimum frequency |
| --- | --- | --- |
| 3 hours | 90 minutes | 0.001 Hz |
| 1 hour | 30 minutes | 0.003 Hz |
| 30 minutes | 15 minutes | 0.01 Hz |
| 15 minutes | 7.5 minutes | 0.02 Hz |

These are supported computational settings, not verified instrument cutoffs.
Two-hour windows fail the current full-bin support check at 0.001 Hz.
Shorter windows require processing the original waveforms again; more overlap
alone gives denser frames while retaining the same averaging duration.

Still pending: network prioritization, instrument usable-frequency validation,
acceptance criteria for all six checks, the three pressure–vertical checks,
paired pressure–vertical cross-spectra/transfer
functions, quality classifications, and combined noise models. Calibration
scores remain uncomputed; a successful PPSD run is not a calibration pass.
