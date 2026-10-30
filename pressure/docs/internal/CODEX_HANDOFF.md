# BRUIT-FM — start here in a fresh Codex chat

**Read `CODEX_NEXT_CHAT.md` first.** The six-check integration in
`docs/six_check_review.md` supersedes the old checkpoint counts below.
The review/scoring/conditional noise-model workflow is implemented, but policy
adoption and independent scientific acceptance remain pending. No production
master, score, waveform or response was changed. The prior local-P report is
`docs/local_p_phase_report.md`.
Latest owner decision: coverage/statistical evidence first, no ratings yet.
Claude's missing-work runner was reviewed with additional epoch/empty-coverage
safeguards; see `docs/reviews/full_data_coverage_integration.md`. No long campaign
was started; the owner awaits a new computer. Historical instructions below
must not trigger repeated calibration batches or premature acceptance policy.
The reviewed `codex/local-p-finish` integration is now on main after the owner
requested cancellation of the separate long batch. Its partial outputs and all
contributor branches are preserved. See `docs/local_p_integration_closeout.md`.
Scores/production masters are unchanged. Do not restart a completed campaign
from this historical note.

**Latest checkpoint:** README progress is complete through **step 49**.
Engineering fixes and optional scientific trial tools are merged into `main`.
**173 tests pass**. The completed 89-station trial has 61 provisional estimates
including C07W's three-event normal-state subset, compared with 46 preserved
baseline gains. It leaves 26 ordinary unresolved cases and two reference-band
exclusions. The original fitter remains the production default; no gains,
responses, masters or scientific scores were corrected.
See [the comparison](../socrates_all_station_comparison.md),
[contributors](../contributions.md) and [CODEX_NEXT_CHAT.md](CODEX_NEXT_CHAT.md).
The remaining sections are historical context, not current execution requests.
User's next priority (14 October 2026): **microseism pressure–vertical check**,
starting with the saved RR52 paired-spectra pilot. Read
[the next-task brief](../microseism_next_chat.md). Establish the source-supported
theory/identifiability before bounded implementation. No new microseism calibration
is implemented yet. Rayleigh troubleshooting and 7D/XO benchmarking are deferred.
Do not rerun completed calibration batches, refit RR52 or apply source gains.

Latest task completed: frozen-gain validation of RR01/RR06/RR48 with RR52 as
control, three flagged-arrival waveform reviews, and an RR05 replacement study
selected using coverage before new errors. RR52 has 14 isolated events (+0.664%
median error); RR06 one; RR01/RR48/RR05 none. No gains or guards changed.
Read [the follow-up](../yv_heldout_arrival_review.md) and CODEX_NEXT_CHAT.md.
Signal/noise and packet-placement research remains open but deferred. Reviewed
cases are development data if windows change; do not rerun completed batches.

This file retains historical project context. Read CODEX_NEXT_CHAT.md first,
then consult this file only as needed. The handoff does not authorize
starting every pending computation or redesigning existing stages.

The available earlier project chats, this chat and the original photographed
WP2 brief have now been reconciled in [project_context.md](../project_context.md).
It records the original multi-network/global noise-model scope, user preferences,
all six checks, limitations and ordered remaining work. Read it before proposing
the next scientific stage. Historical pilot descriptions below remain evidence,
not the current execution status.

## Workspace, inputs and saved state

- Local checkout: `/Users/mohammadamin/Projects/bruit-fm`.
- GitHub: <https://github.com/MohammadAmin-Aminian/bruit-fm>, branch `main`.
  Code and documentation through step 22 are in `220458b`; the compact handoff
  is a subsequent documentation change. Inspect `git log -1` for the latest identifier.
- Environment: conda `bruitfm`, Python 3.13, ObsPy 1.5.1, NumPy/SciPy,
  Matplotlib, pandas and openpyxl. Local interpreter:
  `/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python`.
- Whole archive: `/Users/mohammadamin/Desktop/Data_Bruit_FM/Data`.
  Inventories: `/Users/mohammadamin/Desktop/Data_Bruit_FM/INVENTORIES`.
  Earlier station pilot: `/Users/mohammadamin/Desktop/Data_Bruit_FM/RR52`.
- Shared master: `exports/master_inventory.csv` and `.xlsx`, 16,784
  channel-epoch rows across inventories. There are 227 YV inventory rows;
  172 available channel epochs were processed. Keys are `Seed_id`,
  `channel_start_time`, `channel_end_time`, `inventory_file`.
- Network output: `results/NETWORKS/TEST/<network>/<station>/<epoch>/`;
  master fields identify exact product paths. Original RR52 products are in
  `results/TEST/YV/RR52/`. Backups are under `exports/backups/`.
- Git tracks code, tests, documentation, reference data and curated example plots.
  Raw data, external StationXML inputs, `exports/`, `results/` and full NPZs
  remain local and are excluded from Git. A new chat on this same machine can
  read them; a different machine needs these local files separately.

## Current YV execution state — inspect this before running commands

The archive contains 43 stations, 172 channels/epochs and 45,593 waveform files
(370.08 GB), all matched to inventories. The live `_workflow` summaries were
read for this handoff. [Saved snapshot](../yv_workflow_status.json):

| Stage | Recorded final outcome |
| --- | --- |
| PPSD | 172 complete |
| Base plots | 172 complete |
| Expected pressure | 43 `comparison_only` |
| Nearby comparison | 4 `exploratory_only`; 36 `no_suitable_neighbors`; 2 `insufficient_reference_stations`; 1 `insufficient_common_days` |
| Saved-PPSD Godin | 19 `comparison_only`; 24 `no_comparable_bins` |
| Raw Godin pilot | 39 `comparison_only`; 2 `no_comparable_bins`; 2 processing failures |
| Paired pressure–vertical | 1 `comparison_only` (RR52 network pilot) |
| Rayleigh network | 43 complete: 16 `candidate_comparison_under_assumed_sign`; 25 `unresolved`; 2 `reference_band_unusable` |

The user finished the Rayleigh batch. Its `rayleigh_*` master update ran; all
pre-existing master CSV columns match the pre-run snapshot exactly. No multiplier
was applied to a response, waveform or PPSD, and all six scores remain blank.
RR41/RR57 were excluded by the conservative 0.03–0.07 Hz reference preflight;
RR57's exact frequency-dependent model is borderline inside the cutoff.
[Completed network audit](../yv_rayleigh_completed_summary.json). Do not restart
this completed batch merely because older notes describe the interrupted 7/43 run.

An [initial saved-results review](../yv_rayleigh_saved_results_review.md)
rederived all 481 saved event spectra and 48 candidate gains without production
functions, checked all 16 candidate combined plots and sampled four arrival/
waveform cases. Arithmetic matches to roundoff; independent raw estimation and
full packet/training validation remain pending. Priorities include RR03's window
containing the predicted S arrival, RR07/RR40 fit residuals, coherent general phase
disagreement, RR45/RR53 rotation rejections and depth sensitivity. RR57's flag is
from a conservative preflight bound; its exact model is borderline inside the
cutoff. No gains, scores, statuses or guard were changed.

Subsequent authorized [independent validation](../yv_rayleigh_independent_validation.md)
checked seven raw-count CSD windows and all 48 arrival plots: 20 candidate windows
contain predicted P/S arrivals, and every candidate station has at least one.
46 quiet intervals passed independent continuity screening; RR11/RR25 September
training has fractional-grid boundaries. Rayleigh now opts into a strict native
grid guard before merging; original outputs remain unchanged. Split-tilt and
band/model perturbations are saved; full field acceptance/uncertainty is pending.
34 pressure responses have positive-real Laplace poles; current RESIF reproduces
RR03/RR22 metadata, so do not flip poles or suppress FIR warnings as a fix.
The [equations/code report](../teleseismic_dpg_calibration_report.md) and local
`output/reports/teleseismic_dpg_source_7724bc5.zip` let the user review the original
code while validation continues. Latest suite: 98 tests pass; Ruff passes.

Pilot failures:

- `YV.RR12.00.BDH`: `No complete continuous windows in selected interval`.
- `YV.RR42.00.BDH`: `Requested pilot interval is outside waveform coverage`.

Bounded follow-up reads establish that all 21 selected RR12 windows are flat
at −8,268,016 counts, as are four sampled interior windows in other months.
RR42's sole pressure file covers about one day, verified from its raw header.
[Investigation](../pilot_error_investigation.json). Original pilot/master
statuses remain unchanged; no alternate pilot was run.

These are processing/coverage issues, not failed instrument calibrations.
All six YV calibration score columns are still blank. Re-read the live files
before acting; the snapshot is not a monitor and an additional batch may be active.

## What is implemented

The [README](README.md#completed-steps) lists 22 completed development steps;
[research_log.md](../research_log.md) preserves their numerical checks.

1. `inventory.py`: nested network/station StationXML scanning, full channel
   epochs, duplicate flags, depth estimates, metadata and six score/status fields.
   CSV/Excel synchronization appends missing epochs and preserves existing data.
2. `archive.py`: recursive SDS discovery and bounded MiniSEED reads; real gaps
   and conflicting overlaps are excluded rather than interpolated.
3. `ppsd.py`: full response correction once, channel workers, checkpoints,
   response snapshots and provenance. Pressure spectra use Pa²/Hz; seismic
   spectra use acceleration PSD units. Defaults: 3-hour windows, 50% overlap,
   one-octave smoothing, 1/8-octave steps, lower target 0.001 Hz where supported.
   Ninety-minute spacing is not ninety-minute averaging resolution.
4. `plot_ppsd.py`: PPSD distributions, time spectrograms and all individual PSD
   curves. Seismic distribution/curve plots include Peterson NLNM/NHNM.
5. `pressure_test.py`: Webb bounds and Brown hydroacoustic limits, nominal
   response-band masks, PPSD overlays and diagnostic plots, CSV/JSON and master
   updates. Webb bounds come from supplied instructions; Brown curves are
   approximate digitizations with the µPa²/Hz conversion. See provenance in
   `src/bruitfm/reference_data/pressure_references.json`.
6. `infragravity_test.py`: depth-dependent Godin prediction, matching PPSD-bin
   smoothing, masks, plots and numerical diagnostics. Figure 3 was independently
   rechecked; unsupported bins do not contribute to the result.
7. `godin_pilot.py`: separate raw-pressure 12.42-hour Hann spectra, three windows
   per day, screening high count power at 0.05–0.06 Hz, full response correction,
   retained linear-power means and comparisons. Defaults to seven days near
   mid-recording, not deployment startup.
8. `nearby_test.py`: geodesic distance, relative depth, declared sensor labels,
   common UTC days and distinct-station peer spectra. Preferred 50 km/5% depth
   difference; exploratory up to 100 km/10%. At least 30 qualifying common days,
   12 hours/day coverage and two simultaneous distinct reference stations per
   retained frequency. Settings are proposals, not published acceptance limits.
   Daily dB medians, full-bin band masks, log-frequency interpolation without
   extrapolation and equal-station peer medians are documented in
   [nearby_comparison.md](../nearby_comparison.md). Within-day windows are not
   exactly synchronized. Single-pair results are context only; unknown or
   different sensor labels are exploratory. Reference changes invalidate caches.
9. `network.py`: automatic discovery, explicit selection, failure isolation,
   staged restart, parent-only master writes and saved-PNG viewing. Other
   networks added to the archive/inventories use the same commands.

No new instrument calibration or response modification has been applied.
The first three spectral checks are diagnostics. The paired pressure–vertical spectral foundation is now implemented in
`paired_spectra.py` with network stage `paired`; source-specific checks, tilt
correction and acceptance/classification remain pending. See
[method and commands](../paired_spectra.md).

## Evidence and scientific limits

- Latest local software check: 98 tests pass; Ruff checks/formatting pass.
  Check the latest GitHub CI run after pushing before claiming remote success. A read-only review checked 156 then-available YV NPZs for finite spectra,
  unique times, identity, sample rate and required files. That preceded the
  completed 172-channel run; do not claim all 172 received the same explicit review.
- RR52: 5,249 PSD windows per channel, 20,996 total; four channels at 62.5 Hz.
  The default pressure comparison had 16/19 median levels within provisional
  bounds. Godin had two correlated eligible bins and +16.13 dB median residual.
  Its separate interior-week pilot retained 20/21 windows. These numbers are
  historical RR52 results, not universal YV outcomes.
- Isolated nearby verification: RR43/RR44/RR47/RR48 exploratory; RR52 lacked
  suitable peers in that selected pool. RR44 used three peers over 153 days,
  with -1.99 dB median residual. Source NPZs and scores were preserved.
- Depth is estimated from negative inventory elevation, not verified bathymetry.
  Generic YV `Hydrophone` labels do not independently establish actual DPG/APG/
  hydrophone hardware. The 3 dB response-band heuristic is not a verified cutoff.
- Overlapping windows, smoothed bins and stations are correlated. Spectral
  spread is not a confidence interval; an offset may reflect the environment.
- Existing evalresp FIR-normalization warnings require metadata review. Do not
  silently change responses, declare sensor failure or remove warnings as a fix.
- Numerical/real-data details: [RR52](../rr52_run_summary.json),
  [network pilot](../network_validation_summary.json),
  [nearby pilot](../nearby_validation_summary.json), and method documents in `docs/`.

## What to do next

1. Read live stage summaries and inspect the two pilot failures. Verify raw
   coverage and choose a continuous interior interval or improve interval
   selection if the user requests a fix. Do not use the first deployment week
   or shorten the scientific windows merely to conceal gaps. Do not automatically
   rerun the entire network or rebuild completed PPSDs.
2. The paired pressure–vertical spectral stage has been authorized, implemented
   and verified with synthetic tests and an RR52 interior-week pilot. Network
   integration was initially tested in temporary exports without master writes.
   Live workflow now records one production RR52 paired comparison; this review
   preserved that result. No whole-YV paired batch ran. RR52 retained 167/168 hourly pairs,
   with median squared coherence 0.80214 in 0.1–0.2 Hz. This is an uncorrected
   diagnostic, not a calibration result. Full products remain local; curated
   [validation](../rr52_paired_validation_summary.json) and the plot are tracked.
   The follow-up source review below supersedes the earlier microseism-first
   priority: the user wants their teleseismic calibration method examined.
   Keep scores unchanged.
   Follow-up: the user requested reviewing their thesis/manuscript, Webb/Wayne
   material and earlier calibration packages. The teleseismic method is now the
   next pilot: `rayleigh_calibration.py` is standalone, with **no master or response
   update path**, and is not in network `all`. Eight RR52 mainshocks were examined.
   The initial local-method run gave three amplitude cases with median 0.66449
   and opposite phase. The user's subsequent instruction supersedes treating
   sign uncertainty as a reason to stop: **assume opposite pressure sign (−1)
   for RR52**, explicitly record `assumed_needs_confirmation`, and proceed.
   Do not change source responses/data or select signs independently per event.
   The default now invokes actual **TiSKitPy CleanRotator 2.3.1**. The retained
   `--tilt-method local_eigen` is experimental, `needs_to_be_tested`, with a
   [short method reference](../tilt_methods.md).
   With this default and assumed sign, two of eight mainshocks are candidates;
   Santa Cruz has a predicted overlapping packet, three fail band-mean coherence
   and two have 0.0001 s quiet-training start offsets (phase bound 0.00216° at
   0.06 Hz). The station remains unresolved: fewer than three accepted events.
   Three amplitude-review cases give median 0.66440 (includes the overlap);
   the two candidates give median 0.67044, not a validated station estimate.
   Wide arrival plots and saved manual-window overrides support inspection;
   no automatic picking, timestamp rounding or gap interpolation was introduced.
   [Window/timing guide](../rayleigh_window_review.md),
   [review](../rayleigh_calibration_review.md) and
   [current validation](../rr52_tiskitpy_validation_summary.json) preserve source commits,
   criteria, plots and limitations. The newer packages' original synthetic
   validation is not a previous field calibration. Current ComPy repairs several
   old DPG-code bugs: do not attribute old bugs to the current version.
   Catalogue is local at `results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json`;
   prepared historical disba phase/group/model CSVs and provenance are tracked in
   `src/bruitfm/reference_data/rr52_historical_rayleigh_*`. The historical model
   is not verified bathymetry/local structure or a verbatim PREM profile.
   Production master hashes were unchanged by this review/pilot. Existing raw,
   PPSD, paired and Godin products were preserved.
   Subsequent user request: match paper Figure 5 / thesis Figure 3.13 and add
   magnitude fallback. Every run now writes `paper_figure5_comparison.png` with
   the original semilog 0.005–0.1 Hz axes, labels and colors; broad display
   spectra are distinct from fit/QC spectra. Theory is restricted to the supplied
   phase-velocity range; no extrapolation or source-data correction is fabricated.
   The search starts Mw>7, then Mw≥6.5, then Mw≥6 until `--min-events` QC candidates
   (default three), with `--max-events` now a per-level budget (default eight).
   Preserve depth/phase/coherence/grid/overlap guards and cross-level 24 h
   separation; no repeated attempts. Explicit `--event-id` disables fallback.
   Default 0.02–0.06 Hz remains unchanged: 24 attempted events yield two candidates.
   Separate original-figure **0.03–0.07 Hz** run reaches three candidates after
   four Mw>7 attempts; event-gain median 0.678998, sample scatter 0.004715;
   gain fitted to median ratio 0.674374. This is still provisional, under assumed
   sign, with window/reference/uncertainty review pending and no applied correction.
   [Original-figure guide and commands](../rayleigh_paper_figure.md) and
   [new validation](../rr52_paper_validation_summary.json) give both local run
   paths and curated plots. The Rayleigh solver remains disba 0.7.0's fundamental
   phase/group dispersion table; no new runtime dependency was added.
   Latest request: make calibration ready for all available YV stations, record
   assumptions/checks, append progress points and push. The user then explicitly
   chose to run the batch themselves. **Do not resume the network run on their
   behalf without a new instruction.** The explicit `network rayleigh` stage is
   ready; network `all` does not include it. Defaults: actual TiSKitPy, paper band
   0.03–0.07 Hz, one assumed −1 pressure sign, shared historical RR52 model, and
   bounded >7/≥6.5/≥6 search for three candidates. These are provisional inputs.
   `rayleigh_network.py` binds epochs, isolates unresolved conditions, checkpoints,
   invalidates reuse on source/response/code/input changes, and appends only
   `rayleigh_*` fields in one parent write. Scores/responses are never changed.
   Historically, the agent's run stopped with SIGINT at the user's request after **7/43** epochs:
   RR03/RR07/RR08/RR09 candidate comparisons; RR01/RR05/RR06 unresolved.
   Those checkpoints were preserved. The user subsequently finished **43/43**:
   16 candidate comparisons, 25 unresolved, two unusable reference bands. The
   live report now says `complete` and `master_updated: true`; existing columns
   were verified unchanged, and new `rayleigh_*` diagnostics were appended.
   [Ready commands and outputs](../yv_rayleigh_calibration.md),
   [partial validation](../yv_rayleigh_validation_summary.json), and
   [assumption/checklist flags](../calibration_assumptions.md) are tracked.
   Step 20/21 documents implementation and completed computation, not certified
   field calibration. Master-linked Rayleigh viewing is available. Next review
   saved outcomes and flags; do not rerun or automatically apply candidate gains.
3. Use synchronized raw waveforms: separately saved PPSD magnitudes do not
   contain the phase/cross-spectral information needed. Resolve response epochs,
   gaps, sample rates, physical units, transfer direction and sign conventions
   explicitly. Independently test spectral normalization and known transfers
   before scaling up. Do not assume pressure/vertical channel letters are enough
   for arbitrary networks or imply that this stage establishes calibration.
4. Continue RR52 under the recorded −1 pressure assumption; review arrival
   isolation, timing policy and reference uncertainty before applying gains or
   classifying, and seek hardware sign confirmation for a final correction. A separate screened
   microseism check and local-P check remain pending.
5. Verify actual sensor bands, depth, reference applicability and acceptance
   rules before scoring. Later investigate/correct calibration where justified,
   classify channels and combine accepted noise spectra by depth, season,
   coastal distance and seafloor type. Network priorities remain to be chosen.

The supplied scientific brief asks for pressure and vertical spectra plus
transfer functions; six tests; one master with location/depth and individual
scores; station-level reusable outputs; and final low/high noise statistics
from accepted channels. Scores: 1 pass, 0 fail, intermediate partial, blank
when unperformed/not applicable with a separate status. Proposed categories:
nonworking, valid (all performed checks pass, at least one), invalid (none pass,
at least two performed), unresolved (remaining cases), and rescued after correction.
No classification code or coverage acceptance rule has been implemented.

Method leads from the supplied brief, to verify in primary sources before coding:
Ruan et al., DOI `10.1002/2014JB011162`, especially the microseism pressure–vertical
method; proposed replication network `7D2001`, station `J44A` (availability not
established). Rayleigh methods and OBSIC quality-assessment guidance were also
mentioned. Godin DOI `10.1002/grl.50418`; Doran DOI `10.1029/2019EA000783` is
calibration context, not a 50/100 km neighbor rule. Bibliography is in
[references.bib](../references.bib). Do not invent published thresholds.

## Run commands and implementation constraints

Exact YV commands are in [README](README.md#run-every-completed-stage-for-the-yv-network).
From the checkout, with the `bruitfm` environment:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network check --network YV
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind nearby
```

For headless software verification:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/bruitfm-mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python -m pytest -q
/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python -m ruff check .
/opt/homebrew/Caskroom/miniforge/base/envs/bruitfm/bin/python -m ruff format --check .
```

- Preserve raw waveforms, input StationXML, saved spectra and manually entered
  master results. Use `synchronize_inventory` when appending; fresh exports from
  `save_master_dataset`/the standalone builder now refuse existing masters.
  Cooperating creation/append/result writers share `.master_inventory.lock`.
- Match all four epoch/source keys. Use restricted update helpers and staged
  CSV/Excel backups; preserve extra sheets/styles. New stages need their own
  restricted fields/helper. Do not write production masters during another batch.
  The nearby unfinished-PPSD guard is a workflow check, not a general lock.
- Keep alternative computations in separate products/folders with provenance;
  do not replace the existing 3-hour PPSD configuration.
- Mirror the established module + network-stage + plots + restricted master
  updates + scientific tests + documentation chain. Include cache dependencies
  and empty/unresolved cases. Software failure differs from scientific failure.
- User wants concise explanations, exact commands and low GPT-6.1 usage. Read
  only relevant files; avoid unnecessary full reruns, huge tool output and agents.
- Keep README steps numbered 1–22, not a milestone table. Append subsequent
  implemented steps only when completed. No calendar dates or framing as a
  paper/report in the completed-step descriptions; explain main supervisor points.
- Continue committing/pushing completed authorized work and updating GitHub
  documentation. Do not infer authorization to implement the entire roadmap
  from this handoff request. Do not claim local datasets/results are on GitHub.
- Filesystem writes are permitted in the checkout and temporary folders;
  source data are outside the writable workspace. Use existing Git approval
  handling for commit/push; inspect applicable instructions in the new session.
