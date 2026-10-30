# BRUIT-FM: reconstructed project context

This note reconciles the available earlier BRUIT-FM user messages, this chat,
the original four photographed WP2 instruction pages, the handoff, development
history, scientific method notes and live workflow. Historical/pilot statements
in older documents must be read in that context; the live workflow is authoritative
for execution progress. Raw chat transcripts were not copied into the repository.

## Original objective and scope

The project is WP2: **low and high ocean-bottom noise-level models**. The wider
archive is described as approximately 20 TB, covering 2001–2026 and multiple
networks, seasons, water depths, coastal settings and seafloor types. YV is the
first available network trial, not the whole archive. The original instructions
ask for a priority list of networks and adequate environmental/global coverage.

Deliverables are:

1. An organized, incremental processing system for SDS waveforms and StationXML,
   saving reusable station/channel PPSD NPZs from the instrument-supported
   low-frequency range toward Nyquist, with a computational target of 0.001 Hz.
2. One master with SEED ID, latitude, longitude, depth, complete useful metadata,
   and separate result/score/status fields. The user requested CSV **and** Excel,
   clear labels, no redundant channel/location columns and preserved manual data.
3. Pressure and vertical spectra, coherence and transfer functions, and six
   scientific tests/checks. Only supported frequencies should be interpreted.
4. Defensible acceptance/classification and justified calibration corrections
   where needed, retaining original data and provenance.
5. An accepted-station database and combined/global low/high noise statistics,
   also grouped by depth, coastal distance, seafloor type and season.

The six fields are deliberately labelled as **three spectral diagnostics** and
**three pressure–vertical calibration checks**, following the user's request:

| Group | Check | Current implementation/execution |
|---|---|---|
| Diagnostic 1 | Expected pressure levels: Webb/Brown | Implemented; all 43 YV pressure epochs assessed |
| Diagnostic 2 | Nearby stations at comparable depth | Implemented; all 43 assessed; four exploratory comparisons |
| Diagnostic 3 | Predicted infragravity levels: Godin | Implemented; all 43 assessed; separate raw-window pilot also run |
| Calibration 1 | Earthquake Rayleigh pressure/upward-acceleration ratio | Working method; all 43 YV epochs processed, gains provisional |
| Calibration 2 | Microseism pressure–vertical comparison | Integrated experimental method; 89 stations reviewed, seven repeatable model-conditional estimates; absolute validation/scoring open |
| Calibration 3 | Local-earthquake P-wave pressure–vertical comparison | Reviewed v2 integrated into main: 120-pair staging, three-APG bounded rerun, LT02 conditional with seven events; 63 legacy executed pairs flagged, no automatic score/master write ([phase report](local_p_phase_report.md), [final handoff](local_p_integration_closeout.md)) |

All six score columns remain blank. The requested eventual scoring is 1 for
complete pass, 0 for complete failure, intermediate values for justified partial
pass, and blank for unperformed/not applicable, distinguished by status.
Acceptance rules are not yet adopted. Completed computation is not a passed test.

The read-only six-check integration, categorical review engine, explicit
accepted-channel approval path and station-balanced downstream noise-model
builder are now implemented ([workflow](six_check_review.md)). The default
policy is draft and the real-data accepted pool is empty. Synthetic approvals
test the machinery without approving real instruments. This does not complete
scientific validation or environmental/seasonal/global-model deliverables.

The original proposed categories are nonworking, valid (all performed tests pass,
at least one performed), invalid (none pass, at least two performed), unresolved
(other cases), and rescued after justified correction. Final criteria and coverage
requirements still need agreement. The brief also suggests reference weighting
by expected calibration reliability: APG, hydrophone, DPG. Current nearby code
uses equal-station medians and label screening; that proposed hardware weighting
has not been adopted or validated. Z-score screening is a suggestion, not a rule.

## Confirmed completed YV calculations

The available local YV archive has 43 stations, 172 channel epochs and 45,593
waveform files (370.08 GB). Inventories contain additional unavailable epochs;
these are not fabricated as processed data. The master contains 16,784 epochs
across 29 networks and 1,317 distinct network/station identities in the registered
inventories. Inventory coverage does not establish local waveform availability.

| Product | Live saved outcome |
|---|---|
| PPSD and base plots | 172 complete each |
| Expected pressure | 43 `comparison_only` |
| Nearby comparison | 4 exploratory; 36 no suitable neighbors; 2 insufficient references; 1 insufficient common days |
| Saved-PPSD Godin | 19 comparisons; 24 without comparable bins |
| Raw Godin pilot | 39 comparisons; 2 without bins; RR12/RR42 processing/coverage failures |
| Continuous pressure–vertical diagnostic | One RR52 interior-week comparison, 167 retained hourly pairs |
| Network Rayleigh comparison | **43/43 complete: 16 candidate comparisons, 25 unresolved, 2 unsuitable reference bands** |

The user completed the network Rayleigh run after the agent stopped at 7/43 at
their request. Its final master update ran. All **pre-existing master columns
match the pre-run CSV exactly**; new results/flags are appended as `rayleigh_*`.
The saved workflow is `results/NETWORKS/TEST/_workflow/rayleigh.json`.
No new response calibration, waveform correction or scientific score was applied.
The two unsuitable-band cases are RR41/RR57: the requested 0.03–0.07 Hz band
approaches the adopted water-column resonance margin at their estimated depths.
They are not classified as failed instruments.

RR12's original selected pressure windows and four other sampled interior windows
were flat at −8,268,016 counts. RR42 has approximately one day of pressure data.
These investigations preserve the original failure statuses; neither observation
establishes the physical cause or complete deployment behavior.

Latest recorded software verification: **92 tests passed** and Ruff checks passed.
This is a software check under tested assumptions, not absolute field calibration
certification. The entire 172-channel archive has not received every independent
audit performed on the smaller verification subsets.

## Choices agreed in this chat

- Treat the current Rayleigh implementation as a working method and save its
  limitations explicitly. It recomputes raw-count event spectra; no numerical
  target of 0.67 is built into production fitting.
- Use **actual TiSKitPy CleanRotator 2.3.1** as default. Keep `local_eigen` only
  as an experimental option, marked `needs_to_be_tested`, with a reference.
- Adopt provisional RR52 pressure sign −1 from the user's recollection and
  prior evidence; the network run explicitly extends that assumption across YV.
  It is one supplied convention, not an event-by-event coherence/sign adjustment.
- Use the supplied historical **disba 0.7.0 fundamental Rayleigh phase/group
  table**. Phase velocity enters the water model; group velocity proposes packet
  times. This shared historical model is not independently validated for all YV.
- Reproduce manuscript Figure 5 / thesis Figure 3.13 axes, labels and colors.
  Preserve nominal and provisional-sign phase plots and wide arrival/overlap views.
- Seek three QC candidates, trying Mw>7, then Mw≥6.5, then Mw≥6, with eight
  attempts per level by default. Do not relax QC merely to achieve the count.
- Use the original figure's 0.03–0.07 Hz band for network comparisons. The
  standalone runner's earlier 0.02–0.06 Hz default remains available separately.
- Distinguish the fitted gain of the median ratio from the median of event gains;
  neither event scatter nor agreement with historical values proves accuracy.
- The user prefers to run long network batches themselves. Do not start another
  batch or modify response/source/master results merely from a context review.

## Source review and scientific limits

The existing review records relevant calibration/tilt sections of the thesis,
available manuscript, tilt milestone note and BRUIT-FM material, plus Webb/Zha,
Ruan, Godin, Brown, Doran, TiSKitPy, disba and earlier calibration repositories.
This does not claim every chapter/publication or private Wayne correspondence
was read. Precise reviewed sources/revisions and local document names are in
[the calibration source review](rayleigh_calibration_review.md).

The earlier DPG script has historical filtering/picking/indexing problems;
current ComPy repairs several of them. Do not transfer those old bug claims to
the newer implementation. SeafloorCal and the smaller calibration prototype
provided useful tested components, not previously certified field calibrations.

Remaining concerns include polarity and actual hardware identity; vertical
absolute response; FIR-normalization warnings; depth/model/water assumptions;
sensor geometry; quiet training and tilt stationarity; isolated Rayleigh modes;
clock offsets; pilot QC/window/count heuristics; estimator and selection bias;
scalar gain versus response shape/time constant; and total uncertainty.
[Machine flags and validation checklist](calibration_assumptions.md).

A synthetic test shares the water-model function with its generator. A separate
acoustic-ODE test checks that formula under its assumptions, but neither certifies
the Earth/ocean reference or actual sensor sensitivity. The current RR52 estimate
near 0.67 is consequently a reproducible candidate, not independent proof of the
original calibration. No calibration correction has been applied to saved PPSDs.

## Next work, in order

**Current user priority, 14 October 2026:** move to item 3, the microseism
pressure–vertical check, beginning with the existing RR52 paired-spectra pilot.
Rayleigh follow-up is deferred; no completed calibration is to be rerun.
Read [the fresh-chat task brief](microseism_next_chat.md) for paths, source
questions and limits. This is preparation for a new check, not its completion.
The ordered list below preserves the original roadmap and historical counts.

1. Review the completed Rayleigh network outcomes, using saved plots/spectra,
   flags and metadata; diagnose the 25 unresolved and two reference-band cases.
   The [initial saved-results review](yv_rayleigh_saved_results_review.md) now
   checks all saved spectral arithmetic, candidate combined plots and recorded
   reasons; full packet/training inspection and independent raw estimation remain.
   The [subsequent raw/window review](yv_rayleigh_independent_validation.md) now
   checks seven raw windows, all 48 candidate arrival plots, quiet continuity and
   illustrative sensitivity. Physical mode/response validation and final acceptance
   still remain; timing and positive-pole metadata issues are explicitly recorded.
   Do not rerun the full network or original PPSDs automatically.
2. Independently validate representative Rayleigh measurements with a separate
   implementation, inspect packet isolation, extend event/held-out sampling and
   evaluate window/band/training/model/response uncertainty.
3. Develop the microseism method from primary sources, including the Ruan
   reference and proposed 7D2001/J44A reproduction (data availability unverified).
   Inspect the OBSIC quality-assessment material referenced in the original brief.
4. Develop the local-P method with a suitable local catalogue/picks and verified
   expected ratios; the brief mentions the Auger M2 work as a lead, not a supplied
   numerical model ready to run.
5. Establish sensor frequency limits and scientific acceptance, scoring and
   classification rules; apply justified corrections only in separate, traceable
   products while retaining originals.
6. Prioritize/add other networks, check seasonal/environmental/global coverage,
   and combine accepted spectra into the final low/high noise products.

The earlier 55–60% project estimate was an informal judgment, not a measured
completion fraction. A defensible whole-WP2 percentage requires a work breakdown
including additional networks, validation and final noise-model deliverables.

## Persistent working preferences

Keep computation incremental and memory bounded. Preserve raw SDS, StationXML,
PPSDs, existing manual results and Excel sheets/styles. Use full epoch/source
keys and restricted updates, with backups. Keep long-text metadata intact.
Maintain readable code, Python ≥3.13, scientific tests, Ruff, few necessary
dependencies, exact runnable commands and low usage. Keep GitHub progress as
numbered steps without dates or promotional framing. Continue committing/pushing
completed authorized code/documentation; full raw results and masters stay local.

[Start-here handoff](internal/CODEX_HANDOFF.md) · [Numbered progress](research_log.md) ·
[Run/view commands](yv_rayleigh_calibration.md) ·
[Completed network audit](yv_rayleigh_completed_summary.json).
