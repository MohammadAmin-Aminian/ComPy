> **Note.** This is the development log that served as the project README up to
> version 1.0.0. It records every step in order, with the commands used at the
> time. For installation and current usage see the [README](../BRUITFM_README.md).


# BRUIT-FM

**Station metadata and quality assessment for ocean-bottom noise research.**

[![Checks](https://github.com/MohammadAmin-Aminian/bruit-fm/actions/workflows/checks.yml/badge.svg)](https://github.com/MohammadAmin-Aminian/bruit-fm/actions/workflows/checks.yml)
![Python](https://img.shields.io/badge/Python-3.13%2B-3776AB?logo=python&logoColor=white)
![ObsPy](https://img.shields.io/badge/Metadata-ObsPy%20StationXML-276749)

BRUIT-FM connects StationXML inventories, daily SDS waveforms and saved
ocean-bottom noise spectra. The code has evolved from a metadata exporter
into a restartable workflow for spectral processing, pressure diagnostics
and automatically discovered networks.

## Completed steps

**Inventory and data structure**

1. Read StationXML inventories organized by network or station.
2. Preserve channel metadata epochs and flag duplicates.
3. Build the master table: **16,784 channel-epoch records**, with location,
   estimated depth, sensor and response information.
4. Export matching CSV and Excel master files.
5. Prepare score/status fields for the six quality tests; acceptance decisions remain pending.

**Waveform processing, diagnostics and network workflow**

6. Read four RR52 channels: pressure `BDH`, horizontal `BH1`/`BH2`, and vertical
   `BHZ`; **329 daily files per channel**, covering the available recording period.
7. Calculate response-corrected PPSDs using **3-hour windows, 50% overlap** and
   a 0.001 Hz lower-frequency target. Pressure and seismic spectra use their
   respective physical units; the Nyquist frequency is **31.25 Hz**.
8. Process channels with **four workers** and save **5,249 spectra per channel**
   (**20,996 total**) as reusable NPZ files. Save response/settings records,
   support checkpoint restart, and update the four master CSV/Excel rows.
9. Produce **12 plots**: distributions, time spectrograms and all individual
   spectra. Add NLNM/NHNM references to the seismic distribution and curve plots.
10. Verify saved spectra, plots and data preservation; **62 tests pass**.
    Confirm that the time views represent 3-hour averages at 90-minute spacing,
    and check shorter-window settings for a possible separate calculation.

11. Compare the saved BDH pressure spectra with supplied Webb bounds and
    digitized Brown hydroacoustic curves, within a nominal response band.
    **16 of 19 compared median levels are within the provisional bounds; three
    are above.** Save frequency diagnostics and link them from the master
    CSV/Excel. Status is `comparison_only`; calibration scoring remains pending.
12. Add a **PPSD with color-coded limits** (Webb in magenta, Brown in orange)
    and **test-detail panels** showing median clearance from the bounds and
    the percentage of windows within range. One command displays and saves
    all three pressure-comparison figures and records their paths in CSV/Excel.
13. Implement the **depth-dependent Godin infragravity model** and compare it
    with saved pressure spectra using matching frequency-bin smoothing.
    Save three model/residual plots and diagnostics in CSV/Excel. RR52 has
    only **two fully supported bins**, with a median difference of **+16.13 dB**;
    depth, sensor limits and model applicability still need review.
14. Independently recheck the Godin equation, FFT masks and raw-count PSD.
    Clarify Figure 3 with diagnostic versus eligible residuals, overlapping
    smoothing-bin widths and spectral spread. Numerical results are unchanged.
15. Calculate a separate **12.42-hour Hann-window pilot** from raw pressure
    data in a week near mid-recording. Screen high 0.05–0.06 Hz count spectra,
    retain **20 of 21 windows**, correct the response and average linear power.
    Save three plots, compare with the same week's existing PPSD, and append
    pilot diagnostics to master CSV/Excel. Calibration scoring remains pending.
    [Method and results](godin_pilot.md).
16. Extend the completed workflow to **whole networks discovered from the data
    folder**. Process available channel epochs, preserve existing inventories
    and scores, resume saved work, isolate failures and save all plots and
    diagnostic results in master CSV/Excel. Verify the batch on real pressure
    and seismic data at 50, 62.5 and 100 Hz.
    [Commands in order](network_workflow.md).

17. Compare **nearby pressure stations on common recording days** using saved
    PPSDs. Screen distance, relative depth and declared sensor labels; require
    at least two distinct reference stations and 30 jointly supported days per
    frequency. Save spectral/residual plots, daily heatmaps, selection plots
    and candidate audits; append diagnostics to CSV/Excel without assigning a
    calibration score. Preferred 50 km/5% and exploratory 100 km/10% settings
    remain configurable proposals. [Method and commands](nearby_comparison.md).

18. Calculate **paired pressure and upward vertical-displacement spectra,
    coherence and transfer functions** from synchronized raw waveforms. Validate
    complex response correction, normalization and phase conventions; save
    separate station products and optional restricted master diagnostics. Verify
    an interior RR52 week while preserving existing data and master results.
    Tilt correction and calibration acceptance remain pending.
    [Method, validation and commands](paired_spectra.md).

19. Review the original **teleseismic pressure-gauge calibration method** and
    implement a separate RR52 candidate comparison using predicted Rayleigh
    packets, TiSKitPy CleanRotator quiet-data tilt projection, pressure/acceleration spectra and a
    water-column reference. Independently check grid-search gain conventions
    and preserve phase/arrival review flags. RR52 pressure sign −1 is provisional;
    the retained local tilt method is flagged as needing field tests.
    No instrument correction or score is applied.
    Add the original paper Figure 5 / thesis Figure 3.13 comparison layout and a
    bounded Mw>7 → Mw≥6.5 → Mw≥6 fallback to seek three QC-accepted events.
    [Source review, results and commands](rayleigh_calibration_review.md).

20. Extend the working **Rayleigh calibration comparison to available network
    pressure epochs**. Bind all four response epochs, run TiSKitPy and bounded
    magnitude fallback, checkpoint each station, and append only `rayleigh_*`
    master diagnostics. Preserve existing spectra, scores, responses and extra
    Excel sheets. The user completed all **43 YV pressure epochs**: **16 candidate
    comparisons, 25 unresolved, two unsuitable reference bands**. Gains remain
    provisional and no correction was applied.
    [Completed audit](yv_rayleigh_completed_summary.json) ·
    [Run/view commands](yv_rayleigh_calibration.md).

21. Record **calibration assumptions and remaining validation** as named review
    flags in each result and the master. Document fixed pilot choices, shared
    Earth-model and polarity assumptions, field checks, uncertainty and the
    experimental local tilt method. Verify restricted updates, epoch ambiguity,
    cache invalidation and missing-channel/reference guards; **92 tests pass**.
    [Assumptions and checklist](calibration_assumptions.md).

22. Review **saved Rayleigh events and independent raw spectral estimates**.
    Inspect all 48 candidate arrival plots, screen their quiet intervals, compare
    seven raw-count Welch estimates and test training/model/band sensitivity.
    Flag body-arrival overlap, response conventions and uncertainty. Add a
    calibration-only guard against merging fractional native sample grids;
    preserve the completed run and masters. **98 tests pass**.
    [Equations and code report](teleseismic_dpg_calibration_report.md) ·
    [Independent validation](yv_rayleigh_independent_validation.md).

23. Make **teleseismic tilt correction optional**. Explicit `none` uses only
    pressure and vertical motion, skipping horizontal channels and quiet training.
    Keep both rotation methods available and retain event/response/QC safeguards.
    Add separate network comparison/viewing without master updates. **101 tests
    pass**; original YV results and masters remain unchanged.
    [Run without tilt correction](yv_rayleigh_calibration.md#comparison-without-tilt-correction).

24. Prepare the **five-station YL pilot** with pressure/vertical channels only.
    Accept downloaded daily MiniSEED names without renaming source files, register
    original station XMLs and cache a deployment-specific earthquake catalogue.
    Retain strict timing guards and provisional model/sign assumptions.
    [Five-station commands](yl_five_station_pilot.md).

25. Fix the **YL signed-sensitivity preflight**: accept finite nonzero signed
    sensitivities without changing response gains or polarity. Verify PSD sign
    invariance and complex-response sign preservation. A bounded A02W pressure
    pilot passes; **107 tests pass**. Completed vertical PPSDs remain reusable.
    [Retry instructions](yl_five_station_pilot.md).

26. Separate **event sampling from the minimum station-estimate count**. Default
    Rayleigh sampling now exhausts all eligible catalogue events; it does not stop
    after three successes or eight attempts. Historical limits require explicit
    `bounded` mode. Preserve selection/QC safeguards and all saved runs.
    **109 tests pass**; no real batch was rerun for this change.

27. Add **optional Earth-layer input with disba generation**. Compute and save
    fundamental Rayleigh phase/group dispersion through the requested band;
    record model/curve hashes and solver versions. Retain tabulated input and
    all original products. Compare YL gain conventions and a held-event band
    diagnostic; verified numerical Zha/Lau inputs remain pending. **117 tests pass**.
    [Model options and gain comparison](rayleigh_earth_model_options.md).

28. Prepare an **all-station YL run**: discover the 25 available XML station
    folders and 50 P/Z channels, cache a full-deployment catalogue and select only
    missing/changed PPSDs. Keep original pilot products and calibration masters.
    Add a runner that stops on errors; **119 tests pass**.
    [All-station command](yl_all_station_run.md).

29. Compare the **completed all-station YL run with Zha & Webb**. All 25 epochs
    processed: 22 provisional estimates, three unresolved. In reciprocal
    convention, 21 estimates are 3.4–7.0% below the published coefficients; C07W
    is an outlier requiring investigation. Preserve original results and masters.
    [Station comparison and limits](yl_all_zha_comparison.md).

30. Inspect **YL outliers and unresolved cases using saved results**. C07W's
    discrepancy persists in its accepted magnitude ≥7 subset; B01W's published
    coefficient implies a pressure multiplier beyond the current search limit.
    Record coherence, timing and window rejections for A14W/B01W/C10W. Preserve
    every original fit and score; response/window and method reconciliation remain
    the next targeted checks. [Findings](yl_all_zha_comparison.md).

31. Audit **the repository and completed saved products**. Verify 222 PPSD
    channels, 2,219 Rayleigh spectral products, other pressure diagnostics and
    CSV/Excel agreement while preserving all results. Fix full-interval PPSD
    reuse, pressure-unit aliases, explicit-depth preflight, coordinate-dependent
    cache reuse and auxiliary-fit failure handling. **130 tests pass**;
    physical calibration acceptance remains pending. [Audit](project_audit.md).

32. Add a **read-only matched-event review for YL C07W/C08W/C09W**. Compare
    saved spectra, coherence, phase and exact-epoch nominal responses; retain
    rejected-event labels and links to original arrival/waveform plots. Write
    separate review folders without refitting gains or changing scientific
    results. [Run the viewer](yl_matched_event_review.md).

33. Prepare **isolated all-station 7D and XO runs**. Bind header-identified P/Z
    pairs at equal native rates, fetch complete response epochs, preserve original
    filenames/timestamps and log incompatible-file exclusions. Use separate
    masters/results with checkpoint resume and no applied calibration gains.
    [Commands, channel types and limits](7d_xo_runs.md).

34. Record the **completed 7D run**: 30 PPSDs, 15 calibration jobs, four provisional
    station estimates and 11 unresolved. Retain timing/phase/coherence and event
    sufficiency flags; station-matched published gain benchmarks remain unverified.
    [Saved gains and next checks](7d_xo_runs.md).

35. Record the **completed XO run**: two provisional DPG gains (WD66/WD70),
    four unresolved stations. Identify the EDH hydrophone/APG channel mismatch
    before comparison with Deng et al. (2022); published frequency-dependent
    coefficients do not establish agreement with our scalar estimates.

36. Add an **explicit isolated XO APG runner** for LA33/LA34/LT02 HDH/HHZ.
    Support different native rates through working-copy anti-alias resampling,
    exact response snapshots without PPSD computation, and a recorded 0.1 ppm
    native-rate tolerance. Preserve original samples, responses, masters and saved
    results. [Run command and validation limits](7d_xo_runs.md).

37. Add **automatic all-station XO gauge selection**: APG hardware first, DPG
    otherwise, without EDH fallback. Use a new result root and response snapshots
    without PPSD computation. Validate 1 Hz gain/phase/alias preservation and
    distinguish high-rate HDH DPGs from HDH APGs.
    [All-station runner](7d_xo_runs.md).

38. Review the **completed XO APG and all-gauge runs**. The latest six-station
    run has four provisional gains: LA33/LA34 APGs and WD66/WD70 DPGs. Retain
    LT02/WD48 rejection reasons. APG 1/2 Hz gains agree within 0.00003%; DPG gains
    change by +0.20%/−0.008% while adding accepted events. Numerical stability
    does not establish published agreement. [Saved comparison](calibration_network_summary.md).

39. Check **Ruan's 7D calibration-factor availability**. Confirm the paper
    estimates station scale factors jointly with sediment structure; no verified
    numeric factor table located. Distinguish sediment-delay tables from gains
    and record the 2015 correction to the starting profile c=0.029.
    [Search findings and sources](7d_xo_runs.md).

40. Consolidate **four networks / 89 distinct stations / 46 provisional gains**.
    Recheck saved medians, reciprocals, counts, source hashes and blank scores;
    preserve earlier runs and all original inputs. **143 tests pass**. Update
    numbered progress, research log and handoffs, with explicit external-comparison
    and scientific-validation limits. [Checkpoint](calibration_network_summary.md).

41. Review **Sentinel's engineering changes and Socrates' scientific audit**.
    Preserve QC and spectra at fit boundaries, select matching tilt horizontals,
    protect existing masters and serialize all cooperating master writers.
    Close an additional fresh-export lock race with a failing-before/passing-after
    regression. [Integration verdict](reviews/codex_integration_review.md).

42. Evaluate the integration against **all four completed network workflows**.
    Reproduce 416 saved accepted-event fits and all 46 station estimates; check
    18 raw event windows, a real B01W boundary run and two C07W daily states.
    Preserve original data, responses, results and scores. **155 tests pass**;
    qualify hardware hypotheses, uncertainty and the published YV comparison.
    [Validation record](reviews/codex_integration_validation.json).

43. Add **optional scientific correction trials**: pressure-pole reflection,
    weighted log-gain fitting, C07W vertical-state subsets and diagnostic-only
    gain/corner fitting. Compare each option separately with saved no-tilt
    windows, retaining original calibration defaults and input/QC safeguards.
    [Trial commands and limits](socrates_selected_station_trial.md).

44. Compare the completed **89-station trial** with previous and published
    results. Pole reflection increases YV estimates from 16 to 31; the total is
    61 including C07W's three-event normal-state subset. Weighting worsens YL
    peer published agreement. Check 2,432 protected hashes without modifying
    original products or accepting instrument corrections.
    [Full station comparison](socrates_all_station_comparison.md).

45. Integrate the **reviewed engineering fixes and optional trial tools**.
    Record the contributions of the project owner, Sentinel, Socrates and Codex;
    preserve the original fitter as default and keep new gains provisional.
    Complete **167 software tests**, style checks and the final Git handoff.
    [Contributions and integration history](contributions.md).

46. Prepare **frozen-gain validation on unused earthquakes** for RR01/RR06/RR48,
    with RR52 as control. Exclude prior YV trial events and development probes;
    evaluate every eligible held-out mainshock without refitting the saved gains.
    Report frequency/event errors, response metadata and arrival-isolation flags
    in separate products. **173 tests pass**; physical validation awaits the run.
    [Validation command](yv_heldout_validation.md).

47. Review the completed **four-station frozen-gain study**. RR52 has 14
    isolated usable events with median prediction error +0.664%; RR06 has one,
    RR01/RR48 none. Preserve all 50 attempts per station and retain the new
    DEPAS gains as provisional. [Completed study](yv_heldout_results.md).

48. Inspect **three flagged native RR01/RR06 waveform windows**, retaining
    original QC labels. Select RR05 to replace RR48 using header coverage
    before new errors: 42 covered windows versus nine. Its 50-event frozen-gain
    test yields zero usable events (40 coherence failures, ten input refusals).
    Preserve sources; investigate signal/noise and isolation before acceptance.
    [Waveforms and replacement result](yv_heldout_arrival_review.md).

49. Prepare the **fresh-chat microseism pressure–vertical pilot** at RR52.
    Reuse its saved paired spectra; establish primary-source theory, usable band
    and whether an independent pressure gain is identifiable before implementation.
    Keep its Rayleigh gain fixed, defer unresolved-case troubleshooting and
    preserve all original data/products/masters/scores.
    [Next task and saved pilot](microseism_next_chat.md).

50. Establish the **microseism admittance theory and RR52 identifiability** from
    primary sources, then analyze only the master-linked saved paired spectra.
    Reproduce their aggregate and daily products, verify provenance and nominal
    response support, and show the fixed Rayleigh multiplier without fitting or
    applying it. The P/Z measurement jointly depends on pressure/vertical scale
    and seafloor structure, so it remains diagnostic and unscored.
    [Method and bounded result](rr52_microseism_diagnostic.md).

51. Prepare a **one-station-at-a-time microseism campaign** for the 46 stations
    with baseline Rayleigh estimates. Bind exact P/Z channels and frozen gains
    in a reviewed manifest; reuse RR52's saved spectra and create fresh isolated
    paired products elsewhere. Support the audited XO working-copy resampling,
    but fit no microseism gains and write no masters or scores.
    A sequential wrapper runs all 46 while retaining per-station failure logs.
    [Commands for every station](microseism_baseline_run.md).

52. Compare the user's **completed 46-station microseism campaign** using only
    its explicit log-to-diagnostic paths. Verify all saved artifacts and spectral
    aggregates; report coherence, phase, daily-reference coverage and the effect
    of displaying each frozen gain. The result remains non-identifying and
    unscored. A network-wide sign split (YV near 0°, the other 30 stations near
    180°) must be resolved before pooling signed admittances.
    [Comparison, table and figure](microseism_baseline_comparison.md).

53. Independently implement the **RR52 microseism pressure–vertical check** on saved spectra
    (independent Claude implementation, branch `claude/microseism-rr52`).
    Water-layer theory gives polarity and a structure-free upper bound below
    α/(4H) = 0.1285 Hz; it does not identify an independent gain. Sign −1 is
    supported; g_ref = 0.6787 satisfies the bound, and the best-conditioned
    0.047–0.063 Hz model-conditional range is −1.2…+0.4%. 0.1–0.2 Hz is
    structure-sensitive (Ruan 2014) and cannot calibrate. No gain applied.
    [Method and result](microseism_pressure_vertical_check.md).

54. Run the independent **microseism check at all 43 YV stations** (41 checked;
    RR12/RR42 unresolved). Response/polarity hypotheses explain many P/Z phases;
    RR47 has the opposite relative polarity. Hardware corrections are not
    established. Historical best-band gains are superseded by the frequency
    consistency review in step 56; RR50's discrepancy remains unresolved.
    [All-station results](yv_microseism_check.md).

55. Extend the independent microseism check to **all YV/YL/7D/XO stations** and
    build the four-column comparison: Claude microseism, Rayleigh teleseismic,
    Codex's diagnostic on the same spectra, and published values. Diagnostic
    only. [Commands and columns](microseism_all_networks.md).

56. Integrate Claude's method and Codex's diagnostics while preserving both
    histories. Implement the frozen, YL-designed frequency-consistency rule:
    at least two bands ending at <=0.10 Hz agree within 10%, followed by at
    least three independent passing weeks with <=10% spread. Retain legacy
    best-band results separately. Hash-check saved inputs, reject overlapping
    weeks or differing water/model/channel epochs, and write fresh review
    folders. Completed independent weeks give seven conditional repeatable
    estimates (six stable, one moderate), five outside the YL design set;
    absolute validation remains open. [Completed week evaluation](
    docs/microseism_independent_weeks.md) · [Integration review and next task](
    docs/reviews/microseism_integration.md) ·
    [Frequency/week table](microseism_frequency_table.md).

57. Audit phase in the **exact frozen gain bands** using saved complex spectra:
    22 weeks, 58 bands; 55 within the legacy 10-degree screen and three RR40
    bands approximate. Reproduce exact-epoch complex responses, retain daily
    flags and sparse coverage, and investigate cached water/dry model sensitivity
    without retuning gains. Export 89 full-epoch staging rows with blank scores;
    coherent tilt and absolute response remain unverified.
    [Audit and next bounded pilot](microseism_gain_band_audit.md) ·
    [Table staging export](microseism_table_staging.csv).

58. Implement a bounded **P/H1/H2/Z coupling diagnostic** and test RR52 on two
    days. Retain 46/48 hourly windows, refuse native-grid changes across midnight,
    and reproduce every retained saved P/Z spectrum. Strong low-band horizontal
    association and duplicate horizontal azimuth metadata prevent a negligible-
    tilt claim. Statistical projection is not adopted as a correction; frozen
    gains and scores stay unchanged. [Pilot and remaining constraints](
    docs/rr52_microseism_horizontal_pilot.md).

59. Recover the published RR52 orientation reference and test **saved quiet-trained
    weights on separate June data**, selecting the latest prior training by date
    only. Preserve native channel mapping caveats; do not edit XML azimuths.
    Low-band power falls but high-band coherence degrades, so do not promote an
    automatic correction. Retain conditional gains, blank scores and the existing
    staging table. [Separate-training review and continuation decision](
    docs/rr52_microseism_fixed_projection.md).
60. Implement the **local-earthquake P-wave pressure–vertical check**
    (Claude, branch `claude/local-p-check`). Deng Eq. 5 steady state gives
    conditional gain candidates in shallow water; deep stations allow only the
    time-gated pre-multiple relation, kept diagnostic. Pressure is
    forward-modelled from the vertical and compared through identical
    multitaper windows; an independent finite-difference reference checks both
    relations. The historical pilot is preserved; corrected v2 LT02 gains are
    0.968/0.973 (four/seven events by band, conditional). The LT02 APG clock
    drift is modelled, but timing refinement is not an independent phase test.
    LT02 EDH has very weak recorded activity, not a proven hardware diagnosis. YL is 1 Hz only
    (not applicable). No gain applied, no score.
    [Method](local_p_method.md) · [Data, pilot, commands](local_p_check.md) ·
    [All-station campaign](local_p_campaign_v1.md).

61. Close the **local-P engineering review** with the documented conservative
    noise guard, per-epoch exports, hash-bound inventory aliases, waveform-aware
    resume, safe outputs and parent failure isolation. Compare preserved v1
    results with a three-pair v2 rerun. Stage all 120 pairs, explicitly marking
    63 legacy executed pairs; test the restricted updater only on a private
    master copy. Do not promote deep pre-multiple gains or assign scores.
    [Phase report and safe commands](local_p_phase_report.md) ·
    [Epoch-bound staging](local_p_table_staging.csv).

62. Integrate the reviewed local-P phase into **main** after the owner's long
    batch is cancelled. Preserve all partial products and contributor branches,
    repeat the 332-test suite and saved-input checks, and document the disposition
    of Claude's independent code/science review. No full campaign is restarted,
    master changed or scientific result promoted.
    [Final integration and continuation handoff](local_p_integration_closeout.md).

63. Add read-only **pressure PPSD probability heat maps** for the four available
    networks, using every saved window, native frequency-bin coverage, the classic
    ObsPy PQLX palette and explicit display-only occurrence limits. Preserve
    full-range plots and exact histogram counts; no median curve, gain or raw
    recalculation is required. [Plotting commands](network_ppsd_overlay.md).

64. Integrate **all six checks into one exact-epoch review**: 89 stations,
    118 source identities, preserved manual scores and explicit unresolved joins.
    Add a draft categorical scoring policy, fingerprint-bound scientific
    decisions, an accepted-channel selection audit and a station-balanced
    downstream noise-model builder. Synthetic tests cover the approved path;
    the real run creates no approvals, changes no masters and generates no
    claimed validated noise model. [Review workflow](six_check_review.md).

65. Integrate Claude's **full-data coverage and missing-work workflow** after
    Codex review: explicit archive/epoch/response inventory, reuse, resumable
    isolated PPSDs and three spectral checks, and station-balanced versus
    window-pooled pressure heatmaps. Add exact-epoch reuse and empty-result
    coverage safeguards. No long campaign, new calibration defaults, ratings,
    gains or master updates. [Review and next steps](reviews/full_data_coverage_integration.md).

66. Review pressure-noise references and correct new displays/comparisons:
    exclude unverified handout points attributed to Webb, preserve the original
    Brown transcription and add a traceable approximate full-band one, and
    distinguish environmental regions from sensor/calibration bands. RR52/RR40
    viewers now open a four-panel comparison plus gain-only dB shifts, without
    raw processing, gain refitting or rewriting historical products.
    [Source review and command](pressure_reference_review.md).

67. Review Claude's **YV Rayleigh unresolved-station diagnosis**: preserve the
    16 nominal and 31 experimental estimates, verify response/polarity hypotheses
    and instrument history, and distinguish screening blockers from physical
    causes. Harden the saved-input audit's hashes, exact epoch/P/Z binding and
    safe outputs; reproduce 43 station records without fitting gains. **59
    targeted tests pass**; all 1,104 protected inputs remain unchanged. No raw
    waveform processing, metadata correction, rating or production-default
    change. [Review and next step](reviews/yv_rayleigh_diagnosis_integration.md).

68. Audit pending contributor branches and prepare a **tested deployment
    checkpoint**. Selectively retain the formula review, fix circular phase
    summaries and refuse fractional native-grid tears, while keeping reviewed
    local-P safeguards and scientific defaults. **447 tests pass**; a built wheel
    installs outside the checkout with bundled references and nine safe CLI
    help checks. No data campaign, gain application or scientific approval.
    [Review/dispositions](reviews/deployment_readiness_20261021.md) ·
    [Deployment guide](deployment.md).

[Detailed steps](research_log.md) · [Network commands](network_workflow.md)

## Current state

[Current project handoff to Claude (21 October 2026)](claude_project_handoff_20261021.md):
reviewed main, frozen evidence, next response-confirmation task and execution
boundaries. The latest software checkpoint passed **447 tests** and post-merge
CI; scientific acceptance and the full-data campaign remain separate work.

Preserved calibration workflows cover **YV, YL, 7D and XO: 89 distinct stations
and 46 provisional baseline gains**. The subsequent optional pole/state trial
produces **61 estimates**, including C07W's three-event normal-state subset;
26 stations remain unresolved and two retain reference-band exclusions.
These are numerical QC outcomes, not accepted instrument calibrations.
The six-check checkpoint passed **388 tests**; current coverage integration
verification is recorded in its [review](reviews/full_data_coverage_integration.md).
No gains are applied and
scientific scores remain blank. The owner's separate all-methods batch was
interrupted; its completed base and three shifted weeks remain saved, with
the unfinished -1 week explicitly excluded from completion claims.

The six-check review and downstream model code now exist, but the acceptance
policy remains draft. Scientific approvals, field validation, environmental
coverage/grouping and final global noise-model acceptance are not finished.
Start a read-only review with `bash scripts/run_quality_review.sh`; this is not
a calibration rerun or a request to populate production scores.
The owner now prioritizes coverage/statistical evidence before ratings. Use
`bash scripts/run_full_coverage.sh plan` for the missing-work inventory; postpone
the long execution until the new computation machine is ready.
[Baseline results](calibration_network_summary.md) ·
[Completed trial and publication comparisons](socrates_all_station_comparison.md).

Sentinel's engineering changes, Socrates' review and Codex's optional trial
implementation are merged. Normal calibration commands retain the original
scientific defaults. Weighted fitting remains a sensitivity option: it worsens
YL's systematic published discrepancy. New pole-trial gains still require
response identity, phase/arrival and held-out uncertainty review.
[Contributors and responsibilities](contributions.md).

The completed **YL calibration run covers 25 stations**: 22 provisional gains and
three unresolved cases. Comparing our reciprocal gains with Zha & Webb gives
21 values 3.4–7.0% lower (median −5.4%); C07W differs by +86.8%. Model, frequency
band, event selection and averaging still differ from the publication, so this
does not establish final calibration acceptance. No gains were applied.
[Short comparison and full station table](yl_all_zha_comparison.md).

The full RR52 calculation and its initial pressure/Godin diagnostics are saved.
The YV run completed **172 PPSD and base-plot jobs**, and all 43 pressure channels
were processed through the expected-pressure, nearby and Godin comparisons.
The 43 Godin pilots include **two processing failures**: RR12 has no continuous
windows in its selected interval; RR42's interval is outside waveform coverage.
Bounded reads show all 21 selected RR12 windows are flat; RR42 has only about
one day of pressure coverage. Their original statuses remain preserved.
[Investigation](pilot_error_investigation.json).
The paired pressure–vertical stage is implemented and validated on an interior
RR52 week, with 167 hourly pairs and median microseism squared coherence 0.802.
It remains diagnostic, without tilt correction or calibration scoring.
[Method and results](paired_spectra.md).
The teleseismic RR52 pilot uses TiSKitPy CleanRotator and a documented provisional
−1 pressure sign. Two events pass the pilot checks; one additional coherent case
has a predicted packet overlap. The three amplitude-review cases give median
0.664, close to historical 0.67. Window, timing and model/uncertainty review remain;
no calibration factor is applied.
[Review and comparison](rayleigh_calibration_review.md).
[Open plots and review event windows](rayleigh_window_review.md).
The new [paper Figure 5 comparison](rayleigh_paper_figure.md) preserves its
axes, labels and colors. A separate 0.03–0.07 Hz run finds three candidates with
event-gain median 0.679; its aggregate median-ratio fit is 0.674. The existing
0.02–0.06 Hz default remains unresolved after the bounded magnitude fallback.
These estimates still require window/model/uncertainty review.
The archive contains **43 stations, 172 channels and 45,593 files**, all matched
to inventories. [Workflow snapshot](yv_workflow_status.json); live execution
details are under `results/NETWORKS/TEST/_workflow/`.

An earlier network checkpoint passed **173 automated tests**, covering response
correction, source and master preservation, restart, spectral calculations, plots
and network failure handling.
The network validation used a separate export and an interior week on four real
channels, producing 110 PPSD windows each. All three pressure pilots retained
20 of 21 Hann windows; the existing RR52 pilot arrays were reproduced exactly.
[Validation](network_validation_summary.json).

The nearby-station stage was checked on isolated copies of five completed YV
pressure products. Four support exploratory comparisons; RR52 has no suitable
neighbor in that selected reference pool. RR44 uses three reference stations
over 153 common days. [Nearby verification](nearby_validation_summary.json).

The pressure comparisons remain **diagnostics**, with calibration scores blank.
Response removal uses the response already declared in StationXML; the code
has not estimated or applied a new instrument calibration. Depth estimates,
usable sensor bands and reference applicability still need validation.

The network Rayleigh run is **complete for all 43 available YV pressure epochs**:
16 candidate comparisons, 25 unresolved and two unsuitable reference bands.
The new `rayleigh_*` fields are saved in both masters; all pre-existing CSV
columns were verified unchanged. No response/waveform correction or score was applied.
[Completed audit](yv_rayleigh_completed_summary.json).
Next, review saved outcomes, flagged assumptions and event windows, and define
acceptance before applying gains or assigning scores.
[Reconstructed original scope and remaining work](project_context.md).
[Initial review of saved gains, plots and unresolved cases](yv_rayleigh_saved_results_review.md).

[Detailed steps](research_log.md) · [RR52 processing settings](rr52_methods.md)
· [RR52 saved results](rr52_run_summary.json) · [Network commands](network_workflow.md)
· [Fresh Codex handoff](internal/CODEX_HANDOFF.md)

## Run every completed stage for the YV network

Run in **Terminal**, one command at a time, waiting for each calculation to
finish. YV is the network; these commands select all its available stations
and channels. PPSD/plot stages cover both pressure and seismic channels;
pressure comparisons and the Godin pilot use pressure channels only.

If your PPSD batch is still running, let it finish before running another
command that updates the master. You can then continue from stage 4 below.
Rerun stage 3 if processing was interrupted or channels need retrying.

**0. Enter the checkout and limit numerical threads per worker.**

```bash
cd /Users/mohammadamin/Projects/bruit-fm
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
```

**1. Synchronize StationXML inventories into the existing master.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network inventory
```

Inventory synchronization covers the shared inventory collection. It appends
missing epochs and preserves existing results rather than rebuilding the master.
The following calculations are restricted to YV.

**2. Check YV waveform availability and inventory matches.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network check --network YV
```

**3. Compute/resume full-period PPSDs with four channel workers.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network ppsd --network YV --workers 4
```

Default windows are three hours with 50% overlap. The full StationXML response
is applied once. Existing compatible checkpoints resume automatically.

**4. Save PPSD distributions, time spectrograms and all spectral curves.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network plots --network YV
```

Seismic distribution/curve plots include NLNM and NHNM.

**5. Compare pressure levels with Webb and Brown reference limits.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network pressure --network YV
```

**6. Compare nearby pressure stations at similar depth.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network nearby --network YV
```

Targets are restricted to YV. References can come from any network with ready
pressure products registered in the master; distance/depth/time criteria still
apply. Unsupported comparisons are unresolved rather than failed calibrations.

**7. Compare saved pressure spectra with the Godin infragravity model.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network infragravity --network YV
```

**8. Compute the separate screened full-window Godin pilot.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network pilot --network YV
```

Uses a seven-day interval near the middle of each pressure epoch, with full
12.42-hour Hann windows. This stage reads raw waveforms; stages 4–7 reuse NPZs.

**9. Open saved figures for one station, one type at a time.**

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind ppsd
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind pressure
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind nearby
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind infragravity
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind pilot
```

Replace `RR44` with another processed station. Each view command opens saved
PNGs without repeating calculations or changing the master.

Defaults use `~/Desktop/Data_Bruit_FM/Data`, `~/Desktop/Data_Bruit_FM/INVENTORIES`,
`exports/` and `results/NETWORKS/TEST/`. Results and figure paths are saved in
**both `exports/master_inventory.csv` and `exports/master_inventory.xlsx`**;
backups are under `exports/backups/`. The six calibration scores stay blank
unless explicitly assigned, and existing manual scores are preserved.

Completed steps 1–22 describe the development history. Commands 1–8 above are
the current execution stages; each command includes the earlier supporting
steps automatically. No separate response-removal command is needed.
The provisional Rayleigh comparison is available through the separate command
below. Microseism/local-P checks, applied instrument corrections, quality
classification and combined noise models remain pending.

**10. Run or resume provisional Rayleigh comparisons for YV.**

> **The YV batch is complete; do not rerun this command against the default
> results folder.** Reuse fingerprints include the package code, so after any
> code change this command recomputes every station instead of resuming, using
> the current all-event sampling (`--rayleigh-max-events` applies only with
> `--rayleigh-event-sampling bounded`). It would replace the historical
> `results/NETWORKS/TEST/_workflow/rayleigh.json` entries and the master's
> `rayleigh_*` fields. Use the viewer below for saved results; run any new
> comparison with a new `--results-dir` and `--rayleigh-no-master-update`.

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network rayleigh \
  --network YV \
  --rayleigh-catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --rayleigh-dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --rayleigh-frequency-min 0.03 --rayleigh-frequency-max 0.07 \
  --rayleigh-pressure-sign -1 --rayleigh-min-events 3 --rayleigh-max-events 8
```

After completion, view one station and open the master results:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YV --station RR52 --view-kind rayleigh
open exports/master_inventory.xlsx
```

Change `RR52` to another processed station. Results/flags use `rayleigh_*`
columns; insufficient/reference cases retain blank station gains. This stage
does not rerun earlier computations or apply the estimated multiplier.
[Output details and assumptions](yv_rayleigh_calibration.md).

After reviewing the input check, the stages can be run in sequence. `all` runs
the eight stages above **and then `paired`**, which reads seven days of raw
pressure/vertical waveforms per pressure epoch and appends `paired_*` master fields:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network all --network YV --workers 4
```

Use either the individual commands or `all`. Do not start `all` alongside an
active batch. Stages reuse saved work when fingerprints and output files match;
`all` continues other stages after recorded channel failures and exits nonzero
if any stage has failed/partial results. Inspect `_workflow/<stage>.json` and
retry the relevant command to resolve errors.

## Pressure PPSDs for a whole SDS drive

`scripts/run_fresh_pressure_ppsd.sh` runs inventory, check, PPSD and plots for
the pressure channels of a drive holding `INVENTORIES/` and `SDS/` (or `Data/`),
writing everything to a new folder (default `~/Desktop/BRUIT-FM/fresh_run_<date>`):

```bash
bash scripts/run_fresh_pressure_ppsd.sh /Volumes/DATA_bruitfm inventory
bash scripts/run_fresh_pressure_ppsd.sh /Volumes/DATA_bruitfm check --network 7D
bash scripts/run_fresh_pressure_ppsd.sh /Volumes/DATA_bruitfm all --network YV
```

The network options it uses are also available directly:

* `--pressure-only` restricts the `ppsd` and `plots` stages to pressure channels.
* `--workers max` uses every CPU core but one; the largest channels start first.
* `--skip-uninventoried` sets aside waveform channels without StationXML rows
  (listed in `_workflow/uninventoried.json`) instead of stopping the batch.

For large archives on spinning disks: with `--network`, SDS `YEAR/NET` folders
of other networks are not scanned, and each daily file is indexed from its first
and last MiniSEED records rather than all headers. Filenames with location `00`
match an empty header location when the inventory only has the empty code (the
alias reviewed in `docs/coverage_sources.json`). The 23–24 October run over the
whole drive is described in [the 26 October status](project_status_20261026.md).

## Quick start

For whole-network processing from `~/Desktop/Data_Bruit_FM/Data`, use
`python -m bruitfm.network` in the `bruitfm` environment.
[Run each completed stage in order](network_workflow.md). The workflow
discovers YV and subsequently added networks automatically; `--network YV`
restricts the selection. It updates the existing master without rebuilding
and erasing previous results. Network outputs are under `results/NETWORKS/TEST`.

Requires **Python 3.13 or newer**.

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'

python -m bruitfm.network inventory \
  --inventory-dir /path/to/INVENTORIES \
  --master-dir exports
```

The command creates `master_inventory.csv` and `master_inventory.xlsx` if
needed, or safely appends new inventory epochs to the existing master.
Input files are read without modification. Raw archives and generated exports
stay local and are excluded from Git.

For local convenience, the default inventory folder is
`~/Desktop/Data_Bruit_FM/INVENTORIES`. Use `--inventory-dir` for another location.
Run `python -m bruitfm.network --help` to see the batch options.

### Use from Python or VS Code

```python
from bruitfm.inventory import build_master_dataset, save_master_dataset

master_df = build_master_dataset("/path/to/INVENTORIES")
save_master_dataset(master_df, "exports/new_inventory")

# In an interactive window, evaluate master_df to browse the table.
master_df
```

Importing the module does not scan the archive or write files.
The Python example creates a separate export. Use the network `inventory`
command when updating a master that already contains processing results.
`save_master_dataset` and `python -m bruitfm.inventory` refuse to replace an
existing `master_inventory.csv` or `.xlsx`.

## Calculate PPSDs from RR52

Use the existing master inventory exports to find the StationXML response for
each waveform channel. The data directory may be a station folder such as RR52
or a complete SDS archive; daily filenames must follow
`NET.STA.LOC.CHA.D.YEAR.DAY`.

```bash
conda activate bruitfm
python -m bruitfm.ppsd \
  --data-dir ~/Desktop/Data_Bruit_FM/RR52 --workers 4
```

This processes all matching channel epochs in the supplied folder. To start
with pressure and vertical channels only:

```bash
python -m bruitfm.ppsd \
  --data-dir ~/Desktop/Data_Bruit_FM/RR52 \
  --seed-id YV.RR52.00.BDH --seed-id YV.RR52.00.BHZ
```

Dates are UTC and the end is exclusive. Defaults use three-hour PSD windows,
50% overlap, one-octave smoothing, 1/8-octave period steps, and a target lower
frequency of 0.001 Hz. Histogram amplitude bins span -250 to 100 dB in 1 dB
steps. Reads are approximately one day plus the overlap needed
to avoid losing PSD windows at file/chunk boundaries. Missing samples and
conflicting overlaps are excluded rather than interpolated or zero-filled.
Fractional-sample timestamp changes between daily files are accommodated by
advancing read boundaries. Source timestamps are not changed on disk;
in-memory merging follows ObsPy's sample-grid rules.
Frequency support is validated against the actual FFT length. A supported FFT
frequency is **not** a verified instrument low-frequency cutoff.

Input remains in raw digital counts. ObsPy PPSD corrects the full StationXML
response internally, once. Declared `PA` input units select hydrophone handling
(including the DPG channel), producing pressure PSDs in `Pa^2/Hz`. Motion input
units select normal seismic handling, producing acceleration PSDs in
`m^2/s^4/Hz`. Stored PSD values are in decibels, `10 log10(power density)`,
relative to these physical units. Missing/unsupported responses and sampling-rate
mismatches are reported as processing failures. Do not supply waveforms already
corrected to physical units.

Outputs are organized under:

```text
results/TEST/{network}/{station}/{Seed_id}__{epoch_start}__{inventory_hash}/
  ppsd.npz          Reloadable ObsPy PPSD data, including individual PSD windows
  inventory.xml     Snapshot of the full response for this channel epoch
  processing.json   Parameters, source file paths/sizes/mtimes, inventory hash,
                    processed chunks, warnings, errors, and result summary
```

Checkpoints are saved every ten chunks and at completion. Rerunning the same
command resumes completed chunks; extending a pilot's dates adds spectra
without duplicating existing windows. Settings and source file attributes must
match the checkpoint. Use `--rebuild` to recalculate after changes, or use a
different `--results-dir` to keep an alternative calculation. Source waveform
content is not hashed; the provenance records size and modification time.
`--workers` defaults to one to keep memory/CPU use modest. Independent channel
workers can run in parallel; only their parent process updates the master files.
Do not launch separate processing commands against the same output folders at
the same time.

Matching rows in `exports/master_inventory.csv` and `.xlsx` receive `ppsd_*`
columns: processing status, NPZ/provenance paths, physical units, segment and
source-file counts, parameters, FFT/period-bin frequency coverage, histogram
bounds and out-of-range value counts, processed
dates, requested dates, update time, and errors. Counts and processed dates are
cumulative across resumed runs. `complete` means processing finished without
recorded chunk errors; it does not mean calibration passed or every sample was
available. Other statuses include `partial`, `failed`, `no_data`, and
`insufficient_data`; untouched rows show `not_processed`.

The updater patches only these fields, preserves existing calibration scores,
other Excel sheets and formatting, and backs up both master files under
`exports/backups/` before replacing them. Both outputs must already exist and
contain the same epoch/source keys. The CSV supplies processing metadata;
editing metadata only in Excel does not change the computation. Running the
inventory builder again creates fresh exports, so use a different output
directory if you want to keep existing processing results and manual scores.

Use `--master-dir` and `--results-dir` for other output locations, and
`python -m bruitfm.ppsd --help` for all processing options. Results and backups
stay local and are excluded from Git.

```python
from obspy.signal import PPSD

# Use the paths recorded in the master table / processing.json.
ppsd = PPSD.load_npz("/path/to/ppsd.npz")
```

Use the processing command to add more raw data: it initializes the response
lookup with the correct channel before importing saved spectra. The response
snapshot preserves the calibration used for the saved output. Transfer functions
require paired pressure/vertical cross-spectra and are not computed here.

### Plot saved RR52 results

In VS Code, select the `bruitfm` Python environment, open
`src/bruitfm/plot_ppsd.py`, and choose **Run Python File**. No arguments are needed:
it reads the saved RR52 paths from the master CSV and displays three plots per
channel (`BDH`, `BH1`, `BH2`, `BHZ`), for twelve figures in total. PNG copies are saved in
`results/TEST/YV/RR52/plots/`.

```bash
python -m bruitfm.plot_ppsd
# For a machine without a graphical display:
python -m bruitfm.plot_ppsd --no-show
# Display only one type (repeat --plot-type to select more than one):
python -m bruitfm.plot_ppsd --plot-type spectrogram
```

For a Python/VS Code interactive window:

```python
from bruitfm.plot_ppsd import plot_station

figures = plot_station()
```

Three plot types are saved per channel:

- `__ppsd.png`: the PPSD probability distribution with its coverage panel.
- `__spectrogram.png`: frequency versus UTC window start time, with color showing
  PSD amplitude in dB. Each column represents a three-hour PSD window (50%
  overlap in the default calculation); genuine gaps remain blank.
- `__all_psd_lines.png`: every individual saved PSD curve over frequency, with
  translucent lines and no subsampling or averaging. RR52 has 5,249 per channel.

The seismic distributions and curve plots (`BH1`, `BH2`, `BHZ`) include Peterson's
New Low Noise Model (NLNM) and New High Noise Model (NHNM), using ObsPy's model
data in acceleration PSD units. Models are plotted only over their available
frequencies within the displayed band. Pressure plots (`BDH`) use pressure units.
All plots use the full saved frequency range in Hz. Plotting reads the NPZ files;
processing metadata and calibration scores in the master CSV/Excel are preserved.

### Compare pressure levels with expected ranges

Run from the checkout with the `bruitfm` Python environment selected:

```bash
python -m bruitfm.pressure_test
# Save results without opening a window:
python -m bruitfm.pressure_test --no-show
```

VS Code can also run `src/bruitfm/pressure_test.py` directly. The default uses
RR52's saved pressure NPZ, with no waveform reprocessing. It saves observed
5th/50th/95th percentiles, per-frequency comparisons, JSON provenance and three PNGs
under the channel's `pressure_expected/` directory, and records result links
and diagnostics in both master files.

The three figures show percentile/reference levels, **PPSD with limits**,
and **test details** (median clearance and percentage of windows in range).
Webb limits are **magenta**; Brown limits are **orange** throughout.

The default test status is `comparison_only`; the calibration score stays blank.
Webb bounds are transcribed from the supplied instructions, and Brown Figure 8
is approximately digitized with the µPa²/Hz-to-Pa²/Hz conversion. A nominal
response-based frequency band excludes unsupported smoothing bins. Review
the references and sensor limits before assigning an acceptance result.
See the [comparison details](pressure_comparison.md) for provenance,
RR52 results and optional frequency/score settings.

### Compare with nearby pressure stations

After the network PPSD calculation finishes, run:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network nearby
conda run --no-capture-output -n bruitfm python -m bruitfm.network view --network YV --station RR44 --view-kind nearby
```

The stage reuses saved spectra, supports subsequently added networks, and saves
spectral comparisons, daily residual heatmaps and neighbor-selection plots.
It appends `nearby_*` diagnostics to both CSV and Excel. At least two distinct
reference stations and 30 jointly supported days per frequency are required.
Unsupported comparisons remain unresolved; the second check's score stays blank.
See [criteria, estimator and outputs](nearby_comparison.md).

### Compare with the Godin infragravity model

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.infragravity_test
```

Or run `src/bruitfm/infragravity_test.py` in VS Code with the `bruitfm`
interpreter. It reuses saved pressure data, displays/saves **three plots**
(percentile/model comparison, PPSD with model, and residuals), and records
diagnostics in the master CSV/Excel. Results are beside the NPZ under
`infragravity/`.

The default uses inventory-estimated depth and a nominal sensor band, with a
conservative 0.001–0.02 Hz model comparison range. Only complete smoothing bins
inside both bands contribute. Status is `comparison_only`; the predicted-
infragravity score remains blank. See [model details](infragravity_comparison.md)
for the equation, options, results and limitations.

### Calculate the screened Godin pilot from raw pressure data

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.godin_pilot
```

The default uses a week near the middle of recording coverage. It saves
full-length 12.42-hour Hann spectra, count-based transient screening, linear
pressure means and three plots, alongside the existing three-hour products.
New `infragravity_pilot_*` fields are appended to master CSV/Excel.

View the already saved plots without waveform processing:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.godin_pilot --plot-only
```

Use `--start` and `--days` to select an interval, or `--no-show` for background
processing. The calibration score remains blank. See [pilot details](godin_pilot.md).

## Master table

The first four columns follow the research document:

| Column | Meaning |
| --- | --- |
| `Seed_id` | `network.station.location.channel`, including empty location codes |
| `Lat` | Channel latitude in degrees |
| `Lon` | Channel longitude in degrees |
| `Profondeur` | Estimated depth in metres below sea level, from negative station elevation |

Depth is an **estimate from inventory elevation**, not verified bathymetry.
Nonnegative or missing station elevations leave this estimate blank. Sensor
burial depth is not used as water depth.

The next six columns store the check scores, in document order:

| Group | Result column | Comparison |
| --- | --- | --- |
| Spectral test | `test_pressure_spectrum_vs_expected_noise_levels_score` | Expected pressure-noise levels |
| Spectral test | `test_spectrum_vs_nearby_stations_at_similar_depth_score` | Nearby stations at similar depth |
| Spectral test | `test_pressure_spectrum_vs_predicted_infragravity_levels_score` | Predicted infragravity levels |
| Calibration check | `calibration_check_pressure_vertical_rayleigh_waves_score` | Pressure and vertical motion in earthquake Rayleigh waves |
| Calibration check | `calibration_check_pressure_vertical_microseisms_score` | Pressure and vertical motion in microseisms |
| Calibration check | `calibration_check_pressure_vertical_local_p_waves_score` | Pressure and vertical motion in local earthquake P waves |

The document groups all six under calibration tests. The prefixes distinguish
spectral comparisons from pressure–vertical motion checks; scores are not
instrument correction factors.

Scores use **1 = pass**, **0 = fail**, and intermediate values for partial passes.
They start blank. Each check also has a `_status` column initialized to
`not_tested`, allowing future code to distinguish an unperformed check from a
`not_applicable` check with no score.

Additional populated metadata follows the scores. Identifier components already
in `Seed_id`, globally empty metadata columns, and identical metadata columns
are removed. Scores, statuses, channel epoch dates, sampling rate, response units,
response-stage count, source inventory path and duplicate flag are always
retained because processing depends on their meaning. Repeated epochs are
flagged rather than silently deleted.

When importing CSV, preserve identifiers as text:

```python
import pandas as pd

master_df = pd.read_csv(
    "exports/master_inventory.csv",
    dtype={"Seed_id": "string"},
    low_memory=False,
)
```

## Research roadmap

1. **Inventory:** build and verify the master metadata table — implemented.
2. **Waveforms:** read daily SDS data in bounded time chunks — implemented.
3. **PPSDs:** calculate with ObsPy and save reusable NPZ data, targeting frequencies
   from Nyquist down to approximately 0.001 Hz where instrument/data support it
   — implemented; validate scientific settings on a pilot before scaling up.
4. **Quality checks:** expected-pressure, nearby-station and infragravity comparisons implemented; validate
   reference bounds and usable instrument bands before assigning scores and
   completing the other tests.
5. **Noise models:** combine accepted PPSDs and compare depth, coastal distance,
   seafloor type, sensor type and seasonal coverage.

## Project structure

```text
src/bruitfm/
  inventory.py    StationXML → master table → CSV / Excel
  archive.py      Daily SDS discovery and bounded, gap-preserving reads
  ppsd.py         Response-corrected PPSDs, checkpoints and master-file updates
  plot_ppsd.py    Distributions, time spectrograms and all PSD curves
  pressure_test.py Saved-pressure reference comparisons and master updates
  nearby_test.py  Paired daily pressure spectra, neighbor screening and plots
  infragravity_test.py Depth-dependent Godin comparison and diagnostic plots
  godin_pilot.py  Full-window Hann spectra, transient screening and pilot plots
  network.py     Automatic network discovery, staged batches and saved-plot viewing
  reference_data/ Attributed reference points and digitization provenance
docs/            Research log, methods, verified summary, references and figures
tests/           Automated checks with small synthetic inventories
.github/         GitHub Actions workflow
exports/         Local generated datasets (ignored by Git)
```

## Development

```bash
python -m pytest
ruff check .
ruff format --check .
```

GitHub Actions runs these checks on pushes and pull requests. Tests use synthetic
StationXML fixtures; the research archive is not needed to run them.
