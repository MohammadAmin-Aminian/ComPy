# C07W / C08W / C09W matched-event review

This read-only tool compares saved results for the same earthquakes at the three
stations. It evaluates the supplied, exact-epoch nominal responses in the saved
fit band and plots pressure/acceleration PSDs, measured P/a divided by saved theory,
squared coherence and phase. It does not load raw waveform samples, recompute
spectra, fit new gains, change event acceptance or write masters.

From the checkout, create a fresh review and open it in the browser:

```bash
conda run --no-capture-output -n bruitfm python scripts/compare_yl_matched_events.py --show
```

Default selection includes every event with saved spectra at all three stations,
including QC-rejected/review events, which are labelled. Missing spectra are not
filled. For a smaller large-event review, add `--min-magnitude 7`. For an explicit
shared event, use `--event-id EVENT_ID` (repeat the flag for multiple events).

A new `results/VALIDATION/YL_MATCHED_EVENTS/run_*/` contains:

- `index.html`: plots and links to original arrival/waveform/comparison figures.
- `responses.png` and `responses.csv`: absolute response amplitudes, phase relative
  to C08W, sensitivities, stage gains, units and orientation.
- One five-panel PNG per matched event.
- `events.csv`: saved gains/statuses, normalized ratios, relative UTC bounds,
  predicted first-body-arrival flags and other-arrival IDs.
- `review.json`: input hashes, selected counts, metadata and interpretation limits.

The magnitude ≥7 verification found six shared saved events; only one is a
candidate at all three stations. Original rejected/review statuses remain intact.
All 105 protected workflow/master/response/result/artifact hashes stayed unchanged.
Ruff and direct runtime checks passed. Nominal pressure responses overlap; relative
phase differences at numerical precision are displayed with at least a ±1° range,
so numerical roundoff is not visually magnified into a physical discrepancy.

The same earthquake does not imply identical station windows, source/path effects
or packet isolation. Normalization uses each station's saved provisional theory.
Nominal response agreement does not establish hardware calibration. Review the
linked original arrival plots before interpreting a P/a discrepancy as a sensor
problem. No comparison gain is fitted or applied by this script.
