# Reviewing RR52 event windows and timing

The runner proposes a Rayleigh window from group velocity. It does not claim
that the largest waveform peak is a verified Rayleigh arrival. Review each
event's `*_arrival_review.png` before treating its gain as calibration evidence.

1. Inspect the wide, independently timed pressure and upward-acceleration
   waveforms. Dashed markers are predicted first P/S arrivals (TauP iasp91);
   the group-arrival range comes from the supplied dispersion table across
   0.02–0.06 Hz. The purple region is the proposed analysis window. Gray regions
   identify other predicted Mw≥6 packets, including aftershocks.
2. Compare the observed packet with those predictions and surrounding signal.
   Prefer an isolated packet supported by both channels and an appropriate
   dispersion model. A P pick locates a body wave; it does not alone locate the
   Rayleigh packet. The display is band-limited and Hann-tapered; no fitted
   amplitude gain is applied. RR52 pressure is shown with the provisional −1 sign.
3. Retain, change or exclude a window with a written reason. Review boundaries
   using the waveform and arrival physics, before inspecting the resulting gain;
   do not repeatedly move them to maximize coherence or reproduce 0.67.
4. Rerun the selected event and compare the saved waveform, phase/coherence and
   gain with the original. A manual window does not bypass overlap or coherence
   checks. An amplitude fit alone does not certify an isolated fundamental mode.

The [Santa Cruz overview](figures/rr52/rayleigh_tiskitpy_arrival_review.png)
shows a predicted overlapping packet. Its amplitude-only fit is retained for
review, but it is excluded from the accepted-event count.

## Record an explicit window decision

Each run writes `window_review_template.json` beside `comparison.json`.
Copy that file to a new name such as
`results/RAYLEIGH_REVIEW/rr52_reviewed_windows.json`. Entries are keyed by event ID:

```json
{
  "usc000f1s0": {
    "decision": "exclude",
    "note": "Predicted aftershock packet overlaps the selected Rayleigh packet."
  }
}
```

`pending` retains the proposed window. `exclude` skips the event. `use` requires
both `start_utc` and `end_utc` plus a nonempty `note` explaining the review.
Choose those UTC bounds from the overview; no revised bounds have been invented
for Santa Cruz. A reviewed window must be after the origin, end within 12 hours,
span at most two hours, and contain at least eight 200-second segments (1600 s).
These are current pilot constraints, not universal scientific requirements.

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m bruitfm.rayleigh_calibration \
  --catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --network YV --station RR52 --event-id usc000f1s0 --max-events 1 \
  --window-overrides results/RAYLEIGH_REVIEW/rr52_reviewed_windows.json
```

Every run creates a separate output folder and hashes the review input. The
default tilt method is TiSKitPy; RR52's default pressure sign is provisionally −1.
Use `--pressure-sign 1` only to make an explicit alternate-convention comparison.

## Timing review

The two September events fail the strict sample-alignment guard during quiet
training. BH1 starts **0.0001 s** after the other channels (0.00625 samples at
62.5 Hz). `comparison.json` records each channel's actual start, the offset and
the phase bound `360*f*Δt`: only **0.00216° at 0.06 Hz**. Such an offset cannot
explain the earlier approximately 180° sign difference.

To resolve it, compare MiniSEED record-start timestamps and blockette timing
fields across adjacent files and seismic components; determine whether this is
timestamp quantization or a genuine channel offset. Then assess an explicit,
documented alignment policy against phase/gain uncertainty and synthetic known
transfers. This pilot deliberately leaves the events unresolved until that
policy is justified. It neither rounds timestamps silently nor interpolates
gaps. Clock offsets are separate from the provisional pressure-sign assumption.

The later YV review also found fractional **within-channel** record boundaries
in RR11/RR25 September quiet training. Rayleigh reads now reject them before
ObsPy merge can align records. Exact aligned duplicate samples remain permitted.
This additional guard leaves previous outputs intact; see the
[independent review](yv_rayleigh_independent_validation.md) for offsets and phase bounds.

## Open the saved plots without recomputing

```bash
cd /Users/mohammadamin/Projects/bruit-fm
open docs/figures/rr52/rayleigh_tiskitpy_arrival_review.png
open docs/figures/rr52/rayleigh_tiskitpy_comparison.png
open docs/figures/rr52/rayleigh_tiskitpy_waveform.png
```

The tracked examples show Santa Cruz. Other event plots and spectra are in the
local run folder identified by [the summary](rr52_tiskitpy_validation_summary.json).
