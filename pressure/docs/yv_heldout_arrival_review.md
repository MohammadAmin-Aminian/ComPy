# YV flagged-arrival review and RR48 replacement

Follow-up to the [frozen-gain study](yv_heldout_results.md).
Three bounded native P/Z overviews were read, retaining continuity/grid/epoch
guards. No gain was fitted/applied, original QC label changed or full calibration
batch repeated. New display tool: `scripts/review_yv_heldout_arrivals.py`.
[Selection and provenance](yv_heldout_arrival_review.json) ·
[RR05 summary CSV](yv_heldout_rr05_replacement.csv).

## Waveform observations

| Station/event | Observation | Decision |
|---|---|---|
| RR01, `usb000k1gb`, Mw 6.8 | Raw body-wave energy begins near P; strongest calibration-band P/Z packet is later, around/after the reference group arrival | Retain arrival-review status |
| RR06, `usb000k1gb`, Mw 6.8 | P/S occur early inside the original window; strongest band-limited packet is later | Retain arrival-review status |
| RR06, `usb000kdb4`, Mw 7.1 | S/SKS are near the window start; a smaller band-limited pulse precedes a stronger later packet continuing beyond the original window | Retain arrival-review status |

Plots: [RR01 September](figures/yv_heldout_arrival_review/RR01_usb000k1gb_waveforms.png),
[RR06 September](figures/yv_heldout_arrival_review/RR06_usb000k1gb_waveforms.png),
[RR06 October](figures/yv_heldout_arrival_review/RR06_usb000kdb4_waveforms.png).
Left panels show one-second min/max raw-count envelopes. Right panels show
response-corrected 0.03–0.07 Hz pressure and upward acceleration, with the trial
pole phase reflection in memory and the training pressure sign. No scalar gain
is applied. The 10% Tukey taper and band-limiting are display operations,
distinct from the original 200-second Hann spectral estimator.

The flag is conservative: predicted body arrivals somewhere in a window do not
prove material contamination in the calibration band. These plots support
investigating isolation, but do not prove pure fundamental-mode Rayleigh waves.
Residuals were already viewed; manual window changes on these cases would be
exploratory, not a fresh untouched held-out test. No events were readmitted.

## Replacement selected using coverage before new errors

Four preselected newly estimated stations with long epochs were screened using
waveform headers, with RR48 as comparator. Not every remaining station was ranked.

| Station | Joint epoch end (2013 UTC) | Epoch-valid, P/Z header-covered windows |
|---|---|---:|
| RR05 | November 3, 00:35:38 | **42** |
| RR26 | November 2, 07:24:09 | 41 |
| RR20 | October 30, 19:51:39 | 40 |
| RR16 | October 28, 22:38:16 | 40 |
| RR48 | June 10, 23:28:43 | 9 |

Rule: maximize header-covered, epoch-valid unused windows; break ties by latest
joint epoch end, then alphabetical station ID. RR05's selection was saved
before its new samples/errors were examined. Headers do not certify grid
alignment, coherence or isolation. RR55 was excluded because response identity
remains unclear.

## Completed RR05 test

Run: `results/VALIDATION/YV_HELDOUT/RR05_replacement_20261013T215958Z_fee6a2d2`.
The same 50 previously selected mainshocks, excluded from training, were all
attempted at this new station. **Frozen g = 1.1201435166537377**, from four
pole-only/original-fit training events. No error-based stopping or refitting.

- **Zero QC-usable or isolated events**; no accepted median prediction error.
- **40 coherence-review cases**.
- **Ten input/window refusals:** five outside the epoch, three windows starting
  before origin, two pressure grid offsets (0.0023/0.0024 s).

For Mw 7.1 `usb000kdb4`, median coherence squared is 0.957854 but mean 0.941675.
The unchanged rule requires both ≥0.95. Several middle frequencies are coherent,
while band edges reduce the mean. Other in-epoch Mw≥7 events have means 0.901809
and 0.478228. These are saved-spectrum diagnostics; no bins were removed to
rescue a result. More coverage does not guarantee usable calibration evidence
or show that the saved gain is wrong.

## Preservation, verification and next task

Post-run rechecks matched all protected inputs:

| Check | SHA-256 inputs | Raw-file size/mtime records |
|---|---:|---:|
| Coverage selection | 66 | 2,754 |
| Three waveform reviews | 61 | 1,228 |
| RR05 validation | 62 | 612 |
| Original four-station study | 71 | 2,208 |

These overlap; do not add them as unique files. Raw byte hashes are not claimed.
Waveforms, XMLs, snapshots, masters, original results and scores are unchanged.
Ruff check/format, all three real-window plotter checks and six targeted
held-out tests passed. Production
calibration code/defaults did not change; the previous 173-test checkpoint stands.

Next inspect pressure/vertical signal and noise across saved frequency bins and
how fixed group-arrival windows capture packets. Develop any new window rule
from arrivals/coverage/noise evidence; retain original results and treat reviewed
cases as development data. Validate prospectively on untouched station/event
data. Do not repeatedly select stations, bands or windows until errors improve,
or lower QC to increase acceptance.

Reproduce the three plots in a new folder:

```bash
PYTHONPATH=src:. conda run --no-capture-output -n bruitfm python \
  scripts/review_yv_heldout_arrivals.py \
  --validation-dir results/VALIDATION/YV_HELDOUT/run_20261013T213126Z_87e24fe6
```

RR05 is already finished. View saved plots without recomputation:

```bash
open results/VALIDATION/YV_HELDOUT/RR05_replacement_20261013T215958Z_fee6a2d2/heldout_summary.png
open docs/figures/yv_heldout_arrival_review/*.png
```
