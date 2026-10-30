# Microseism gain-band phase and model audit

## Outcome

The saved-spectrum follow-up is complete. It checks all seven frozen-rule
reportable stations, not a newly selected favorable subset: **22 passing weeks,
58 eligible gain bands**. Fifty-five band/week cases meet the existing 10-degree
phase screen; three RR40 cases are approximate (10–20 degrees). No case exceeds
20 degrees in the weekly aggregate. This is a descriptive audit, **not a new
acceptance rule or an absolute calibration**. The seven v1 estimates and all
original selections are unchanged.

| Station | Weeks | Weekly bands within 10 degrees | Largest weekly RMS | Daily bands within 10 degrees | Microseism vs fixed Rayleigh |
|---|---:|---:|---:|---:|---:|
| YV.RR40 | 4 | 10/13 | 11.89 degrees | 48/91 | +0.80% |
| YV.RR52 | 3 | 9/9 | 4.36 degrees | 40/63 | +0.97% |
| YL.C05W | 3 | 6/6 | 7.61 degrees | 17/28 | -1.73% |
| YL.C09W | 3 | 8/8 | 6.57 degrees | 36/44 | +1.57% |
| 7D.J45A | 3 | 8/8 | 3.59 degrees | 43/56 | -4.14% |
| 7D.J54A | 3 | 7/7 | 4.18 degrees | 33/47 | -5.54% |
| 7D.J67A | 3 | 7/7 | 7.38 degrees | 28/42 | no baseline |

Daily fractions are not pass scores or independent statistical trials. Among
371 day/band cases, 245 are consistent, 27 approximate, five inconsistent,
76 have fewer than ten coherent bins and 18 have only one hourly window.
Single-window coherence is not interpreted as a successful check.

The three weekly RR40 flags are:

- Shift +6 weeks: 0.02710–0.03443 Hz, RMS 11.89 degrees;
  0.03443–0.04375 Hz, RMS 10.46 degrees.
- Shift +12 weeks: 0.02710–0.03443 Hz, RMS 10.66 degrees.

Exact edges are in the [machine-readable evidence](microseism_gain_band_audit.json).
All five inconsistent daily cases are RR52: 14 June, 23/28 July and 2/8 September
2013, with RMS approximately 20.7–27.6 degrees in individual gain bands.
Strong weekly agreement therefore does not mean every day is well described.
Do not remove these days retrospectively to improve the frozen estimates.

## What was checked

The explicit input is [microseism_weeks.json](microseism_weeks.json), produced
from the completed 24-run review. Source station selections and band estimates
are reproduced, station/week identities and non-overlap are checked, and saved
spectral means, transfer and coherence arithmetic are verified.

Phase is calculated from the original complex NPZ spectra, not rounded CSV
phases. Each original gain interval is half-open, `low <= f < high`, below the
first water resonance. The existing whole-band hypothesis is prescribed; it
is **not refitted per band**. Alternatives are reported diagnostically only.
Pressure and vertical pole-reflection terms enter with opposite signs in
`H = D/P`. Nominal and reversed labels represent relative phase conventions,
not independently verified hardware polarity.

The legacy descriptive settings are retained: coherence squared >=0.5,
at least ten selected frequency bins, RMS <=10 degrees for consistency and
10–20 degrees for approximate agreement. These are not newly justified scoring
thresholds. Daily results are formed by averaging the original complex hourly
cross spectra and auto spectra first, never averaging phase or single-window
coherence values. Adjacent frequency bins and days are not independent.

Full complex pressure and displacement responses reproduce from the exact-epoch
XML snapshots with zero numerical discrepancy on the saved grid. All seven
vertical channels have metadata dip -90 degrees. This verifies the recorded
processing convention, not the physical orientation, clock or absolute gain.

## 7D model/band discrepancy

On exactly the eligible bins, changing the cached historical mode-0 curve from
water-loaded to dry changes the band-median H1 estimates by at most 0.32%
(J45A), 0.59% (J54A) and 0.28% (J67A). The corresponding H2 changes are at most
0.38%, 0.56% and 0.41%. Thus that particular model choice does not explain the
4–6% station-median offset from Rayleigh. This is not an uncertainty bound on
all Earth models, modal mixtures, depth, response or water assumptions.

Rayleigh lies inside the model H1–H2 band envelope in only 2/8 J45A and 1/7 J54A
band/week cases. However, the *pooled* frozen envelopes still contain those
references. Neither comparison is a confidence interval or a statistical
rejection of Rayleigh. The cached structure-free upper limits are reported
separately. No model is inverted to force agreement.

Coverage also matters: J54A's +2 week has 73 retained hours; J67A's base/+1 weeks
have 98/74. Some days contribute only one hour. The frozen rule has no newly
added coverage threshold, so its original labels remain, with these limitations
explicitly attached. Physical causes of the offset are unresolved.

## Tilt and vertical-response identifiability

Writing nominal traces as `P_nom = P_true/gP` and `Z_nom = Z_true/gZ`, an ideal
mode has `H_nom = (gP/gZ) * eta`. A pressure-coherent vertical contaminant adds
another complex term. A real contribution can change amplitude while retaining
near-zero or near-180-degree phase and high coherence. A regression test
explicitly demonstrates this blind spot. Phase consistency cannot certify a
pressure gain or exclude coherent tilt.

[Ruan et al. (2014)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JB011162)
establish the structure/mode dependence of the microseism displacement/pressure
ratio. [Bell, Forsyth & Ruan (2015)](https://doi.org/10.1785/0120140054)
use horizontal channels to assess/remove tilt and describe its temporal
variation. These support a bounded P/H1/H2/Z investigation, not automatic
pressure-based subtraction from Z: that could remove the calibration signal.

A scalar error in vertical sensitivity is also shared by pressure–vertical
Rayleigh calibration; it is not automatically an explanation for the difference
between the methods. Absolute vertical response still needs an independent
constraint, and frequency/time-dependent errors remain possible.

Header-only preflight confirms local BDH/BH1/BH2/BHZ daily files on the same
62.5 Hz grid for RR52 on 10 June 2013 and RR40 on 8 June 2013. This is not a
gap/quality audit or a processed P/H/Z product. The configured 7D master has
only P/Z rows for these stations; horizontal availability outside that configured
input set has not been established. No data were downloaded.

## Table integration and continuation

**Completed continuation:** the [bounded RR52 horizontal pilot](rr52_microseism_horizontal_pilot.md)
now checks two days. It identifies strong low-band horizontal-associated motion,
ambiguous horizontal azimuth metadata and two midnight grid refusals. It does
not identify a physical tilt correction or validate an absolute gain. The
original next-pilot plan below is retained as context.

[microseism_table_staging.csv](microseism_table_staging.csv) contains **89 unique
full pressure-channel epoch keys**, including inventory source. It is a staging
export only; no production master is changed. Status counts are:

| Research status | Rows |
|---|---:|
| `conditional_repeatable_unvalidated` | 7 |
| `insufficient_passing_weeks` | 35 |
| `no_passing_weeks` | 42 |
| `unresolved_pairing` | 5 |

All 89 existing-schema microseism score fields are blank. Only the seven
reportable conditional estimates have a numeric gain. Phase flags, pairing
reasons, provenance and unresolved tilt/vertical-response assumptions remain
separate fields. These statuses are not adopted pass/fail classifications.
The five pairing-unresolved rows are a more specific subset of the original
47 no-passing-week rows, not a change to v1 results.

Next is one bounded P/H1/H2/Z pilot: RR52's base week, with 10 June and flagged
14 June retained for comparison, and RR40's base interval as a control. Freeze
the diagnostic procedure before inspecting new amplitudes. Verify full response,
timing and gap handling and distinguish physical horizontal Rayleigh motion from
tilt leakage; high horizontal coherence alone does not identify tilt. Do not
reselect weeks, fit toward Rayleigh, overwrite spectra, or apply a gain.
Station-specific model/response evidence and final scoring criteria remain
necessary before claiming an accepted absolute calibration. C07W time-state
analysis remains separate; its retracted single-week claim stays retracted.

## Reproduce the read-only audit

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m scripts.audit_microseism_gain_bands \
  --review docs/microseism_weeks.json
```

This is an audit, not a calibration batch. It writes a fresh folder under
`results/REVIEWS/MICROSEISM_GAIN_BANDS/`; omit `--docs` to preserve the curated
snapshot. The current full audit is
`run_20261015T144203Z_90633bfe/audit.json`; its SHA-256 is in the curated evidence.
It verifies **2,466 protected files** before and after, with zero changed hashes.
No raw waveform samples, new modal solution or calibration batch were processed;
no gain, response, master, original result or score was changed.
