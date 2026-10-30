# Isolated 7D and XO calibration runs

The initial runners use downloaded pressure/vertical pairs, full FDSN response
metadata, 3-hour PPSDs, base plots and all eligible Rayleigh events without tilt
correction. New XO gauge runners use exact response snapshots without PPSD.
Outputs and CSV/Excel masters stay under each new run root; the original master,
waveforms, XMLs, YV/YL results and scores are preserved. No gain is applied.

```bash
bash scripts/run_7d_all.sh
bash scripts/run_xo_all.sh
```

Run each command from the checkout. Checkpoints support resume after interruption.
The preparation scripts create header-named symlink views because many downloaded
filenames say location `00` while the actual traces have blank location. XO also
has SDS filenames with an extra `.mseed` suffix. Originals are not renamed or
rewritten, and timestamps/sample rates are not adjusted.

Mixed-rate files and files inconsistent with fetched metadata are explicitly
excluded from the prepared view, with reasons in `preparation.json`. This is
whole-file exclusion, not correction or recovery of their valid subsections.
The initial 7D preparation logged 81 excluded files across the inspected candidate
channels; XO logged zero. Original excluded files remain untouched. Future changes
to sources/scope require a fresh results directory rather than overwriting a run.

## 7D

`results/NETWORKS/7D_ALL_NO_TILT/`: 15 stations / 30 selected P/Z epochs.
Eleven stations use native LDH/LHZ at 1 Hz. J35A/J36A/J43A/J44A use BDH/BHZ at
50 Hz. The original inventory lacked full response stages, so the runner fetched
complete responses into its separate input folder.

The user completed all 30 PPSDs and all 15 calibration station jobs. Four stations
have provisional pressure multipliers, each from three accepted events:

| Station | Pressure multiplier g | Reciprocal sensitivity multiplier 1/g |
|---|---:|---:|
| J30A | 0.980613 | 1.019770 |
| J38A | 0.779660 | 1.282611 |
| J45A | 0.828663 | 1.206764 |
| J54A | 1.079624 | 0.926249 |

The other 11 stations remain unresolved. J29A/J55A/J67A have two accepted events,
below the existing three-event estimate minimum. Several other stations have
native timing offsets near 0.029 s; phase/coherence, coverage and gain-boundary
rejections also occur. No timing tolerance or gain bound was relaxed.

No verified station-matched published coefficient table was located in the
sources checked. Ruan et al. (2014), DOI 10.1002/2014JB011162, studies Cascadia
D/P ratios and sediment structure; Doran et al. (2019), DOI 10.1029/2019EA000783,
reports SCB station calibrations. Neither establishes a direct numerical gain
benchmark for these J-series results from the material reviewed. Do not substitute
Zha's YL coefficients or report a percentage agreement without matching references.

Follow-up source search: Ruan (2014) **does estimate a per-station DPG amplitude
scale factor**, jointly with sediment structure; the limitation is locating the
numeric fitted factors, not their computation. The publisher text and
[author-hosted PDF](https://www.princeton.edu/~youyir/papers/Ruan_etal-Rayleighwave-Admittance-sediment-JGR-2014.pdf)
do not provide a verified station-factor table. The author's publication page,
related Cascadia studies and indexed public GitHub/Zenodo/Figshare/MGDS records
did not yield the original factors. This search is not proof that no public copy
exists. Request the fitted factors, scaling direction and nominal responses from
the authors before a quantitative comparison.

The [author-uploaded Bell, Ruan & Forsyth (2015) follow-up](https://www.researchgate.net/publication/281845214_Shear_Velocity_Structure_of_Abyssal_Plain_Sediments_in_Cascadia),
DOI 10.1785/0220150101, also corrects the 2014 **starting** sediment profile
parameter: c=0.029, not 0.29. This is not a pressure gain or a change to the
2014 fitted inversion parameter c. Preserve this correction if constructing a
Cascadia model; no saved model/gain was changed during this search.

## XO

`results/NETWORKS/XO_ALL_NO_TILT/`: six stations / 12 selected P/Z epochs.
LA33/LA34/LT02 use EDH hydrophones and HHZ at 100 Hz. Their HDH absolute pressure
gauges run at 125 Hz and are not included in this equal-rate calibration. Thus
these results must not be labelled APG calibrations. WD48/WD66/WD70 use LDH DPGs
and LHZ at 1 Hz. The fetched response epochs correspond to AACSE 2018–2019.

The resumed XO run is complete: WD66 has pressure multiplier 0.9076893658
(sensitivity multiplier 1.1016984859; 10 accepted events); WD70 has 0.9299621720
(1.0753125559; 16 accepted events). LA33/LA34/LT02/WD48 have zero accepted events.
LA33/LA34 chiefly fail coherence, with 8/6 gain-boundary rejections. LT02 chiefly
fails coherence. WD48 chiefly lacks continuous usable event windows. Window,
mainshock-separation and timing guards also exclude events. These are not final
instrument-failure classifications.

### Separate APG run

The original preparation explicitly excluded HDH. This was an equal-rate software
choice, not evidence that APGs cannot be calibrated. The new explicit APG runner
selects **LA33/LA34/LT02 HDH (Paroscientific APG, nominal 125 Hz) / HHZ (100 Hz)**:

```bash
bash scripts/run_xo_apg.sh
```

It writes only `results/NETWORKS/XO_APG_NO_TILT/` and opens the saved calibration
plots when complete. Preparation binds full exact response epochs in `_responses`
and a separate master. Full-archive PPSD computation is not required; the saved
response snapshots are checked against the live exact epochs at fit frequencies.
No PPSD status or spectrum is fabricated. The completed EDH/DPG results remain
unchanged.

To run **all downloaded XO stations**, automatically using APG where the response
identifies APG hardware and DPG otherwise:

```bash
bash scripts/run_xo_all_gauges.sh
```

This separate runner writes `results/NETWORKS/XO_ALL_GAUGES_NO_TILT/`, skips PPSD,
uses all eligible events, and opens the calibration plots when complete. Current
downloaded stations: LA33/LA34/LT02 use HDH APG / HHZ; WD48/WD66/WD70 use LDH DPG /
LHZ. WD66/WD70's 100 Hz HDH channels are also DPGs; channel code alone does not
identify APG hardware. EDH hydrophones are not selected. Native equal-rate DPG
pairs are preferred, and all working copies use 1 Hz (32-second padding), within
the tested low-frequency passband. The three-station APG runner remains separate.
For a later expanded download requiring a fresh root, set `BRUIT_XO_RUN_ROOT`
before running; existing scope/provenance mismatches remain errors.

`--rayleigh-event-resample-rate 2` explicitly permits different native rates with
no tilt. Actual native UTC timestamps are used to grid continuous observations
at the nominal native rate, then polyphase anti-alias filtering (Kaiser beta 12)
and Lanczos interpolation (a=12) place working copies on a common 2 Hz UTC grid.
Padding is 16 seconds on both sides. Complex native responses are divided out
once in the analysis band after resampling counts; original inventories retain
their native rates. The fit band must lie below one tenth of the working rate.
Synthetic analytic P/a signals test gain and phase preservation and suppression
of high-frequency aliasing. Real gaps, conflicting overlaps and response-epoch
crossings are rejected. Default processing still refuses unequal rates.

Some APG MiniSEED records encode 125 Hz as 125.00000762939453 Hz. The runner uses
an explicit **0.1 ppm maximum native-rate tolerance**, records original header
rates in `preparation.json`, and uses actual timestamps during gridding. It does
not rewrite raw headers, samples or timestamps. Arbitrary rate differences remain
errors. The real preparation includes all three stations, six epochs and no
excluded files. A bounded first-event P/Z check produced response-corrected
spectra at all three stations; this was not a full calibration batch.

**UTC gridding is not logger clock correction.** LT02's separate APG/seismometer
logger timing remains a scientific review requirement. Polarity, borrowed
dispersion, nominal responses and scalar-gain assumptions remain provisional.

### Published XO comparison

[Deng et al. (2022), DOI 10.1029/2021GL096952](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GL096952)
uses teleseismic body waves and separate frequency bands. The text reports
WD48 coefficient about 0.95 and WD70 coefficients from 0.91 to 1.71; WD66's
band-specific values are plotted in Figure 4. LT02 is an APG validation with
pressure/seismic timing synchronization. The original EDH test is not a matched
APG comparison. The coefficient direction and exact band-specific values have
not yet been verified; supplementary download was blocked. Do not claim agreement
or calculate percentage errors from unmatched scalar/range values.

## Scientific settings and next checks

Both initial runs use an explicit provisional pressure sign +1, band 0.03–0.07 Hz,
the historical RR52 dispersion table, and unchanged QC/grid/gain safeguards.
This borrowed model is not verified for Cascadia or Alaska. Scalar estimates remain
provisional, and EDH hydrophone applicability at this band requires particular
review. Check station-specific pressure hardware, nominal response shape, polarity,
physical sensor band, depth/model and packet isolation before accepting gains.

The separate APG and all-gauge runs are now complete. The latest has four
provisional estimates and two unresolved cases; 143 tests pass. See the
[station table, stability checks and numerical provenance](calibration_network_summary.md).

Next: review saved APG event timing and unresolved coverage;
review 7D timing/phase and event sufficiency without altering timestamps; reconcile
published method/response/model conventions; coordinate core changes with Sentinel.
