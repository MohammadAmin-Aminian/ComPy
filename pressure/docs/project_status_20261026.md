# BRUIT-FM project status, 26 October 2026

Continuation guide after the move to the new computer and the first pressure
PPSD run over the whole data drive. It extends
[the 22 October status](project_status_20261022.md); the rules listed there are
unchanged. Nothing here adopts a calibration, threshold, rating or response
correction.

## New computer

* Apple M5 Pro, 15 CPU cores, 48 GB. User account `mohammad-amin` (the old one
  was `mohammadamin`), so saved products from the old computer, which record
  absolute paths, were **not** reused. Everything was recomputed.
* Code: `~/Desktop/BRUIT-FM/bruit-fm`; environment: Homebrew Miniforge,
  conda env `bruitfm` (Python 3.13); `conda init zsh` done.
* Data: USB drive `/Volumes/DATA_bruitfm` (18 TB, spinning disk), with
  `INVENTORIES/` (182 StationXML files, 29 networks) and `SDS/<year>/<net>/...`
  for 2001–2026 (about 14.5 TB). The Claude Code binary needed its own Full Disk
  Access entry to read the drive.

## Run of 23–24 October (pressure channels only)

Outputs: `~/Desktop/BRUIT-FM/fresh_run_2026-10-23/` (see its `README.md` and
`report/`). Settings unchanged: 3 h windows, 50% overlap, 0.001 Hz to Nyquist,
histogram −250 to 100 dB, gaps skipped, full response removed on the PSD.

| StationXML source | Networks |
|---|---|
| Drive `INVENTORIES/` | 25 networks |
| EarthScope FDSN download | 7D (the drive's `7D.xml` has no response stages) |
| FDSN downloads via the IRIS federator | 4N AU IM KP X9 XE YD YF YI Z4 ZZ, plus missing 7B BDH, OO BDO, ZM EDH/LDH |

Download URLs, times and SHA-256 sums are in each run's `inventories/SOURCE.txt`.

| Result | Count |
|---|---:|
| Channel epochs | 917 (32 networks, 694 stations) |
| Complete | 878 |
| Failed (all metadata) | 31 |
| Too little data / no data | 6 / 2 |
| PSD windows | 2,808,989 (about 175,600 channel-days) |

Run time 14 h 21 min, limited by the USB disk (80–120 MB/s).

## Code changes in this branch

* `network`: `--pressure-only`, `--workers max`, `--skip-uninventoried`;
  largest channels scheduled first; one BLAS thread per worker; SDS scan skips
  unselected `YEAR/NET` folders and avoids one stat call per file.
* `archive`: the reviewed `00` → `""` location alias; files indexed from their
  first and last records (29× faster; on 144 real files the start always
  matched; 12 1T files with drifting clocks ended up to 6 ms earlier, covered by
  a 1 s margin when selecting files to read); empty (0-byte) files skipped.
* `ppsd.select_epoch`: identical repeated StationXML station blocks count as one
  epoch (Z3 repeats each block 2–3 times); differing copies are still refused.

## Open issues (owner decisions)

### Failed channels by cause

| Class | Channels | Mismatch | Possible fix |
|---|---|---|---|
| A. Rate rounding in headers | 7D FN09A, FS20B, G26B HDH; 7D J43A, J65A, M01A BDH; OO AXBA1, AXCC1 HDH; YD.OBS05 BDH; Z4 CP02, CP10 HDH (11) | 0.06–0.16 ppm in all or some files, e.g. 125.0000076 Hz | A ≤1 ppm rate tolerance (a policy change: J43A and XO LA33 were blocked on purpose) |
| B. Small real rate offset | 7D.FN19A HDH (1) | 1.95 ppm, 124.99976 Hz in every file | A 2 ppm tolerance, or leave out |
| C. Wrong rate in StationXML | 1T BDH ×15 (62.0 vs 62.5 Hz); 2A.43 BDH/HDH (32/128 vs 32.768/131.072 Hz); XE.EE02 BDH (50 vs 200 Hz) (18) | 0.8–300% | Corrected StationXML from the operators |
| D. Wrong response | 4N HDH ×8; KP.KPDR HDH (9) | see below | Corrected StationXML |

0.06 ppm is about 5 µs per day, far below one sample interval; class A is a
header representation of the nominal rate rather than a different clock rate.

### 4N response (diagnosed 26 October)

The waveforms are normal int32 counts. In the FDSN StationXML, stage 1
(HTI-94-SSQ hydrophone, Pa → V) has gain −1.655 × 10⁸, which is physically
impossible: the sensor is about −165 dB re 1 V/µPa, i.e. 5.3 × 10⁻³ V/Pa. The
value looks like "−165.5 dB" entered as a number. The passband response is
therefore 1.3 × 10¹⁴ counts/Pa; the stated sensitivity (1.6 × 10⁶) is quoted at
100 Hz, above the 50 Hz Nyquist, where the FIR filters have removed the signal.
A stage-1 gain of 10^(−165.5/20) × 10⁶ V/Pa would raise all 4N PSDs by about
210 dB, to about 0 dB re 1 Pa²/Hz at 1 Hz, in line with other networks. Not
applied; the operator should confirm.

### Other

* **Constant-data windows**: 135 channels contain 3 h windows whose PSD is
  below the −250 dB floor in every bin (zero power, about −3077 dB): dead or
  zero-filled data. `bruitfm.screening` flags them; plots now show them as gaps
  and scale to real data, and report medians leave them out. The saved PPSDs
  keep every window. The QC stages do not screen yet; deciding to do so is part
  of the QC step.
* **Not computable**: 1T DH3 and UX3 (about 85 channels) have no StationXML at
  any FDSN centre.

## Reports (26 October)

Copies of the three reports, with their build scripts and run logs, are in
[`docs/reports/`](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/README.md):

* the whole-drive pressure PPSD run;
* the YV pilot (six stations, all QC tests and calibrations; RR52 reproduced
  exactly, RR09 g = 1.1937, RR40 g = 0.7426);
* the report for Wayne Crawford: 7D.J38A through every test and calibration
  (Rayleigh g = 0.7805 ± 0.016 from 16 events; microseism 0.780–0.794; local P
  from one event), and the PPSDs of 869 channels grouped by sensor type, depth,
  ocean basin and network. The J38A Godin tests used a hand-set 4 mHz lower
  band edge, because the 3 dB rule leaves this DPG no bins.

## All-station QC (26 October, provisional)

The pressure, nearby and Godin tests ran on all **846 distinct pressure epochs**
(the earlier "878 complete" counted 32 repeated Z3 rows), each in the band where
its response is within 20 dB of the sensitivity. Scores are provisional, set by
Claude at the owner's request and not yet reviewed: Test 1 in-range fraction
1 at ≥ 0.8, 0 at ≤ 0.4; Test 2 |median difference| 1 at ≤ 3 dB, 0 at ≥ 10 dB
(≥ 2 references); Test 3 |Godin residual| 1 at ≤ 10 dB, 0 at ≥ 20 dB; linear
between. Result: 456 good, 171 questionable, 189 bad, 5 nonworking, 8 4N
metadata errors, 17 7D LDI (internal, not ocean pressure). Whole-network level
offsets suggest StationXML gain problems at 2A and YD (a "dummy" DPG response,
50–70 dB low), Z4 (~55 dB low) and X9 (~24 dB high). Details:
[`reports/2026-10-26_all_station_qc/`](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-26_all_station_qc).

## Next

1. Resolve the issues above (4N response, screening rule, metadata decisions).
2. Run the `pressure`, `nearby` and `infragravity` stages on all runs.
3. Rayleigh, microseism and local-P calibration checks on every network with
   pressure and vertical channels, then ratings and the noise model.
