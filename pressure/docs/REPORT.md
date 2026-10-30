# BRUIT-FM v1.2.0: results of the all-network campaign

**Ocean-bottom pressure-gauge noise, quality control and calibration across 30 networks**

Mohammad-Amin Aminian · Institut de Physique du Globe de Paris (IPGP) · October 2026

> **Status.** All numbers are outputs of BRUIT-FM v1.2.0 with *provisional*
> thresholds. No gain or response correction has been applied to any data. The
> thresholds, metadata corrections and station acceptance still need to be
> agreed before the results are used to build noise models.
>
> **Gain-corrected version.** Everything in this report is also redone after
> multiplying each pressure channel by its calibration gain (265 channels
> corrected), in
> [`reports/2026-10-28_gain_corrected/`](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_gain_corrected/README.md).
> This report itself shows the uncorrected results.

---

## Contents

1. [Summary](#1-summary)
2. [Data](#2-data)
3. [Processing](#3-processing)
4. [Pressure noise across the networks](#4-pressure-noise-across-the-networks)
5. [Quality-control tests](#5-quality-control-tests)
6. [Calibration checks](#6-calibration-checks)
7. [Worked example: one station through every check](#7-worked-example-one-station-through-every-check)
8. [Provisional noise model](#8-provisional-noise-model)
9. [Provisional decisions: policy, metadata and acceptance](#9-provisional-decisions-policy-metadata-and-acceptance)
10. [Limitations and next steps](#10-limitations-and-next-steps)
11. [Reproducing these results](#11-reproducing-these-results)
12. [References](#12-references)

---

## 1. Summary

- **Pressure PPSDs.** 846 unique pressure channel epochs from 30 networks have
  response-corrected PPSDs. They cover the Pacific, Atlantic, Indian, Southern
  and Arctic oceans and the Mediterranean, at water depths from 52 m to 5788 m
  (median 2779 m).
- **Seismometer PPSDs.** 1,345 seismometer channel epochs (457 vertical, 888
  horizontal) at 378 of the pressure stations also have PPSDs, so pressure and
  ground motion can be compared station by station.
- **Three spectral QC tests on every pressure channel.** **456 good, 171
  questionable, 189 bad**. Thirty more could not be classed: 5 are mostly dead
  data, 8 have StationXML that cannot be tested and 17 measure pressure inside
  the instrument.
- **Most "bad" channels look like metadata problems, not broken instruments.**
  Whole networks fail together, offset by tens of dB from comparable stations.
  Examples are 2A and YD (a "dummy" DPG response, 50–70 dB low), Z4 (about
  55 dB low) and X9 (about 24 dB high).
- **Three calibration checks on every pressure channel with a vertical
  seismometer.** Rayleigh waves, microseism admittance and local P waves were
  checked. 215 channels have at least one scored estimate. Where Rayleigh and
  microseism can both be scored (43 channels), their median difference is
  +0.1 dB and 34 agree within 3 dB. At the reference station YV.RR52 the network
  run gives g = 0.682 (published 0.67). At 7D.J38A, the full-event Rayleigh and
  microseism estimates agree within 1 % (0.781 and 0.780–0.794).
- **Final class after all six checks.** **439 good, 188 questionable, 189
  bad.**
- **Provisional decisions and noise model.** Under a frozen provisional policy,
  396 channels at 311 stations are accepted. Five networks are flagged for
  metadata: 2A, YD and Z4 sit 47–59 dB low, X9 17 dB high, and 4N cannot be
  evaluated (§9). The station-balanced noise model of the accepted stations is
  given by ocean basin, sensor type, water depth and local season (§8).
- **Observed infragravity levels lie above the Godin et al. (2013) model.** The
  excess is about +10 to +20 dB in the Pacific and Indian oceans, consistent
  with Webb (1998).

## 2. Data

| Item | Value |
|---|---|
| Archive | 18 TB SDS archive of ocean-bottom deployments |
| Metadata | StationXML from the archive. Downloaded from EarthScope/FDSN when the archive copy lacked response stages (7D, 4N, AU, IM, X9, XE, YD, YF, YI, Z4, ZZ) |
| Pressure channels | 846 unique epochs: 357 differential pressure gauges (DPG), 50 absolute pressure gauges (APG), 422 hydrophones, 17 internal pressure (LDI) |
| Seismometer channels | 1,345 epochs (457 vertical, 888 horizontal) at 378 pressure stations: Trillium Compact, Trillium 240, Güralp, geophones and others |
| Networks (30) | 1L 1T 2A 2D 4N 5F 7B 7D 7K 8J AU IM X3 X9 XE XI XW YD YF YI YO YR YV Z3 Z4 ZM ZS ZU ZX ZZ |

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/network/basin_map.png?raw=true" width="90%"
       alt="Pressure stations coloured by ocean basin">
  <br><em>Figure 1. Pressure stations by ocean basin.</em>
</p>

## 3. Processing

**PPSDs.** Each channel epoch is read day by day from the SDS archive and
processed with ObsPy's `PPSD`:

- 3-hour windows with 50 % overlap;
- a lower frequency limit of 1 mHz;
- the full StationXML response removed from the spectrum, giving Pa²/Hz for
  pressure and (m/s²)²/Hz for seismometers.

Windows of constant (dead) data are counted and reported. Each result keeps the
exact response used (`inventory.xml`) and a provenance file with the code hash
and settings.

**Usable band.** A channel is tested only where its response lies within 20 dB
of its sensitivity. Outside that band the spectrum is dominated by the response
correction, not the ocean. A hydrophone with no response at 1–20 mHz is
therefore *not scored* on the infragravity test rather than failed.

**Scale.** The campaign ran on one workstation with every CPU in use. The
seismometer PPSDs alone took about 12 hours. Every stage is checkpointed, and
the run resumed cleanly after a reboot and several disconnections of the data
drive.

## 4. Pressure noise across the networks

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/network/by_ocean_basin.png?raw=true" width="100%"
       alt="Pressure PSDs grouped by ocean basin">
  <br><em>Figure 2. Median pressure PSD of every channel, grouped by ocean basin.
  Black: group median. Red bars: Webb (1998) Fig. 2 ranges (≤ 0.2 Hz). Green band:
  Brown et al. (2014) range (≥ 0.2 Hz).</em>
</p>

The groups differ in the way Webb (1998) describes:

- the Pacific has the strongest microseisms;
- the Mediterranean is the quietest basin;
- the noise notch between the infragravity band and the microseism peak
  (about 0.03–0.1 Hz) deepens with water depth (Figure 3).

The few curves that rise steeply above 1 Hz or sit tens of dB below the group
are the metadata cases discussed in §9 (step 2).

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/network/by_depth_medians.png?raw=true" width="85%"
       alt="Median pressure PSD by water depth">
  <br><em>Figure 3. Median pressure PSD by water-depth class. The infragravity
  band (< 20 mHz) falls with depth as surface gravity waves decay.</em>
</p>

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/network/by_sensor_type_medians.png?raw=true" width="85%"
       alt="Median pressure PSD by sensor type">
  <br><em>Figure 4. Median pressure PSD by sensor type (DPG, APG, hydrophone).</em>
</p>

## 5. Quality-control tests

Each pressure channel gets three independent scores between 0 and 1. Between
the two bounds the score is linear.

| Test | Metric | Score 1 | Score 0 | Not scored when |
|---|---|---|---|---|
| **1 · Expected levels** | Fraction of frequency bins whose median lies within Webb (1998) (≤ 0.2 Hz) or Brown et al. (2014) (≥ 0.2 Hz) | ≥ 0.8 | ≤ 0.4 | < 3 bins in the usable band |
| **2 · Nearby stations** | Median of daily differences from the median of nearby stations (≤ 100 km, depth within 10 %, ≥ 30 common days) | ≤ 3 dB | ≥ 10 dB | < 2 reference stations |
| **3 · Infragravity** | Median residual to Godin et al. (2013) Eq. 6 at the station's depth, 1–20 mHz | ≤ 10 dB | ≥ 20 dB | < 3 bins in the usable band |

The class of each channel is set from its scores:

| Class | Rule |
|---|---|
| **good** | mean score ≥ 0.75 and no test at 0 |
| **bad** | mean score < 0.4 |
| **questionable** | anything in between |
| **nonworking** | ≥ 50 % of windows are dead data |

### Results

| Test | Pass | Partial | Fail | Not scored |
|---|---:|---:|---:|---:|
| 1 · Expected levels | 369 | 269 | 183 | 25 |
| 2 · Nearby stations | 216 | 37 | 26 | 567 |
| 3 · Infragravity (Godin) | 331 | 112 | 97 | 306 |

Two thirds of the channels cannot be scored on Test 2. Most deployments have no
second station within 100 km, at a similar depth and with 30 overlapping days.

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_all_station_qc/report/img/metrics.png?raw=true" width="100%"
       alt="Distributions of the three QC metrics">
  <br><em>Figure 5. Distributions of the three QC metrics. Green dashes: score 1
  inside. Red dashes: score 0 outside. The pile at −40 dB in Test 3 is 2A, YD and
  Z4, which sit about 50 dB below the model.</em>
</p>

| Sensor type | good | questionable | bad | nonworking | metadata error | not ocean pressure |
|---|---:|---:|---:|---:|---:|---:|
| APG (absolute) | 8 | 41 | 1 | 0 | 0 | 0 |
| DPG | 224 | 64 | 65 | 4 | 0 | 0 |
| Hydrophone | 224 | 66 | 123 | 1 | 8 | 0 |
| Internal pressure (LDI) | 0 | 0 | 0 | 0 | 0 | 17 |
| **Total** | **456** | **171** | **189** | **5** | **8** | **17** |

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_all_station_qc/report/img/class_map.png?raw=true" width="90%"
       alt="Map of channels coloured by QC class">
  <br><em>Figure 6. One dot per channel, coloured by class. Red clusters are
  mostly whole networks: 2A (East Pacific Rise), YD (Ionian Sea), YR (Oregon
  margin), 1T (Mayotte) and Z4 (off central Chile).</em>
</p>

<details>
<summary><b>Results by network</b> (click to expand)</summary>

| Network | Channels | Good | Questionable | Bad | Other | Mean score | Note |
|---|---:|---:|---:|---:|---:|---:|---|
| 1L | 6 | 5 | 1 | 0 | 0 | 0.95 | |
| 1T | 92 | 25 | 4 | 63 | 0 | 0.32 | Microseism peak ~34 dB low, 5 Hz floor ~20 dB high: response shape suspect |
| 2A | 32 | 0 | 0 | 32 | 0 | 0.01 | ~50 dB low at 0.15 Hz; StationXML sensor is a "dummy" LDEO DPG |
| 2D | 2 | 0 | 2 | 0 | 0 | 0.50 | |
| 4N | 8 | 0 | 0 | 0 | 8 | | Known StationXML defect (stage-1 gain −1.655e8) |
| 5F | 3 | 3 | 0 | 0 | 0 | 0.96 | |
| 7B | 30 | 11 | 19 | 0 | 0 | 0.72 | |
| 7D | 222 | 130 | 61 | 11 | 20 | 0.78 | 17 LDI internal-pressure channels |
| 7K | 12 | 12 | 0 | 0 | 0 | 0.99 | |
| 8J | 41 | 33 | 5 | 3 | 0 | 0.87 | |
| AU | 3 | 2 | 0 | 1 | 0 | 0.67 | |
| IM | 18 | 15 | 2 | 1 | 0 | 0.87 | |
| X3 | 7 | 7 | 0 | 0 | 0 | 0.87 | |
| X9 | 6 | 0 | 0 | 6 | 0 | 0.00 | SIO DPG ~24 dB high across the band: gain suspect |
| XE | 29 | 23 | 3 | 2 | 1 | 0.82 | |
| XI | 26 | 26 | 0 | 0 | 0 | 1.00 | |
| XW | 51 | 12 | 29 | 10 | 0 | 0.54 | |
| YD | 6 | 0 | 0 | 6 | 0 | 0.00 | Same "dummy" LDEO DPG response; ~70 dB low at 0.15 Hz |
| YF | 9 | 9 | 0 | 0 | 0 | 1.00 | |
| YI | 5 | 5 | 0 | 0 | 0 | 1.00 | |
| YO | 30 | 12 | 16 | 2 | 0 | 0.71 | |
| YR | 50 | 3 | 15 | 32 | 0 | 0.31 | Shallow (~1.4 km) Oregon margin; mostly fails Godin and expected levels |
| YV | 44 | 38 | 3 | 2 | 1 | 0.86 | |
| Z3 | 22 | 15 | 5 | 2 | 0 | 0.76 | |
| Z4 | 7 | 0 | 3 | 4 | 0 | 0.21 | LDEO APG ~55 dB low at 0.15 Hz: response gain suspect |
| ZM | 34 | 23 | 2 | 9 | 0 | 0.72 | |
| ZS | 4 | 3 | 1 | 0 | 0 | 0.90 | |
| ZU | 26 | 25 | 0 | 1 | 0 | 0.96 | |
| ZX | 15 | 13 | 0 | 2 | 0 | 0.87 | |
| ZZ | 6 | 6 | 0 | 0 | 0 | 1.00 | |

</details>

## 6. Calibration checks

Each check estimates a pressure multiplier *g*. A value of 1 means the
StationXML response agrees with the method; 0.5 means the response overstates
the pressure by a factor of 2. Every check needs a co-located vertical
seismometer, and each was run on every pressure channel that has one. **No *g*
is applied to the data.**

| Method | Physics | Band | Station estimate needs |
|---|---|---|---|
| **Rayleigh waves** | Pressure/vertical-acceleration ratio of teleseismic M ≥ 5.5 Rayleigh waves against a dispersion model (historical RR52 dispersion for every network). TiSKitPy tilt cleaning where both horizontals exist | 30–70 mHz | ≥ 3 accepted events (bounded sampling, at most 8) |
| **Microseisms** | Pressure/vertical admittance of coherent microseism energy through a water layer of the station's depth (one week of hourly pairs) | best coherent sub-band, mostly 0.03–0.1 Hz | squared coherence ≥ 0.5 |
| **Local P waves** | Pressure/vertical ratio of the P arrival of local earthquakes (USGS, M ≥ 2.5, ≤ 320 km) | 1–16 Hz, by octave | ≥ 3 events in a band |

The provisional score uses e = |20 log₁₀ g|. It is 1 for e ≤ 1.5 dB (within
about ±19 %) and 0 for e ≥ 6 dB (a factor of 2), linear in between. Where both
pressure signs were possible, the sign with more accepted events is kept.

### Results

| Method | Channels with an estimate | Pass (≤ 1.5 dB) | Partial | Fail (≥ 6 dB) | Not scored |
|---|---:|---:|---:|---:|---:|
| Rayleigh waves | 54 | 20 | 26 | 8 | 792 |
| Microseisms | 309 (167 with coh² ≥ 0.5) | 52 | 70 | 45 | 679 |
| Local P waves | 54 | 23 | 18 | 13 | 792 |

In all, 215 of the 846 pressure channels have at least one scored calibration.

**The methods agree.** For the 43 channels with both a Rayleigh and a scored
microseism estimate, the median difference is **+0.1 dB**. 34 of them agree
within 3 dB and 41 within 6 dB (Figure 9). The Rayleigh estimates cluster
around g ≈ 1.0, with 10th–90th percentiles of 0.72–1.46. At the YV reference
stations the network run gives:

| Station | Rayleigh g | Microseism g (coh²) | Local P g | Published |
|---|---:|---:|---:|---:|
| YV.RR52 (DPG) | 0.682 | 0.628 (0.82) | — | 0.67 |
| YV.RR40 (DPG) | 0.774 | 0.921 (0.66) | 0.52 | 0.74 |
| YV.RR09 (hydrophone) | 1.206 | 1.353 (0.73) | — | — |

RR52 matches the published 0.67 within 2 %. The single-station pilot gave
0.6787 with the same method.

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_seismic_calibration/report/img/calibration/gain_histograms.png?raw=true" width="100%"
       alt="Distributions of the calibration multipliers">
  <br><em>Figure 8. Calibration multipliers in dB (0 = StationXML response).
  Green dashes: score 1 inside. Red dashes: score 0 outside. Values beyond ±12 dB
  are piled at the ends. 104 of the 127 microseism values in the piles come from
  low-coherence weeks (coh² < 0.5) and are not scored.</em>
</p>

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_seismic_calibration/report/img/calibration/method_comparison.png?raw=true" width="100%"
       alt="Calibration methods compared channel by channel">
  <br><em>Figure 9. The three methods compared channel by channel. Black line:
  equal gains.</em>
</p>

**Why most channels have no calibration.** Every method needs a usable vertical
seismometer recorded alongside the pressure channel:

- **Geophones only.** XW, YR, X3 and part of 1T carry only short-period geophones
  (channel codes `EL?`/`DH?`). These cannot record 30–70 mHz Rayleigh waves and
  were not part of the seismic PPSD run.
- **Different sample rates.** Examples are 7D BDH at 40 Hz against BHZ at 50 Hz,
  and YO BDH at 20 Hz against HHZ at 100 Hz. The methods refuse to resample
  implicitly. The same sensor's 1 Hz channel (for example 7D LDH with LHZ) is
  calibrated instead.
- **Missing metadata or data.** These are verticals without a saved response,
  stations with several overlapping vertical epochs, and weeks without
  synchronised data.
- **Too few events.** 202 channels ran Rayleigh but had fewer than 3 accepted
  events.

**Final class.** Combining the three QC tests and the three calibration scores
gives **439 good, 188 questionable and 189 bad**, plus 5 nonworking, 8 with a
metadata error and 17 that are not ocean pressure. Compared with the QC-only
classes, 29 channels moved from good to questionable and 12 from questionable
to good.

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_seismic_calibration/report/img/calibration/final_class_map.png?raw=true" width="90%"
       alt="Map of final classes">
  <br><em>Figure 10. Final class of every pressure channel after the QC tests and
  calibration checks.</em>
</p>

### Seismometer noise

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_seismic_calibration/report/img/seismic/by_ocean_basin_vertical_medians.png?raw=true" width="85%"
       alt="Vertical seismometer noise by ocean basin">
  <br><em>Figure 11. Median vertical-acceleration PSD by ocean basin, with the
  Peterson (1993) NLNM and NHNM (dashed). 1,345 seismometer channel epochs (457
  vertical, 888 horizontal) at 378 stations in 20 networks have PPSDs.</em>
</p>

The Indian Ocean verticals are the noisiest below 0.05 Hz, by about 20 dB over
the other basins. The Atlantic verticals are the quietest in the 0.02–0.06 Hz
noise notch.

**Tilt correction** was tested on two pilot stations (YV.RR52 and 2A.28) and then
not run network-wide, by decision of the project owner. At RR52, ComPy cleaning
in 6-hour windows removes about 14 dB of vertical noise at 2–10 mHz. The method
is available as `bruitfm-tilt`.

The full calibration report, with every figure, is in
[`reports/2026-10-28_seismic_calibration/`](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-28_seismic_calibration)
([PDF](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_seismic_calibration/report/BRUIT-FM_seismic_calibration_report.pdf)).
The completed table, with one row per pressure epoch and every score and
metric, is
[`table/pressure_qc_calibration_table.csv`](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_seismic_calibration/table/pressure_qc_calibration_table.csv).

## 7. Worked example: one station through every check

This section takes one well-recorded station, **7D.J38A**, through each QC test
and each calibration method, with one figure per check. J38A is a Cascadia
Initiative SIO differential pressure gauge (DPG) at 2734 m on the Juan de Fuca
plate. It recorded from 23 November 2011 to 17 May 2012, giving 2,816
three-hour PSD windows with no dead windows. Its final class is **good**.

J38A's local P result rests on a single event. The local P example therefore
uses **YV.RR40** (INSU DPG, 4780 m, RHUM-RUM), where six events give a
station estimate.

| Check | Station | Result | Verdict |
|---|---|---|---|
| QC 1 · Expected levels | J38A | 62.5 % of 40 bins within the limits; all misses above Brown at 0.22–0.74 Hz | partial |
| QC 2 · Nearby stations | J38A | 6 neighbours, 176 days: median +0.6 dB, RMS 1.8 dB | pass |
| QC 3 · Infragravity | J38A | +14.6 dB above Godin Eq. 6 (11 bins, 6–14 mHz) | partial |
| Calibration 1 · Rayleigh | J38A | g = 0.7805 ± 0.016 from 16 events | consistent with microseism |
| Calibration 2 · Microseism | J38A | g = 0.780–0.794, coh² 0.975 | contains the Rayleigh value |
| Calibration 3 · Local P | RR40 | g = 0.75 at 2–4 Hz (6 events); published 0.74 | consistent at 2–4 Hz |

### QC test 1: expected pressure levels

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/j38a_test1_levels.png?raw=true" width="85%"
       alt="J38A pressure PSD compared with the Webb 1998 and Brown 2014 limits">
</p>

*Figure 12. 7D.J38A BDH pressure PSD (2,816 windows): median (blue line) and
5th–95th percentiles (blue band). It is shown against the Webb (1998) Fig. 2
limits at ≤ 0.2 Hz (red bars) and the digitised Brown et al. (2014) range at
≥ 0.2 Hz (green). Black circles: median inside the limits. Orange triangles:
median above. Grey: outside the −3 dB response band (34 mHz–7.7 Hz), not
tested. Lower panel: instrument response relative to its reference. 25 of 40
bins (62.5 %) lie within. All 15 misses are above Brown on the Pacific
microseism slope (0.22–0.74 Hz), where Brown's in-water hydrophone references
are expected to be exceeded by a seafloor gauge under North Pacific swell.*

### QC test 2: nearby stations

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/j38a_test2_nearby.png?raw=true" width="85%"
       alt="J38A compared with the median of its six neighbours">
</p>

*Figure 13. Top: J38A daily pressure spectrum against the equal-weight median of
six neighbouring Cascadia gauges within 100 km and 10 % water depth. Bottom:
median daily difference (line), its 5th–95th spread (band) and the RMS of the
daily differences (dashed). Over 176 common days the median difference is
**+0.6 dB** and the RMS **1.8 dB**. Between 0.01 and 1 Hz they are +1.0 dB and
1.9 dB, well inside the 3 dB pass bound.*

### QC test 3: infragravity waves (Godin model)

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/j38a_test3_godin.png?raw=true" width="85%"
       alt="J38A infragravity band compared with the Godin 2013 model">
</p>

*Figure 14. J38A pressure PSD at 1–20 mHz (blue: median and 5th–95th
percentiles) against Godin et al. (2013) Eq. 6 at the station's depth of
2734 m. Orange: model unsmoothed. Dashed orange: model averaged over the PSD
bins. Dotted purple: Eq. 6 without the tanh(kH)^−3/2 factor. Blue points: the
11 bins compared (6–14 mHz). Grey: below the 4 mHz band edge. The data lie
**+14.6 dB** above the model, a partial score. The plateau of about 45 dB near
5 mHz matches Webb's (1998) Pacific levels. This suggests that the excess
comes from the energetic North Pacific long-wave climate rather than the
gauge.*

### Calibration 1: Rayleigh waves

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/j38a_cal1_rayleigh.png?raw=true" width="85%"
       alt="J38A Rayleigh-wave pressure/vertical calibration">
</p>

*Figure 15. Rayleigh-wave calibration of J38A LDH against LHZ (1 Hz) at
30–70 mHz, with pressure sign +1 and every eligible M ≥ 5.5 USGS event.
Top: pressure/vertical-acceleration coherence of the 16 accepted events (green)
and their median. Bottom: measured, theoretical and corrected P/a ratios. The
station estimate is the median of the per-event multipliers, **g = 0.7805 ±
0.016**. The earlier 7D run gave 0.7797 from 3 events. The title value (0.79)
is the fit to the stacked median ratio.*

### Calibration 2: microseisms

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/img/j38a_cal2_microseism.png?raw=true" width="85%"
       alt="J38A microseism admittance calibration">
</p>

*Figure 16. Microseism calibration of J38A over one week (15–22 February 2012;
118 of 168 one-hour windows kept). The panels show pressure/vertical coherence,
phase against the expected polarity, admittance, the implied multiplier and the
apparent velocity of the water-layer model. In the best coherent sub-band
(36–49 mHz, coh² 0.975) the mode-0 model brackets the multiplier at **g =
0.780–0.794**. This bracket contains the independent Rayleigh estimate of
0.7805.*

### Calibration 3: local P waves

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_yv_pilot/report/img/cal_localp_RR40.png?raw=true" width="70%"
       alt="RR40 local P-wave calibration by frequency band">
</p>

*Figure 17. Local P-wave calibration of YV.RR40 BDH against BHZ, using local
earthquakes recorded during the RHUM-RUM deployment. Top: pressure multiplier
per frequency bin for each accepted event (colours). Black bars: the band median
and its spread, labelled with the number of events. Middle: phase of observed
minus predicted. Bottom: squared coherence, above 0.8 throughout. At 2–4 Hz six
events give **g = 0.75**. This agrees with the pilot's Rayleigh estimate (0.743), its
microseism bracket (0.713–0.751) and the published value (0.74). Above 4 Hz
the estimate falls to about 0.44 at both YV DPGs, which still needs explaining,
so local P is used as a diagnostic.*

The full J38A report, with every figure, is in
[`reports/2026-10-26_j38a_report/`](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-26_j38a_report)
([PDF](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-26_j38a_report/BRUIT-FM_pressure_QC_report_2026-10-26.pdf)).
The RR40 example comes from the
[YV pilot report](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-26_yv_pilot).

## 8. Provisional noise model

The noise model combines the PPSDs of the **396 provisionally accepted
channels** (§9) into station-balanced percentile curves. These cover 311
stations and 1.44 million three-hour windows. It is **provisional**: the station
list comes from the provisional policy of §9, not from an agreed acceptance, and
the same command reruns once the stations are approved.

**Method** (`bruitfm-group-noise-model`):

1. Keep one channel per station and sensor type, the one with the widest usable
   band.
2. Use only PSD bins inside the usable response band (within 20 dB of the
   sensitivity).
3. Drop dead windows: those more than 30 dB below the channel's median window
   (0.4 % of windows).
4. Pool each station's windows into one histogram, then weight stations equally
   at each frequency.
5. Report the 5th, 10th, 50th, 90th and 95th percentiles, with the number of
   stations at each frequency.

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_all.png?raw=true" width="70%"
       alt="Provisional pressure noise model of the 311 accepted stations">
  <br><em>Figure 18a. Provisional pressure noise model of the 311 accepted
  stations: median (line), 10–90 % (dark band) and 5–95 % (light band), each
  station weighted equally. Red bars: Webb (1998) Fig. 2 ranges. Green: Brown et
  al. (2014) range at ≥ 0.2 Hz. The step near 15 Hz is where the 40–50 Hz
  channels stop and fewer, higher-rate stations remain.</em>
</p>

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_all_heatmap.png?raw=true" width="85%"
       alt="Heat map of the station-balanced pressure PSD density">
  <br><em>Figure 18b. The same model as a heat map: the station-balanced
  probability of each 1 dB level at each frequency (colour, % per dB), with the
  median (solid) and the 5th and 95th percentiles (dashed). The bright ridge is
  the typical deep-water spectrum. The faint upper branch in the notch
  (0.02–0.1 Hz, up to 50 dB) comes from shallow sites, where swell reaches the
  seafloor.</em>
</p>

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_basin_lines.png?raw=true" width="100%"
       alt="Median PSD of every accepted station, by ocean basin">
  <br><em>Figure 18c. Every accepted PPSD as a line: the median PSD of each
  station over its whole record (blue), inside its usable band, with the group
  median (black), by ocean basin. The lines show which stations make up each
  percentile band. For example, the Pacific curves that stay high through the
  notch come from the 30 accepted Pacific stations shallower than 1 km, 19 of
  them 7D gauges on the Cascadia shelf.</em>
</p>

Each grouping has the same three kinds of figure: percentile bands, a heat map
and all station lines.

| Grouping | Percentile bands | Heat map | All station lines |
|---|---|---|---|
| All stations | [bands](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_all.png) | [heat map](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_all_heatmap.png) | [lines](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_all_lines.png) |
| Ocean basin | [bands](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_basin.png) | [heat map](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_basin_heatmap.png) | [lines](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_basin_lines.png) |
| Water depth | [bands](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_depth_bin.png) | [heat map](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_depth_bin_heatmap.png) | [lines](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_depth_bin_lines.png) |
| Sensor type | [bands](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_type.png) | [heat map](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_type_heatmap.png) | [lines](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_type_lines.png) |
| Local season | [bands](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_season.png) | [heat map](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_season_heatmap.png) | [lines](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_season_lines.png) |

**Median level (10th–90th percentile) in dB re 1 Pa²/Hz** at five frequencies:

| Group | Stations | 5 mHz (infragravity) | 0.05 Hz (notch) | 0.2 Hz (microseism) | 1 Hz | 5 Hz |
|---|---:|---|---|---|---|---|
| **All** | 311 | 36 (26, 48) | −10 (−22, 20) | 22 (4, 36) | −6 (−16, 0) | −36 (−46, −24) |
| Pacific | 128 | 38 (32, 54) | −14 (−24, 48) | 28 (14, 40) | −4 (−12, 2) | −38 (−46, −20) |
| Atlantic | 89 | 30 (22, 40) | −8 (−14, 6) | 14 (4, 28) | −10 (−20, 0) | −34 (−44, −26) |
| Indian | 62 | 32 (26, 40) | −14 (−24, 2) | 22 (4, 32) | −8 (−16, −2) | −40 (−50, −28) |
| Mediterranean | 11 | – | 2 (−10, 6) | 6 (−4, 24) | −6 (−16, 2) | −28 (−36, −24) |
| < 1 km depth | 29 | 50 (40, 62) | 48 (0, 60) | 14 (2, 30) | −6 (−14, 0) | −34 (−42, −14) |
| > 5 km depth | 43 | 34 (26, 38) | −10 (−20, 0) | 34 (10, 42) | −8 (−18, 0) | −42 (−48, −26) |
| DPG | 129 | 38 (32, 50) | −14 (−24, 18) | 28 (10, 38) | −4 (−12, 2) | −38 (−46, −28) |
| Hydrophone | 174 | 30 (22, 38) | −8 (−20, 4) | 16 (2, 30) | −10 (−18, −2) | −36 (−46, −24) |
| APG | 8 | 52 (44, 62) | 50 (28, 60) | 10 (−4, 26) | −8 (−14, 0) | −14 (−16, −14) |
| Local winter | 261 | 36 (30, 52) | −10 (−22, 20) | 26 (12, 38) | −6 (−14, 2) | −36 (−46, −26) |
| Local summer | 270 | 34 (24, 44) | −10 (−22, 22) | 16 (0, 34) | −8 (−18, 0) | −36 (−46, −24) |

– means fewer than 3 stations cover that frequency. The full curves for every
group are in the CSV files. APGs sit about 60 dB above DPGs and hydrophones in
the notch, which matches their known high noise floor there. APGs should
therefore be modelled separately.

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_depth_bin_medians.png?raw=true" width="85%"
       alt="Median pressure noise by water depth">
  <br><em>Figure 19. Median pressure noise by water depth. The infragravity
  cut-off moves to lower frequency as the water deepens, because surface
  gravity waves shorter than about twice the depth do not reach the seafloor.
  Stations shallower than 1 km have no noise notch: swell and wind waves reach
  the gauge at 0.05–0.1 Hz. The deepest stations have the strongest microseism
  peak.</em>
</p>

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_basin.png?raw=true" width="100%"
       alt="Provisional noise model by ocean basin">
  <br><em>Figure 20. Provisional noise model by ocean basin. The Pacific has the
  strongest microseisms and the widest notch spread, because its 90th
  percentile includes shallow sites. The Mediterranean has the weakest
  microseisms. The Southern Ocean stations are hydrophones without
  low-frequency response, so only frequencies above 0.5 Hz are covered.</em>
</p>

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model/noise_model_season_medians.png?raw=true" width="85%"
       alt="Median pressure noise by local season">
  <br><em>Figure 21. Median pressure noise by local season. Winter is DJF in the
  north and JJA in the south. The microseism peak is about 10 dB stronger in
  winter than in summer (26 against 16 dB at 0.2 Hz), and the infragravity band
  about 2 dB stronger. Above 1 Hz the seasons agree within 2 dB.</em>
</p>

Figures by sensor type and every CSV are in
[`reports/2026-10-28_provisional_decisions/noise_model/`](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-28_provisional_decisions/noise_model).
A first version built from all 352 *good* stations, before the acceptance rules
of §9, is kept in
[`reports/2026-10-28_noise_model_provisional/`](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-28_noise_model_provisional).
Its medians agree with this one within 2 dB in every cell of the table, except
the shallow (< 1 km) notch, which is 4 dB lower there.

## 9. Provisional decisions: policy, metadata and acceptance

Three decisions need the project supervisor: freezing the thresholds, dealing
with suspect metadata and accepting stations. They are taken here
**provisionally**, so that the noise model can be built now. Each is written
down so it can be reviewed and changed, and none applies a gain.

### Step 1: policy

All the rules of §3–§6, and the acceptance and metadata rules below, are frozen
as
[`policy_v1_provisional.json`](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/policy_v1_provisional.json)
(version 1.0-provisional). They are unchanged from the values set before the
results were seen.

### Step 2: metadata flags

For every channel, the median PSD is compared with the median of the good
stations in the same ocean basin, at 0.1–2 Hz inside the channel's usable band.
A network is flagged when its median offset is at least 15 dB over at least 3
channels.

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/metadata_offsets.png?raw=true" width="100%"
       alt="Level offset of each channel from comparable good gauges, by network">
  <br><em>Figure 22. Level offset of each pressure channel from the good-station
  median of its basin, 0.1–2 Hz. Dots: channels coloured by final class. Bars:
  network median, purple and starred when flagged. Grey: within ±15 dB. Four
  networks sit far outside as a block; 1T's bad channels scatter between −10
  and −25 dB with a network median of −10 dB.</em>
</p>

| Network | Sensor | Channels | Network offset | Multiplier that would match comparable gauges | Independent check | Flag |
|---|---|---:|---:|---:|---|---|
| Z4 | DPG | 7 | −59 dB | ×900 | Godin residual −43 dB | **metadata** |
| 2A | DPG | 32 | −47 dB | ×233 | Godin residual −43 dB; StationXML sensor is a "dummy" LDEO DPG | **metadata** |
| YD | DPG | 6 | −47 dB | ×220 | Godin residual −42 dB; same "dummy" response | **metadata** |
| X9 | DPG | 6 | +17 dB | ×0.135 | **Rayleigh g = 0.14 and microseism g = 0.12 on the same gauges** | **metadata** |
| 4N | Hydrophone | 8 | – | – | Stage-1 gain −1.655e8; cannot be evaluated | **metadata** |
| 1T | Hydrophone | 92 | −10 dB | ×3.2 | Response shape (microseism low, 5 Hz floor high) | not flagged; 59 bad |
| X3 | Hydrophone | 7 | +14 dB | ×0.21 | all 7 good | below the flag; worth a look |

The multiplier is the pressure factor g that would bring the channel to the
level of comparable gauges. **No gain is applied.** The offset method agrees
with the calibration checks: the implied multiplier is within 6 dB of the
Rayleigh g on 41 of 54 channels and of the microseism g on 118 of 167. For X9
the two methods agree within 15 %, which is strong evidence of a gain error of
about 7 in its StationXML. Each flagged network still needs checking against
the operator's calibration sheets.

### Step 3: acceptance

| Decision | Channels | Rule |
|---|---:|---|
| **accepted** | **396** (311 stations) | final class *good* and ≥ 7 days of PSD windows, in a network without a metadata flag |
| review | 228 | final class *questionable* (185), or *good* with < 7 days of data (43) |
| rejected (metadata) | 59 | network flagged in step 2 (2A, Z4, X9, YD) or StationXML that cannot be evaluated (4N) |
| rejected | 163 | final class *bad* (141), internal pressure (17) or nonworking (5) |

The 7-day minimum keeps very short records from carrying the same weight in the
station-balanced noise model as year-long deployments. Each accepted channel
carries its accepted band (the usable band of the QC tests) and time range (the
whole epoch).

<p align="center">
  <img src="https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/acceptance_map.png?raw=true" width="90%"
       alt="Map of provisional decisions">
  <br><em>Figure 23. Provisional decision for every pressure channel. Accepted
  stations cover every ocean basin. The rejected clusters are mostly whole
  networks: the Cascadia margin (YR, and the 7D internal-pressure channels), 1T
  at Mayotte, and the metadata networks 2A, YD, Z4, X9 and 4N.</em>
</p>

The decision table, with one row per pressure epoch, its decision and reason,
the accepted band and time range, and the metadata offset, is
[`acceptance.csv`](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/acceptance.csv)
(also [.xlsx](https://github.com/MohammadAmin-Aminian/bruit-fm/blob/v1.3.0/docs/reports/2026-10-28_provisional_decisions/acceptance.xlsx)).

## 10. Limitations and next steps

**Limitations**

- **Thresholds, metadata flags and acceptance are provisional.** They were
  chosen before looking at the results (§9) and have not yet been reviewed.
- **The references are generic.** Webb (1998) and Brown et al. (2014) are
  global ranges, and the Godin model was fitted to New Zealand data.
  Basin-specific references would sharpen Tests 1 and 3.
- **Test 2 coverage is limited.** Only a third of the channels have a usable
  neighbour.
- **Calibration methods have limits.** The Rayleigh and microseism methods need a
  co-located vertical seismometer at a usable sample rate. Local P needs
  several well-recorded local events. Some stations have neither.
- **Depth comes from the StationXML elevation** and has not been checked against
  bathymetry.

**Next steps**

1. Review the provisional policy (§9, step 1) using the distributions in §5
   and §6, and adopt it or change it.
2. Check the StationXML of the flagged networks (2A, YD, Z4, X9, 4N), and of 1T
   and X3, against the operators' calibration sheets.
3. Review the 228 channels marked *review* and confirm the accepted list.
4. Rerun `run.sh` in `docs/reports/2026-10-28_provisional_decisions/` with the
   agreed policy. It rebuilds the flags, the acceptance list and the noise
   model, which is then labelled `approved`.
5. Extend the calibration coverage. Add the geophone channels (`EL?`/`DH?`) for
   the local P check, and decide on explicit resampling where pressure and
   vertical rates differ.

## 11. Reproducing these results

The commands in the [README](../BRUITFM_README.md#quick-start) reproduce every stage.
The dated reports in [`docs/reports/`](reports/) contain the exact scripts and
logs of each run:

| Report | Contents |
|---|---|
| [2026-10-23 fresh pressure PPSD](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-23_fresh_pressure_ppsd) | Pressure PPSDs of the whole drive, failures by cause |
| [2026-10-26 YV pilot](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-26_yv_pilot) | Six YV stations through every test and calibration; RR52 reproduces the published gain |
| [2026-10-26 J38A report](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-26_j38a_report) | 7D.J38A in full, and the pressure PPSDs of all networks by group |
| [2026-10-26 all-station QC](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-26_all_station_qc) | The three QC tests on all 846 channels; scored table |
| [2026-10-28 seismic PPSDs and calibrations](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-28_seismic_calibration) | Seismometer PPSDs; three calibration checks on every pressure channel; completed table |
| [2026-10-28 provisional decisions](https://github.com/MohammadAmin-Aminian/bruit-fm/tree/v1.3.0/docs/reports/2026-10-28_provisional_decisions) | Provisional policy, metadata flags, acceptance list and the noise model of the accepted stations |

## 12. References

- Brown, D., Ceranna, L., Prior, M., Mialle, P. & Le Bras, R. J. (2014). The IDC
  seismic, hydroacoustic and infrasound global low and high noise models.
  *Pure and Applied Geophysics*, 171, 361–375. doi:10.1007/s00024-012-0573-6
- Godin, O. A., Zabotin, N. A., Sheehan, A. F., Yang, Z. & Collins, J. A. (2013).
  Power spectra of infragravity waves in a deep ocean. *Geophysical Research
  Letters*, 40, 2159–2165. doi:10.1002/grl.50418
- Webb, S. C. (1998). Broadband seismology and noise under the ocean. *Reviews of
  Geophysics*, 36(1), 105–142.
- Peterson, J. (1993). Observations and modeling of seismic background noise.
  *USGS Open-File Report* 93-322.
- Further method references (Rayleigh, microseism and local-P calibration, tilt
  correction) are in [references.bib](references.bib).
