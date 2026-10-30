# Local-P campaign v1: all available stations

**Historical v1 results, not the current protocol.** The independent review
found the documented noise-window guard was not enforced. Keep these results
unchanged for reproducibility, not as validated calibration. See the
[v2 phase report](local_p_phase_report.md) and
[current staging table](local_p_table_staging.csv). Its three-APG rerun replaces
only those rows; the other 63 executed pairs remain labelled legacy protocol.
Polarity dominance in v1 is not an overall-confidence classification. Statements
below about hardware failure or response errors are hypotheses, not proven
physical diagnoses. The old table-regeneration command must use fresh
`results/LOCAL_P/` destinations; it must not overwrite these historical docs.

Full run of the local-P check (calibration check 3) on every high-rate
pressure/vertical pair in the local archive. Method:
[local_p_method.md](local_p_method.md); commands and pilot:
[local_p_check.md](local_p_check.md). Diagnostic outcomes only: no gain is
applied, no master is written, no score is assigned.

## Provenance

| Item | Value |
|---|---|
| Campaign | `results/LOCAL_P/campaigns/all_v1/` (run by the user), code `efa4c3e`, finished 2026-10-15 15:36 UTC |
| APG timing rerun | `results/LOCAL_P/campaigns/apg_timing_v2/`, code `ff45286`, XO LT02/LA33/LA34 HDH only; supersedes those three rows (‡) |
| Export | `results/LOCAL_P/campaigns/all_v1_export/` (120 master-epoch rows) |
| Catalogues | USGS ComCat, `results/LOCAL_P/catalogs/` (hashes in [local_p_check.md](local_p_check.md)) |
| Table | [local_p_campaign_v1.csv](local_p_campaign_v1.csv), regenerate with `scripts/summarize_local_p_campaign.py` (§5) |
| Inputs | All inventory/catalogue/waveform hashes unchanged in both runs |

## 1. Outcomes

| Status | Pairs | Meaning |
|---|---|---|
| `conditional_gain_candidate` | 1 | XO.LT02 APG, steady state, ≥ 3 events |
| `diagnostic_only` | 17 | pre-multiple (deep water) results |
| `unresolved` | 22 | events processed, none accepted |
| `no_suitable_events` | 26 | YV: no ComCat M ≥ 2.5 event within 300 km in the record |
| `not_applicable_sample_rate` | 39 | YL (all 25), XO WD48 and the 1 Hz LDH pairs |
| `no_full_pressure_response` | 15 | 7D BXH |

66 pairs ran, none failed. Rejection reasons over all events processed in
`all_v1`: lag refinement (the APG search later replaced) 219, vertical onset not measured
163, pressure onset not measured 156, onset mismatch 65, no valid band 33,
missing waveform day 21 (7D archive covers only 34–74 days in late 2011).

## 2. Findings

1. **Polarity is the most consistent result.** All 8 YV stations with results
   are −1 relative to the nominal response, agreeing with the Rayleigh assumed
   YV sign (−1) and the YV microseism near-0° phase family, at 1–16 Hz with a
   different wavefield. XO pressure gauges (DPG, EDH, LT02 APG) and 7D QEP
   stations J28A/J30A/J38A are +1, agreeing with their assumed sign.
2. **7D lead:** J35A and J36A (SIO ABALONES, 50 Hz) are −1 (1–2 events; J36A
   7/7 bands) while QEP stations are +1. Their BDH metadata sign should be
   checked before signed 7D analysis. Few events: a lead, not a conclusion.
3. **XO SIO DPGs** (WD70 41 events, WD66 10): ≈0.77–0.78 at 1–4 Hz, rising to
   ≈0.86–0.95 above 4 Hz. Rayleigh/microseism on the 1 Hz LDH channel give
   0.91–0.93; Deng et al. found WD66 coefficients falling to 0.79 by 0.07 Hz.
   Different band and channel: context, not validation.
4. **YV DPGs** (RR36, RR40, RR50): 1–2 Hz values (0.64, 0.67, 0.51) sit near
   or below the low-frequency gains (0.65, 0.74, 0.72), then fall to ≈0.4 above
   4 Hz with many indeterminate phases. This suggests the nominal YV pressure
   response above ~2 Hz does not match the instruments. RR53 (0.06–0.10 from
   one event) is an outlier.
5. **XO HTI-90U hydrophones (LA33/LA34 EDH):** 0.29 → 1.66 from 1–2 to
   8–16 Hz with many indeterminate phases: the nominal low-frequency response
   (corner) is very likely wrong. **LT02 EDH is dead** in 2018.
6. **XO.LT02 APG**, the only steady-state site: 0.972 (0.5–1 Hz, n = 5) and
   1.021 (1–2 Hz, n = 9), consistent (scalar 0.999), model interval ±10–17%.
   Expected ≈1 for an absolute gauge; still conditional (few events, shared
   vertical response).
7. **APG timing** ([local_p_check.md](local_p_check.md) §2.1): LT02 drifts
   −0.483 s/day (reproduces Deng's offsets within 0.3 s); LA33/LA34 are
   synchronized but unresolved because their 2–10 Hz onset SNR is low. The
   single LA34 polarity in `all_v1` was a lag-search artefact (withdrawn).
8. **Coverage limits:** only one shallow station exists, so only one
   conditional gain candidate is possible; 7D has 33 catalogue events and
   34–74 archive days; most YV stations have no local events.

## 3. Station table

Multipliers g (pressure multiplier, `g·p_nominal = p_model`), median over
events, with event count. ‡ = row from `apg_timing_v2`. Rayleigh g:
0.03–0.07 Hz baseline; † = different (1 Hz LDH) channel. Microseism g:
frozen-rule conditional value (*) or legacy weekly median. Only pairs that
ran are listed; all 120 rows are in the CSV.

| Pressure channel | Status | Depth (m) | Accepted/processed | g 0.5-1 Hz (n) | g 1-2 Hz (n) | g 2-4 Hz (n) | g 4-8 Hz (n) | g 8-16 Hz (n) | Polarity | Rayleigh g | Microseism g | Main rejection |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7D.J28A..BDH | diagnostic_only | 2867 | 2/6 |  | 0.81 (1) | 0.77 (1) | 1.06 (2) |  | +1 |  |  | missing waveform file |
| 7D.J29A..BDH | unresolved | 2849 | 0/5 |  |  |  |  |  |  |  | 1.002 | pressure vertical onset mismatch |
| 7D.J30A..BDH | diagnostic_only | 2824 | 1/3 |  | 1.15 (1) | 1.02 (1) | 1.10 (1) | 1.08 (1) | +1 | 0.981† | 0.943 | vertical onset not measured |
| 7D.J35A..BDH | diagnostic_only | 2662 | 1/8 |  |  |  |  | 0.87 (1) | −1 |  | 0.889 | vertical onset not measured |
| 7D.J36A..BDH | diagnostic_only | 2821 | 2/6 |  | 0.73 (1) | 0.96 (2) | 0.96 (2) | 1.00 (2) | −1 |  |  | pressure onset not measured |
| 7D.J38A..BDH | diagnostic_only | 2734 | 1/3 |  |  | 0.59 (1) | 0.94 (1) | 0.93 (1) | +1 | 0.780† | 0.729 | vertical onset not measured |
| 7D.J43A..BDH | unresolved | 2654 | 0/6 |  |  |  |  |  |  |  | 1.074 | vertical onset not measured |
| 7D.J44A..BDH | unresolved | 2724 | 0/5 |  |  |  |  |  |  |  |  | pressure vertical onset mismatch |
| 7D.J45A..BDH | unresolved | 2757 | 0/3 |  |  |  |  |  |  | 0.829† | 0.794* | pressure vertical onset mismatch |
| 7D.J52A..BDH | unresolved | 2613 | 0/2 |  |  |  |  |  |  |  |  | missing waveform file |
| 7D.J53A..BDH | unresolved | 2686 | 0/2 |  |  |  |  |  |  |  |  | missing waveform file |
| 7D.J54A..BDH | unresolved | 2662 | 0/1 |  |  |  |  |  |  | 1.080† | 1.020* | missing waveform file |
| 7D.J55A..BDH | unresolved | 2800 | 0/1 |  |  |  |  |  |  |  | 0.795 | missing waveform file |
| 7D.J63A..BDH | unresolved | 2882 | 0/1 |  |  |  |  |  |  |  |  | missing waveform file |
| 7D.J67A..BDH | unresolved | 2612 | 0/2 |  |  |  |  |  |  |  | 0.973* | missing waveform file |
| XO.LA33..EDH | diagnostic_only | 4478 | 32/80 |  | 0.29 (9) | 0.66 (27) | 1.16 (32) | 1.64 (29) | +1 | 0.934† | 0.886 | pressure onset not measured |
| XO.LA33..HDH ‡ | unresolved | 4478 | 0/80 |  |  |  |  |  |  | 0.934 | 0.886 | pressure onset not measured |
| XO.LA34..EDH | diagnostic_only | 4612 | 27/80 |  | 0.30 (12) | 0.58 (21) | 1.23 (26) | 1.68 (22) | +1 | 0.963† | 0.926 | vertical onset not measured |
| XO.LA34..HDH ‡ | unresolved | 4612 | 0/80 |  |  |  |  |  |  | 0.963 | 0.926 | pressure onset not measured |
| XO.LT02..EDH | unresolved | 152 | 0/80 |  |  |  |  |  |  |  |  | pressure onset not measured |
| XO.LT02..HDH ‡ | conditional_gain_candidate | 152 | 9/80 | 0.97 (5) | 1.02 (9) |  |  |  | +1 |  |  | lag refinement failed |
| XO.WD66..HDH | diagnostic_only | 4941 | 10/80 |  | 0.78 (5) | 0.78 (8) | 0.86 (8) | 0.95 (10) | +1 | 0.909† | 0.929 | pressure onset not measured |
| XO.WD70..HDH | diagnostic_only | 3253 | 41/80 |  | 0.77 (15) | 0.77 (32) | 0.88 (22) | 0.89 (27) | +1 | 0.930† | 0.930 | pressure vertical onset mismatch |
| YV.RR01.00.BDH | no_suitable_events | 4298 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR03.00.BDH | no_suitable_events | 4340 | 0/0 |  |  |  |  |  |  | 1.124 |  |  |
| YV.RR05.00.BDH | no_suitable_events | 4092 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR06.00.BDH | no_suitable_events | 4216 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR07.00.BDH | no_suitable_events | 4370 | 0/0 |  |  |  |  |  |  | 1.079 |  |  |
| YV.RR08.00.BDH | no_suitable_events | 4190 | 0/0 |  |  |  |  |  |  | 1.160 | 1.210 |  |
| YV.RR09.00.BDH | unresolved | 2976 | 0/11 |  |  |  |  |  |  | 1.194 | 1.209 | vertical onset not measured |
| YV.RR10.00.BDH | diagnostic_only | 2310 | 1/13 |  |  | 1.09 (1) | 1.18 (1) | 1.67 (1) | −1 |  | 1.160 | vertical onset not measured |
| YV.RR11.00.BDH | diagnostic_only | 3473 | 2/12 |  | 1.15 (2) | 1.37 (1) |  |  | −1 | 1.083 | 1.099 | vertical onset not measured |
| YV.RR12.00.BDH | unresolved | 3185 | 0/2 |  |  |  |  |  |  |  |  | gap or missing data |
| YV.RR13.00.BDH | no_suitable_events | 4130 | 0/0 |  |  |  |  |  |  | 1.224 |  |  |
| YV.RR14.00.BDH | unresolved | 3420 | 0/1 |  |  |  |  |  |  | 1.096 | 1.101 | vertical onset not measured |
| YV.RR16.00.BDH | no_suitable_events | 4426 | 0/0 |  |  |  |  |  |  |  | 1.197 |  |
| YV.RR17.00.BDH | no_suitable_events | 2205 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR18.00.BDH | no_suitable_events | 4743 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR19.00.BDH | no_suitable_events | 4901 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR20.00.BDH | no_suitable_events | 4820 | 0/0 |  |  |  |  |  |  |  | 1.392 |  |
| YV.RR21.00.BDH | no_suitable_events | 4782 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR22.00.BDH | no_suitable_events | 4920 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR25.00.BDH | no_suitable_events | 4759 | 0/0 |  |  |  |  |  |  | 1.163 |  |  |
| YV.RR26.00.BDH | no_suitable_events | 4259 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR28.00.BDH | no_suitable_events | 4550 | 0/0 |  |  |  |  |  |  | 0.820 | 0.792 |  |
| YV.RR29.00.BDH | no_suitable_events | 4829 | 0/0 |  |  |  |  |  |  | 0.730 | 0.732 |  |
| YV.RR30.00.BDH | no_suitable_events | 5140 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR31.00.BDH | no_suitable_events | 2709 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR34.00.BDH | no_suitable_events | 4265 | 0/0 |  |  |  |  |  |  | 0.585 | 0.571 |  |
| YV.RR36.00.BDH | diagnostic_only | 3560 | 10/15 |  | 0.64 (5) | 0.51 (8) | 0.43 (9) | 0.39 (10) | −1 | 0.654 | 0.606 | pressure vertical onset mismatch |
| YV.RR38.00.BDH | diagnostic_only | 4560 | 1/5 |  |  | 0.78 (1) | 0.70 (1) | 0.40 (1) | −1 | 0.636 | 0.596 | vertical onset not measured |
| YV.RR40.00.BDH | diagnostic_only | 4780 | 6/9 |  | 0.67 (3) | 0.76 (6) | 0.44 (6) | 0.43 (6) | −1 | 0.743 | 0.749* | s too close |
| YV.RR41.00.BDH | unresolved | 5430 | 0/3 |  |  |  |  |  |  |  |  | vertical onset not measured |
| YV.RR42.00.BDH | no_suitable_events | 4771 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR43.00.BDH | unresolved | 4264 | 0/1 |  |  |  |  |  |  |  |  | s too close |
| YV.RR44.00.BDH | unresolved | 4548 | 0/1 |  |  |  |  |  |  |  |  | no valid band |
| YV.RR45.00.BDH | no_suitable_events | 2822 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR46.00.BDH | unresolved | 3640 | 0/1 |  |  |  |  |  |  |  |  | vertical onset not measured |
| YV.RR47.00.BDH | unresolved | 4582 | 0/1 |  |  |  |  |  |  |  |  | s too close |
| YV.RR48.00.BDH | unresolved | 4830 | 0/1 |  |  |  |  |  |  |  |  | pressure vertical onset mismatch |
| YV.RR50.00.BDH | diagnostic_only | 4118 | 3/26 |  | 0.51 (1) | 0.60 (3) | 0.45 (2) | 0.42 (2) | −1 | 0.717 | 0.817 | vertical onset not measured |
| YV.RR52.00.BDH | diagnostic_only | 2918 | 1/22 |  |  |  | 0.56 (1) | 0.44 (1) | −1 | 0.679 | 0.685* | vertical onset not measured |
| YV.RR53.00.BDH | diagnostic_only | 2940 | 1/13 |  |  | 0.10 (1) | 0.06 (1) |  | −1 |  |  | vertical onset not measured |
| YV.RR55.00.BDH | no_suitable_events | 4462 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR56.00.BDH | no_suitable_events | 4230 | 0/0 |  |  |  |  |  |  |  |  |  |
| YV.RR57.00.BDH | no_suitable_events | 5200 | 0/0 |  |  |  |  |  |  |  |  |  |

## 4. Caveats

- Deep-water rows use the time-gated pre-multiple relation and are diagnostic
  only. Many entries rest on 1–3 events.
- All methods share the nominal vertical response; agreement is not proof of
  absolute calibration.
- Local P (0.5–16 Hz) and the Rayleigh/microseism bands (≤ 0.1 Hz) differ by
  10–300×, and some rows compare different channels.

## 5. Reproduce the table

```bash
conda run --no-capture-output -n bruitfm python scripts/summarize_local_p_campaign.py \
  --campaign-dir results/LOCAL_P/campaigns/all_v1 \
  --supersede results/LOCAL_P/campaigns/apg_timing_v2 \
  --csv results/LOCAL_P/campaigns/all_v1_summary.csv \
  --markdown results/LOCAL_P/campaigns/all_v1_table.md
```

## 6. Open decisions

Pre-multiple relation as a calibration reference; treatment of frequency-
dependent multipliers; scoring policy; and whether deep APGs with verified
timing may skip the 2–10 Hz pressure-onset rule (see the handoff).
