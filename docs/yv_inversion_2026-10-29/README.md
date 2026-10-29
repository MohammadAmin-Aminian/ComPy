# Vs inversion of the YV broadband compliance (2026-10-29)

Bayesian shear-velocity inversion of the eight YV (RHUM-RUM) broadband compliance functions computed with ComPy 2.2 on
2026-10-28 (`docs/yv_compliance_2026-10-28`). The method, the audit of the earlier inversion code and the validation are
in [`docs/INVERSION.md`](../INVERSION.md).

## Inputs and choices

- **Compliance**: gain from ComPy's own Rayleigh-wave calibration (matches Aminian et al. 2025, Table 1). RR36, RR38, RR40
  and RR50 use the stable epoch before their pressure/vertical sensitivity step.
- **Covariance**:
  - Bendat & Piersol random error;
  - covariance of the monthly stacks divided by the number of months, shrunk halfway to the diagonal;
  - a 5 % pressure-gain error shared by all frequencies;
  - a 1 % floor.
- **Model** (Zha & Webb 2016):
  - Hamilton sediment of thickness 0–1000 m;
  - crustal layers of 100, 200, 400, 800 and 1600 m, then one more down to the CRUST1.0 Moho at the station;
  - a mantle half-space;
  - Vp = 1.85 Vs and Christensen & Shaw density;
  - uniform Vs priors plus a second-difference roughness prior whose scale is sampled.
- **Sampler**: 8 DRAM chains × 80 000 iterations (30 000 burn-in), started from prior draws, plus a prior-only run. The
  fixed-Vp and free-Vp/Vs variants use 8-temperature parallel tempering.

## Main results

Vs in m/s as median (5–95 %). The last column gives the depths above which 50 % and 90 % of the lowest-frequency
sensitivity lie; deeper Vs reflects the prior and the layer geometry.

| Station | Record | Band [mHz] | χ²/N | R-hat | ESS | Sediment [m] | Vs 0.5 km | Vs 1 km | Vs 2 km | 50/90 % sensitivity [km] |
|---|---|---|---|---|---|---|---|---|---|---|
| RR28 | whole | 7.3-13.4 | 3.7/6 | 1.009 | 896 | 330 (240-398) | 2340 (1616-3860) | 3040 (2398-3929) | 3771 (3330-4156) | 0.4/6.6 |
| RR29 | whole | 7.3-12.2 | 3.6/5 | 1.007 | 885 | 87 (13-138) | 3267 (2502-4342) | 3591 (2940-4342) | 3942 (3437-4347) | 2.1/8.2 |
| RR34 | whole | 7.3-13.4 | 2.7/6 | 1.029 | 305 | 204 (37-303) | 2315 (1698-3688) | 2940 (2206-3628) | 3428 (2929-3809) | 0.9/6.9 |
| RR36 | stable epoch | 7.3-15.9 | 3.5/8 | 1.015 | 643 | 75 (9-148) | 2358 (1855-3160) | 2857 (2268-3271) | 3474 (3074-3755) | 1.3/6.9 |
| RR38 | stable epoch | 8.5-14.6 | 21.6/6 | 1.008 | 1489 | 251 (226-274) | 4386 (3726-4651) | 4545 (4245-4675) | 4601 (4439-4686) | 0.2/5.2 |
| RR40 | stable epoch | 8.5-13.4 | 2.5/5 | 1.010 | 612 | 155 (27-215) | 3043 (2220-4240) | 3411 (2788-4150) | 3788 (3434-4168) | 0.8/6.0 |
| RR50 | stable epoch | 7.3-13.4 | 2.3/6 | 1.011 | 1082 | 75 (12-111) | 3468 (2676-4408) | 3805 (3208-4461) | 4147 (3784-4501) | 2.0/8.0 |
| RR52 | whole | 7.3-17.1 | 3.2/9 | 1.009 | 798 | 42 (5-78) | 2918 (2529-3645) | 3349 (2993-3728) | 3826 (3607-4030) | 1.9/7.1 |

## Findings

- **Convergence and fit.** All main runs converged (R-hat ≤ 1.03). Seven fit to χ² below N, so the covariance is slightly
  conservative.
- **RR52** agrees with the GJI 2025 description: 2.1–2.5 km/s in the top 300 m, 3.8 km/s at 1.5–3.1 km and 4.3 km/s
  below.
- **RR38** does not fit with Vp/Vs 1.85 (χ² 21.6 for 6) or with Brocher Vp.
  - With a free crustal Vp/Vs it moves to the upper bound (2.33; 90 %: 2.15–2.39) and still fits worse than the other
    stations (χ² 9.9).
  - With Vp fixed at CRUST1.0 it fits (χ² 2.3), but this implies Vp/Vs ≈ 5 in the upper crust.
  - Either the crust has a very high Poisson ratio (serpentinization) or the low-frequency compliance is biased low. RR38
    also has the largest monthly scatter (5–13 %). Treat it with caution.
- **Vp/Vs.** Compliance barely constrains it: free in [1.7, 2.4], the posteriors span most of that range at seven
  stations (medians 2.0–2.15).
- **Fixed Vp** at CRUST1.0 (the assumption stated in Aminian et al. 2025) makes the posteriors multimodal. Plain chains
  reach R-hat 2.3; tempering is required. It fits RR50 (χ² 19 for 6, not converged) and RR52 (χ² 8.5 for 9) poorly.
- **Gain, record, Moho.** The gain choice (`gains.csv` vs ComPy), whole record vs stable epoch, and a ±1 km Moho change
  Vs by less than the posterior width.
- **Sediment.** The sediment model matters most in the top kilometre (variant `sedscale`). The slow layer is 330 m at
  RR28 (CRUST1.0 260 m), 87 m at RR29 (GlobSed 110 m), 204 m at RR34 (CRUST1.0 410 m), and 40–250 m at the ridge
  stations. At the ridge stations it is more likely porous extrusives than sediment (Zha & Webb 2016).

## Files

- `scripts/`:
  - `a_invert.py`: all variants and the prior-only run;
  - `b_occam.py`: Occam inversion from CRUST1.0;
  - `c_figures.py`: figures and summary tables;
  - `d_tables.py`: variant comparison;
  - `e_report.py`: HTML report;
  - `run_all.sh`, `run_fixedvp.sh`, `run_vpvs.sh`: runners.
- `results/<STA>/<STA>_<variant>.json`: parameter quantiles, R-hat/ESS, χ², sampler settings and the CRUST1.0
  column. The draws (`.npz`, about 3 MB each) are kept in `~/Desktop/BRUIT-FM/inversion_yv_2026-10-29/results`.
- `results/tables/`: `inversion_summary.{csv,json}`, `main.md`, `variants.{md,json}`.
- `figures/`:
  - `<STA>_inversion.png`: data fit, Vs posterior, sensitivity, sediment;
  - `vs_profiles.png`;
  - `compliance_fit.png`;
  - `sensitivity.png`;
  - `sediment.png`.

Reproduce (about 30 minutes on 8 cores with Numba): `scripts/run_all.sh`, `run_fixedvp.sh`, `run_vpvs.sh`, then `b_occam.py
<STA>`, `c_figures.py`, `d_tables.py`, `e_report.py`. CRUST1.0 must be unpacked in
`compliance_yv_2026-10-28/literature/crust1`.
