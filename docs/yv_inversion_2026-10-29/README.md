# Vs inversion of the YV broadband compliance (2026-10-29)

Bayesian shear-velocity inversion of the eight YV (RHUM-RUM) broadband compliance functions computed with ComPy 2.2 on
2026-10-28 (`docs/yv_compliance_2026-10-28`). The method, the audit of the earlier inversion code and the validation are
in [`docs/INVERSION.md`](../INVERSION.md).

## Second round: gain sampled, epochs joined, paper-style figures

- **Gain.** The pressure gain is sampled as a free factor within ±5 % of the Rayleigh value, one per epoch
  (`gain_bounds`), instead of a 5 % covariance term.
- **After-jump data.** For RR36, RR38, RR40 and RR50, the data after the sensitivity jump (`g_after_epoch.py`) are
  stacked with their own Rayleigh gain and inverted together with the stable epoch.
- **Final model.** Crustal Vp/Vs is free in [1.7, 2.4]; Z/P estimator; tempering with 8 replicas.
- **Bracket.** `final_h2` repeats the final run with 1/(P/Z) (`h_estimators.py`).
- **Figures.** `i_paper_figures.py` writes the figures in the style of Aminian et al. (2025):
  - `paper_fig9_heatmaps.png`;
  - `paper_fig_misfit.png`;
  - `paper_fig8_compliance.png`;
  - `paper_fig10_profiles.png`;
  - `<STA>_paper.png`.

| Station | χ²/N (H1) | χ²/N (H2) | R-hat | Vp/Vs | Gain used → sampled | Sediment [m] | Vs 0.5 km | Vs 1 km | Vs 2 km | Vs 3 km |
|---|---|---|---|---|---|---|---|---|---|---|
| RR28 | 2.9/6 | 17.8/6 | 1.012 | 2.15 (1.81–2.38) | data: 0.780 → 0.777 | 311 (138–398) | 2157 (1450–3760) | 2789 (2150–3725) | 3421 (2910–3931) | 3421 (2910–3931) |
| RR29 | 2.9/5 | 12.4/5 | 1.005 | 2.15 (1.79–2.38) | data: 0.693 → 0.689 | 91 (15–139) | 3186 (2333–4252) | 3426 (2709–4217) | 3678 (3106–4218) | 3678 (3106–4218) |
| RR34 | 2.4/6 | 33.9/6 | 1.008 | 2.08 (1.75–2.37) | data: 0.573 → 0.574 | 219 (50–309) | 2258 (1723–3727) | 2811 (2203–3560) | 3185 (2757–3686) | 3185 (2757–3686) |
| RR36 | 6.8/16 | 39.3/16 | 1.017 | 2.12 (1.80–2.37) | stable: 0.642 → 0.656; after: 1.097 → 1.071 | 72 (10–143) | 2236 (1756–3320) | 2578 (2082–3094) | 3090 (2740–3594) | 3090 (2740–3594) |
| RR38 | 33.0/12 | 19.6/12 | 1.015 | 2.20 (1.90–2.38) | stable: 0.596 → 0.570; after: 0.638 → 0.667 | 224 (138–251) | 3062 (2307–4369) | 3722 (3025–4402) | 4058 (3697–4471) | 4058 (3697–4471) |
| RR40 | 4.2/10 | 31.5/10 | 1.009 | 2.01 (1.73–2.36) | stable: 0.739 → 0.717; after: 0.786 → 0.810 | 164 (40–219) | 2964 (2092–4165) | 3269 (2581–3978) | 3556 (3078–3994) | 3556 (3078–3994) |
| RR50 | 5.7/12 | 85.1/12 | 1.005 | 2.11 (1.82–2.37) | stable: 0.681 → 0.685; after: 1.161 → 1.148 | 80 (14–112) | 3389 (2502–4346) | 3672 (2980–4363) | 3964 (3495–4429) | 3964 (3495–4429) |
| RR52 | 3.3/9 | 11.4/9 | 1.007 | 2.04 (1.74–2.36) | data: 0.666 → 0.664 | 46 (6–81) | 2823 (2332–3632) | 3193 (2779–3672) | 3590 (3283–4013) | 3590 (3283–4013) |

Vs in m/s, median (5–95 %); gains are the absolute values (Rayleigh gain → posterior median).

### Findings

- **Gain resolution.** One epoch cannot pin the absolute gain, which trades off with the Vs level. Two epochs resolve
  their relative gain:

  | Station | After/stable gain: inversion | Rayleigh |
  |---|---|---|
  | RR50 | 1.677 | 1.705 |
  | RR36 | 1.635 | 1.708 (1 event) |
  | RR40 | 1.127 | 1.063 |

  The RR50 and RR36 jumps are gain changes that the Rayleigh calibration corrects.
- **RR38.** Its after-jump compliance differs in shape (+27 % at the lowest, +7 % at the highest frequency), and it
  misfits (χ² 33/12). 1/(P/Z) fits it better than Z/P, so its pressure channel is suspect. Leave RR38 out of the
  interpretation.
- **Estimator.** 1/(P/Z) misfits the other seven stations (χ² 11–85), so their noise is on the vertical and Z/P is
  right.
- **Vp/Vs and Vs level.** With Vp/Vs free, Vs at 1–2 km is 3–11 % lower than with Vp/Vs 1.85, and its intervals
  include that uncertainty.
- **Errors in the publication figure code** (`docs/INVERSION.md`, F1–F2):
  - the heat-map depths are 5/3 too deep relative to their labels;
  - the serpentinization formula disagrees with the Fig. 10 caption.

## First round

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
  - `a_invert.py`: all variants and the prior-only runs (`main`, `final`);
  - `g_after_epoch.py`: compliance after the sensitivity jump;
  - `h_estimators.py`: Z/P and 1/(P/Z) estimators;
  - `i_paper_figures.py`: figures in the style of the paper;
  - `j_final_tables.py`: tables of the final runs;
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
