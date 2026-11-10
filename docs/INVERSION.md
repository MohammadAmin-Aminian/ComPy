# Shear-velocity inversion in ComPy 3

This note describes the default inversion, `compy.inversion.auto`: what seafloor compliance measures, the errors, the model and its priors, the sampler, the outputs and their limits. The tests behind each choice are in [inversion_overhaul_2026-11/REPORT.md](inversion_overhaul_2026-11/REPORT.md). The inversions of ComPy 2.3 to 2.8 are described in [INVERSION_2x.md](INVERSION_2x.md); they still run, as checks.

## 1. Run it

One command on a stage C result:

```bash
compy-invert RR52_compliance_compy.npz --lat -20.47 --lon 68.11 --crust1 ~/crust1 --fetch-crust1
```

The same in Python:

```python
import numpy as np
from compy.inversion import auto, dram

d = np.load("RR52_compliance_compy.npz")
keep = d["cov_usable"]  # bins with a healthy estimate
data = auto.compliance_data(
    d["cov_freqs"][keep],
    np.interp(d["cov_freqs"][keep], d["freqs"], d["compliance"]),
    float(d["depth"]),
    cov_rel=d["cov_rel"][np.ix_(keep, keep)],
)
column = dram.crust1_column(-20.47, 68.11, "crust1")
model = auto.model_from_crust1(column)
result = auto.invert(data, model, gain_sigma=auto.GAIN_SIGMA["compy"])
result.save("RR52_auto")  # .npz, .json, _profile.csv
print(result.quality())
```

In the pipeline, stage C saves the covariance (`cov_rel`, `cov_freqs`, `cov_usable`) and stage E (`pipeline/e_invert.py`) runs this inversion by default. A compliance curve without its windows (a CSV with `freq_hz`, `compliance`, `uncertainty`) also works; its covariance then takes the measured structure from the uncertainties (`errors.default_covariance`).

Nothing is tuned per station. A run takes 1.5 to 10 minutes on six cores.

## 2. What compliance measures

Normalized compliance is `eta = k u_z / tau_zz` at the seafloor under an infragravity wave of wavenumber `k`, with `omega^2 = g k tanh(kH)` (Crawford 2004, eqs 1 and 4). Three properties decide what an inversion can and cannot say.

1. **It measures stiffness, not velocity.** For a half-space under slow loading, `eta = (lambda + 2 mu) / (2 mu (lambda + mu)) = (1 - nu) / mu` (Crawford 2004, eq. 5). Vs follows from the stiffness only with a density and a Vp/Vs. Between Vp/Vs 1.73 and 2.4 the same compliance gives Vs 11 % apart; a 10 % density error is 5 % in Vs. In relative terms compliance is as sensitive to Vp as to Vs: summed over a normal oceanic crust, `d ln eta / d ln Vs` is about -0.8 to -1.0 and `d ln eta / d ln Vp` about -1.0 to -1.2. The usual statement that it is 2 to 9 times more sensitive to Vs (Crawford 1994, eq. 2.21) is per m/s.
2. **A thin soft layer acts through its P modulus.** A layer much thinner than `1/k`, bonded to stiffer rock, is squeezed uniaxially and adds `k h / (rho Vp^2)` to compliance whatever its Vs. Under 3000 m of water, 300 m of sediment adds 60 to 80 % to the compliance of bare rock at 4 to 6 mHz, and halving the sediment Vs changes that by 1 to 3 %. So sediment thickness and sediment P modulus trade directly with deep crustal Vs, and the sediment Vs matters only above about `kh = 0.15`.
3. **Depth follows wavelength.** Sensitivity lies between `lambda / 2 pi` and `lambda / 4` (Crawford 2004). A deep station with a band of 7 to 17 mHz constrains broad averages between about 0.7 and 6 km below the seafloor, not the top few hundred metres and not individual thin layers.

The forward model is the minor-vector propagator of Crawford et al. (1991) with the full equations of motion. It agrees with a global-matrix solution in 50-digit arithmetic to 4e-7 over 60 random models (`scripts/e01_forward_checks.py`). A velocity gradient needs about 16 sublayers for 0.3 %.

The measured ratio is corrected for the attraction of the wave and the free-air gradient: `eta = k (omega^2 |Z/P| + a_w/p) / (omega^2 + 2g/R)` with `a_w/p = (2 pi G / g) cosh(kH) exp(-kH)` (Crawford et al. 1998, as quoted by Zha et al. 2014, eq. 2; re-derived in the report). The correction raises compliance, by 18 % at 7 mHz under 2900 m of water.

## 3. Data and errors

Measured on 172 stations of the 2026 archive (`scripts/e02_error_model.py`):

| Quantity | Median |
|---|---|
| Bendat & Piersol error of the stack | 0.15 % |
| Scatter when windows are resampled (bootstrap) | 0.23 %, with correlation 0.48 between neighbouring bins |
| The same in 30-day blocks | 0.38 %; 59 % of it a common level, 11 % a tilt |
| First half of the record against the second | level 0.7 %, tilt 0.8 % |
| Two independent pressure calibrations of the same gauge | 2-3 % apart; 16-84 % range -5.4 to +2.6 % |

The bin-to-bin shape of a stacked curve is known to a few tenths of a per cent. Its level and its slope are not. The covariance therefore has four parts (`auto.compliance_data`):

- **Stack**: the bootstrap covariance of the stack from its own windows, plus a slow level and a slow tilt from 30-day blocks (`errors.stack_covariance`).
- **Tilt**: an error linear in log frequency across the band, for the low-frequency response of the pressure gauge: 4 % edge to edge when the response is unchecked, 2 % when it was corrected with the tide (`compy.response_check`), more where the tide and the phase estimates of that correction disagree.
- **Smooth**: 0.5 % with a correlation length of 0.5 in ln f, for what a 1-D elastic model and the estimator leave out.
- **Gain**: not in the covariance. The pressure gain is a sampled factor with a Gaussian prior, 5 % for a calibrated gauge and 20 % for a nominal one.

There is no independent "floor" per bin. The 1 % floor of ComPy 2.x was five times the real scatter, which is why its fits had chi-square/N near 0.4, and it was too strict on slope errors.

Band. Stage C chooses the band from coherence. Two further limits are automatic: bins above 0.8 of the cut-off frequency `f_c = sqrt(g / 2 pi H)` are dropped, and bins at the band edges whose bootstrap scatter exceeds 2.5 times that of the quiet half of the band are dropped (`errors.usable_bins`). The second catches the collapse of the estimate where the pressure signal no longer reaches the seafloor, which coherence alone does not: at shelf station 7D.J09B it removes everything above 34 mHz.

Estimator. H1 = |S_pz|/S_pp is right when the incoherent noise is on the vertical. Pressure noise would make H1 rise with the pressure level; it falls by 1.5 % from the quietest to the loudest fifth of windows instead (median), so H1 is kept and H2, which is 9 % higher, is not used as a bracket.

## 4. Model and priors

`auto.SeafloorModel`, built for a station by `auto.model_from_crust1`:

| Part | Parameters | Prior | Source |
|---|---|---|---|
| Soft sediment | thickness | uniform 0-2000 m, or Gaussian from a sediment map if given | Zha & Webb (2016) use 0-2000 m |
| | Vs scale on the Hamilton profile | lognormal, sigma 0.3, bounds 0.4-2.5 | Hamilton (1976, 1979) as coded in ComPy 2.3 |
| | P-modulus factor on Hamilton Vp | lognormal, sigma 0.15 | see point 2 of section 2 |
| | density | Nafe-Drake curve of Vp (1.64 g/cc at the seafloor), not a constant 2.0 | Brocher (2005, eq. 1), via a published implementation |
| Crust | Vs of each layer; layers 100 m doubling in thickness | ln Vs Gaussian about the reference, sigma 0.3, correlation 0.7 between neighbours; bounds 300-4800 m/s | layers: Zha & Webb (2016) |
| | Vp | `oceanic_vp(Vs)` times a sampled factor, sigma 7 %, one for the top 1500 m and one below | mafic line of Brocher (2005, eq. 8), Vp/Vs 1.82-1.85 |
| | density | Carlson & Raskin (1984) of Vp, times a sampled factor, sigma 3 % | as in Crawford (1994, eq. 6.1) |
| Moho | depth below basement | Gaussian about the reference, sigma 30 % (at least 1.5 km) | CRUST1.0 |
| Mantle | Vs | uniform 3900-4900 m/s, above the deepest crustal Vs | Moho rule requested by the author; Shen et al. (2013) |
| Pressure gain | factor | Gaussian about 1, sigma 5 % or 20 % | section 3 |

The reference is the CRUST1.0 column of the station's 1-degree cell, or any layered model the user gives (`SeafloorModel(reference=...)`), for example from a refraction profile. Where CRUST1.0 has 500 m of sediment or more, its sediment layers join the reference and sigma becomes 0.5.

Why this prior. The earlier inversion used a flat prior on each layer with a hierarchical penalty on second differences. That penalty extrapolates the shallow gradient into the layers the data do not see and pushes them to the bounds. A Gaussian prior about a reference is proper: where the data say nothing the result returns to the reference with the stated width, and the figure shows it.

The reference is a prior, not a starting model. Chains start from independent draws of the prior, so the result does not depend on where they start; R-hat compares chains that began far apart.

## 5. Sampling

Delayed-rejection adaptive Metropolis with parallel tempering (`compy.inversion.dram`; Haario et al. 2001, 2006; ladder 1.6^0 to 1.6^5), eight chains. Sampling repeats with twice the iterations, up to three rounds, until every parameter has rank-normalized R-hat at most 1.05 and bulk effective sample size at least 400 (Vehtari et al. 2021). The prior is drawn directly, not by a chain. A seed makes a run repeat exactly.

## 6. Outputs

- `<stem>_profile.csv`: Vs against depth (5, 16, 50, 84, 95 %), the prior, the stiffness `M = mu / (1 - nu)`, the width ratio and a `resolved` flag.
- `<stem>.json`: every parameter with its prior, fit, convergence, quality label with reasons, and `resolution`.
- `<stem>.npz`: all draws, the prior draws and the data.
- `<stem>.png`: data and fit, residuals, Vs with prior and reference, and what the data constrain (`compy.plotting.auto.station_figure`).

**Resolution.** For each parameter, `resolution` is the width of the posterior over the width of the prior (16-84 %, of ln Vs for velocities). A depth is marked resolved when the ratio of the layer that holds it is at most 0.6 and the depth is above the level that holds 90 % of the sensitivity of the lowest frequency. Typical deep-ocean values: sediment thickness 0.02-0.25; top 700 m of crust 0.55-1.0; 0.7-6 km 0.15-0.7; the Vp, density and gain factors 1.0. Read unresolved depths as the prior.

**Quality.** *good*: converged, chi-square within the 95 % point for N data, the gain, Vp and density factors within 2 sigma of their priors, some depth resolved, no resolved layer at a hard bound. *usable*: weaker limits (R-hat 1.1, twice the chi-square point, 3 sigma). *poor*: anything else. A pulled factor usually means a wrong calibration: at 7D.J35A (pressure gain 4.4) the inversion asks for a gain factor of 0.83, Vp 25 % above the rock relation and Vs of 4.7 km/s.

## 7. Limits

- **Absolute Vs rests on the Vp/Vs and density priors.** Their uncertainty is in the intervals, but a site whose Vp/Vs differs from the mafic line by more than the 7 % prior (thick serpentinite, a melt body) will be biased. Give a measured Vp profile as the reference where one exists; this is what Crawford (1994, 2004) did.
- **Sediment thickness assumes the sediment P modulus** (Hamilton Vp, Nafe-Drake density, 15 %). It is a compliance thickness of soft material, not a stratigraphic one.
- **The top few hundred metres of basement are not resolved** at deep stations, and thin layers are never resolved one by one.
- **The slope of the curve carries the deep information**, and its error allowance (2-4 %) is a judgement supported by section 3, not a measurement at each station. A 3 % slope error moves Vs at 2-5 km by about 6 %.
- **1-D.** Ridge axes, slopes and seamounts are not 1-D (Crawford et al. 1998; Hulme et al. 2005).
- **Not obtained for this work:** Crawford et al. (1991, 1998), Ruan et al. (2014), Bell et al. (2015), Doran & Laske (2019) and Brocher (2005) in the original. Their content is used only as quoted by sources that were read; the Ruan sediment relation is not used at all.

## References

Read in full or in the sections used: Aminian et al. (2025) GJI 242, ggaf253. Crawford (1994) PhD thesis, UCSD. Crawford (2004) GJI 157, 1130-1145. Crawford et al. (2002) EPSL. Hulme et al. (2003, 2005) GJI. Janiszewski et al. (2019) GJI 217. Mosher et al. (2021) GJI 227. Zha (2014) PhD thesis, Columbia (contains Zha et al. 2014 and Zha & Webb 2016).

Used through a secondary source, named where used: Brocher (2005) BSSA 95 (R package hydrogeo). Carlson & Raskin (1984) (Crawford 1994). Christensen & Shaw (1970) (Zha 2014). Crawford et al. (1991, 1998) JGR (Crawford 1994, 2004; Zha 2014). Haario et al. (2001, 2006) (ComPy 2.3; Zha 2014). Hamilton (1976, 1979) (ComPy 2.3; Zha 2014). Shen et al. (2013) (Zha 2014). Vehtari et al. (2021) (ComPy 2.3).
