# Alternative inversions

ComPy 2.8 adds inversions that are separate from the main Bayesian workflow in
`compy.inversion.dram`, to test how much a Vs model depends on the method. They have their
own forward model and do not call `dram`, `propagator` or `occam`. Stage E stays the
main result: these are checks.

| Module | What it does |
|---|---|
| `compy.inversion.impedance` | Forward model by surface-impedance recursion |
| `compy.inversion.literature_mcmc` | Monte Carlo inversion with the parameterizations of Zha & Webb (2016) and Zheng et al. (2024) |
| `compy.inversion.minimum_structure` | Minimum-structure inversion after Crawford (2004) |
| `pipeline/e_alternative_inversions.py` | All three on one station's stage C result |

## The Moho rule

Vs must increase across the Moho. Without this rule, the last crustal layer was faster
than the mantle in 7 to 81 % of the stage E draws at eight YV stations. It is a hard
prior: a model that breaks it is rejected, as in the Bayesian Monte Carlo inversions
of Shen et al. (2013), and mantle Vs in `literature_mcmc` goes up to 4.9 km/s, their
limit.

| Where | How |
|---|---|
| `CrustalModel(require_moho_jump=True)` | Rejected in the target. Off by default, so published scripts repeat |
| Every variant of `pipeline/e_invert.py` | On |
| `ZhaWebb2016Model`, `Zheng2024Model` | Always on |
| `invert_minimum_structure` | On by default: every layer stays at least 1 m/s below the fixed half-space |
| `invert_occam(moho_jump=True)` | Off by default; same bound |

**Stage E results therefore change compared with ComPy 2.7.**

## Forward model

Depth z is positive downward and fields go as exp(i k x). In each layer,

```
nu_a = sqrt(k^2 - omega^2/Vp^2),  nu_b = sqrt(k^2 - omega^2/Vs^2)
gamma = mu (2 k^2 - omega^2/Vs^2)

mode          u_x     u_z     t_xz               t_zz
P, decaying   i k     -nu_a   -2 i mu k nu_a     gamma
P, growing    i k     +nu_a   +2 i mu k nu_a     gamma
S, decaying   +nu_b   i k     -gamma             -2 i mu k nu_b
S, growing    -nu_b   i k     -gamma             +2 i mu k nu_b
```

The impedance Z, with traction = Z displacement, starts from the decaying modes of the
half-space and is carried to the seafloor. Decaying amplitudes refer to the layer top
and growing ones to the layer base, so no exponential grows. A pressure p gives
t_xz = 0 and t_zz = -p, so `eta = k u_z / p = -k (Z^-1)_22`.

| Check | Largest relative difference |
|---|---|
| Exact homogeneous half-space (in the tests) | below 1e-10; measured 1e-13 |
| Matrix-exponential integration of the elastic equations (in the tests) | below 1e-7; measured 5e-10 on these models |
| `compy.inversion.dram` forward, 300 random models | 6e-6 |
| 60-digit calculation, the 3 models tried, picked by how much the two double-precision codes disagree | new 1.3e-6, `dram` 4.9e-6 |

The errors at the 1e-6 level are rounding: the P and S modes are nearly degenerate when
the water wave is much slower than Vs. They are about a thousand times smaller than the
data errors.

## Monte Carlo with published parameterizations

The likelihood is `exp(-(chi2 + penalty) / 2)` with independent errors. There are 10
chains from prior draws. The proposal is an adaptive Gaussian, frozen after burn-in. The
defaults are 400,000 iterations, 200,000 burn-in and thinning by 40. The papers use far
fewer, but synthetic tests here had not converged at 50,000. The misfit filter of each
paper is then applied.

| | `ZhaWebb2016Model` | `Zheng2024Model` |
|---|---|---|
| Sediment thickness | 0-2000 m | 0-2000 m |
| Sediment Vs | Fixed, Hamilton (1979) | Top and base sampled, linear |
| Sediment Vp | Fixed, Hamilton (1979) | 1.16 Vs + 1.36 km/s (Castagna et al. 1985), Vp ≤ 4 km/s, 1.7 ≤ Vp/Vs ≤ 15 |
| Sediment density | 2000 kg/m³ | 2000 kg/m³ |
| Crust | 100 m layers doubling to the Moho, each Vs sampled (100-4500 m/s) | 4 cubic B-splines to 5 km, increasing with depth, ≤ 4500 m/s |
| Crustal Vp | 1.85 Vs | Brocher (2005) |
| Crustal density | Fixed by the starting Vp: 1850 + 0.165 Vp | 2800 kg/m³ |
| Mantle Vs | 3800-4900 m/s, above the last crustal Vs | 3800-4900 m/s, above the last crustal Vs |
| Mantle Vp, density | 1.85 Vs; from the starting Vp | Brocher (2005); 3000 kg/m³ |
| Smoothing | Squared second differences of crustal Vs, weight 1 per (km/s)² | None |
| Retained | chi2/N < max(2 min, min + 0.25) | rms < 2 min, or min + 0.5 if min < 0.5 |

Where these differ from the papers:

- **Zha & Webb:**
  - Sediment Vs is Hamilton (1979), the alternative Zha tests, because the coefficients of
    Ruan et al. (2014) were not available.
  - Thickness is continuous, because a 50 m grid made chains stick in neighbouring cells
    (R-hat 1.3 on a synthetic test).
  - Zha's mantle range of 3800-4200 m/s suits a hot back-arc and is replaced.
  - The paper gives no crustal Vs limits or roughness weight, so those are defaults here.
    The weight moved Vs at 1 km by up to about 15 % when varied from 0 to 10.
- **Zheng et al.:**
  - The sampler is Metropolis, not No-U-Turn.
  - There are no Ps delays or GlobSed prior.
  - The published model stops at 5 km. Here the last spline value continues to the Moho,
    above a mantle half-space.
- **Both:** the Moho and starting Vp come from the caller, CRUST1.0 in the pipeline.

The data errors in the pipeline are the standard deviation of the monthly curves, the
statistic Zheng et al. compute on daily curves. It is larger than a standard error.

Not included:

- Ball et al. (2014): the full text was not available.
- Mosher et al. (2021): it uses a neural network, not Monte Carlo.
- Doran & Laske (2021): it follows Ball and Zha.

## Minimum-structure inversion

Crawford (2004, equations 10 to 14) with the Constable et al. (1987) iteration:

```
m = ((W J)^T W J + mu D^T D)^-1 (W J)^T W (d - F(m0) + J m0)
```

The model m is Vs, D holds second differences and W = diag(1/sigma). Layers start at
50 m and grow by a factor of 1.1. Vp, density and the half-space are fixed, and the
target is chi2 = N. Each iteration keeps the best fit until the target is reached,
then the smoothest fitting model.

Differences from the paper:

- Each solve is bounded, because clipping an unbounded solution stalled on real data.
  Vs lies between 50 m/s and Vp/√2, and the Moho rule applies.
- The starting Vp/Vs of 1.87 is from Crawford (1994).

`invert_occam` remains for ln(Vs) smoothing, and the two can differ by a median of 3-29 %.

## Eight YV stations

Same compliance in every column; RR36, RR38, RR40 and RR50 use their stable epoch. Each
cell gives Vs at 2 km below the seafloor [m/s] and chi2 / N. Stage E uses its full
covariance, the Monte Carlo setups the monthly standard deviation, and minimum structure the
marginal errors of stage E, so chi2 values compare only within a column. Bold marks
chi2 > N.

| Station | Stage E, free Vp/Vs | Stage E, Brocher Vp | Stage E, CRUST1.0 Vp | Zha & Webb | Zheng et al. | Minimum structure |
|---|---|---|---|---|---|---|
| RR28 | 3369, 3.7 / 6 | 3962, **7.3** / 6 | 3072, 3.7 / 6 | 3788, 1.3 / 6 | 3738, 2.3 / 6 | 2498, 2.1 / 6 |
| RR29 | 3598, 3.8 / 5 | 4194, **6.1** / 5 | 3527, **5.1** / 5 | 3892, 1.0 / 5 | 4005, **6.9** / 5 | 4259, **293** / 5 |
| RR34 | 3166, 3.0 / 6 | 3606, 3.5 / 6 | 4650, **70** / 6 | 3251, 1.2 / 6 | 3337, 1.8 / 6 | 4445, **390** / 6 |
| RR36 | 3145, 3.3 / 8 | 3761, 5.2 / 8 | 2739, 3.5 / 8 | 3439, 1.3 / 8 | 3826, 4.1 / 8 | 2318, 1.2 / 8 |
| RR38 | 4498, **10.5** / 6 | 4624, **32** / 6 | 3437, 2.5 / 6 | 4299, **10.1** / 6 | 4188, **6.3** / 6 | 2849, 3.6 / 6 |
| RR40 | 3668, 2.3 / 5 | 4006, 2.1 / 5 | 3552, 2.3 / 5 | 3980, 0.9 / 5 | 3935, 1.4 / 5 | 2375, 5.0 / 5 |
| RR50 | 3949, 2.5 / 6 | 4424, 2.4 / 6 | 4241, **71** / 6 | 4263, 1.1 / 6 | 4073, **7.7** / 6 | 4137, **227** / 6 |
| RR52 | 3607, 3.2 / 9 | 4086, 3.0 / 9 | 3838, **9.9** / 9 | 4031, 1.3 / 9 | 4234, 3.2 / 9 | 3954, **13** / 9 |

Vs at 2 km changes by hundreds of m/s, and at some stations by more than 1000 m/s, with
the Vp assumption and the method, often with fits of similar quality. A fixed CRUST1.0 Vp
cannot fit RR34 and RR50 in stage E or RR29, RR34, RR50 and RR52 in the minimum-structure
inversion, and at RR34, RR50 and RR52 stage E pushes the gain factor to its lower bound
of 0.95. Mantle Vs is weakly constrained by these data.

## Use

```python
from compy.inversion import (
    Zheng2024Model,
    ZhaWebb2016Model,
    sample_literature_model,
    crawford_thin_layers,
    invert_minimum_structure,
    crust1_column,
    crust1_model,
)

column = crust1_column(lat, lon, "path/to/crust1")
result = sample_literature_model(
    Zheng2024Model(column["moho_depth"]), depth, freqs, eta, sigma, processes=8
)
print(result.summary()["max_rhat"], result.summary()["moho_jump"])

grid, edges = crawford_thin_layers(crust1_model(column))
smooth = invert_minimum_structure(eta, freqs, depth, grid, sigma)
print(smooth.chi2, smooth.reached_target)
```

From a pipeline run:

```bash
python e_alternative_inversions.py RR52 minimum_structure zha_webb_2016 zheng_2024
python e_invert.py RR52 final final_brocher final_crust1
```

## References

- Ball, J. S., et al. (2014). G-cubed 15, 5051.
- Brocher, T. M. (2005). BSSA 95, 2081.
- Castagna, J. P., et al. (1985). Geophysics 50, 571.
- Constable, S. C., Parker, R. L. & Constable, C. G. (1987). Geophysics 52, 289.
- Crawford, W. C. (1994). PhD thesis, UC San Diego.
- Crawford, W. C. (2004). GJI 157, 1130.
- Doran, A. K. & Laske, G. (2021). GJI, doi:10.1093/gji/ggaa470.
- Haario, H., Saksman, E. & Tamminen, J. (2001). Bernoulli 7, 223.
- Hamilton, E. L. (1979). JASA 66, 1093.
- Mosher, S. G., et al. (2021). GJI 227, 1879.
- Ruan, Y., Forsyth, D. W. & Bell, S. W. (2014). JGR 119, 6357.
- Shen, W., et al. (2013). GJI 192, 807.
- Zha, Y. (2015). PhD thesis, Columbia University.
- Zha, Y. & Webb, S. C. (2016). JGR 121, doi:10.1002/2015JB012688.
- Zheng, M., et al. (2024). JGR 129, e2024JB029862.
