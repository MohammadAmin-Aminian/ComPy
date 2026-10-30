# Teleseismic differential-pressure-gauge calibration: code and equations

**Purpose:** estimate a scalar pressure-amplitude multiplier using earthquake
Rayleigh waves and a pressure/upward-acceleration water-column reference.
These are provisional estimates; no correction has been applied to source
responses, waveforms or saved PPSDs. This report describes the production code
at revision `7724bc5`, before further validation changes.

## 1. Code supplied and how to run it

- Full station implementation: [rayleigh_calibration.py](../../compy_pressure/rayleigh_calibration.py).
- Network wrapper: [rayleigh_network.py](../../compy_pressure/rayleigh_network.py).
- Response and archive helpers: [paired_spectra.py](../../compy_pressure/paired_spectra.py),
  [archive.py](../../compy_pressure/archive.py), [ppsd.py](../../compy_pressure/ppsd.py).
- Scientific tests: [test_rayleigh_calibration.py](../tests/test_rayleigh_calibration.py).
- Prepared dispersion/model: `src/bruitfm/reference_data/rr52_historical_rayleigh_*`.
- Source archive: `output/reports/teleseismic_dpg_source_7724bc5.zip`, with package
  sources, reference tables, tests, project dependencies and this report.
  External waveforms, StationXML, masters and the cached event catalogue remain local.

The following standalone command creates a **separate RR52 comparison folder**
and does not update the master. It recomputes this station; use the viewer below
if you only want yesterday's saved figures. The explicit band avoids the
standalone runner's earlier 0.02-0.06 Hz default.

```bash
cd /Users/mohammadamin/Projects/bruit-fm
conda run --no-capture-output -n bruitfm python -m bruitfm.rayleigh_calibration \
  --network YV --station RR52 \
  --catalog results/RAYLEIGH_REVIEW/catalog_usgs_m5p5_2012_2013.json \
  --dispersion-file src/bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv \
  --frequency-min 0.03 --frequency-max 0.07 \
  --pressure-sign -1 --tilt-method tiskitpy \
  --min-events 3 --max-events 8
```

View the existing result without computing:

```bash
conda run --no-capture-output -n bruitfm python -m bruitfm.network view \
  --network YV --station RR52 --view-kind rayleigh
```

## 2. Physical reference

Let $P$ be positive compressional pressure [Pa], $a$ upward acceleration [m/s²],
$H$ water depth [m], $\rho$ water density [kg/m³], $\alpha$ acoustic sound speed
[m/s], and $c_p(f)$ Rayleigh **phase** speed [m/s]. Define

$$\omega=2\pi f,\qquad r=\omega\sqrt{\alpha^{-2}-c_p(f)^{-2}}.$$

The adopted reference is

$$T(f)=\frac{P(f)}{a(f)}=\rho\frac{\tan(rH)}{r}
=\rho H\frac{\tan(rH)}{rH}.$$

It has units kg/m² and approaches $\rho H$ as $f\to0$. When $r$ is imaginary,
the equivalent expression uses $\tanh$ instead of $\tan$; at $r=0$ use the limit.

Under the adopted homogeneous acoustic model, with depth coordinate $z$ downward,
$P''+r^2P=0$, an ideal pressure-release boundary $P(0)=0$, and floor condition
$P'(H)=\rho a_{up}$. Therefore $P(z)=B\sin(rz)$ and dividing floor pressure by
floor acceleration gives the equation above. This explains the adopted equation;
it does not independently validate its ocean/solid-Earth assumptions.

The code excludes real $rH\ge0.9\pi/2$, before the first tangent pole. Defaults
are $\rho=1028$ kg/m³ and $\alpha=1500$ m/s. Depth is estimated from negative
**station elevation** in StationXML. No measured bathymetry is substituted.
The historical velocity model must be validated for each station/path.

The physical reference follows the documented Zha & Webb (2016) Appendix-A
method, DOI [10.1002/2015JB012688](https://doi.org/10.1002/2015JB012688), as reviewed
in [the source review](rayleigh_calibration_review.md). It is not a full DPG leak/
capillary response model: this stage fits only a scalar within the selected band.
The [publisher's Appendix A](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012688)
was rechecked for this report; the adopted pressure-release/acoustic reference
and scalar-gain limitation are consistent with that method.

## 3. Measured spectra, response correction and tilt

Each synchronized pressure/BH1/BH2/BHZ count segment is demeaned and multiplied
by a periodic Hann taper $w_n$. With NumPy's negative-exponent forward DFT,

$$P_k(f)=\frac{\mathcal F\{w(C_p-\bar C_p)\}}{R_p(f)},\qquad
V_{j,k}(f)=\frac{\mathcal F\{w(C_j-\bar C_j)\}}{R_j^{VEL}(f)}.$$

$R_p$ is the full complex counts/Pa response; $R_j^{VEL}$ is counts/(m/s).
Upward-positive Z is enforced from its dip. Responses are divided **once**.

Actual **TiSKitPy CleanRotator 2.3.1** estimates a three-component vertical
projection $\boldsymbol b$ from a preceding two-hour quiet interval. Its public
`apply()` is used to recover the projection weights. On event coefficients,

$$V_{up,k}=b_1V_{1,k}+b_2V_{2,k}+b_ZV_{Z,k},\qquad
A_k=i\omega V_{up,k}.$$

Quiet copies are converted to velocity and anti-aliased to at most 1 Hz;
event spectra retain their native rate. Constructor band: 0.001-0.01 Hz;
internal optimizer band: 0.001-0.005 Hz; rotation limit: 10 degrees. Orthogonal
sensor axes and stable tilt are assumptions. Pressure is **not** regressed out
of Z because the Rayleigh pressure signal is the calibration target.
[TiSKitPy's official documentation](https://tiskitpy.readthedocs.io/latest/classes/clean_rotator.html)
describes the vertical-noise-minimizing rotation and marks `tfs()` unvalidated;
BRUIT-FM uses `apply()`, leaving horizontal axes unrotated.

For the retained positive-frequency bins, with sampling rate $f_s$,
$q=2/(f_s\sum_n w_n^2)$. Average the powers and complex cross-power over $K$
segments **before** forming ratios:

$$S_{pp}=q\langle|P_k|^2\rangle,\quad
S_{aa}=q\langle|A_k|^2\rangle,\quad
S_{pa}=q\langle sP_kA_k^*\rangle.$$

Here $s=-1$ is the user's provisional pressure-sign assumption for YV. It changes
cross-phase, not the amplitude ratio. Nominal phase is retained separately.
The factor 2 applies to interior one-sided bins, not DC/Nyquist; neither endpoint
is in the current fit band.

$$R(f)=\sqrt{S_{pp}/S_{aa}},\quad H_1(f)=S_{pa}/S_{aa},\quad
\gamma^2(f)=\frac{|S_{pa}|^2}{S_{pp}S_{aa}},\quad
\phi(f)=\arg S_{pa}.$$

The production gain uses $R$, an unsigned PSD amplitude ratio. $|H_1|$ is saved
as an alternative estimator, and phase is a separate QC check. The two estimators
have different noise biases; their agreement alone is not absolute calibration.

## 4. Gain objective and interpretation

For accepted frequency bins $i$, the exact objective is

$$g=\arg\min_{g>0}\sum_i(gR_i-T_i)^2
=\frac{\sum_iR_iT_i}{\sum_iR_i^2}.$$

The analytical value is checked against the grid 0.1-3.0 at spacing 0.001.
The saved gain is the analytical value, not restricted to three decimals.
Boundary grid minima are unresolved. The objective is unweighted and fits
$gR$ to $T$; fitting $GT$ to $R$ and then inverting is a different loss.

If $K_{nom}$ is nominal sensitivity in counts/Pa, the candidate update would be
$P_{candidate}=gP_{nom}$ and $K_{candidate}=K_{nom}/g$ for amplitude. A separate
signed pressure multiplier is $sg$; polarity still requires confirmation.
For RR52, median event gain is 0.678998, implying amplitude sensitivity about
1.473 times nominal. This is not applied. The plotted gain of the median ratio,
0.674374, is a distinct aggregate estimate. Neither is a hard-coded target.

## 5. Event selection and current limitations

The prepared table was generated with **disba 0.7.0** fundamental-mode phase
and group dispersion. Runtime reads this table; it does not solve dispersion
again. Group speed predicts $t_R=t_0+D/c_g(0.04\,Hz)$ for a 30-minute window;
phase speed enters $T$. TauP `iasp91` supplies separate P/S arrival predictions.
The predicted window is a candidate, not an automatic proof of a Rayleigh mode.
The current event selector has no explicit minimum epicentral-distance filter;
some selected large earthquakes are regional. Magnitude alone does not establish
a teleseismic or isolated fundamental-mode wavefield.

(Historical pilot policy, now `--event-sampling bounded`; the current default
`all` evaluates every eligible event at all three levels.)
Search magnitude >7, then >=6.5, then >=6; at most eight attempts per level and
stop after three QC candidates. Source depth <150 km and 24-hour event separation
are explicit pilot choices. Quiet screening tries 6/12/24/48 hours before origin,
using the documented catalogue-based coda heuristic. It can miss transients.

Nine nonoverlapping 200-second segments in a 30-minute window give 0.005 Hz bins
(at least eight segments are required). Candidate rules: mean and median squared
coherence >=0.95; at least three bins with coherence >=0.95 and phase within
+/-20 degrees; no predicted overlapping >=6 earthquake packet; non-boundary gain.
The three-event count and thresholds are proposals, not instrument standards.

**Current outcome:** 43 YV stations processed; 16 provisional candidate estimates,
25 unresolved and two conservative reference-band exclusions. Existing outputs
remain preserved. The saved-spectrum arithmetic audit passed, but an independent
raw spectral estimate, packet/quiet-window reviews and total uncertainty are
still required. RR03's reviewed window contains predicted S energy timing;
RR57 is borderline under the conservative resonance preflight. Depth changes of
+/-5% can move gains much more than their event scatter. Actual gauge type,
polarity, vertical absolute response, FIR metadata, depth, Earth/water model,
tilt stationarity and frequency-dependent residuals need verification.

## 6. Exact core production functions

These excerpts are copied verbatim from the production module. They are review
excerpts, not the entire standalone program; the source archive contains all
package dependencies needed by the module. The full spectral/event workflow is
in `corrected_fft()`, `event_spectra()`, `estimate_tilt()` and `run_calibration()`.

```python
def real_array(value):
    if np.ma.isMaskedArray(value) and np.any(np.ma.getmaskarray(value)):
        raise ValueError("Masked values are not measurements")
    array = np.asarray(value)
    if array.dtype.kind not in "iuf" or not np.isfinite(array).all():
        raise ValueError("Expected finite real values")
    return array.astype(float)
```

```python
def water_ratio(f, phase_velocity, depth, density, sound_speed):
    """Zha/Webb A5: P/upward-a, SI units; reject the first acoustic pole."""
    f, c = np.broadcast_arrays(real_array(f), real_array(phase_velocity))
    if any(not np.isfinite(v) or v <= 0 for v in (depth, density, sound_speed)):
        raise ValueError("Water parameters must be positive and finite")
    if (
        not np.isfinite(f).all()
        or not np.isfinite(c).all()
        or (f < 0).any()
        or (c <= 0).any()
    ):
        raise ValueError("Invalid frequency or phase velocity")
    q2 = sound_speed**-2 - c**-2
    x = 2 * np.pi * f * depth * np.sqrt(abs(q2))
    if np.any((q2 >= 0) & (x >= 0.9 * np.pi / 2)):
        raise ValueError("Band approaches first water-column resonance")
    result = np.ones_like(x)
    small = x < 1e-5
    result[small] = 1 + np.sign(q2[small]) * x[small] ** 2 / 3 + 2 * x[small] ** 4 / 15
    real = ~small & (q2 >= 0)
    imaginary = ~small & (q2 < 0)
    result[real] = np.tan(x[real]) / x[real]
    result[imaginary] = np.tanh(x[imaginary]) / x[imaginary]
    return density * depth * result
```

```python
def fit_multiplier(observed, expected):
    """Fit g*observed=expected; g multiplies pressure, 1/g multiplies counts/Pa."""
    observed, expected = real_array(observed), real_array(expected)
    if observed.shape != expected.shape or observed.ndim != 1 or len(observed) < 3:
        raise ValueError("At least three corresponding frequency bins required")
    if (
        not (np.isfinite(observed).all() and np.isfinite(expected).all())
        or np.any(observed <= 0)
        or np.any(expected <= 0)
    ):
        raise ValueError("Fit ratios must be finite and positive")
    # Normalize first so SI scaling cannot overflow the objective.
    scale = np.max(expected)
    r, t = observed / scale, expected / scale
    g = float(np.dot(r, t) / np.dot(r, r))
    grid = np.linspace(0.1, 3.0, 2901)
    loss = np.mean((grid[:, None] * r - t) ** 2, axis=1)
    index = int(np.argmin(loss))
    if index in {0, len(grid) - 1}:
        raise ValueError(
            "Grid minimum lies on a boundary; expand justified search bounds"
        )
    if abs(g - grid[index]) > (grid[1] - grid[0]) / 2 + 1e-10:
        raise ValueError("Analytical and grid-search solutions disagree")
    return (
        {
            "pressure_multiplier": g,
            "sensitivity_counts_per_pa_multiplier": 1 / g,
            "grid_pressure_multiplier": float(grid[index]),
            "relative_rms": float(np.linalg.norm(g * r - t) / np.linalg.norm(t)),
        },
        grid,
        loss,
    )
```

[Assumptions/checklist](calibration_assumptions.md) · [Initial saved-results audit](yv_rayleigh_saved_results_review.md).

## 7. Subsequent validation and one reader safeguard

The source archive above freezes the code used for yesterday's calculation.
The current checkout additionally rejects within-channel fractional sample-grid
offsets before merging calibration records. This prevents silent alignment of
overlapping records and records a timing-review reason. The pressure formula,
gain objective, QC thresholds and original saved results have not been changed.

The follow-up found 20 of 48 candidate windows contain predicted P/S arrival
times; seven independent raw-count Welch estimates match the saved spectra;
46 quiet intervals passed independent continuity screening, with two flagged
fractional-grid boundaries. Pressure metadata at 34 stations contain positive-real
Laplace poles. The RESIF service reproduces the local RR03/RR22 poles and RR52
DPG metadata; matching archive metadata does not confirm physical correctness.
No poles, FIR coefficients, source samples or responses were changed.
[Independent validation report](yv_rayleigh_independent_validation.md).
