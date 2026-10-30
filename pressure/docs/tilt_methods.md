# Tilt correction used by the Rayleigh pilot

## Optional: no tilt correction

Use `--tilt-method none` for the standalone calibration or
`--rayleigh-tilt-method none` for the network stage. This uses only the pressure
channel and the matched vertical channel. It response-corrects vertical velocity,
makes it upward-positive from inventory dip, and forms acceleration with `i*2πf`.
It does not require BH1/BH2, a preceding quiet interval or TiSKitPy. No rotation
is estimated or applied. Event-window timing, response, coherence, phase and
gain safeguards still apply.

Results explicitly record `tilt_correction_applied: false`, `method: none` and
null quiet-training/tilt measurements; tilt-only review flags are omitted.
Each standalone run saves a new folder and leaves the master unchanged.
The existing default remains `tiskitpy` for compatibility with earlier commands;
select `none` explicitly for the user's unrotated teleseismic method.

See [YV commands](yv_rayleigh_calibration.md#comparison-without-tilt-correction)
for a separate network comparison and viewing without master updates.

## Default: TiSKitPy CleanRotator

The default is the actual `tiskitpy.CleanRotator` class, pinned to **2.3.1** in
the optional `rayleigh` extra. There is no silent fallback to the local method.
Install it in the BRUIT-FM environment with:

```bash
python -m pip install -e '.[rayleigh]'
```

A preceding two-hour quiet interval is screened against the supplied earthquake
catalogue. Temporary copies of BH1/BH2/BHZ are response-corrected into the same
velocity units; Z is made upward-positive using inventory dip. A tapered
deconvolution prefilter is 0.0003/0.0005/0.03/0.05 Hz. Only these quiet training
copies are anti-aliased to at most 1 Hz using `scipy.signal.resample_poly`.
The earthquake spectra retain the original sample rate.

`CleanRotator` receives `filt_band=(0.001, 0.01)`, `H_over_Z=1` and sensor
components `("1", "2")`. Version 2.3.1's internal optimizer also applies its
0.001–0.005 Hz filter; both bands are recorded. Its internal earthquake download
is disabled because the supplied catalogue has already screened the interval.
This screening cannot guarantee that every transient has been excluded.

The fitted rotation is applied to three unit-basis traces with the public
`apply(horiz_too=False, rot_limit=10, set_dtype=False)` method. The resulting
vertical coefficients are used on the response-corrected event Fourier
coefficients. This reproduces the library's linear vertical projection while
preserving horizontal traces and floating-point precision. We do not use its
`tfs()` method, which the documentation identifies as unvalidated.

RR52 horizontal axes are assumed orthogonal in the sensor frame; both inventory
azimuths are zero, so their geographic orientation is unverified. A small fitted
angle and lower quiet power do not establish removal of all target-band noise.
Pressure is not subtracted from Z: it is the signal being measured for calibration.

Reference: Wayne Crawford, *TiSKitPy*, version 2.3.1,
[CleanRotator documentation](https://tiskitpy.readthedocs.io/latest/classes/clean_rotator.html)
and [source](https://github.com/WayneCrawford/tiskitpy).

## Experimental local method

`--tilt-method local_eigen` retains BRUIT-FM's earlier implementation, introduced
in commit `df0b449`. Every use saves `validation_status: "needs_to_be_tested"`.
It is not the default.

The method builds a real covariance matrix from the three response-corrected
quiet velocity components in 0.001–0.01 Hz. Its smallest-eigenvalue eigenvector
minimizes projected vertical power; the axis is oriented upward and limited to
10° from the nominal vertical. This is a local covariance method, not a copy of
the complete CleanRotator processing chain. Synthetic checks exist; comparative
field testing of contamination, preprocessing, angle limits and event gains
remains pending.

Implementation reference: `quiet_axis()` and `tilt_direction()` in
[rayleigh_calibration.py](../../compy_pressure/rayleigh_calibration.py). Mathematical
background: Jolliffe & Cadima (2016), *Principal component analysis: a review and
recent developments*, [doi:10.1098/rsta.2015.0202](https://doi.org/10.1098/rsta.2015.0202).
That reference explains covariance/eigenvector methods; it does not validate
this OBS correction.
