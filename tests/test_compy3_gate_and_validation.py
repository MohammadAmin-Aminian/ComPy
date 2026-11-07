"""ComPy 3 audit: the stacked estimator is unchanged, the new gate removes its bias.

The reference below is the 2.8.0 body of ``stacked_compliance`` applied to
precomputed spectra, frozen here so that any change to the default path is
caught bit for bit. The archive test recomputes every published compliance
in ``docs/results`` from its stored stack.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy import signal

from compy import anomaly
from compy.compliance import (
    GJI2025,
    SelectionCriteria,
    _compliance,
    stacked_compliance_from_spectra,
)
from compy.numerics import (
    FREE_AIR_GRADIENT,
    compliance_uncertainty,
    wavenumber,
    welch_averages,
)

ROOT = Path(__file__).resolve().parents[1]
FS, N, NSEG, DEPTH = 1.25, 2625, 1024, 4000.0


def _reference_stack(
    f, windows, depth, gain, band, criteria, nseg, n_samples, output_band
):
    """ComPy 2.8.0 ``stacked_compliance`` loop, on (start, zz, pp, zp, coh) tuples."""
    f_min, f_max = band
    rejected = {"coherence": 0, "pressure_psd": 0, "acceleration_psd": 0}
    zz_sum = pp_sum = zp_sum = coh_sum = None
    curves, starts = [], []
    for start, zz, pp, zp, coh in windows:
        sel = (f >= f_min) & (f <= f_max)
        pp, zp = pp * gain**2, zp * gain
        reason = criteria.check(f, zz, pp, coh, sel)
        if reason is not None:
            rejected[reason] += 1
            continue
        if zz_sum is None:
            zz_sum, pp_sum = np.zeros_like(zz), np.zeros_like(pp)
            zp_sum, coh_sum = np.zeros_like(zp), np.zeros_like(coh)
        zz_sum += zz
        pp_sum += pp
        zp_sum += zp
        coh_sum += coh
        ratio = np.zeros_like(f)
        np.divide(np.abs(zp), pp, out=ratio, where=pp > 0)
        curves.append(_compliance(f, depth, ratio, pp, True))
        starts.append(start)
    n = len(curves)
    szz, spp, spz = zz_sum / n, pp_sum / n, zp_sum / n
    coherence = np.zeros_like(szz)
    np.divide(np.abs(spz), np.sqrt(szz * spp), out=coherence, where=szz * spp > 0)
    coherence = np.clip(coherence, 0, 1)
    ratio = np.zeros_like(szz)
    np.divide(np.abs(spz), spp, out=ratio, where=spp > 0)
    eta = _compliance(f, depth, ratio, spp, True)
    averages = n * welch_averages(n_samples, nseg)
    lo, hi = output_band if output_band is not None else band
    keep = (f >= lo) & (f <= hi)
    omega2 = (2 * np.pi * f[keep]) ** 2
    measured = (
        wavenumber(2 * np.pi * f[keep], depth)
        * omega2
        * ratio[keep]
        / (omega2 + FREE_AIR_GRADIENT)
    )
    curves = np.asarray(curves)
    return dict(
        freqs=f[keep], compliance=eta[keep],
        uncertainty=compliance_uncertainty(measured, coherence[keep], averages),
        coherence=coherence[keep], phase_deg=np.degrees(np.angle(spz[keep])),
        window_mean=curves[:, keep].mean(axis=0),
        window_std=curves[:, keep].std(axis=0, ddof=1 if n > 1 else 0),
        n_windows=n, n_averages=averages, window_starts=starts, rejected=rejected,
        mean_coherence_all=coh_sum / n, S_zz=szz, S_pp=spp, S_pz=spz,
    )  # fmt: skip


def _windows(n_windows=60, coherence=0.9, seed=0):
    """Per-window spectra of displacement-like Z = h P + noise, varying sea state."""
    rng = np.random.default_rng(seed)
    b, a = signal.butter(2, [0.004, 0.04], btype="band", fs=FS)
    kw = dict(
        fs=FS, nperseg=NSEG, noverlap=NSEG // 2, window="hann", detrend="constant"
    )
    out = []
    for i in range(n_windows):
        scale = np.exp(0.6 * rng.standard_normal())
        p = 30.0 * scale * signal.lfilter(b, a, rng.standard_normal(N))
        z = (
            2e-11 * p
            + 2e-11 * 30 * np.sqrt(1 / coherence**2 - 1) * rng.standard_normal(N) * 0.2
        )
        f, zz = signal.welch(z, **kw)
        _, pp = signal.welch(p, **kw)
        _, zp = signal.csd(p, z, **kw)
        coh = np.clip(np.abs(zp) / np.sqrt(zz * pp), 0, 1)
        out.append((float(i * 2100), zz, pp, zp, coh))
    return f, out


@pytest.mark.parametrize("output_band", [None, (0.002, 0.04)])
def test_default_gate_reproduces_the_280_estimator_bit_for_bit(output_band):
    f, windows = _windows()
    band = (0.006, 0.02)
    sel = (f >= band[0]) & (f <= band[1])
    gate = float(
        np.median([np.median(w[4][sel]) for w in windows])
    )  # half the windows pass
    criteria = SelectionCriteria(gate, None, None, "test")
    ref = _reference_stack(f, windows, DEPTH, 0.7, band, criteria, NSEG, N, output_band)
    new = stacked_compliance_from_spectra(f, windows, DEPTH, 0.7, band, criteria, NSEG,
                                          n_samples=N, output_band=output_band)  # fmt: skip
    assert new.gate == "inclusive"
    assert 0 < new.n_windows < len(windows)
    for key in ("freqs", "compliance", "uncertainty", "coherence", "phase_deg",
                "window_mean", "window_std", "mean_coherence_all"):  # fmt: skip
        assert np.array_equal(getattr(new, key), ref[key]), key
    for key in ("S_zz", "S_pp", "S_pz"):
        assert np.array_equal(new.stack[key], ref[key]), key
    assert new.n_windows == ref["n_windows"] and new.n_averages == ref["n_averages"]
    assert new.rejected == ref["rejected"] and new.window_starts == ref["window_starts"]
    assert np.all(new.bin_counts == new.n_windows)


def test_exclusive_gate_matches_inclusive_when_every_window_passes():
    f, windows = _windows(coherence=0.99)
    criteria = SelectionCriteria(0.1, None, None, "all")
    inc = stacked_compliance_from_spectra(
        f, windows, DEPTH, 1.0, (0.006, 0.02), criteria, n_samples=N
    )
    exc = stacked_compliance_from_spectra(f, windows, DEPTH, 1.0, (0.006, 0.02), criteria,
                                          n_samples=N, gate="exclusive")  # fmt: skip
    np.testing.assert_allclose(exc.compliance, inc.compliance, rtol=1e-12)
    np.testing.assert_allclose(exc.uncertainty, inc.uncertainty, rtol=1e-12)
    assert exc.gate == "exclusive" and exc.n_windows == inc.n_windows


def test_unknown_rejection_reason_is_counted():
    class Odd:
        def check(self, f, zz, pp, coh, band):
            return "storm" if pp[band].mean() > np.median(pp[band]) * 1.5 else None

    f, windows = _windows()
    result = stacked_compliance_from_spectra(
        f, windows, DEPTH, 1.0, (0.006, 0.02), Odd(), n_samples=N
    )
    assert result.rejected.get("storm", 0) + result.n_windows == len(windows)


def test_windows_with_another_frequency_axis_raise():
    from compy.compliance import _ComplianceStack

    f, windows = _windows(n_windows=2)
    stack = _ComplianceStack(DEPTH, 1.0, (0.006, 0.02), SelectionCriteria(0.0, None, None),
                             "inclusive", True)  # fmt: skip
    stack.add(f, *windows[0][1:], 0, N)
    with pytest.raises(ValueError, match="sampling rate or nseg"):
        stack.add(2 * f, *windows[1][1:], 1, N)
    with pytest.raises(ValueError, match="gate"):
        stacked_compliance_from_spectra(
            f, windows, DEPTH, 1.0, (0.006, 0.02), n_samples=N, gate="median"
        )


def test_inclusive_gate_biases_h1_high_and_exclusive_gate_does_not():
    # Steady sea state, coherence below the gate: selection on noise.
    inc = anomaly.simulate_gate_bias(0.85, seastate_sigma=0.0, n_windows=1500, seed=1)
    exc = anomaly.simulate_gate_bias(
        0.85, seastate_sigma=0.0, n_windows=1500, seed=1, exclusive=True
    )
    assert 0.25 < inc["acceptance"] < 0.6
    assert inc["bias_pct"] > 2.0
    assert abs(exc["bias_pct"]) < 0.7
    # Variable sea state: selection mostly on signal strength, small bias.
    var = anomaly.simulate_gate_bias(0.85, seastate_sigma=1.0, n_windows=1500, seed=2)
    assert abs(var["bias_pct"]) < 0.5


ARCHIVE = sorted(ROOT.glob("docs/results/*/results/**/*_compliance_*.npz"))


@pytest.mark.skipif(not ARCHIVE, reason="archived results not present")
def test_every_published_compliance_is_reproduced_from_its_stored_stack():
    for path in ARCHIVE:
        with np.load(path) as d:
            f, spp, spz, szz = d["all_freqs"], d["S_pp"], d["S_pz"], d["S_zz"]
            ratio = np.zeros_like(spp)
            np.divide(np.abs(spz), spp, out=ratio, where=spp > 0)
            eta = _compliance(f, float(d["depth"]), ratio, spp, True)
            keep = (f >= d["freqs"][0] - 1e-12) & (f <= d["freqs"][-1] + 1e-12)
            np.testing.assert_allclose(
                eta[keep], d["compliance"], rtol=1e-13, err_msg=str(path)
            )
            coh = np.abs(spz[keep]) / np.sqrt(szz[keep] * spp[keep])
            np.testing.assert_allclose(
                coh, d["coherence"], rtol=1e-12, err_msg=str(path)
            )
    assert len(ARCHIVE) >= 60


def test_stacked_compliance_on_streams_matches_the_array_path():
    obspy = pytest.importorskip("obspy")
    from compy.compliance import spectra, stacked_compliance

    rng = np.random.default_rng(5)
    b, a = signal.butter(2, [0.004, 0.04], btype="band", fs=FS)
    windows = []
    for i in range(12):
        p = 30.0 * signal.lfilter(b, a, rng.standard_normal(N))
        z = 2e-11 * p + 1e-11 * rng.standard_normal(N)
        t0 = obspy.UTCDateTime(2013, 1, 1) + i * N / FS
        windows.append(obspy.Stream([
            obspy.Trace(z, header=dict(channel="BHZ", sampling_rate=FS, starttime=t0)),
            obspy.Trace(p, header=dict(channel="BDH", sampling_rate=FS, starttime=t0)),
        ]))  # fmt: skip
    criteria = SelectionCriteria(0.5, None, None)
    streams = stacked_compliance(windows, DEPTH, 0.8, (0.006, 0.02), criteria, NSEG)
    arrays = [(w[0].stats.starttime, *spectra(w, NSEG)[1:]) for w in windows]
    f = spectra(windows[0], NSEG)[0]
    ref = _reference_stack(
        f, arrays, DEPTH, 0.8, (0.006, 0.02), criteria, NSEG, N, None
    )
    assert np.array_equal(streams.compliance, ref["compliance"])
    assert np.array_equal(streams.uncertainty, ref["uncertainty"])
