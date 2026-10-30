"""The pre-2.0 compliance estimator dropped the wave-attraction term."""

import numpy as np
import pytest

from compy_numerics import legacy_wave_attraction_correction, wavenumber


def _estimators(f, depth, ratio, gamma):
    """Old (pre-2.0) and new (2.0) compliance from the same spectra."""
    w2 = (2 * np.pi * f) ** 2
    k = wavenumber(2 * np.pi * f, depth)
    kh = k * depth
    pa = 2 * np.pi * (6.6743e-11 / 9.8) * (1 + np.exp(-2 * kh)) / 2
    old = k * gamma * w2 / (w2 + 3.07e-6) * (pa + ratio)  # compy.py before 2.0
    new = k * gamma * (w2 * ratio + pa) / (w2 + 3.07e-6)  # compy_processing 2.0
    return old, new


@pytest.mark.parametrize("depth", [2918.0, 4118.0])
def test_correction_converts_legacy_to_current_estimator(depth):
    f = np.geomspace(0.006, 0.02, 12)
    ratio = 1e-7 * (1 + 3 * (f - 0.006) / 0.014)  # displacement / pressure
    gamma = np.linspace(0.8, 0.97, 12)
    old, new = _estimators(f, depth, ratio, gamma)
    corrected = old + legacy_wave_attraction_correction(f, depth, gamma)
    np.testing.assert_allclose(corrected, new, rtol=1e-12)
    rel = (new - old) / old
    assert np.all(rel > 0) and np.all(np.diff(rel) < 0)  # largest at low frequency


def test_correction_validates_inputs():
    with pytest.raises(ValueError):
        legacy_wave_attraction_correction([0.0, 0.01], 3000)
    with pytest.raises(ValueError):
        legacy_wave_attraction_correction([0.01], 3000, coherence=1.2)
