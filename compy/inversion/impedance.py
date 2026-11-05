"""Normalized compliance by a surface-impedance recursion.

A forward model written independently of ``compy.inversion.propagator``, so that
the two can check each other. Depth z is positive downward. In each elastic layer,
potentials proportional to exp(i k x) give four modes:

    nu_a = sqrt(k**2 - omega**2 / Vp**2),  nu_b = sqrt(k**2 - omega**2 / Vs**2)
    gamma = mu (2 k**2 - omega**2 / Vs**2)

    mode          u_x     u_z     t_xz               t_zz
    P, z decays   i k     -nu_a   -2 i mu k nu_a     gamma
    P, z grows    i k     +nu_a   +2 i mu k nu_a     gamma
    S, z decays   +nu_b   i k     -gamma             -2 i mu k nu_b
    S, z grows    -nu_b   i k     -gamma             +2 i mu k nu_b

The impedance Z maps displacement to traction, (t_xz, t_zz) = Z (u_x, u_z). It
starts from the decaying modes of the half-space and is carried up one layer at a
time. Decaying amplitudes are referred to the layer top and growing amplitudes to
its base, so no exponential grows and thick layers stay stable. At the seafloor a
pressure p gives t_xz = 0 and t_zz = -p, hence

    eta = k u_z / p = -k (Z^-1)_22

The model is elastic and 1-D, and the water enters only through k.

Accuracy, against a 60-digit calculation of 300 random layered models: the largest
relative error is about 1e-6, which is rounding from the near degeneracy of the P
and S modes when the water wave is much slower than Vs. ``exact_halfspace_compliance``
and ``matrix_exponential_compliance`` are references for tests, not for sampling.
"""

import numpy as np
from scipy.linalg import expm

from compy.numerics import GRAVITY, positive

try:
    from numba import njit

    _jit = njit(cache=True)
except ImportError:  # Numba is optional; the same code runs as plain Python.

    def _jit(function):
        return function


def gravity_wave_wavenumber(omega, depth, gravity=GRAVITY):
    """Wavenumber [rad/m] of the surface gravity wave, by Newton iterations.

    Solves omega**2 = g k tanh(k H) from the shallow-water value k = omega / sqrt(g H).

    Args:
        omega (array): Positive angular frequencies [rad/s].
        depth (float): Water depth [m].
        gravity (float): Gravitational acceleration [m/s**2].

    Returns:
        ndarray: Wavenumbers with the shape of ``omega``.
    """
    depth, gravity = positive(depth, "water depth"), positive(gravity, "gravity")
    omega = np.asarray(omega, float)
    k = omega / np.sqrt(gravity * depth)
    for _ in range(60):
        step = (gravity * k * np.tanh(k * depth) - omega**2) / (
            gravity * np.tanh(k * depth) + gravity * k * depth / np.cosh(k * depth) ** 2
        )
        k = k - step
        if np.all(np.abs(step) <= 1e-15 * k):
            break
    return k


@_jit
def _multiply(a, b):
    return (
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    )


@_jit
def _inverse(a):
    det = a[0] * a[3] - a[1] * a[2]
    return (a[3] / det, -a[1] / det, -a[2] / det, a[0] / det)


@_jit
def _modes(k, omega, rho, vp, vs):
    """Displacement and traction of the decaying and growing modes, as 2 x 2 tuples."""
    mu = rho * vs * vs
    kb2 = (omega / vs) ** 2
    nu_a = np.sqrt(complex(k * k - (omega / vp) ** 2, 0.0))
    nu_b = np.sqrt(complex(k * k - kb2, 0.0))
    gamma = complex(mu * (2.0 * k * k - kb2), 0.0)
    ik = complex(0.0, k)
    shear = 2.0 * mu * ik
    return (
        (ik, nu_b, -nu_a, ik),
        (-shear * nu_a, -gamma, gamma, -shear * nu_b),
        (ik, -nu_b, nu_a, ik),
        (shear * nu_a, -gamma, gamma, shear * nu_b),
        nu_a,
        nu_b,
    )


@_jit
def _impedance_compliance(k, omega, thickness, rho, vp, vs):
    """Compliance at each (k, omega); the last row is the half-space."""
    out = np.empty(k.shape[0])
    last = thickness.shape[0] - 1
    for j in range(k.shape[0]):
        kj, om = k[j], omega[j]
        u_d, t_d, _, _, _, _ = _modes(kj, om, rho[last], vp[last], vs[last])
        z = _multiply(t_d, _inverse(u_d))
        for i in range(last - 1, -1, -1):
            u_d, t_d, u_g, t_g, nu_a, nu_b = _modes(kj, om, rho[i], vp[i], vs[i])
            ea, eb = np.exp(-nu_a * thickness[i]), np.exp(-nu_b * thickness[i])
            # Radiation condition at the layer base: t = Z u gives growing per decaying amplitude.
            ud_base = (u_d[0] * ea, u_d[1] * eb, u_d[2] * ea, u_d[3] * eb)
            td_base = (t_d[0] * ea, t_d[1] * eb, t_d[2] * ea, t_d[3] * eb)
            z_ud, z_ug = _multiply(z, ud_base), _multiply(z, u_g)
            decaying = (
                td_base[0] - z_ud[0],
                td_base[1] - z_ud[1],
                td_base[2] - z_ud[2],
                td_base[3] - z_ud[3],
            )
            growing = (
                t_g[0] - z_ug[0],
                t_g[1] - z_ug[1],
                t_g[2] - z_ug[2],
                t_g[3] - z_ug[3],
            )
            r = _multiply(_inverse(growing), decaying)
            r = (-r[0] * ea, -r[1] * ea, -r[2] * eb, -r[3] * eb)  # referred to the top
            gu, gt = _multiply(u_g, r), _multiply(t_g, r)
            u_top = (u_d[0] + gu[0], u_d[1] + gu[1], u_d[2] + gu[2], u_d[3] + gu[3])
            t_top = (t_d[0] + gt[0], t_d[1] + gt[1], t_d[2] + gt[2], t_d[3] + gt[3])
            z = _multiply(t_top, _inverse(u_top))
        out[j] = -kj * _inverse(z)[3].real
    return out


class ImpedanceForward:
    """Compliance of layered models at a fixed water depth and fixed frequencies.

    Args:
        depth (float): Water depth [m].
        freqs (array): Positive frequencies [Hz].
        gravity (float): Gravitational acceleration [m/s**2].

    Calling the instance with arrays ``thickness, rho, vp, vs`` (last row the
    half-space, its thickness ignored) returns compliance [1/Pa]. There are no
    argument checks in the call, because samplers make millions of them.
    """

    def __init__(self, depth, freqs, gravity=GRAVITY):
        self.depth = positive(depth, "water depth")
        self.freqs = np.asarray(freqs, float)
        if self.freqs.ndim != 1 or not self.freqs.size or np.any(~(self.freqs > 0)):
            raise ValueError("frequencies must be a nonempty vector of positive values")
        self.omega = 2 * np.pi * self.freqs
        self.k = gravity_wave_wavenumber(self.omega, self.depth, gravity)

    def __call__(self, thickness, rho, vp, vs):
        return _impedance_compliance(self.k, self.omega, thickness, rho, vp, vs)


def compliance_by_impedance(depth, freqs, model, gravity=GRAVITY):
    """Normalized compliance [1/Pa] of ``model`` by the surface-impedance recursion.

    Args:
        depth (float): Water depth [m].
        freqs (array): Positive frequencies [Hz].
        model (array): Rows ``[thickness, rho, vp, vs]``; the last row is the half-space.
        gravity (float): Gravitational acceleration [m/s**2].

    Returns:
        ndarray: Compliance per frequency.
    """
    m = np.ascontiguousarray(model, dtype=float)
    if m.ndim != 2 or m.shape[1] != 4 or len(m) < 1 or np.any(~np.isfinite(m)):
        raise ValueError("model must have finite columns [thickness, rho, vp, vs]")
    if np.any(m[:, 1:] <= 0) or np.any(m[:-1, 0] <= 0):
        raise ValueError("density, velocities and finite thicknesses must be positive")
    return ImpedanceForward(depth, freqs, gravity)(
        m[:, 0].copy(), m[:, 1].copy(), m[:, 2].copy(), m[:, 3].copy()
    )


def exact_halfspace_compliance(k, omega, rho, vp, vs):
    """Exact compliance of a homogeneous half-space (Rayleigh-function form).

    With c = omega / k, r = sqrt(1 - c**2 / Vp**2) and s = sqrt(1 - c**2 / Vs**2):
    eta = (c**2 / Vs**2) r / (mu [4 r s - (2 - c**2 / Vs**2)**2]). For c -> 0 this is
    the static (lambda + 2 mu) / (2 mu (lambda + mu)).
    """
    c = np.asarray(omega, float) / np.asarray(k, float)
    r, s = np.sqrt(1 - (c / vp) ** 2), np.sqrt(1 - (c / vs) ** 2)
    ratio = (c / vs) ** 2
    return ratio * r / (rho * vs**2 * (4 * r * s - (2 - ratio) ** 2))


def _stress_displacement_system(k, omega, rho, vp, vs):
    """Matrix A of dy/dz = A y for y = (u_x, u_z, t_xz, t_zz), from the elastic equations."""
    mu = rho * vs * vs
    modulus = rho * vp * vp
    lam = modulus - 2 * mu
    i = 1j * k
    return np.array(
        [
            [0, -i, 1 / mu, 0],
            [-i * lam / modulus, 0, 0, 1 / modulus],
            [
                -rho * omega**2 + k**2 * (modulus - lam**2 / modulus),
                0,
                0,
                -i * lam / modulus,
            ],
            [0, -rho * omega**2, -i, 0],
        ],
        complex,
    )


def matrix_exponential_compliance(k, omega, thickness, rho, vp, vs):
    """Compliance at one (k, omega) by integrating the elastic equations with expm.

    It shares no mode formulas with the recursion: it starts from the two downward-
    decaying eigenvectors of the half-space system and propagates them up. It is slow
    and is meant as a test reference.
    """
    w, v = np.linalg.eig(_stress_displacement_system(k, omega, rho[-1], vp[-1], vs[-1]))
    y = v[:, np.real(w) < 0]
    for i in range(len(thickness) - 2, -1, -1):
        y = (
            expm(
                -_stress_displacement_system(k, omega, rho[i], vp[i], vs[i])
                * thickness[i]
            )
            @ y
        )
        y = y / np.max(np.abs(y), axis=0)
    amplitudes = np.linalg.solve(y[2:, :], np.array([0.0, -1.0], complex))
    return float(k * (y[:2, :] @ amplitudes)[1].real)
