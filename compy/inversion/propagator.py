"""Forward model: normalized compliance of a layered elastic seafloor.

``raydep``, ``argdtray``, ``gravd`` and ``dtanh`` keep the names of Wayne
Crawford's MATLAB compliance code they were ported from (minor-vector
P-SV propagator; Crawford et al. 1991). ``normalized_compliance`` evaluates a
model at a set of frequencies and ``compliance_kernels`` gives its Vs
sensitivity.
"""

import numpy as np

from compy.inversion.models import _validate_model
from compy.numerics import positive, wavenumber


def gravd(W, h):
    """Gravity-wave wavenumber in rad/m; scalar input retains legacy 1D output."""
    return np.atleast_1d(wavenumber(W, h))


# stable hyperbolic tangent
def dtanh(x):
    """Stable hyperbolic tangent, including x=0 and negative arguments."""
    return np.tanh(x)


def argdtray(wd, h):
    """Cosine and scaled-sine propagator terms of one layer (Crawford's argdtray).

    The phase is ``wd * sqrt(|h|)``; the wave propagates for ``h <= 0`` and is
    evanescent (hyperbolic terms) for ``h > 0``.
    """

    hh = np.sqrt(abs(h))  # % magnitude of wavenumber/freq
    th = wd * hh  # % number of waves (or e-foldings) in layer in radians
    if th >= 1.5e-14:
        if h <= 0:  # % propagating wave
            c = np.cos(th)
            s = -np.sin(th) / hh
        else:  # % evenescent wave
            d = np.exp(th)
            c = 0.5 * (d + 1 / d)
            s = -0.5 * (d - 1 / d) / hh
    else:
        c = 1
        s = -wd

    return c, s


def raydep(P, om, d, ro, vp2, vs2):
    """
        %RAYDEP	Propagator matrix sol'n for P-SV waves, minor vector method
    %   u   =  horizontal velocity AT TOP OF EACH LAYER
    %   v   =  vertical velocity AT TOP OF EACH LAYER
    %   sigzx = horizontal stress AT TOP OF EACH LAYER
    %   sigzz = vertical stress AT TOP OF EACH LAYER
    %   {  Normalized compliance = -k*v/(omega*sigzz)  }
    % [v u sigzz sigzx] = raydep(p,omega,d,ro,vp2,vs2)
    %   p     = slowness (s/m) of surface wave
    %   omega = angular frequency (radians/sec) of surface wave
    %   d     = thicknesses of the model layers (meters?)
    %   rho   = density of the layer (kg/m^3) (= gm/cc * 1000)
    %   vp2   = compressional velocity squared (m/s)^2
    %   vs2   = shear velocity squared (m/s)^2

    % W.C. Crawford 1-5-89
    %   x = stress-displacement vector:
    %		( vertical velocity,  horiz velocity, vert stress, horiz stress )
    %	y = minor vector matrix
    """

    mu = ro * vs2
    n = len(d)
    ist = n - 1
    psq = P * P
    r2 = 2 * mu[ist] * P
    # % R and S are the "Wavenumbers" of compress and shear waves in botlayer
    # % RoW and SoW are divided by ang freq
    RoW = np.sqrt(psq - 1 / vp2[ist])
    SoW = np.sqrt(psq - 1 / vs2[ist])
    ym = np.zeros((ist + 1, 5))
    i = ist
    y = np.zeros((5,))
    x = np.zeros((i + 1, 4))
    y[3 - 1] = RoW
    y[4 - 1] = -SoW
    y[1 - 1] = (RoW * SoW - psq) / ro[i]
    y[2 - 1] = r2 * y[1 - 1] + P
    y[5 - 1] = ro[i] - r2 * (P + y[2 - 1])
    ym[i, :] = y
    # %*****PROPAGATE UP LAYERS*********
    while i > 0:
        i = i - 1
        ha = psq - 1 / vp2[i]
        ca, sa = argdtray(om * d[i], ha)
        hb = psq - 1 / vs2[i]
        cb, sb = argdtray(om * d[i], hb)
        hbs = hb * sb
        has = ha * sa
        r1 = 1 / ro[i]
        r2 = 2 * mu[i] * P
        b1 = r2 * y[1 - 1] - y[2 - 1]
        g3 = (y[5 - 1] + r2 * (y[2 - 1] - b1)) * r1
        g1 = b1 + P * g3
        g2 = ro[i] * y[1 - 1] - P * (g1 + b1)
        e1 = cb * g2 - hbs * y[3 - 1]
        e2 = -sb * g2 + cb * y[3 - 1]
        e3 = cb * y[4 - 1] + hbs * g3
        e4 = sb * y[4 - 1] + cb * g3
        y[3 - 1] = ca * e2 - has * e4
        y[4 - 1] = sa * e1 + ca * e3
        g3 = ca * e4 - sa * e2
        b1 = g1 - P * g3
        y[1 - 1] = (ca * e1 + has * e3 + P * (g1 + b1)) * r1
        y[2 - 1] = r2 * y[1 - 1] - b1
        y[5 - 1] = ro[i] * g3 - r2 * (y[2 - 1] - b1)
        ym[i, :] = y

    y[5 - 1] / np.sqrt(y[1 - 1] * y[1 - 1] + y[2 - 1] * y[2 - 1])
    ynorm = 1 / y[3 - 1]
    y[1 - 1 : 4] = np.array([0, -ynorm, 0, 0])
    # %*****PROPAGATE BACK DOWN LAYERS*********
    while i <= ist:
        x[i, 1 - 1] = (
            -ym[i, 2 - 1] * y[1 - 1] - ym[i, 3 - 1] * y[2 - 1] + ym[i, 1 - 1] * y[4 - 1]
        )
        x[i, 2 - 1] = (
            -ym[i, 4 - 1] * y[1 - 1] + ym[i, 2 - 1] * y[2 - 1] - ym[i, 1 - 1] * y[3 - 1]
        )
        x[i, 3 - 1] = (
            -ym[i, 5 - 1] * y[2 - 1] - ym[i, 2 - 1] * y[3 - 1] - ym[i, 4 - 1] * y[4 - 1]
        )
        x[i, 4 - 1] = (
            ym[i, 5 - 1] * y[1 - 1] - ym[i, 3 - 1] * y[3 - 1] + ym[i, 2 - 1] * y[4 - 1]
        )
        ls = i
        if i >= 2 - 1:
            sum = np.hypot(x[i, 0], x[i, 1])
            1 / vs2[i]
            if sum < 1e-4:
                break

        ha = psq - 1 / vp2[i]
        ca, sa = argdtray(om * d[i], ha)
        hb = psq - 1 / vs2[i]
        cb, sb = argdtray(om * d[i], hb)
        hbs = hb * sb
        has = ha * sa
        r2 = 2 * P * mu[i]
        e2 = r2 * y[2 - 1] - y[3 - 1]
        e3 = ro[i] * y[2 - 1] - P * e2
        e4 = r2 * y[1 - 1] - y[4 - 1]
        e1 = ro[i] * y[1 - 1] - P * e4
        e6 = ca * e2 - sa * e1
        e8 = cb * e4 - sb * e3
        y[1 - 1] = (ca * e1 - has * e2 + P * e8) / ro[i]
        y[2 - 1] = (cb * e3 - hbs * e4 + P * e6) / ro[i]
        y[3 - 1] = r2 * y[2 - 1] - e6
        y[4 - 1] = r2 * y[1 - 1] - e8
        i = i + 1
    #
    # if x(1,3) == 0
    #  error('vertical surface stress = 0 in DETRAY');
    # end
    ist = ls
    v = x[:, 1 - 1]
    u = x[:, 2 - 1]
    zz = x[:, 3 - 1]
    zx = x[:, 4 - 1]
    return v, u, zz, zx


def normalized_compliance(depth, freq, model):
    """calculate normalized compliance for a model and water depth

    model=[thick(m) rho(kg/m^3) vp(m/s) vs(m/s)]; last row is half-space
    """

    model = _validate_model(model)
    depth = positive(depth, "water depth")
    freq = np.asarray(freq, float)
    if (
        freq.ndim != 1
        or not len(freq)
        or np.any(~np.isfinite(freq))
        or np.any(freq <= 0)
    ):
        raise ValueError("frequencies must be a nonempty positive finite vector")
    thick = model[:, 0]
    rho = model[:, 1]
    vpsq = model[:, 2] * model[:, 2]
    vssq = model[:, 3] * model[:, 3]
    omega = 2 * np.pi * freq
    k = gravd(omega, depth)
    p = k / omega
    ncomp = np.zeros((len(p)))

    for i in np.arange((len(p))):
        v, u, sigzz, sigzx = raydep(p[i], omega[i], thick, rho, vpsq, vssq)
        ncomp[i] = -k[i] * v[1 - 1] / (omega[i] * sigzz[1 - 1])

    #   u   =  horizontal velocity AT TOP OF EACH LAYER
    #   v   =  vertical velocity AT TOP OF EACH LAYER
    #   sigzx = horizontal stress AT TOP OF EACH LAYER
    #   sigzz = vertical stress AT TOP OF EACH LAYER

    if not np.all(np.isfinite(ncomp)):
        raise ValueError("model is outside the finite real-valued propagator domain")
    return ncomp


def resample_model(model, dz, max_depth=None):
    """Split a layered model into slices of at most ``dz`` (m) down to ``max_depth``.

    Slice edges include every model interface, so the resampled model is the
    same Earth and has the same compliance. The half-space is sliced as well
    when ``max_depth`` exceeds the finite thickness, and a half-space with its
    properties is kept below. Returns the sliced model and slice mid-depths (m).
    """
    model = _validate_model(model)
    dz = positive(dz, "dz")
    total = model[:-1, 0].sum()
    max_depth = total if max_depth is None else positive(max_depth, "max_depth")
    interfaces = np.cumsum(model[:-1, 0])
    grid = np.arange(0, max_depth + dz / 2, dz)
    edges = np.unique(
        np.concatenate(
            [grid[grid < max_depth], interfaces[interfaces < max_depth], [max_depth]]
        )
    )
    edges = edges[np.concatenate([[True], np.diff(edges) > 1e-9 * max(dz, 1.0)])]
    mid = 0.5 * (edges[:-1] + edges[1:])
    rows = model[np.searchsorted(np.cumsum(model[:-1, 0]), mid, side="right")].copy()
    rows[:, 0] = np.diff(edges)
    return np.vstack([rows, model[-1]]), mid


def compliance_kernels(depth, freq, model, dz=50.0, max_depth=None, rel=0.01):
    """Numerical Vs sensitivity kernels of normalized compliance.

    K[i, j] = d ln(eta(f_i)) / d ln(Vs) per km of depth, for slice j of
    thickness ``dz`` (Vp and density held fixed, as in the inversion). Summing
    ``K * dz / 1000`` over depth gives the response to a uniform relative Vs
    change. Computed by finite differences, the approach of Doran & Laske
    (2020); depth sensitivity peaks between about lambda/2pi and lambda/4 of
    the water wavelength (Crawford 2004). Kernels are mostly negative (a
    stiffer slice lowers compliance). With Vp fixed, raising Vs also lowers
    lambda, so thin shallow slices can have small positive kernels at long
    periods; the Vs sensitivity of such layers is then weak.

    Returns ``(z_mid_m, K)`` with ``K`` of shape ``(len(freq), n_slices)``.
    """
    rel = positive(rel, "rel")
    fine, mid = resample_model(model, dz, max_depth)
    base = np.log(normalized_compliance(depth, freq, fine))
    kernels = np.empty((len(base), len(mid)))
    for j in range(len(mid)):
        trial = fine.copy()
        trial[j, 3] *= 1 + rel
        if trial[j, 3] >= trial[j, 2]:
            raise ValueError("perturbed Vs reaches Vp; reduce rel")
        kernels[:, j] = (np.log(normalized_compliance(depth, freq, trial)) - base) / (
            rel * fine[j, 0] / 1000
        )
    return mid, kernels
