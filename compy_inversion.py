"""Bayesian shear-velocity inversion of normalized seafloor compliance.

The workflow has four parts, each usable on its own:

1. ``ComplianceData`` holds measured normalized compliance and a full data
   covariance. Its terms are the random (Bendat & Piersol) error, the
   month-to-month scatter of the estimate, a pressure-gain error common to
   all frequencies (Aminian et al. 2025, section 4.1.5), and a relative floor
   for forward-model error.
2. ``CrustalModel`` is the parameterization of Zha & Webb (2016, section 3.1).
   A sediment layer of variable thickness has a prescribed Vs(z) profile.
   Below it, crustal layers 100 m thick double in thickness down to a Moho
   depth taken from CRUST1.0 (``crust1_column``). A mantle half-space closes
   the model. Vp and density follow Vs.
3. ``sample`` runs independent chains of a delayed-rejection adaptive
   Metropolis sampler (DRAM; Haario et al. 2001, 2006). The smoothing weight
   of the crustal roughness prior is itself sampled (hierarchical Bayes;
   Malinverno & Briggs 2004). Adaptation stops at the end of burn-in, so the
   retained draws come from a fixed Metropolis kernel.
4. ``Posterior`` provides Vs-depth quantiles, convergence diagnostics
   (Vehtari et al. 2021) and the data fit. Comparing a posterior with a
   prior-only run (``prior_only=True``) shows which depths the data resolve.

The forward model is the minor-vector P-SV propagator of ``inv_compy.raydep``
(Crawford et al. 1991), evaluated at the seafloor only. Numba compiles it
when installed; otherwise the same code runs as plain Python.

All quantities are SI: m, m/s, kg/m**3, Pa, Hz.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from compy_numerics import FREE_AIR_GRADIENT, GRAVITATIONAL_CONSTANT, GRAVITY
from compy_numerics import integer, positive, wavenumber

try:
    from numba import njit

    _jit = njit(cache=True)
except ImportError:  # Numba is optional; the kernel is plain Python.

    def _jit(function):
        return function


# --------------------------------------------------------------------------
# Forward model
# --------------------------------------------------------------------------


@_jit
def _cos_sin(wd, h):
    """Layer propagator terms of ``inv_compy.argdtray``."""
    hh = math.sqrt(abs(h))
    th = wd * hh
    if th >= 1.5e-14:
        if h <= 0.0:  # propagating wave
            return math.cos(th), -math.sin(th) / hh
        if th > 700.0:
            return math.inf, -math.inf
        e = math.exp(th)  # evanescent wave
        return 0.5 * (e + 1.0 / e), -0.5 * (e - 1.0 / e) / hh
    return 1.0, -wd


@_jit
def _eta_kernel(k, omega, thick, rho, vp, vs):
    """Normalized compliance -k v / (omega sigma_zz) at the top of the stack.

    The upward pass of ``inv_compy.raydep``. With the surface vector
    (0, -1/y3, 0, 0), the downward pass gives v = 1, sigma_zx = 0 and
    sigma_zz = y5 / y3 at the seafloor, so no downward pass is needed.
    """
    out = np.empty(k.shape[0])
    b = thick.shape[0] - 1
    for j in range(k.shape[0]):
        om = omega[j]
        p = k[j] / om
        psq = p * p
        ra = psq - 1.0 / (vp[b] * vp[b])
        rb = psq - 1.0 / (vs[b] * vs[b])
        if ra <= 0.0 or rb <= 0.0:  # half-space not evanescent
            out[j] = np.nan
            continue
        row, sow = math.sqrt(ra), math.sqrt(rb)
        r2 = 2.0 * rho[b] * vs[b] * vs[b] * p
        y1 = (row * sow - psq) / rho[b]
        y2 = r2 * y1 + p
        y3 = row
        y4 = -sow
        y5 = rho[b] - r2 * (p + y2)
        for i in range(b - 1, -1, -1):
            ha = psq - 1.0 / (vp[i] * vp[i])
            ca, sa = _cos_sin(om * thick[i], ha)
            hb = psq - 1.0 / (vs[i] * vs[i])
            cb, sb = _cos_sin(om * thick[i], hb)
            hbs = hb * sb
            has = ha * sa
            r1 = 1.0 / rho[i]
            r2 = 2.0 * rho[i] * vs[i] * vs[i] * p
            b1 = r2 * y1 - y2
            g3 = (y5 + r2 * (y2 - b1)) * r1
            g1 = b1 + p * g3
            g2 = rho[i] * y1 - p * (g1 + b1)
            e1 = cb * g2 - hbs * y3
            e2 = -sb * g2 + cb * y3
            e3 = cb * y4 + hbs * g3
            e4 = sb * y4 + cb * g3
            y3 = ca * e2 - has * e4
            y4 = sa * e1 + ca * e3
            g3 = ca * e4 - sa * e2
            b1 = g1 - p * g3
            y1 = (ca * e1 + has * e3 + p * (g1 + b1)) * r1
            y2 = r2 * y1 - b1
            y5 = rho[i] * g3 - r2 * (y2 - b1)
        out[j] = -k[j] * y3 / (om * y5)
    return out


class Forward:
    """Normalized compliance of layered models at fixed depth and frequencies.

    Args:
        depth (float): Water depth [m].
        freqs (array): Positive frequencies [Hz].

    Calling the instance with a model ``[thickness, rho, vp, vs]`` (last row
    the half-space, its thickness ignored) returns compliance [1/Pa], or
    raises ``ValueError`` if the propagator is not finite.
    """

    def __init__(self, depth, freqs):
        self.depth = positive(depth, "water depth")
        self.freqs = np.asarray(freqs, float)
        if self.freqs.ndim != 1 or not self.freqs.size or np.any(~(self.freqs > 0)):
            raise ValueError("frequencies must be a nonempty vector of positive values")
        self.omega = 2 * np.pi * self.freqs
        self.k = np.asarray(wavenumber(self.omega, self.depth), float)

    def __call__(self, model):
        m = np.ascontiguousarray(model, dtype=float)
        if m.ndim != 2 or m.shape[1] != 4 or len(m) < 1:
            raise ValueError("model must have columns [thickness, rho, vp, vs]")
        eta = _eta_kernel(self.k, self.omega, m[:, 0].copy(), m[:, 1].copy(),
                          m[:, 2].copy(), m[:, 3].copy())  # fmt: skip
        if not np.all(np.isfinite(eta)):
            raise ValueError("model is outside the finite real-valued propagator domain")
        return eta


def forward(depth, freqs, model):
    """Normalized compliance [1/Pa] of ``model`` (see ``Forward``)."""
    return Forward(depth, freqs)(model)


def attraction_compliance(freqs, depth):
    """Part of gravity-corrected compliance that does not scale with gain.

    ``compy_processing.gravity_corrected_compliance`` adds the wave
    attraction k (a_w/p) / (omega**2 + FREE_AIR_GRADIENT), which does not
    depend on the pressure calibration. A relative gain error d therefore
    changes compliance by -d (eta - this term).
    """
    f = np.asarray(freqs, float)
    k = wavenumber(2 * np.pi * f, depth)
    kh = k * depth
    pa = 2 * np.pi * (GRAVITATIONAL_CONSTANT / GRAVITY) * (1 + np.exp(-2 * kh)) / 2
    return k * pa / ((2 * np.pi * f) ** 2 + FREE_AIR_GRADIENT)


# --------------------------------------------------------------------------
# Rock-physics relations
# --------------------------------------------------------------------------


def hamilton_vs(z):
    """Vs [m/s] of terrigenous silt-clay at depth ``z`` [m] below the seafloor.

    Hamilton (1976, 1979): 116 + 4.65 z (0-36 m), 237 + 1.28 z (36-120 m),
    322 + 0.58 z (120-650 m). The last segment is extended below 650 m, as by
    Zha & Webb (2016, Figure 8).
    """
    z = np.asarray(z, float)
    return np.where(z < 36, 116 + 4.65 * z, np.where(z < 120, 237 + 1.28 * z, 322 + 0.58 * z))


def hamilton_vp(z):
    """Vp [m/s] of terrigenous sediment at depth ``z`` [m] (Hamilton 1979).

    Vp = 1.511 + 1.304 D - 0.741 D**2 + 0.257 D**3 km/s with D in km.
    """
    d = np.asarray(z, float) / 1000.0
    return 1000.0 * (1.511 + 1.304 * d - 0.741 * d**2 + 0.257 * d**3)


def brocher_vp(vs):
    """Vp [m/s] from Vs [m/s] (Brocher 2005, eq. 9; Vs up to 4.5 km/s)."""
    km = np.asarray(vs, float) / 1000.0
    return 1000.0 * (0.9409 + 2.0947 * km - 0.8206 * km**2 + 0.2683 * km**3 - 0.0251 * km**4)


def christensen_shaw_density(vp):
    """Density [kg/m**3] from Vp [m/s]: 1.85 + 0.165 Vp (g/cc, km/s).

    Christensen & Shaw (1970) for basalts, as used by Zha & Webb (2016).
    """
    return 1850.0 + 0.165 * np.asarray(vp, float)


# --------------------------------------------------------------------------
# CRUST1.0
# --------------------------------------------------------------------------

CRUST1_LAYERS = ("water", "ice", "upper_sediments", "middle_sediments", "lower_sediments",
                 "upper_crust", "middle_crust", "lower_crust", "mantle")  # fmt: skip


def crust1_column(lat, lon, directory):
    """CRUST1.0 column of the 1x1 degree cell containing (lat, lon).

    Args:
        lat, lon (float): Degrees; longitude in [-180, 360).
        directory (path): Folder with ``crust1.bnds``, ``crust1.vp``,
            ``crust1.vs`` and ``crust1.rho`` (Laske et al. 2013,
            https://igppweb.ucsd.edu/~gabi/crust1.html).

    Returns:
        dict: ``cell`` centre, ``water_depth`` [m], ``layers`` (name,
        thickness [m], vp, vs [m/s], rho [kg/m**3]) below the water for
        layers of nonzero thickness, the ``mantle`` properties, and
        ``moho_depth`` and ``sediment_thickness`` below the seafloor [m].
    """
    directory = Path(directory)
    lon = (float(lon) + 180.0) % 360.0 - 180.0
    row = min(max(int(math.floor(90.0 - float(lat))), 0), 179) * 360
    row += min(max(int(math.floor(lon + 180.0)), 0), 359)
    values = {}
    for key in ("bnds", "vp", "vs", "rho"):
        with open(directory / f"crust1.{key}") as handle:
            for i, line in enumerate(handle):
                if i == row:
                    values[key] = np.array(line.split(), float)
                    break
    thickness = -np.diff(values["bnds"]) * 1000.0
    layers = [
        (CRUST1_LAYERS[i], float(thickness[i]), 1000 * values["vp"][i],
         1000 * values["vs"][i], 1000 * values["rho"][i])
        for i in range(2, 8) if thickness[i] > 0
    ]  # fmt: skip
    return dict(
        cell=(89.5 - row // 360, -179.5 + row % 360),
        water_depth=float(thickness[0]),
        layers=layers,
        mantle=(1000 * values["vp"][8], 1000 * values["vs"][8], 1000 * values["rho"][8]),
        moho_depth=float(sum(layer[1] for layer in layers)),
        sediment_thickness=float(sum(t for name, t, *_ in layers if "sediment" in name)),
    )


def crust1_model(column):
    """Layered ``[thickness, rho, vp, vs]`` model of a ``crust1_column``."""
    rows = [[t, rho, vp, vs] for _, t, vp, vs, rho in column["layers"]]
    vp, vs, rho = column["mantle"]
    return np.array(rows + [[0.0, rho, vp, vs]])


# --------------------------------------------------------------------------
# Data and covariance
# --------------------------------------------------------------------------


@dataclass
class ComplianceData:
    """Measured normalized compliance with its full covariance.

    Build instances with ``ComplianceData.build`` or ``from_npz``.

    Attributes:
        freqs (array): Frequencies [Hz].
        eta (array): Normalized compliance [1/Pa].
        depth (float): Water depth [m].
        covariance (array): Data covariance [1/Pa**2].
        components (dict): The named terms that sum to ``covariance``.
        station (str): Label.
    """

    freqs: np.ndarray
    eta: np.ndarray
    depth: float
    covariance: np.ndarray
    components: dict = field(default_factory=dict)
    station: str = ""

    def __post_init__(self):
        self.freqs = np.asarray(self.freqs, float)
        self.eta = np.asarray(self.eta, float)
        self.covariance = np.asarray(self.covariance, float)
        n = len(self.eta)
        if self.freqs.shape != (n,) or n == 0 or np.any(~(self.freqs > 0)):
            raise ValueError("freqs and eta must be equal-length vectors of positive frequencies")
        if np.any(~np.isfinite(self.eta)) or np.any(self.eta <= 0):
            raise ValueError("eta must be finite and positive")
        if self.covariance.shape != (n, n) or not np.allclose(self.covariance, self.covariance.T):
            raise ValueError("covariance must be a symmetric N x N matrix")
        positive(self.depth, "water depth")
        try:
            chol = np.linalg.cholesky(self.covariance)
        except np.linalg.LinAlgError as exc:
            raise ValueError("covariance must be positive definite") from exc
        self._whitener = np.linalg.inv(chol)

    @classmethod
    def build(cls, freqs, eta, depth, *, random=None, subsets=None, shrinkage=0.5,
              gain_rel=0.05, gain_vector=None, floor_rel=0.01, station=""):  # fmt: skip
        """Assemble the covariance from its physical sources.

        Args:
            freqs, eta, depth: As in the class.
            random (array): Random error per frequency [1/Pa], e.g. the
                Bendat & Piersol error of the stacked estimate.
            subsets (array): Estimates from independent subsets of the
                record, shape (n_subsets, n_freqs), e.g. monthly stacks. Their
                covariance divided by n_subsets is the error of the full
                stack under non-stationary conditions; it is shrunk towards
                its diagonal by ``shrinkage`` because few subsets are
                available.
            gain_rel (float): Relative standard error of the pressure gain,
                0.05 in Aminian et al. (2025, Table 1). It is fully correlated
                across frequency.
            gain_vector (array): Change of compliance per unit relative gain
                change; defaults to eta minus ``attraction_compliance``.
            floor_rel (float): Independent relative error per frequency for
                what the 1-D forward model and the estimator neglect.
            station (str): Label.
        """
        freqs, eta = np.asarray(freqs, float), np.asarray(eta, float)
        n = len(eta)
        terms = {}
        if random is not None:
            terms["random"] = np.diag(np.broadcast_to(np.asarray(random, float), (n,)) ** 2)
        if subsets is not None:
            s = np.atleast_2d(np.asarray(subsets, float))
            if s.shape[1] != n or s.shape[0] < 3:
                raise ValueError("subsets must have shape (n_subsets >= 3, n_freqs)")
            if not 0 <= shrinkage <= 1:
                raise ValueError("shrinkage must lie in [0, 1]")
            c = np.cov(s, rowvar=False) / s.shape[0]
            terms["subsets"] = (1 - shrinkage) * c + shrinkage * np.diag(np.diag(c))
        if gain_rel:
            v = eta - attraction_compliance(freqs, depth) if gain_vector is None else gain_vector
            v = np.asarray(v, float)
            terms["gain"] = gain_rel**2 * np.outer(v, v)
        if floor_rel:
            terms["floor"] = np.diag((floor_rel * eta) ** 2)
        if not terms:
            raise ValueError("at least one covariance term is required")
        return cls(freqs, eta, depth, sum(terms.values()), terms, station)

    @classmethod
    def from_npz(cls, path, *, band=None, subsets=None, **kwargs):
        """Load a ComPy 2.2 compliance file (``stacked_compliance`` output).

        Uses the keys ``freqs``, ``compliance``, ``uncertainty``, ``band``
        and ``depth``. Frequencies outside ``band`` (default: the stored band)
        are dropped; ``subsets`` must already be restricted to that band.
        Other keywords go to ``build``.
        """
        d = np.load(path)
        lo, hi = d["band"] if band is None else band
        keep = (d["freqs"] >= lo - 1e-12) & (d["freqs"] <= hi + 1e-12)
        return cls.build(d["freqs"][keep], d["compliance"][keep], float(d["depth"]),
                         random=d["uncertainty"][keep], subsets=subsets, **kwargs)  # fmt: skip

    @property
    def sigma(self):
        """Marginal standard deviation per frequency [1/Pa]."""
        return np.sqrt(np.diag(self.covariance))

    def chi2(self, prediction):
        """Generalized chi-square of a prediction (Mahalanobis distance**2)."""
        r = self._whitener @ (self.eta - np.asarray(prediction, float))
        return float(r @ r)


# --------------------------------------------------------------------------
# Parameterization
# --------------------------------------------------------------------------


@dataclass
class CrustalModel:
    """Sediment, doubling crustal layers and a mantle half-space.

    The free parameters, in order, are the sediment thickness (if
    ``sediment_max`` > 0), a sediment Vs scale (if its bounds differ), one Vs
    per crustal layer, the mantle Vs and a Moho shift (if its bounds differ).

    Attributes:
        moho_depth (float): Moho depth below the seafloor [m], e.g.
            ``crust1_column(...)["moho_depth"]``.
        first_layer (float): Thickness of the top crustal layer [m].
        growth (float): Ratio between successive crustal layers.
        min_last (float): Smallest thickness left for the deepest crustal
            layer, which absorbs sediment and Moho changes [m].
        sediment_max (float): Upper bound of the sediment thickness [m];
            0 removes the sediment.
        sediment_scale (tuple): Bounds on a multiplier of the sediment
            profile; (1, 1) fixes the profile, as in Zha & Webb (2016).
        sediment_vs, sediment_vp (callable): Profiles of depth below the
            seafloor [m -> m/s]; Hamilton (1976, 1979) by default.
        sediment_density (float): Constant sediment density [kg/m**3].
        sediment_sublayers (int): Equal sublayers representing the sediment.
        vs_crust, vs_mantle (tuple): Uniform Vs bounds [m/s].
        moho_shift (tuple): Bounds on a Moho depth change [m]; (-1000, 1000)
            reproduces the freedom allowed by Aminian et al. (2025).
        vp_vs (float or None): Crust and mantle Vp/Vs; None uses
            ``brocher_vp``. Density follows ``christensen_shaw_density``.
        roughness_order (int): Order of the Vs differences between crustal
            layers penalized by the smoothing prior.
        sigma_rough (tuple): Bounds [m/s] of the roughness scale; the
            smoothing weight alpha = 1/sigma**2 has a Jeffreys prior there.
    """

    moho_depth: float
    first_layer: float = 100.0
    growth: float = 2.0
    min_last: float = 500.0
    sediment_max: float = 1000.0
    sediment_scale: tuple = (1.0, 1.0)
    sediment_vs: object = hamilton_vs
    sediment_vp: object = hamilton_vp
    sediment_density: float = 2000.0
    sediment_sublayers: int = 8
    vs_crust: tuple = (500.0, 4700.0)
    vs_mantle: tuple = (3800.0, 4700.0)
    moho_shift: tuple = (0.0, 0.0)
    vp_vs: float | None = 1.85
    roughness_order: int = 2
    sigma_rough: tuple = (20.0, 3000.0)

    def __post_init__(self):
        positive(self.moho_depth, "moho_depth")
        positive(self.first_layer, "first_layer")
        integer(self.sediment_sublayers, "sediment_sublayers")
        if self.growth < 1 or self.sediment_max < 0:
            raise ValueError("growth must be >= 1 and sediment_max >= 0")
        for name in ("sediment_scale", "vs_crust", "vs_mantle", "moho_shift", "sigma_rough"):
            lo, hi = getattr(self, name)
            if not lo <= hi:
                raise ValueError(f"{name} must be (lower, upper)")
        shallowest = self.moho_depth + self.moho_shift[0] - self.sediment_max - self.min_last
        fixed, edge = [], 0.0
        while edge + self.first_layer * self.growth ** len(fixed) <= shallowest:
            fixed.append(self.first_layer * self.growth ** len(fixed))
            edge += fixed[-1]
        if shallowest <= 0:
            raise ValueError("Moho too shallow for the sediment and min_last bounds")
        self.fixed = np.array(fixed)
        self.n_crust = len(fixed) + 1
        names, bounds = [], []
        if self.sediment_max > 0:
            names.append("sediment_thickness")
            bounds.append((0.0, float(self.sediment_max)))
        if self.sediment_max > 0 and self.sediment_scale[0] < self.sediment_scale[1]:
            names.append("sediment_vs_scale")
            bounds.append(tuple(map(float, self.sediment_scale)))
        names += [f"vs_crust_{i + 1}" for i in range(self.n_crust)] + ["vs_mantle"]
        bounds += [tuple(map(float, self.vs_crust))] * self.n_crust
        bounds.append(tuple(map(float, self.vs_mantle)))
        if self.moho_shift[0] < self.moho_shift[1]:
            names.append("moho_shift")
            bounds.append(tuple(map(float, self.moho_shift)))
        self.names = names
        self.bounds = np.array(bounds, float)
        self._index = {name: i for i, name in enumerate(names)}
        if self.roughness_order >= self.n_crust:
            raise ValueError("roughness_order must be smaller than the number of crustal layers")

    @property
    def n_rough(self):
        """Number of differences in the roughness prior."""
        return self.n_crust - self.roughness_order

    def _get(self, theta, name, default):
        i = self._index.get(name)
        return default if i is None else theta[i]

    def crust_vs(self, theta):
        i = self._index["vs_crust_1"]
        return np.asarray(theta[i:i + self.n_crust], float)

    def roughness(self, theta):
        """Sum of squared ``roughness_order`` differences of crustal Vs [m**2/s**2]."""
        return float(np.sum(np.diff(self.crust_vs(theta), n=self.roughness_order) ** 2))

    def elastic(self, vs):
        """(Vp, rho) of crust or mantle rows for Vs [m/s]."""
        vp = brocher_vp(vs) if self.vp_vs is None else self.vp_vs * np.asarray(vs, float)
        return vp, christensen_shaw_density(vp)

    def interfaces(self, theta):
        """Sediment base and crustal layer bottoms below the seafloor [m]."""
        h = self._get(theta, "sediment_thickness", 0.0)
        moho = self.moho_depth + self._get(theta, "moho_shift", 0.0)
        return h, np.append(h + np.cumsum(self.fixed), moho)

    def layers(self, theta):
        """Layered model ``[thickness, rho, vp, vs]`` (last row half-space)."""
        theta = np.asarray(theta, float)
        h, bottoms = self.interfaces(theta)
        rows = []
        if h > 1.0:
            edges = np.linspace(0.0, h, self.sediment_sublayers + 1)
            mid = 0.5 * (edges[:-1] + edges[1:])
            vs = self._get(theta, "sediment_vs_scale", 1.0) * self.sediment_vs(mid)
            vp = np.maximum(self.sediment_vp(mid), brocher_vp(vs))
            for t, a, b in zip(np.diff(edges), vp, vs):
                rows.append((t, self.sediment_density, a, b))
        vs = self.crust_vs(theta)
        vp, rho = self.elastic(vs)
        for t, a, b, r in zip(np.diff(np.append(h, bottoms)), vp, vs, rho):
            rows.append((t, r, a, b))
        vsm = theta[self._index["vs_mantle"]]
        vpm, rhom = self.elastic(vsm)
        rows.append((0.0, float(rhom), float(vpm), float(vsm)))
        return np.array(rows, float)

    def vs_at(self, theta, z):
        """Vs [m/s] of the model at depths ``z`` [m] below the seafloor."""
        model = self.layers(theta)
        tops = np.concatenate([[0.0], np.cumsum(model[:-1, 0])])
        return model[np.searchsorted(tops, np.asarray(z, float), side="right") - 1, 3]

    def in_bounds(self, theta):
        return bool(np.all(theta >= self.bounds[:, 0]) and np.all(theta <= self.bounds[:, 1]))

    def random(self, rng):
        """A draw from the uniform part of the prior."""
        return rng.uniform(self.bounds[:, 0], self.bounds[:, 1])

    def settings(self):
        out = {k: v for k, v in asdict(self).items() if not callable(v)}
        out.update(sediment_vs=getattr(self.sediment_vs, "__name__", "custom"),
                   sediment_vp=getattr(self.sediment_vp, "__name__", "custom"),
                   names=self.names, fixed_layers=self.fixed.tolist())  # fmt: skip
        return out


# --------------------------------------------------------------------------
# Sampler
# --------------------------------------------------------------------------


def truncated_gamma(rng, shape, rate, lower, upper):
    """Draw from Gamma(shape, rate) restricted to [lower, upper]."""
    if rate <= 0:  # density proportional to x**(shape - 1)
        u = rng.random()
        return (lower**shape + u * (upper**shape - lower**shape)) ** (1 / shape)
    for _ in range(64):
        x = rng.gamma(shape, 1 / rate)
        if lower <= x <= upper:
            return float(x)
    from scipy.stats import gamma

    dist = gamma(shape, scale=1 / rate)
    lo, hi = dist.cdf(lower), dist.cdf(upper)
    if hi - lo < 1e-300:
        return float(lower if dist.mean() < lower else upper)
    return float(np.clip(dist.ppf(lo + rng.random() * (hi - lo)), lower, upper))


class _Target:
    """Log-likelihood, roughness and prediction of one parameter vector."""

    def __init__(self, data, model, prior_only):
        self.data, self.model, self.prior_only = data, model, prior_only
        self.forward = Forward(data.depth, data.freqs)

    def __call__(self, theta):
        if not self.model.in_bounds(theta):
            return None
        try:
            pred = self.forward(self.model.layers(theta))
        except ValueError:
            return None
        loglik = 0.0 if self.prior_only else -0.5 * self.data.chi2(pred)
        return loglik, self.model.roughness(theta), pred


def dram_chain(target, theta0, bounds, *, n_iter, n_burn, seed, thin=1, n_rough=0,
               alpha_bounds=(1.0, 1.0), alpha0=None, adapt_start=500, adapt_every=50,
               target_acceptance=0.234, delayed_rejection=True, dr_scale=0.1,
               initial_scale=0.02):  # fmt: skip
    """One chain of delayed-rejection adaptive Metropolis with Gibbs smoothing.

    The target density is exp(loglik(theta) - alpha R(theta) / 2) inside
    ``bounds``. With ``n_rough`` > 0, alpha has the prior
    alpha**(n_rough/2 - 1) on ``alpha_bounds`` (a Jeffreys prior times the
    normalization of the Gaussian roughness prior) and is drawn from its
    Gamma conditional after every Metropolis step.

    The first-stage proposal is Gaussian with covariance lambda**2 (2.38**2/d)
    C, where C is the covariance of the second half of the chain so far
    (Haario et al. 2001; discarding the first half forgets the transient of an
    over-dispersed start) and lambda is tuned by Robbins-Monro towards
    ``target_acceptance`` (Andrieu & Thoms 2008). Both are updated every
    ``adapt_every`` iterations between ``adapt_start`` and ``n_burn``, then
    frozen.
    A rejected first stage is followed by a second, ``dr_scale`` times
    smaller proposal accepted with the delayed-rejection probability of
    Tierney & Mira (1999), which keeps detailed balance (Haario et al. 2006).

    Args:
        target (callable): theta -> (loglik, roughness, prediction) or None
            outside the support.
        theta0 (array): Starting point inside the support.
        bounds (array): (d, 2) lower and upper bounds, used for the initial
            proposal scale.

    Returns:
        dict of retained draws after ``n_burn``, every ``thin`` iterations:
        ``theta``, ``alpha``, ``loglik``, ``roughness``, ``pred``; and the
        acceptance rates of both stages after burn-in.
    """
    rng = np.random.default_rng(seed)
    x = np.array(theta0, float)
    d = len(x)
    width = bounds[:, 1] - bounds[:, 0]
    width = np.where(width > 0, width, 1.0)
    cov = np.diag((initial_scale * width) ** 2)
    chol = np.linalg.cholesky(cov)
    chol_inv = np.linalg.inv(chol)
    eps = np.diag((1e-4 * width) ** 2)
    sd = 2.38**2 / d
    log_scale = 0.0
    history = np.empty((max(n_burn, 1), d))
    window_taken = window_tried = 0
    state = target(x)
    if state is None:
        raise ValueError("theta0 lies outside the support of the target")
    loglik, rough, pred = state
    lo_a, hi_a = alpha_bounds
    alpha = float(alpha0 if alpha0 is not None else math.sqrt(lo_a * hi_a))

    def log_post(ll, r):
        return ll - 0.5 * alpha * r

    lp = log_post(loglik, rough)
    keep = []
    tried = np.zeros(2)
    taken = np.zeros(2)
    for it in range(n_iter):
        z = chol @ rng.standard_normal(d)
        y1 = x + z
        s1 = target(y1)
        a1 = 0.0
        if s1 is not None:
            lp1 = log_post(s1[0], s1[1])
            a1 = 1.0 if lp1 >= lp else math.exp(lp1 - lp)
        stage = None
        if rng.random() < a1:
            stage, new = 0, (y1, s1, lp1)
        elif delayed_rejection:
            y2 = x + dr_scale * (chol @ rng.standard_normal(d))
            s2 = target(y2)
            if s2 is not None:
                lp2 = log_post(s2[0], s2[1])
                a21 = 0.0 if s1 is None else min(1.0, math.exp(min(lp1 - lp2, 0.0)))
                if a21 < 1.0:
                    u2 = chol_inv @ (y1 - y2)
                    ux = chol_inv @ (y1 - x)
                    log_a2 = (lp2 - lp - 0.5 * (u2 @ u2) + 0.5 * (ux @ ux)
                              + math.log1p(-a21) - math.log1p(-a1))  # fmt: skip
                    if math.log(max(rng.random(), 1e-300)) < log_a2:
                        stage, new = 1, (y2, s2, lp2)
        if it >= n_burn:
            tried[0] += 1
            tried[1] += stage != 0 and delayed_rejection
            if stage is not None:
                taken[stage] += 1
        if stage is not None:
            x, (loglik, rough, pred), lp = new[0], new[1], new[2]
        if n_rough:
            alpha = truncated_gamma(rng, 0.5 * n_rough, 0.5 * rough, lo_a, hi_a)
            lp = log_post(loglik, rough)
        if it < n_burn:  # adaptive proposal, frozen after burn-in
            history[it] = x
            window_tried += 1
            window_taken += stage == 0
            if it >= adapt_start and (it - adapt_start) % adapt_every == 0:
                rate = window_taken / window_tried
                log_scale += (1 + (it - adapt_start) / adapt_every) ** -0.6 * (rate - target_acceptance) * 2
                window_taken = window_tried = 0
                recent = history[it // 2:it + 1]
                cov = math.exp(2 * log_scale) * sd * (np.cov(recent, rowvar=False) + eps)
                chol = np.linalg.cholesky(cov)
                chol_inv = np.linalg.inv(chol)
        elif (it - n_burn) % thin == 0:
            keep.append((x.copy(), alpha, loglik, rough, pred))
    if not keep:
        raise ValueError("n_iter must exceed n_burn")
    with np.errstate(invalid="ignore", divide="ignore"):
        rates = taken / tried
    return dict(
        theta=np.array([k[0] for k in keep]),
        alpha=np.array([k[1] for k in keep]),
        loglik=np.array([k[2] for k in keep]),
        roughness=np.array([k[3] for k in keep]),
        pred=np.array([k[4] for k in keep]),
        acceptance=float(rates[0]),
        acceptance_dr=float(rates[1]),
        proposal_cov=cov,
    )


def _chain_job(job):
    data, model, prior_only, seed, theta0, kwargs = job
    target = _Target(data, model, prior_only)
    lo, hi = model.sigma_rough
    return dram_chain(target, theta0, model.bounds, seed=seed, n_rough=model.n_rough,
                      alpha_bounds=(1 / hi**2, 1 / lo**2), **kwargs)  # fmt: skip


def sample(data, model, *, n_chains=8, n_iter=60000, n_burn=20000, thin=10, seed=0,
           prior_only=False, processes=None, starts=None, **kwargs):  # fmt: skip
    """Run independent DRAM chains and return a ``Posterior``.

    Chains start from independent draws of the uniform prior (over-dispersed
    starts, so that R-hat is meaningful) unless ``starts`` gives one
    parameter vector per chain. Draws come from ``_Target``: the Gaussian
    likelihood of ``data`` (or none with ``prior_only``) and the uniform and
    roughness priors of ``model``. ``processes`` sets the worker count
    (1 runs serially). Remaining keywords go to ``dram_chain``.
    """
    rng = np.random.default_rng(seed)
    target = _Target(data, model, prior_only)
    if starts is None:
        starts = []
        while len(starts) < n_chains:
            theta = model.random(rng)
            if target(theta) is not None:
                starts.append(theta)
    if len(starts) != n_chains:
        raise ValueError("give one start per chain")
    seeds = rng.integers(2**31, size=n_chains)
    settings = dict(n_iter=n_iter, n_burn=n_burn, thin=thin, **kwargs)
    jobs = [(data, model, prior_only, int(s), np.asarray(t, float), settings)
            for s, t in zip(seeds, starts)]  # fmt: skip
    if processes == 1 or n_chains == 1:
        chains = [_chain_job(j) for j in jobs]
    else:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=processes) as pool:
            chains = list(pool.map(_chain_job, jobs))
    stack = {key: np.stack([c[key] for c in chains])
             for key in ("theta", "alpha", "loglik", "roughness", "pred")}  # fmt: skip
    return Posterior(
        data=data, model=model, prior_only=prior_only, **stack,
        acceptance=np.array([c["acceptance"] for c in chains]),
        acceptance_dr=np.array([c["acceptance_dr"] for c in chains]),
        settings=dict(n_chains=n_chains, seed=seed, **settings),
    )  # fmt: skip


# --------------------------------------------------------------------------
# Posterior
# --------------------------------------------------------------------------


@dataclass
class Posterior:
    """Retained draws of ``sample``: arrays shaped (chains, draws, ...)."""

    data: ComplianceData
    model: CrustalModel
    prior_only: bool
    theta: np.ndarray
    alpha: np.ndarray
    loglik: np.ndarray
    roughness: np.ndarray
    pred: np.ndarray
    acceptance: np.ndarray
    acceptance_dr: np.ndarray
    settings: dict

    @property
    def names(self):
        return self.model.names

    def flat(self, max_draws=None, seed=0):
        """Pooled parameter draws, optionally a random subset."""
        theta = self.theta.reshape(-1, self.theta.shape[-1])
        if max_draws is not None and len(theta) > max_draws:
            theta = theta[np.random.default_rng(seed).choice(len(theta), max_draws, replace=False)]
        return theta

    def chi2(self):
        """Generalized chi-square of every draw (chains, draws)."""
        return np.array([[self.data.chi2(p) for p in chain] for chain in self.pred])

    def diagnostics(self):
        """Rank-normalized split R-hat and bulk ESS per parameter and for log alpha."""
        from compy_diagnostics import rank_rhat_ess

        out = {}
        series = {name: self.theta[:, :, i] for i, name in enumerate(self.names)}
        series["log_alpha"] = np.log(self.alpha)
        for name, x in series.items():
            if np.ptp(x) == 0:
                out[name] = dict(rhat=None, ess_bulk=None)
                continue
            r, ess = rank_rhat_ess(x)
            out[name] = dict(rhat=float(r), ess_bulk=float(ess))
        return out

    def vs_draws(self, z, max_draws=4000, seed=0):
        """Vs [m/s] of pooled draws at depths ``z`` [m], shape (draws, len(z))."""
        return np.array([self.model.vs_at(t, z) for t in self.flat(max_draws, seed)])

    def vs_quantiles(self, z, q=(0.05, 0.16, 0.5, 0.84, 0.95), max_draws=4000):
        """Pointwise Vs quantiles at depths ``z``, shape (len(q), len(z))."""
        return np.quantile(self.vs_draws(z, max_draws), q, axis=0)

    def summary(self):
        """JSON-ready summary: parameter quantiles, diagnostics and fit."""
        flat = self.theta.reshape(-1, self.theta.shape[-1])
        diag = self.diagnostics()
        params = {}
        for i, name in enumerate(self.names):
            q = np.quantile(flat[:, i], [0.05, 0.16, 0.5, 0.84, 0.95])
            params[name] = dict(zip(("p05", "p16", "median", "p84", "p95"), map(float, q)),
                                **diag[name])  # fmt: skip
        sigma = 1 / np.sqrt(self.alpha.ravel())
        chi2 = self.chi2().ravel()
        return dict(
            station=self.data.station,
            prior_only=self.prior_only,
            n_data=len(self.data.eta),
            parameters=params,
            sigma_rough=dict(zip(("p05", "median", "p95"),
                                 map(float, np.quantile(sigma, [0.05, 0.5, 0.95]))),
                             **diag["log_alpha"]),  # fmt: skip
            chi2=dict(zip(("p05", "median", "p95"), map(float, np.quantile(chi2, [0.05, 0.5, 0.95])))),
            max_rhat=max((v["rhat"] for v in diag.values() if v["rhat"] is not None), default=None),
            min_ess=min((v["ess_bulk"] for v in diag.values() if v["ess_bulk"] is not None),
                        default=None),  # fmt: skip
            acceptance=self.acceptance.tolist(),
            acceptance_dr=self.acceptance_dr.tolist(),
            settings=self.settings,
            model=self.model.settings(),
        )

    def save(self, path):
        """Write ``<path>.npz`` (draws and data) and ``<path>.json`` (summary)."""
        path = Path(path)
        np.savez_compressed(path.with_suffix(".npz"), theta=self.theta, alpha=self.alpha,
                            loglik=self.loglik, roughness=self.roughness, pred=self.pred,
                            freqs=self.data.freqs, eta=self.data.eta,
                            covariance=self.data.covariance, depth=self.data.depth,
                            names=np.array(self.names))  # fmt: skip
        path.with_suffix(".json").write_text(json.dumps(self.summary(), indent=1, default=float))


def resolution_ratio(posterior, prior, z, q=(0.16, 0.84)):
    """Posterior over prior width of the Vs interval ``q`` at depths ``z``.

    Values well below 1 mark depths constrained by the data; values near 1
    mark depths where the posterior only repeats the prior.
    """
    post = np.quantile(posterior.vs_draws(z), q, axis=0)
    pri = np.quantile(prior.vs_draws(z), q, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (post[1] - post[0]) / (pri[1] - pri[0])
