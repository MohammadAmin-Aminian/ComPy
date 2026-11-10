"""Automatic Bayesian inversion of seafloor compliance (ComPy 3 default).

One call, ``invert``, takes a compliance curve with its covariance and a
reference Earth model at the station and returns shear velocity against depth
with credible intervals, the depth range the data resolve, and a quality label.
Nothing is tuned by hand.

What differs from ``compy.inversion.dram.CrustalModel`` and why
(docs/inversion_overhaul_2026-11/REPORT.md gives the tests):

* **Proper prior about a reference model.** ln Vs of each crustal layer is
  Gaussian about the reference (CRUST1.0 at the station, or a model given by the
  user), with standard deviation ``vs_sigma`` and correlation
  ``vs_correlation`` between neighbouring layers. Where the data say nothing,
  the posterior returns to the reference with the prior width. The earlier
  flat prior with second-difference smoothing extrapolated the shallow
  gradient into the unresolved deep layers and pushed them to their bounds.
* **Vp and density are uncertain, not fixed.** Compliance measures stiffness.
  Under a half-space it is (1 - nu)/mu, so Vs follows only once Vp/Vs and
  density are given. Vp follows Vs by ``oceanic_vp`` with a sampled factor for
  the upper and for the lower crust, and density follows Vp with a sampled
  factor. Their uncertainty is carried into the Vs intervals.
* **Sediment enters through its P modulus as well as its Vs.** A soft layer
  much thinner than the wavelength adds k h/(rho Vp**2) to compliance whatever
  its Vs. Thickness, a Vs scale and a P-modulus factor are free; density
  follows Vp (Nafe-Drake) instead of a constant 2000 kg/m**3.
* **Errors as measured.** ``compliance_data`` builds the covariance from the
  bootstrap covariance of the stack (``compy.inversion.errors``) plus a tilt
  term for the pressure-response shape and a smooth model-error term. The
  pressure gain is sampled with a Gaussian prior.
* **Resolution is reported.** ``InversionResult.profile`` marks the depths where
  the posterior is narrower than the prior and within reach of the data.

The sampler, the forward model and the convergence diagnostics are those of
``compy.inversion.dram``. All quantities are SI.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from compy.inversion import dram
from compy.inversion.errors import default_covariance
from compy.numerics import positive

# --------------------------------------------------------------------------
# Rock relations
# --------------------------------------------------------------------------


def oceanic_vp(vs):
    """Vp [m/s] expected for Vs [m/s] in igneous oceanic crust.

    Brocher (2005): the mafic line, Vs = 2.88 + 0.52 (Vp - 5.25) km/s for
    5.25 <= Vp <= 7.25 km/s (his eq. 8, for mafic rocks and gabbros), continued
    above with its end ratio Vp/Vs = 1.85; his regression fit (eq. 9) below
    Vs = 1.9 km/s, where porous and cracked rock has a higher Vp/Vs; a straight
    line between the two. The coefficients were taken from a published
    implementation (R package hydrogeo), not from the paper itself. The
    relation is a prior mean: ``SeafloorModel`` samples a factor on it.
    """
    km = np.asarray(vs, float) / 1000.0
    brocher = 0.9409 + 2.0947 * km - 0.8206 * km**2 + 0.2683 * km**3 - 0.0251 * km**4
    v19 = 0.9409 + 2.0947 * 1.9 - 0.8206 * 1.9**2 + 0.2683 * 1.9**3 - 0.0251 * 1.9**4
    bridge = v19 + (5.25 - v19) * (km - 1.9) / (2.88 - 1.9)
    mafic = 5.25 + (km - 2.88) / 0.52
    high = 7.25 * km / 3.92
    return 1000.0 * np.where(km < 1.9, brocher,
                             np.where(km < 2.88, bridge, np.where(km < 3.92, mafic, high)))  # fmt: skip


def nafe_drake_density(vp):
    """Density [kg/m**3] from Vp [m/s]: the Nafe-Drake curve (Brocher 2005, eq. 1).

    For sediments and rocks with 1.5 <= Vp <= 8.5 km/s, mafic rocks excepted.
    Coefficients from a published implementation (R package hydrogeo).
    """
    v = np.clip(np.asarray(vp, float) / 1000.0, 1.5, 8.5)
    return 1000.0 * (
        1.6612 * v - 0.4721 * v**2 + 0.0671 * v**3 - 0.0043 * v**4 + 0.000106 * v**5
    )


def carlson_raskin_density(vp):
    """Density [kg/m**3] of oceanic crust from Vp [m/s] (Carlson & Raskin 1984).

    rho = 3.50 - 3.79/Vp below 6.65 km/s and 3.81 - 5.99/Vp above, as used for
    compliance by Crawford (1994, eq. 6.1). The two pieces differ by 0.7 % at
    6.65 km/s; the relation is kept as published.
    """
    v = np.maximum(np.asarray(vp, float) / 1000.0, 1.0)
    return 1000.0 * np.where(v < 6.65, 3.50 - 3.79 / v, 3.81 - 5.99 / v)


#: Vp [m/s] at which the Carlson & Raskin and the Nafe-Drake curves cross.
DENSITY_SWITCH_VP = 2938.0


def crust_density(vp):
    """Density [kg/m**3] of crustal rock from Vp [m/s].

    Carlson & Raskin (1984) for oceanic crust above 2.94 km/s. Below it that
    curve falls under the Nafe-Drake curve (and under water density near
    1.5 km/s), so the Nafe-Drake value is used there. The two curves cross at
    2.94 km/s, so the relation is continuous.
    """
    vp = np.asarray(vp, float)
    return np.where(
        vp < DENSITY_SWITCH_VP, nafe_drake_density(vp), carlson_raskin_density(vp)
    )


def plane_strain_modulus(layers):
    """M = mu/(1 - nu) = 2 rho Vs**2 (1 - Vs**2/Vp**2) [Pa] of model rows.

    The stiffness that quasi-static compliance of a half-space measures:
    eta = 1/M (Crawford 2004, eq. 5).
    """
    layers = np.asarray(layers, float)
    rho, vp, vs = layers[..., 1], layers[..., 2], layers[..., 3]
    return 2.0 * rho * vs**2 * (1.0 - (vs / vp) ** 2)


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------
@dataclass
class SeafloorModel:
    """Sediment, crustal layers and mantle with a prior about a reference model.

    Parameters, in order: sediment thickness, sediment Vs scale, sediment
    P-modulus factor (the three only if ``sediment_max`` > 0), one Vs per
    crustal layer, mantle Vs, Moho depth below basement, ln of the Vp factor
    of the upper and of the lower crust, ln of the crustal density factor.

    The crustal layers form a fixed grid below basement that reaches the
    deepest Moho the prior allows. The Moho cuts this grid: layers above it
    are crust, the mantle half-space lies below, and the Vs of grid layers
    under the Moho has no effect (it keeps its prior). So the crustal
    thickness is estimated together with Vs instead of being held at the
    reference value.

    Attributes:
        reference (array): Layered model ``[thickness, rho, vp, vs]`` of the
            crystalline crust (no sediment) with the mantle half-space as its
            last row, e.g. ``reference_from_crust1``. It sets the Moho depth
            below basement, the prior mean of crustal Vs, and the mantle.
        first_layer, growth (float): Top crustal layer thickness [m] and the
            ratio between successive layers; 100 m doubling follows Zha & Webb
            (2016) and matches how resolution falls with depth.
        moho_prior (tuple or None): (mean, sigma) [m] of the Gaussian prior on
            the Moho depth below basement. None: the reference crustal
            thickness, with sigma the larger of 30 % of it and 1500 m.
        moho_bounds (tuple or None): Hard bounds [m]. None: from the larger
            of 1000 m and a quarter of the mean, to 1.6 times the mean.
        sediment_max (float): Upper bound of sediment thickness [m]; 0 for none.
        sediment_prior (tuple or None): (mean, sigma) [m] of a Gaussian prior
            on sediment thickness, e.g. from a sediment map; None is uniform.
        sediment_vs, sediment_vp (callable): Reference profiles of depth [m];
            Hamilton (1976, 1979) by default.
        sediment_vs_sigma (float): Sigma of ln(Vs scale).
        sediment_modulus_sigma (float): Sigma of ln(P-modulus factor).
        sediment_sublayers (int): Sublayers, thinner towards the seafloor.
        vs_bounds, vs_mantle (tuple): Hard bounds [m/s].
        vs_sigma (float): Prior sigma of ln(Vs/reference) of crustal layers.
        vs_correlation (float): Prior correlation of neighbouring layers.
        vp_sigma (float): Sigma of ln(Vp factor) on ``oceanic_vp``.
        upper_crust (float): Depth below basement [m] that separates the two
            Vp factors.
        density_sigma (float): Sigma of ln(density factor) on ``crust_density``.
        require_moho_jump (bool): Reject models whose mantle Vs is not above
            the deepest crustal Vs.
    """

    reference: np.ndarray
    first_layer: float = 100.0
    growth: float = 2.0
    moho_prior: tuple | None = None
    moho_bounds: tuple | None = None
    sediment_max: float = 2000.0
    sediment_prior: tuple | None = None
    sediment_vs: object = dram.hamilton_vs
    sediment_vp: object = dram.hamilton_vp
    sediment_vs_sigma: float = 0.3
    sediment_modulus_sigma: float = 0.15
    sediment_sublayers: int = 12
    vs_bounds: tuple = (300.0, 4800.0)
    vs_mantle: tuple = (3900.0, 4900.0)
    vs_sigma: float = 0.3
    vs_correlation: float = 0.7
    vp_sigma: float = 0.07
    upper_crust: float = 1500.0
    density_sigma: float = 0.03
    require_moho_jump: bool = True

    #: Interface with ``dram.sample``: no hierarchical roughness term.
    n_rough = 0
    sigma_rough = (1.0, 1.0)

    def __post_init__(self):
        ref = np.asarray(self.reference, float)
        if ref.ndim != 2 or ref.shape[1] != 4 or len(ref) < 2:
            raise ValueError(
                "reference must have rows [thickness, rho, vp, vs] and a half-space"
            )
        self.reference = ref
        self.moho_depth = float(ref[:-1, 0].sum())  # of the reference
        positive(self.moho_depth, "reference crustal thickness")
        if self.moho_prior is None:
            self.moho_prior = (self.moho_depth, max(0.3 * self.moho_depth, 1500.0))
        if self.moho_bounds is None:
            mean = self.moho_prior[0]
            self.moho_bounds = (max(1000.0, 0.25 * mean), 1.6 * mean)
        self.moho_prior = tuple(map(float, self.moho_prior))
        self.moho_bounds = tuple(map(float, self.moho_bounds))
        if not 0 < self.moho_bounds[0] < self.moho_bounds[1]:
            raise ValueError(
                "moho_bounds must be (lower, upper) with 0 < lower < upper"
            )
        for name in ("first_layer", "vs_sigma", "vp_sigma", "density_sigma",
                     "sediment_vs_sigma", "sediment_modulus_sigma"):  # fmt: skip
            positive(getattr(self, name), name)
        if not 0 <= self.vs_correlation < 1:
            raise ValueError("vs_correlation must lie in [0, 1)")
        if self.growth < 1 or self.sediment_max < 0:
            raise ValueError("growth must be >= 1 and sediment_max >= 0")
        deepest = self.moho_bounds[1]
        grid, edge = [], 0.0
        while edge < deepest - 1e-9:
            t = self.first_layer * self.growth ** len(grid)
            if edge + 1.5 * t >= deepest:  # the last layer takes what is left
                t = deepest - edge
            grid.append(t)
            edge += t
        if len(grid) < 2:
            raise ValueError("the Moho bounds leave fewer than two crustal layers")
        self.fixed = np.array(grid)
        self.n_crust = len(grid)
        self._tops = np.concatenate([[0.0], np.cumsum(self.fixed)[:-1]])
        # Reference Vs of each grid layer: thickness-weighted mean of the reference
        # crust over its depth range (the deepest crustal value below the reference Moho).
        edges = np.concatenate([[0.0], np.cumsum(self.fixed)])
        self.reference_vs = np.array([
            dram._interval_mean(ref[:-1, 0], ref[:-1, 3], ref[-2, 3], a, b)
            for a, b in zip(edges[:-1], edges[1:])
        ])  # fmt: skip
        self._upper = 0.5 * (edges[:-1] + edges[1:]) < self.upper_crust
        self.mantle_ratio = float(ref[-1, 2] / ref[-1, 3])
        self.mantle_density = float(ref[-1, 1])
        names, bounds = [], []
        if self.sediment_max > 0:
            names += ["sediment_thickness", "sediment_vs_scale", "sediment_modulus"]
            bounds += [(0.0, float(self.sediment_max)), (0.4, 2.5), (0.6, 1.6)]
        names += [f"vs_crust_{i + 1}" for i in range(self.n_crust)] + ["vs_mantle"]
        bounds += [tuple(map(float, self.vs_bounds))] * self.n_crust
        bounds.append(tuple(map(float, self.vs_mantle)))
        names += ["moho_depth", "vp_upper", "vp_lower", "density"]
        bounds += [self.moho_bounds, (-0.25, 0.25), (-0.25, 0.25), (-0.12, 0.12)]
        self.names = names
        self.bounds = np.array(bounds, float)
        self._index = {name: i for i, name in enumerate(names)}
        c = self.vs_correlation
        corr = c ** np.abs(
            np.subtract.outer(np.arange(self.n_crust), np.arange(self.n_crust))
        )
        self._precision = np.linalg.inv(corr) / self.vs_sigma**2
        self._chol = np.linalg.cholesky(corr) * self.vs_sigma

    # -- interface used by dram.sample -------------------------------------
    def roughness(self, theta):
        """No roughness term: the prior is in ``log_prior``."""
        return 0.0

    def _get(self, theta, name, default):
        i = self._index.get(name)
        return default if i is None else theta[i]

    def crust_vs(self, theta):
        """Crustal layer shear velocities in ``theta`` [m/s]."""
        i = self._index["vs_crust_1"]
        return np.asarray(theta[i : i + self.n_crust], float)

    def interfaces(self, theta):
        """Sediment base and crustal layer bottoms below the seafloor [m] (to the Moho)."""
        h = self._get(theta, "sediment_thickness", 0.0)
        moho = theta[self._index["moho_depth"]]
        bottoms = np.cumsum(self.fixed)
        return h, h + np.append(bottoms[bottoms < moho], moho)

    def sediment_rows(self, theta):
        """Sediment sublayers ``[thickness, rho, vp, vs]`` (empty below 1 m)."""
        h = self._get(theta, "sediment_thickness", 0.0)
        if h <= 1.0:
            return np.empty((0, 4))
        edges = (
            h * (np.arange(self.sediment_sublayers + 1) / self.sediment_sublayers) ** 2
        )
        mid = 0.5 * (edges[:-1] + edges[1:])
        vs = self._get(theta, "sediment_vs_scale", 1.0) * self.sediment_vs(mid)
        vp0 = np.asarray(self.sediment_vp(mid), float)
        rho = nafe_drake_density(vp0)
        vp = vp0 * math.sqrt(self._get(theta, "sediment_modulus", 1.0))
        vp = np.maximum(
            vp, dram.brocher_vp(vs)
        )  # stiff sediment is not slower than rock of its Vs
        return np.column_stack([np.diff(edges), rho, vp, vs])

    def layers(self, theta):
        """Layered model ``[thickness, rho, vp, vs]`` (last row half-space)."""
        theta = np.asarray(theta, float)
        vs = self.crust_vs(theta)
        factor = np.where(self._upper, theta[self._index["vp_upper"]],
                          theta[self._index["vp_lower"]])  # fmt: skip
        vp = oceanic_vp(vs) * np.exp(factor)
        rho = crust_density(vp) * math.exp(theta[self._index["density"]])
        moho = theta[self._index["moho_depth"]]
        n = max(
            int(np.searchsorted(self._tops, moho - 1.0, side="left")), 1
        )  # layers above
        thick = self.fixed[:n].copy()
        thick[-1] = moho - self._tops[n - 1]
        vsm = theta[self._index["vs_mantle"]]
        crust = np.column_stack([thick, rho[:n], vp[:n], vs[:n]])
        mantle = np.array([[0.0, self.mantle_density, self.mantle_ratio * vsm, vsm]])
        return np.vstack([self.sediment_rows(theta), crust, mantle])

    def vs_at(self, theta, z):
        """Vs [m/s] of the model at depths ``z`` [m] below the seafloor."""
        model = self.layers(theta)
        tops = np.concatenate([[0.0], np.cumsum(model[:-1, 0])])
        return model[np.searchsorted(tops, np.asarray(z, float), side="right") - 1, 3]

    def column_at(self, theta, z, column):
        """Any model column (1 rho, 2 Vp, 3 Vs) or 'M' at depths ``z`` [m]."""
        model = self.layers(theta)
        tops = np.concatenate([[0.0], np.cumsum(model[:-1, 0])])
        rows = model[np.searchsorted(tops, np.asarray(z, float), side="right") - 1]
        return plane_strain_modulus(rows) if column == "M" else rows[:, column]

    def in_bounds(self, theta):
        """True if every parameter lies inside its hard bounds."""
        return bool(
            np.all(theta >= self.bounds[:, 0]) and np.all(theta <= self.bounds[:, 1])
        )

    def random(self, rng):
        """A draw from the prior (see ``prior_draw``)."""
        return self.prior_draw(rng)

    # -- prior ---------------------------------------------------------------
    def log_prior(self, theta):
        """Log prior density inside the bounds (up to a constant)."""
        lp = 0.0
        ix = self._index
        if "sediment_thickness" in ix:
            if self.sediment_prior is not None:
                mean, sigma = self.sediment_prior
                lp -= 0.5 * ((theta[ix["sediment_thickness"]] - mean) / sigma) ** 2
            s = math.log(theta[ix["sediment_vs_scale"]])
            m = math.log(theta[ix["sediment_modulus"]])
            lp -= 0.5 * (s / self.sediment_vs_sigma) ** 2 + s
            lp -= 0.5 * (m / self.sediment_modulus_sigma) ** 2 + m
        d = np.log(self.crust_vs(theta) / self.reference_vs)
        lp -= 0.5 * float(d @ self._precision @ d) + float(
            np.sum(np.log(self.crust_vs(theta)))
        )
        lp -= (
            0.5
            * ((theta[ix["moho_depth"]] - self.moho_prior[0]) / self.moho_prior[1]) ** 2
        )
        lp -= 0.5 * (theta[ix["vp_upper"]] / self.vp_sigma) ** 2
        lp -= 0.5 * (theta[ix["vp_lower"]] / self.vp_sigma) ** 2
        lp -= 0.5 * (theta[ix["density"]] / self.density_sigma) ** 2
        return float(lp)

    def valid(self, theta):
        """Inside the bounds, Poisson ratio >= 0 everywhere, and the Moho rule."""
        if not self.in_bounds(theta):
            return False
        layers = self.layers(theta)
        if np.any(layers[:, 2] < math.sqrt(2) * layers[:, 3]):
            return False
        return not (self.require_moho_jump and layers[-1, 3] <= layers[-2, 3])

    def prior_draw(self, rng, max_tries=100000):
        """One exact draw from the prior, including its hard constraints."""
        ix = self._index
        lo, hi = self.bounds[:, 0], self.bounds[:, 1]
        for _ in range(max_tries):
            theta = np.empty(len(self.names))
            if "sediment_thickness" in ix:
                if self.sediment_prior is None:
                    theta[ix["sediment_thickness"]] = rng.uniform(
                        0.0, self.sediment_max
                    )
                else:
                    theta[ix["sediment_thickness"]] = rng.normal(*self.sediment_prior)
                theta[ix["sediment_vs_scale"]] = math.exp(
                    rng.normal(0, self.sediment_vs_sigma)
                )
                theta[ix["sediment_modulus"]] = math.exp(
                    rng.normal(0, self.sediment_modulus_sigma)
                )
            i = ix["vs_crust_1"]
            theta[i : i + self.n_crust] = self.reference_vs * np.exp(
                self._chol @ rng.standard_normal(self.n_crust)
            )
            theta[ix["vs_mantle"]] = rng.uniform(*self.vs_mantle)
            theta[ix["moho_depth"]] = rng.normal(*self.moho_prior)
            theta[ix["vp_upper"]] = rng.normal(0, self.vp_sigma)
            theta[ix["vp_lower"]] = rng.normal(0, self.vp_sigma)
            theta[ix["density"]] = rng.normal(0, self.density_sigma)
            if np.all(theta >= lo) and np.all(theta <= hi) and self.valid(theta):
                return theta
        raise ValueError("no valid prior draw; check the reference model and bounds")

    def reference_theta(self, sediment=None):
        """Parameter vector of the reference model (the centre of the prior).

        ``sediment`` is the sediment thickness [m]; default the prior mean, or
        a tenth of ``sediment_max`` for a uniform prior.
        """
        theta = np.zeros(len(self.names))
        ix = self._index
        if "sediment_thickness" in ix:
            if sediment is None:
                sediment = (self.sediment_prior[0] if self.sediment_prior is not None
                            else 0.1 * self.sediment_max)  # fmt: skip
            theta[ix["sediment_thickness"]] = sediment
            theta[ix["sediment_vs_scale"]] = theta[ix["sediment_modulus"]] = 1.0
        i = ix["vs_crust_1"]
        theta[i : i + self.n_crust] = self.reference_vs
        theta[ix["vs_mantle"]] = self.reference[-1, 3]
        theta[ix["moho_depth"]] = self.moho_prior[0]
        theta = np.clip(theta, self.bounds[:, 0], self.bounds[:, 1])
        if self.require_moho_jump:  # the reference itself may break the rule
            top = max(theta[i : i + self.n_crust].max() + 100.0, theta[ix["vs_mantle"]])
            theta[ix["vs_mantle"]] = min(top, self.vs_mantle[1])
            theta[i : i + self.n_crust] = np.minimum(theta[i : i + self.n_crust],
                                                     theta[ix["vs_mantle"]] - 100.0)  # fmt: skip
        return theta

    def settings(self):
        """JSON-ready dictionary of the model settings."""
        out = {
            k: v
            for k, v in asdict(self).items()
            if not callable(v) and k != "reference"
        }
        out.update(reference=self.reference.tolist(), names=self.names,
                   fixed_layers=self.fixed.tolist(), reference_vs=self.reference_vs.tolist(),
                   moho_depth=self.moho_depth, model="SeafloorModel",
                   sediment_vs=getattr(self.sediment_vs, "__name__", "custom"),
                   sediment_vp=getattr(self.sediment_vp, "__name__", "custom"))  # fmt: skip
        return out


#: CRUST1.0 sediment thickness [m] from which the sediment layers become part of
#: the reference model instead of being left to the soft-sediment unit.
THICK_SEDIMENT = 500.0


def reference_from_crust1(column, sediment=None):
    """Reference model from a ``dram.crust1_column``.

    Args:
        column (dict): CRUST1.0 column at the station.
        sediment (bool or None): Keep the CRUST1.0 sediment layers as the top
            of the reference. None keeps them when they are thicker than
            ``THICK_SEDIMENT``: on shelves, slopes and fans the section below
            the seafloor is kilometres of consolidated sediment, which the
            igneous crust of the reference cannot stand for. Thin open-ocean
            sediment is left to the soft-sediment unit of ``SeafloorModel``.

    Returns:
        ndarray: rows ``[thickness, rho, vp, vs]`` with the mantle half-space last.
    """
    if sediment is None:
        sediment = column.get("sediment_thickness", 0.0) >= THICK_SEDIMENT
    rows = [[t, rho, vp, vs] for name, t, vp, vs, rho in column["layers"]
            if sediment or "sediment" not in name]  # fmt: skip
    if not rows:
        raise ValueError("the CRUST1.0 column has no crust")
    vp, vs, rho = column["mantle"]
    return np.array(rows + [[0.0, rho, vp, vs]], float)


def model_from_crust1(column, **overrides):
    """``SeafloorModel`` of a station from its CRUST1.0 column, with the default rules.

    Where CRUST1.0 has at least ``THICK_SEDIMENT`` of sediment, its sediment
    layers join the reference and the prior on ln Vs is widened to 0.5, because
    a 1-degree sediment model is a weaker guide than normal oceanic crust is.
    Elsewhere the reference is the crystalline crust with sigma 0.3. Keyword
    arguments override any ``SeafloorModel`` field.
    """
    thick = column.get("sediment_thickness", 0.0) >= THICK_SEDIMENT
    settings = dict(vs_sigma=0.5 if thick else 0.3)
    settings.update(overrides)
    return SeafloorModel(reference_from_crust1(column), **settings)


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
#: Gaussian prior sigma of the pressure-gain factor by calibration source.
#: Two independent calibrations of 73 archive stations agree to 2-3 % (robust)
#: with a 16-84 % range of -5.4 to +2.6 %, so 5 % for a calibrated gauge.
GAIN_SIGMA = {"compy": 0.05, "bruitfm": 0.05, "calibrated": 0.05, "nominal": 0.2}
#: Tilt of the curve (relative, band edge to band edge) left by the pressure
#: response when it has not been checked against the tide, and when it has.
RESPONSE_TILT = {"unchecked": 0.04, "corrected": 0.02}
#: Smooth model error (2-D structure, estimator; relative) and its length in ln f.
MODEL_ERROR = dict(smooth_rel=0.005, smooth_length=0.5)


def compliance_data(freqs, eta, depth, *, cov_rel=None, sigma=None, band=None,
                    max_fraction_fc=0.8, response="unchecked", tilt_rel=None,
                    smooth_rel=MODEL_ERROR["smooth_rel"],
                    smooth_length=MODEL_ERROR["smooth_length"], station=""):  # fmt: skip
    """``dram.ComplianceData`` with the default error model of ComPy 3.

    Args:
        freqs, eta (array): Frequencies [Hz] and normalized compliance [1/Pa].
        depth (float): Water depth [m].
        cov_rel (array): Relative covariance of the stack over ``freqs``
            (``compy.inversion.errors.stack_covariance``).
        sigma (array): One-sigma errors [1/Pa], used through
            ``errors.default_covariance`` when ``cov_rel`` is not given.
        band (tuple): Frequencies kept [Hz]; default all given.
        max_fraction_fc (float): Upper frequency limit as a fraction of the
            infragravity cut-off f_c = sqrt(g/(2 pi H)); bins above it are
            dropped because the pressure signal vanishes there.
        response (str): "unchecked" or "corrected" pressure response; sets the
            tilt term unless ``tilt_rel`` is given.
        smooth_rel, smooth_length: Smooth model-error term.
        station (str): Label.

    The gain is not in the covariance: ``invert`` samples it.
    """
    f, eta = np.asarray(freqs, float), np.asarray(eta, float)
    depth = positive(depth, "water depth")
    if cov_rel is None:
        if sigma is None:
            raise ValueError("give cov_rel or sigma")
        cov_rel = default_covariance(f, eta, np.asarray(sigma, float))
    cov_rel = np.asarray(cov_rel, float)
    fc = math.sqrt(dram.GRAVITY / (2 * math.pi * depth))
    lo, hi = (f.min(), f.max()) if band is None else band
    keep = (f >= lo - 1e-12) & (f <= min(hi, max_fraction_fc * fc) + 1e-12)
    if keep.sum() < 3:
        raise ValueError(f"fewer than 3 frequencies between {lo:.4f} Hz and "
                         f"{max_fraction_fc:.2f} f_c = {max_fraction_fc * fc:.4f} Hz")  # fmt: skip
    f, eta = f[keep], eta[keep]
    cov = cov_rel[np.ix_(keep, keep)] * np.outer(eta, eta)
    tilt = RESPONSE_TILT[response] if tilt_rel is None else tilt_rel
    data = dram.ComplianceData.build(f, eta, depth, random_cov=cov, gain_rel=0.0, floor_rel=0.0,
                                     tilt_rel=tilt, smooth_rel=smooth_rel,
                                     smooth_length=smooth_length, station=station)  # fmt: skip
    return data


# --------------------------------------------------------------------------
# Inversion
# --------------------------------------------------------------------------
#: Sampler settings of ``invert`` (first round; iterations double each round).
SAMPLER = dict(n_chains=8, n_iter=60000, n_burn=20000, thin=10)
#: Parallel-tempering ladder.
TEMPERATURES = tuple(float(v) for v in 1.6 ** np.arange(6))
#: Convergence required of every parameter.
CONVERGENCE = dict(rhat=1.05, ess=400.0)


@dataclass
class InversionResult:
    """Posterior of ``invert`` with the prior and the summaries built from them."""

    posterior: dram.Posterior
    prior_theta: np.ndarray
    rounds: int
    converged: bool
    gain_sigma: float
    info: dict = field(default_factory=dict)

    @property
    def model(self):
        """The ``SeafloorModel`` inverted."""
        return self.posterior.model

    @property
    def data(self):
        """The ``ComplianceData`` inverted."""
        return self.posterior.data

    def draws(self, z, column=3, prior=False, max_draws=4000, seed=0):
        """Model column (1 rho, 2 Vp, 3 Vs, 'M') of draws at depths ``z``."""
        theta = self.prior_theta if prior else self.posterior.flat(max_draws, seed)
        n = len(self.model.names)
        return np.array(
            [self.model.column_at(t[:n], z, column) for t in theta[:max_draws]]
        )

    def reach(self, fraction=0.9):
        """Depth [m] above which ``fraction`` of the sensitivity of the lowest frequency lies.

        From the Vs kernel (Vp/Vs fixed) of the posterior median model.
        """
        theta = np.median(self.posterior.flat(), axis=0)[: len(self.model.names)]
        layers = self.model.layers(theta)
        zmax = 2.0 * self.model.moho_bounds[1]
        z, kern = dram.vs_kernels(self.data.depth, self.data.freqs[:1], layers,
                                  dz=min(100.0, zmax / 200), max_depth=zmax)  # fmt: skip
        return dram.sensitivity_depth(z, kern[0], fraction)

    def parameter_resolution(self):
        """Posterior over prior width of every parameter (16-84 % interval).

        For velocities the width is that of ln Vs. Values near 1 mean the data
        add nothing to the prior; small values mean the data constrain the
        parameter. A parameter-wise measure is used instead of a pointwise one
        because the prior of Vs at a fixed depth mixes sediment, crust and
        mantle draws and is far wider than any of them.
        """
        flat = self.posterior.flat()
        out = {}
        for i, name in enumerate(self.posterior.names):
            a = np.quantile(flat[:, i], [0.16, 0.84])
            b = np.quantile(self.prior_theta[:, i], [0.16, 0.84])
            if name.startswith("vs_"):
                a, b = np.log(a), np.log(b)
            out[name] = (
                float((a[1] - a[0]) / (b[1] - b[0])) if b[1] > b[0] else float("nan")
            )
        return out

    def profile(self, z=None, q=(0.05, 0.16, 0.5, 0.84, 0.95), resolved_ratio=0.6):
        """Vs against depth with intervals, the prior, and the resolved range.

        Each depth takes the resolution of the parameter that holds it in the
        posterior median model: the sediment Vs scale above the median sediment
        base, the Vs of the crustal layer below it, the mantle Vs below the
        median Moho. A depth counts as resolved when that ratio
        (``parameter_resolution``) is at most ``resolved_ratio`` and the depth
        is above ``reach(0.9)``. At 0.6 the data have cut the prior standard
        deviation of ln Vs by at least 40 %. Thin layers near the seafloor are
        seldom resolved one by one even where their average is: use
        ``width_ratio`` as a continuous measure.

        Returns:
            dict of arrays over ``z``: ``z``, ``vs`` (len(q), len(z)),
            ``vs_prior``, ``modulus`` (quantiles of M = mu/(1 - nu) [Pa]),
            ``width_ratio``, ``resolved`` (bool), and the scalars ``reach``,
            ``sediment`` and ``moho`` (posterior medians, the Moho below the
            seafloor) and ``q``.
        """
        model = self.model
        if z is None:
            z = np.unique(np.concatenate([np.arange(0.0, 500.0, 10.0),
                                          np.arange(500.0, model.moho_bounds[1] + 2000.0, 50.0)]))  # fmt: skip
        z = np.asarray(z, float)
        post, pri = self.draws(z), self.draws(z, prior=True)
        qs, qp = np.quantile(post, q, axis=0), np.quantile(pri, q, axis=0)
        res = self.parameter_resolution()
        med = np.median(self.posterior.flat(), axis=0)
        h = (
            float(med[model._index["sediment_thickness"]])
            if "sediment_thickness" in model._index
            else 0.0
        )
        moho = h + float(med[model._index["moho_depth"]])
        layer = np.clip(
            np.searchsorted(model._tops, z - h, side="right") - 1, 0, model.n_crust - 1
        )
        ratio = np.array([res[f"vs_crust_{k + 1}"] for k in layer])
        ratio = np.where(z >= moho, res["vs_mantle"], ratio)
        if h > 0:
            ratio = np.where(z < h, res["sediment_vs_scale"], ratio)
        reach = self.reach(0.9)
        return dict(z=z, vs=qs, vs_prior=qp, modulus=np.quantile(self.draws(z, "M"), q, axis=0),
                    width_ratio=ratio, resolved=(ratio <= resolved_ratio) & (z <= reach),
                    reach=reach, sediment=h, moho=moho, q=np.asarray(q))  # fmt: skip

    def quality(self):
        """Label ("good", "usable" or "poor") with the reasons.

        good: converged; median chi-square at most the 95 % point of a
        chi-square with N degrees of freedom; the gain, the Vp factors and the
        density factor within 2 sigma of their prior means; some depth
        resolved. usable: R-hat <= 1.1, chi-square at most twice that point,
        those factors within 3 sigma, some depth resolved. poor: anything else.

        A factor pulled far from its prior means the data cannot be fitted by
        a plausible Earth with the given calibration: at 7D.J35A, with a
        pressure gain of 4.4, the inversion needs a gain factor of 0.83, Vp
        25 % above the rock relation and crustal Vs at 4.7 km/s.
        """
        from scipy.stats import chi2 as chi2_dist

        s = self.posterior.summary()
        n = s["n_data"]
        limit = float(chi2_dist.ppf(0.95, n))
        chi2 = s["chi2"]["median"]
        par = s["parameters"]
        gain = par.get("gain", {}).get("median", 1.0)
        model = self.model
        pulls = {"gain": abs(gain - 1.0) / self.gain_sigma}
        for name, sigma in (("vp_upper", model.vp_sigma), ("vp_lower", model.vp_sigma),
                            ("density", model.density_sigma)):  # fmt: skip
            pulls[name] = abs(par[name]["median"]) / sigma
        worst = max(pulls, key=pulls.get)
        resolved = bool(np.any(self.profile()["resolved"]))
        res_ratio = self.parameter_resolution()
        reasons = []
        if not self.converged:
            reasons.append(
                f"not converged (max R-hat {s['max_rhat']:.3f}, min ESS {s['min_ess']:.0f})"
            )
        if chi2 > limit:
            reasons.append(f"misfit: chi2 {chi2:.1f} above {limit:.1f} for {n} data")
        for name, pull in pulls.items():
            if pull > 2:
                what = (
                    f"factor {gain:.3f}"
                    if name == "gain"
                    else f"ln factor {par[name]['median']:+.3f}"
                )
                reasons.append(
                    f"{name} pulled {pull:.1f} sigma from its prior ({what})"
                )
        if not resolved:
            reasons.append("no depth resolved beyond the prior")
        lo, hi = model.vs_bounds
        at_bound = [name for name in model.names if name.startswith("vs_crust_")
                    and res_ratio.get(name, 1.0) <= 0.6
                    and not 1.05 * lo < par[name]["median"] < 0.98 * hi]  # fmt: skip
        if at_bound:
            reasons.append("resolved Vs at a hard bound: " + ", ".join(at_bound))
        if not reasons:
            label = "good"
        elif (
            s["max_rhat"] <= 1.1
            and chi2 <= 2 * limit
            and resolved
            and pulls[worst] <= 3
        ):
            label = "usable"
        else:
            label = "poor"
        return dict(label=label, reasons=reasons, chi2=chi2, chi2_limit=limit, n_data=n,
                    max_rhat=s["max_rhat"], min_ess=s["min_ess"], gain=gain,
                    gain_pull=pulls["gain"], pulls=pulls)  # fmt: skip

    def summary(self):
        """JSON-ready record: posterior summary, quality, resolution, prior."""
        s = self.posterior.summary()
        prof = self.profile()
        res = prof["z"][prof["resolved"]]
        s.update(
            quality=self.quality(), rounds=self.rounds, converged=self.converged,
            gain_sigma=self.gain_sigma, reach_m=prof["reach"],
            resolved_from_m=float(res.min()) if len(res) else None,
            resolved_to_m=float(res.max()) if len(res) else None,
            resolution=self.parameter_resolution(),
            prior={name: dict(zip(("p05", "p16", "median", "p84", "p95"),
                                  map(float, np.quantile(self.prior_theta[:, i],
                                                         [0.05, 0.16, 0.5, 0.84, 0.95]))))
                   for i, name in enumerate(self.posterior.names)},
            covariance_terms={k: float(np.sqrt(np.median(np.diag(v)) )
                                       / np.median(self.data.eta))
                              for k, v in self.data.components.items()},
            info=self.info,
        )  # fmt: skip
        return s

    def save(self, path):
        """Write ``<path>.npz`` (draws, prior draws, data), ``.json`` and ``_profile.csv``."""
        path = Path(path)
        self.posterior.save(path)
        with np.load(path.with_suffix(".npz")) as d:
            content = {k: d[k] for k in d.files}
        np.savez_compressed(
            path.with_suffix(".npz"), prior_theta=self.prior_theta, **content
        )
        path.with_suffix(".json").write_text(
            json.dumps(self.summary(), indent=1, default=float)
        )
        prof = self.profile()
        cols = ["depth_m"] + [f"vs_p{round(100 * v):02d}" for v in prof["q"]]
        cols += [f"prior_p{round(100 * v):02d}" for v in prof["q"]]
        cols += [f"modulus_p{round(100 * v):02d}_pa" for v in prof["q"]] + [
            "width_ratio",
            "resolved",
        ]
        table = np.column_stack([prof["z"], prof["vs"].T, prof["vs_prior"].T, prof["modulus"].T,
                                 prof["width_ratio"], prof["resolved"].astype(float)])  # fmt: skip
        np.savetxt(path.with_name(path.name + "_profile.csv"), table, delimiter=",",
                   header=",".join(cols), comments="", fmt="%.6g")  # fmt: skip


def invert(data, model, *, gain_sigma=0.05, seed=0, processes=None, max_rounds=3,
           n_prior=4000, temperatures=TEMPERATURES, convergence=None, info=None,
           **sampler):  # fmt: skip
    """Invert compliance for a ``SeafloorModel`` and summarize the result.

    Chains start from draws of the prior (the first at the reference model), so
    the result does not depend on a starting model: R-hat compares chains that
    began far apart. Sampling repeats with twice the iterations, up to
    ``max_rounds`` times, until every parameter has R-hat and bulk effective
    sample size within ``convergence``. The prior is drawn directly
    (``model.prior_draw``), not by a chain.

    Args:
        data (dram.ComplianceData): e.g. from ``compliance_data``.
        model (SeafloorModel): Parameterization and prior.
        gain_sigma (float): Prior sigma of the pressure-gain factor (see
            ``GAIN_SIGMA``); the factor is bounded at 4 sigma.
        seed (int): Makes the run reproducible.
        processes (int): Worker processes (1 runs serially).
        max_rounds (int): Rounds of sampling at most.
        n_prior (int): Prior draws kept for the comparison.
        temperatures (tuple): Parallel-tempering ladder, starting at 1.
        convergence (dict): ``rhat`` and ``ess`` limits.
        **sampler: Overrides of ``SAMPLER`` (n_chains, n_iter, n_burn, thin).

    Returns:
        InversionResult
    """
    conv = dict(CONVERGENCE, **(convergence or {}))
    settings = dict(SAMPLER, **sampler)
    gain_sigma = positive(gain_sigma, "gain_sigma")
    bounds = (max(1 - 4 * gain_sigma, 0.2), 1 + 4 * gain_sigma)
    rng = np.random.default_rng(seed)
    prior = []
    while len(prior) < n_prior:
        theta = model.prior_draw(rng)
        g = rng.normal(1.0, gain_sigma)
        if bounds[0] <= g <= bounds[1]:
            prior.append(np.append(theta, g))
    prior = np.array(prior)
    target = dram._Target(data, model, False, bounds, gain_sigma)
    starts = [np.append(model.reference_theta(), 1.0)]
    if target(starts[0]) is None:
        starts = []
    for theta in prior:
        if len(starts) >= settings["n_chains"]:
            break
        if target(theta) is not None:
            starts.append(theta)
    if len(starts) < settings["n_chains"]:
        raise ValueError(
            "the forward model fails on the prior draws; check the reference model"
        )
    post, converged, rounds = None, False, 0
    for rounds in range(1, max_rounds + 1):
        post = dram.sample(data, model, seed=seed + rounds - 1, processes=processes, starts=starts,
                           gain_bounds=bounds, gain_sigma=gain_sigma,
                           temperatures=temperatures, **settings)  # fmt: skip
        diag = post.diagnostics()
        rhat = max(v["rhat"] for v in diag.values() if v["rhat"] is not None)
        ess = min(v["ess_bulk"] for v in diag.values() if v["ess_bulk"] is not None)
        converged = rhat <= conv["rhat"] and ess >= conv["ess"]
        if converged:
            break
        settings["n_iter"], settings["n_burn"] = (
            2 * settings["n_iter"],
            2 * settings["n_burn"],
        )
    return InversionResult(post, prior, rounds, converged, gain_sigma, dict(info or {}))


def load_result(path):
    """Read an ``InversionResult`` written by ``InversionResult.save``.

    ``path`` is the common stem (with or without ``.npz``). The model is
    rebuilt from the settings in the JSON file; custom sediment profiles are
    not stored, so a result made with them must be given the same callables
    again through ``result.model``.
    """
    from dataclasses import fields

    path = Path(path)
    path = path.with_suffix("") if path.suffix in (".npz", ".json") else path
    record = json.loads(path.with_name(path.name + ".json").read_text())
    settings = record["model"]
    allowed = {f.name for f in fields(SeafloorModel)} - {
        "reference",
        "sediment_vs",
        "sediment_vp",
    }
    kwargs = {k: (tuple(v) if isinstance(v, list) else v)
              for k, v in settings.items() if k in allowed}  # fmt: skip
    model = SeafloorModel(np.array(settings["reference"], float), **kwargs)
    with np.load(path.with_name(path.name + ".npz")) as d:
        data = dram.ComplianceData(d["freqs"], d["eta"], float(d["depth"]), d["covariance"],
                                   station=record.get("station", ""))  # fmt: skip
        names = [str(n) for n in d["names"]]
        post = dram.Posterior(
            data=data, model=model, prior_only=False, theta=d["theta"], alpha=d["alpha"],
            loglik=d["loglik"], roughness=d["roughness"], pred=d["pred"],
            acceptance=np.array(record["acceptance"]),
            acceptance_dr=np.array(record["acceptance_dr"]), settings=record["settings"],
            trace_chi2=d["trace_chi2"], names_extra=names[len(model.names):],
        )  # fmt: skip
        prior = d["prior_theta"]
    return InversionResult(post, prior, record["rounds"], record["converged"],
                           record["gain_sigma"], record.get("info", {}))  # fmt: skip
