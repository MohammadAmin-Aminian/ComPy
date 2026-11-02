"""Layered starting and reference models for the compliance inversion.

Models are arrays of layers ``[thickness, density, vp, vs]`` in SI units.
Includes generic sediment-over-crust templates, CRUST1.0 starting models,
the hard-coded RHUM-RUM station templates of the paper, and the empirical
Vp/density relations (Poisson ratio, Gardner, Brocher 2005).
"""

import numpy as np

from compy.numerics import integer, positive


def _depth_profile(model, depth_grid):
    """Sample a layered model on a fixed depth grid; last row is half-space."""
    interfaces = np.cumsum(model[:-1, 0])
    return model[np.searchsorted(interfaces, depth_grid, side="right"), 3]


def _validate_model(model):
    model = np.asarray(model, dtype=float)
    if model.ndim != 2 or model.shape[1] != 4 or len(model) < 2:
        raise ValueError(
            "model must have at least two rows and columns [thickness, rho, Vp, Vs]"
        )
    if (
        np.any(~np.isfinite(model))
        or np.any(model[:, 0] < 0)
        or np.any(model[:, 1:] <= 0)
    ):
        raise ValueError(
            "model must be finite, with nonnegative thickness and positive rho, Vp, Vs"
        )
    if np.any(model[:-1, 0] <= 0) or np.any(model[:, 2] <= model[:, 3]):
        raise ValueError(
            "finite layers must have positive thickness and Vp must exceed Vs"
        )
    return model


def exponential_model(iteration, first_layer=200, n_layer=15, power_factor=1.15):
    """Generic exponential-thickness crustal model without sediments."""
    return sediment_crust_model(
        iteration, first_layer, n_layer, power_factor, n_sediment_layer=0
    )


def sediment_crust_model(
    iteration,
    first_layer=200,
    n_layer=10,
    power_factor=1.17,
    sediment_thickness=80,
    n_sediment_layer=1,
):
    """Generic sediment plus graded-crust starting model (SI units).

    Use ``crust1_starting_model`` for a station's CRUST1.0 column. The last
    row is the elastic half-space, represented by zero thickness.
    """
    iteration = integer(iteration, "iteration")
    n_layer = integer(n_layer, "n_layer")
    n_sediment_layer = integer(n_sediment_layer, "n_sediment_layer", minimum=0)
    first_layer = positive(first_layer, "first_layer")
    power_factor = positive(power_factor, "power_factor")
    sediment = np.empty((0, 4))
    if n_sediment_layer:
        sediment_thickness = positive(sediment_thickness, "sediment_thickness")
        weights = np.geomspace(1, sediment_thickness, n_sediment_layer)
        vs = np.geomspace(340, 600, n_sediment_layer)
        sediment = np.column_stack(
            [
                sediment_thickness * weights / weights.sum(),
                gardner_density(vs),
                vp_from_vs(vs),
                vs,
            ]
        )
    vs = np.linspace(2700, 4300, n_layer)
    crust = np.column_stack(
        [
            first_layer * power_factor ** np.arange(n_layer),
            gardner_density(vs),
            vp_from_vs(vs),
            vs,
        ]
    )
    model = np.vstack([sediment, crust, [0, 3340, 8120, 4510]])
    _validate_model(model)
    chain = np.zeros((*model.shape, iteration))
    chain[:, :, 0] = model
    grid = np.arange(max(1, int(np.ceil(model[:-1, 0].sum()))))
    prior = _depth_profile(model, grid)[None, :, None]
    return chain, prior, np.zeros_like(prior)


def crust1_starting_model(column, sublayers=1, mantle_layer=10000.0):
    """Starting model ``[thickness, rho, vp, vs]`` from CRUST1.0 at a station.

    Replaces the hard-coded ``rhumrum_starting_model`` templates. Read the column with
    ``compy.inversion.crust1_column(lat, lon, crust1_dir)``, which picks the
    1x1 degree CRUST1.0 cell containing the station.

    Args:
        column (dict): ``crust1_column`` at the station.
        sublayers (int): Equal sublayers per CRUST1.0 crystalline crustal
            layer (1, 2, 3 and 4 give the 3, 6, 9 and 12 crustal layers of
            the former templates). Sediment layers stay whole.
        mantle_layer (float): Thickness of a finite mantle layer above the
            half-space [m], as in the former templates; 0 omits it.

    Returns:
        ndarray: One row per layer; the last row is the mantle half-space
        (thickness 0, ignored by the propagator).
    """
    sublayers = integer(sublayers, "sublayers")
    rows = []
    for name, thickness, vp, vs, rho in column["layers"]:
        n = 1 if "sediment" in name else sublayers
        rows += [[thickness / n, rho, vp, vs]] * n
    vp, vs, rho = column["mantle"]
    if mantle_layer:
        rows.append([positive(mantle_layer, "mantle_layer"), rho, vp, vs])
    rows.append([0.0, rho, vp, vs])
    return _validate_model(np.array(rows, float))


def rhumrum_starting_model(iteration, n_layer=3, sta="RR38"):
    """Hard-coded RHUM-RUM starting templates (deprecated).

    Most templates do not match CRUST1.0 at their station; use
    ``crust1_starting_model`` instead.
    """
    import warnings

    warnings.warn(
        "rhumrum_starting_model templates do not match CRUST1.0 at most stations; use "
        "crust1_starting_model(compy.inversion.crust1_column(lat, lon, dir))",
        DeprecationWarning,
        stacklevel=2,
    )
    iteration = integer(iteration, "iteration")
    if n_layer not in (3, 6, 9, 12):
        raise ValueError("n_layer must be 3, 6, 9, or 12")
    if sta not in (
        "RR28",
        "RR29",
        "RR34",
        "RR36",
        "RR38",
        "RR40",
        "RR50",
        "RR52",
    ) and n_layer not in (3, 6):
        raise ValueError("generic station template supports only 3 or 6 layers")
    dep = 0
    if sta == "RR28":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            # Thickness, Density, Vp, Vs
            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [260.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

    elif sta == "RR29":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

    elif sta == "RR34":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                    [10000.0, 3230.0, 7530.0, 4370.0],
                ]
            )

    elif sta == "RR36":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

    elif sta == "RR38":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

    elif sta == "RR40":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [10.0, 1820.0, 1750.0, 340.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                ]
            )

    elif sta == "RR50":
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [700.0, 2550.0, 5000.0, 2700.0],
                    [1540.0, 2850.0, 6500.0, 3700.0],
                    [4750.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [350.0, 2550.0, 5000.0, 2700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [770.0, 2850.0, 6500.0, 3700.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [2375.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [233.0, 2550.0, 5000.0, 2700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [513.0, 2850.0, 6500.0, 3700.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [1583.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [175.0, 2550.0, 5000.0, 2700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [385.0, 2850.0, 6500.0, 3700.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [1187.0, 3050.0, 7100.0, 4050.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3100.0, 7530.0, 4190.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )
    elif sta == "RR52":
        sediment_thickness = 0
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [1170.0, 2400.0, 4700.0, 2540.0],
                    [1170.0, 2760.0, 6300.0, 3590.0],
                    [4560.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 3, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [585.0, 2400.0, 4700.0, 2540.0],
                    [585.0, 2400.0, 4700.0, 2540.0],
                    [585.0, 2760.0, 6300.0, 3590.0],
                    [585.0, 2760.0, 6300.0, 3590.0],
                    [2280.0, 2960.0, 6900.0, 3940.0],
                    [2280.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 9:
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [390.0, 2400.0, 4700.0, 2540.0],
                    [390.0, 2400.0, 4700.0, 2540.0],
                    [390.0, 2400.0, 4700.0, 2540.0],
                    [390.0, 2760.0, 6300.0, 3590.0],
                    [390.0, 2760.0, 6300.0, 3590.0],
                    [390.0, 2760.0, 6300.0, 3590.0],
                    [1520.0, 2960.0, 6900.0, 3940.0],
                    [1520.0, 2960.0, 6900.0, 3940.0],
                    [1520.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        else:
            n_layer = 12
            starting_model = np.zeros([n_layer + 3, 4, iteration])
            starting_model[:, :, 0] = np.array(
                [
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2400.0, 4700.0, 2540.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [292.0, 2760.0, 6300.0, 3590.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [1140.0, 2960.0, 6900.0, 3940.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )
    # elif sta == "A422A":
    else:
        # sediment_thickness = 0
        if n_layer == 3:
            starting_model = np.zeros([n_layer + 4, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [2000.0, 2400.0, 4700.0, 2000.0],
                    [2000.0, 2760.0, 6300.0, 2500.0],
                    [3000.0, 2760.0, 6300.0, 3200.0],
                    [3000.0, 2960.0, 6900.0, 3800.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

        elif n_layer == 6:
            starting_model = np.zeros([n_layer + 5, 4, iteration])

            starting_model[:, :, 0] = np.array(
                [
                    [1000.0, 2400.0, 4700.0, 2000.0],
                    [1000.0, 2400.0, 4700.0, 2000.0],
                    [1000.0, 2760.0, 6300.0, 2500.0],
                    [1000.0, 2760.0, 6300.0, 2500.0],
                    [1500.0, 2960.0, 6900.0, 3200.0],
                    [1500.0, 2960.0, 6900.0, 3200.0],
                    [1500.0, 2960.0, 6900.0, 3800.0],
                    [1500.0, 2960.0, 6900.0, 3800.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3190.0, 7730.0, 4310.0],
                    [10000.0, 3230.0, 7830.0, 4370.0],
                ]
            )

    vs0 = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    for i in range(0, len(starting_model)):
        vs0[
            0,
            int(np.sum(starting_model[:, :, 0][:, 0][0:i])) : int(
                np.sum(starting_model[:, :, 0][:, 0][0 : i + 1])
            ),
            0,
        ] = starting_model[:, 3][i][0]

    vsi = np.zeros([1, int(np.sum(starting_model[:, 0, 0])), 1])
    return (starting_model, vs0, vsi)


def brocher_elastic(vs):
    """Vp and density from Vs with Brocher's (2005) empirical crustal relations.

    Vp (eq. 9, valid for Vs up to about 4.5 km/s) and the Nafe-Drake density
    fit (eq. 1). Input and output in SI units: Vs and Vp in m/s, density in
    kg/m**3. Intended for ``elastic_from_vs`` in the samplers.
    """
    km = np.asarray(vs, dtype=float) / 1000.0
    vp = 0.9409 + 2.0947 * km - 0.8206 * km**2 + 0.2683 * km**3 - 0.0251 * km**4
    rho = (
        1.6612 * vp
        - 0.4721 * vp**2
        + 0.0671 * vp**3
        - 0.0043 * vp**4
        + 0.000106 * vp**5
    )
    return vp * 1000.0, rho * 1000.0


def vp_from_vs(vs, p=0.25):
    """Vp in m/s from Vs in m/s and a stable isotropic Poisson ratio."""
    vs = np.asarray(vs, float)
    if (
        not np.isfinite(p)
        or not -1 < p < 0.5
        or np.any(~np.isfinite(vs))
        or np.any(vs <= 0)
    ):
        raise ValueError(
            "Vs must be finite and positive; Poisson ratio must lie in (-1, .5)"
        )
    return vs * np.sqrt((1 - p) / (0.5 - p))


def gardner_density(vs, p=0.25):
    """Gardner density in kg/m³ for Vp expressed in m/s.

    Empirical relation; it is not automatically applied by the v2 sampler.
    """
    return 1000 * 0.31 * vp_from_vs(vs, p) ** 0.25


def station_starting_model_mean(sta):
    """Mean starting Vs profile [m/s] of the station group of ``sta`` on the 100 heat-map rows."""
    start_model = np.zeros(100)
    if sta == "RR28" or sta == "RR29" or sta == "RR34":
        start_model[0:2] = 340
        start_model[2:7] = 2540
        start_model[7:16] = 3590
        start_model[16:44] = 3940
        start_model[44:100] = 4310

    elif (
        sta == "RR52"
        or sta == "RR50"
        or sta == "RR40"
        or sta == "RR38"
        or sta == "RR36"
    ):
        start_model[0:5] = 2540
        start_model[5:14] = 3590
        start_model[14:42] = 3940
        start_model[42:100] = 4310

    return start_model


def reference_models(p=0.25):
    # Possion Ratio
    """Five published Vp profiles (e.g. SWIR-64-EW) on 100 depth rows, converted to Vs with Poisson ratio ``p``."""
    ref_models = np.zeros([5, 100])

    # SWIR-64-EW_profile
    ref_models[0][0:2] = 3500
    ref_models[0][2:3] = 4070
    ref_models[0][3:5] = 4730
    ref_models[0][5:7] = 5200
    ref_models[0][7:10] = 5900
    ref_models[0][10:13] = 6200
    ref_models[0][13:15] = 6600
    ref_models[0][15:17] = 6900
    ref_models[0][17:20] = 7120
    ref_models[0][20:23] = 7230
    ref_models[0][23:39] = 7600
    ref_models[0][39:100] = 7780

    # SWIR-64-NS_profile.txt
    ref_models[1][0:1] = 2938
    ref_models[1][1:4] = 4214
    ref_models[1][4:7] = 4994
    ref_models[1][7:10] = 5881
    ref_models[1][10:13] = 6414
    ref_models[1][13:19] = 7018
    ref_models[1][19:23] = 7374
    ref_models[1][23:28] = 7589
    ref_models[1][28:35] = 7770
    ref_models[1][35:100] = 7846

    # SWIR65-66E.txt
    ref_models[2][0:12] = 3680
    ref_models[2][12:34] = 6540
    ref_models[2][34:100] = 6910

    # BestFit
    ref_models[3][0:1] = 3500
    ref_models[3][1:2] = 3730
    ref_models[3][2:3] = 3910
    ref_models[3][3:5] = 4060
    ref_models[3][5:7] = 4230
    ref_models[3][7:8] = 4420
    ref_models[3][8:19] = 4630
    ref_models[3][9:10] = 4870
    ref_models[3][10:11] = 5080
    ref_models[3][11:12] = 5280
    ref_models[3][12:13] = 5460
    ref_models[3][13:14] = 5610

    ref_models[3][14:15] = 5760
    ref_models[3][15:16] = 5910
    ref_models[3][16:17] = 6060
    ref_models[3][17:18] = 6210
    ref_models[3][18:20] = 6360
    ref_models[3][20:21] = 6520
    ref_models[3][21:22] = 6670
    ref_models[3][22:23] = 6850

    ref_models[3][23:24] = 7030
    ref_models[3][24:25] = 7200
    ref_models[3][25:27] = 7370
    ref_models[3][27:28] = 7510
    ref_models[3][28:35] = 7640
    ref_models[3][35:100] = 8000

    # Atlantis Bank M uller 1997
    ref_models[4][0:1] = 5900
    ref_models[4][1:4] = 6000
    ref_models[4][4:7] = 6200
    ref_models[4][7:10] = 6300
    ref_models[4][10:13] = 6414
    ref_models[4][13:19] = 6500
    ref_models[4][19:23] = 6600
    ref_models[4][23:28] = 6789
    ref_models[4][28:30] = 6800
    ref_models[4][30:100] = 8000

    return ref_models / np.sqrt((1 - p) / (0.5 - p))


def station_starting_model(sta):
    """Starting Vs profile [m/s] of RHUM-RUM station ``sta`` on the 100 depth rows of the paper's heat maps."""
    start_model = np.zeros(100)
    if sta == "RR28":
        start_model[0:1] = 340
        start_model[1:6] = 2700
        start_model[6:16] = 3700
        start_model[16:45] = 4050
        start_model[45:99] = 4510

    elif sta == "RR29":
        start_model[0:1] = 340
        start_model[1:5] = 2700
        start_model[5:14] = 3700
        start_model[14:44] = 4050
        start_model[44:99] = 4510

    elif sta == "RR34":
        start_model[0:2] = 340
        start_model[2:7] = 2700
        start_model[7:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4500

    elif sta == "RR36":
        start_model[0:4] = 2700
        start_model[4:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4370

    elif sta == "RR38":
        start_model[0:4] = 2700
        start_model[4:15] = 3700
        start_model[15:42] = 4050
        start_model[42:99] = 4190

    elif sta == "RR40":
        start_model[0:4] = 2700
        start_model[4:13] = 3700
        start_model[13:42] = 4050
        start_model[42:100] = 4360

    elif sta == "RR50":
        start_model[0:7] = 2540
        start_model[7:14] = 3590
        start_model[14:42] = 3940
        start_model[42:100] = 4190

    elif sta == "RR52":
        start_model[0:7] = 2660
        start_model[7:14] = 3670
        start_model[14:42] = 4022
        start_model[42:100] = 4330

    else:
        # elif  sta == "A422A":
        start_model[0:13] = 2200
        start_model[13:25] = 2700
        start_model[25:42] = 3500
        start_model[42:60] = 3800
        start_model[60:100] = 4310
    return start_model
