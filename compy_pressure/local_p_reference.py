"""Independent finite-difference reference for local-P validation tests.

A 1-D staggered-grid acoustic simulation of a water layer (free surface on top)
over a fluid half-space, driven by a pulse inside the half-space. A plane wave
with horizontal slowness q reduces exactly to 1-D with vertical modulus K/eta^2
in each layer (eta = sqrt(1 - q^2 c^2)). The code never evaluates the closed
forms of local_p_physics, so agreement is a check rather than a tautology.
Pressure is compression-positive; returned vertical velocity is upward.
"""

import numpy as np


def simulate_plane_wave(
    *,
    depth,
    water_speed=1500.0,
    water_density=1028.0,
    solid_speed=2500.0,
    solid_density=2000.0,
    q=0.0,
    dz=2.5,
    duration=20.0,
    peak_frequency=3.0,
    source_below_seafloor=600.0,
    sponge=800.0,
    courant=0.5,
):
    """Return t, seafloor pressure (Pa), upward velocity (m/s), first-arrival time."""
    eta_w = np.sqrt(1 - (q * water_speed) ** 2)
    eta_s = np.sqrt(1 - (q * solid_speed) ** 2)
    if not (np.isfinite(eta_w) and np.isfinite(eta_s)):
        raise ValueError("Slowness makes a layer evanescent")
    n_water = round(depth / dz)
    if not np.isclose(n_water * dz, depth):
        raise ValueError("Depth must be a multiple of dz")
    n_source = n_water + round(source_below_seafloor / dz)
    n_sponge = round(sponge / dz)
    n = n_source + round(200 / dz) + n_sponge + 1  # pressure nodes, z down from 0
    z = np.arange(n) * dz
    # Pressure node properties: water for z < depth, solid at and below it.
    in_water = z < depth - 1e-9
    modulus = np.where(
        in_water,
        water_density * water_speed**2 / eta_w**2,
        solid_density * solid_speed**2 / eta_s**2,
    )
    # Velocity nodes sit between pressure nodes i and i+1.
    rho_v = np.where(z[:-1] + dz / 2 < depth, water_density, solid_density)
    vertical_speed = max(water_speed / eta_w, solid_speed / eta_s)
    dt = courant * dz / vertical_speed
    steps = int(np.ceil(duration / dt))
    damping = np.zeros(n)
    sponge_start = n - n_sponge
    ramp = (np.arange(n) - sponge_start).clip(0) / max(n_sponge, 1)
    damping = np.exp(-0.15 * ramp**2)

    p = np.zeros(n)
    v = np.zeros(n - 1)  # positive downward
    t0 = 1.5 / peak_frequency
    record_p = np.empty(steps)
    record_v = np.empty(steps)
    for step in range(steps):
        time = step * dt
        v -= dt / (rho_v * dz) * (p[1:] - p[:-1])
        v *= damping[:-1]
        p[1:-1] -= dt * modulus[1:-1] / dz * (v[1:] - v[:-1])
        arg = (np.pi * peak_frequency * (time - t0)) ** 2
        p[n_source] += (1 - 2 * arg) * np.exp(-arg)  # Ricker pressure source
        p[0] = 0.0  # free surface
        p *= damping
        record_p[step] = p[n_water]
        # Normal velocity is continuous; average the two adjacent velocity nodes.
        record_v[step] = -0.5 * (v[n_water - 1] + v[n_water])
    t = np.arange(steps) * dt
    arrival = t0 + source_below_seafloor / (solid_speed / eta_s)
    return t, record_p, record_v, arrival
