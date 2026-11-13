"""E1: checks of the forward model that do not use ComPy's own references.

1. A global-matrix solution written here, in 50-digit arithmetic (mpmath):
   unknown up- and down-going P and SV amplitudes in every layer, continuity of
   displacement and traction at interfaces, decaying waves in the half-space,
   t_xz = 0 and t_zz = -1 at the seafloor. Compared with ``dram.Forward`` and
   ``impedance.compliance_by_impedance`` on random models, including slow sediment.
2. Quasi-static half-space: eta -> (1 - nu)/mu.
3. Convergence of a linear-gradient layer under finer slicing.
4. What the data constrain: models with the same M = mu/(1 - nu) profile but
   different Vp/Vs and density.

Run: PYTHONPATH=<worktree> python e01_forward_checks.py <outdir>
"""

import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

import compy
from compy.inversion.dram import Forward
from compy.inversion.impedance import compliance_by_impedance
from compy.numerics import wavenumber

mp.mp.dps = 50


def global_matrix_compliance(depth, freqs, model):
    """Normalized compliance from one linear system for all layer amplitudes."""
    model = np.asarray(model, float)
    n = len(model)
    out = []
    for f in freqs:
        om = mp.mpf(2) * mp.pi * mp.mpf(float(f))
        k = mp.mpf(float(wavenumber(2 * np.pi * f, depth)))
        # unknowns: finite layer j has (A_dn, A_up, B_dn, B_up); half-space (A_dn, B_dn)
        nun = 4 * (n - 1) + 2
        a = mp.zeros(nun, nun)
        b = mp.zeros(nun, 1)

        def fields(j, z_rel, which):
            """Rows (u_x, u_z, t_xz, t_zz) of mode ``which`` of layer j at depth z_rel
            below the layer top. Potentials phi = exp(-/+ r z), psi = exp(-/+ s z),
            u_x = d phi/dx - d psi/dz, u_z = d phi/dz + d psi/dx, common exp(ikx)."""
            h, rho, vp, vs = (mp.mpf(float(v)) for v in model[j])
            mu = rho * vs * vs
            lam = rho * vp * vp - 2 * mu
            r = mp.sqrt(k * k - om * om / (vp * vp))
            s = mp.sqrt(k * k - om * om / (vs * vs))
            ik = mp.mpc(0, 1) * k
            if which in ("Pd", "Pu"):
                sg = -1 if which == "Pd" else 1
                # reference depth: top for decaying, bottom for growing (no overflow)
                zr = z_rel if which == "Pd" else z_rel - h
                e = mp.exp(sg * r * zr)
                ux, uz = ik * e, sg * r * e
                dux, duz = sg * r * ux, sg * r * uz
            else:
                sg = -1 if which == "Sd" else 1
                zr = z_rel if which == "Sd" else z_rel - h
                e = mp.exp(sg * s * zr)
                ux, uz = -sg * s * e, ik * e
                dux, duz = sg * s * ux, sg * s * uz
            txz = mu * (dux + ik * uz)
            tzz = lam * (ik * ux + duz) + 2 * mu * duz
            return ux, uz, txz, tzz

        def cols(j):
            if j == n - 1:
                return [(4 * (n - 1), "Pd"), (4 * (n - 1) + 1, "Sd")]
            return [(4 * j + i, w) for i, w in enumerate(("Pd", "Pu", "Sd", "Su"))]

        row = 0
        for c, w in cols(0):  # seafloor: t_xz = 0, t_zz = -1
            fl = fields(0, mp.mpf(0), w)
            a[row, c], a[row + 1, c] = fl[2], fl[3]
        b[row + 1] = -1
        row += 2
        for j in range(n - 1):  # interfaces: continuity of the four fields
            hj = mp.mpf(float(model[j, 0]))
            for c, w in cols(j):
                fl = fields(j, hj, w)
                for i in range(4):
                    a[row + i, c] = fl[i]
            for c, w in cols(j + 1):
                fl = fields(j + 1, mp.mpf(0), w)
                for i in range(4):
                    a[row + i, c] = -fl[i]
            row += 4
        x = mp.lu_solve(a, b)
        uz = sum(x[c] * fields(0, mp.mpf(0), w)[1] for c, w in cols(0))
        # eta = k u_z / tau_zz with tau_zz = -1 and z down: u_z down is positive
        out.append(float(mp.re(k * uz)))
    return np.array(out)


def random_model(rng, sediment):
    rows = []
    if sediment:
        for _ in range(rng.integers(1, 4)):
            vs = rng.uniform(80, 900)
            rows.append(
                [
                    rng.uniform(20, 500),
                    rng.uniform(1500, 2100),
                    rng.uniform(1500, 2400),
                    vs,
                ]
            )
    z = 0
    for _ in range(rng.integers(2, 7)):
        vs = rng.uniform(1200, 4200)
        rows.append(
            [
                rng.uniform(100, 3000),
                rng.uniform(2300, 3000),
                vs * rng.uniform(1.6, 2.6),
                vs,
            ]
        )
    vs = rng.uniform(4000, 4700)
    rows.append([0, 3300, 1.8 * vs, vs])
    return np.array(rows)


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    rec = {"compy": compy.__file__}
    worst = {"dram": 0.0, "impedance": 0.0}
    sign_ok = True
    n_models = 60
    for i in range(n_models):
        depth = float(rng.choice([150, 500, 1500, 3000, 5000]))
        fc = np.sqrt(9.79329 / (2 * np.pi * depth))
        freqs = np.geomspace(0.004, min(0.9 * fc, 0.06), 5)
        m = random_model(rng, sediment=i % 2 == 0)
        ref = np.abs(global_matrix_compliance(depth, freqs, m))
        sign_ok &= bool(np.all(global_matrix_compliance(depth, freqs[:1], m) > 0))
        worst["dram"] = max(
            worst["dram"], float(np.max(np.abs(Forward(depth, freqs)(m) / ref - 1)))
        )
        worst["impedance"] = max(
            worst["impedance"],
            float(np.max(np.abs(compliance_by_impedance(depth, freqs, m) / ref - 1))),
        )
    rec["random_models"] = dict(
        n=n_models, max_rel_diff=worst, reference_positive=sign_ok
    )

    # 2. quasi-static half-space
    depth, freqs = 4000.0, np.array([0.003, 0.006, 0.012])
    hs = {}
    for vs, r in [(3000.0, 1.75), (1500.0, 2.5), (300.0, 6.0)]:
        rho = 2500.0
        m = np.array([[0, rho, r * vs, vs]])
        mu = rho * vs**2
        nu = (r**2 - 2) / (2 * (r**2 - 1))
        hs[f"vs{vs:.0f}_r{r}"] = (Forward(depth, freqs)(m) / ((1 - nu) / mu)).tolist()
    rec["halfspace_over_static"] = hs

    # 3. gradient discretisation: 2 km layer, Vs 1000 -> 3500 linearly, over a half-space
    depth, freqs = 3000.0, np.geomspace(0.005, 0.02, 8)

    def gradient(nsub):
        z = (np.arange(nsub) + 0.5) / nsub
        vs = 1000 + 2500 * z
        rows = [[2000 / nsub, 2400 + 300 * zz, 1.9 * v, v] for zz, v in zip(z, vs)]
        return np.array(rows + [[0, 2900, 6800, 3800]])

    ref = Forward(depth, freqs)(gradient(1024))
    rec["gradient_convergence"] = {
        str(n): float(np.max(np.abs(Forward(depth, freqs)(gradient(n)) / ref - 1)))
        for n in (2, 4, 8, 16, 32, 64, 128)
    }

    # 4. same M = mu/(1-nu) profile, different Vp/Vs and density
    base = np.array(
        [
            [300, 2000, 1800, 400],
            [500, 2400, 4200, 2000],
            [1500, 2700, 5800, 3100],
            [4500, 2950, 6900, 3850],
            [0, 3300, 8000, 4500],
        ],
        float,
    )

    def big_m(m):
        mu = m[:, 1] * m[:, 3] ** 2
        return 2 * mu * (1 - (m[:, 3] / m[:, 2]) ** 2)  # = mu / (1 - nu)

    def with_ratio(m, r, rho_scale):
        """Model with Vp/Vs = r and density scaled, keeping M of every layer."""
        new = m.copy()
        new[:, 1] *= rho_scale
        mu = big_m(m) / (2 * (1 - 1 / r**2))
        new[:, 3] = np.sqrt(mu / new[:, 1])
        new[:, 2] = r * new[:, 3]
        return new

    eta0 = Forward(depth, freqs)(base)
    same_m = {}
    for r in (1.7, 1.85, 2.1, 2.4, 3.0):
        for s in (0.9, 1.0, 1.1):
            alt = with_ratio(base, r, s)
            same_m[f"r{r}_rho{s}"] = dict(
                max_eta_change=float(
                    np.max(np.abs(Forward(depth, freqs)(alt) / eta0 - 1))
                ),
                vs_change_crust=float(alt[2, 3] / base[2, 3] - 1),
            )
    rec["same_M_different_vpvs_rho"] = same_m
    # and the converse: same Vs, different Vp/Vs
    conv = {}
    for r in (1.7, 1.85, 2.1, 2.4):
        alt = base.copy()
        alt[1:, 2] = r * alt[1:, 3]
        conv[f"r{r}"] = float(np.max(np.abs(Forward(depth, freqs)(alt) / eta0 - 1)))
    rec["same_vs_different_vpvs_eta_change"] = conv
    (out / "e01_forward_checks.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
