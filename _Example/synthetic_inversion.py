"""Small offline numerical example: python _Example/synthetic_inversion.py."""

import numpy as np
import inv_compy as inv


def main():
    model = np.array(
        [
            [100.0, 2200.0, 3000.0, 1500.0],
            [1000.0, 2800.0, 6000.0, 3500.0],
            [0.0, 3300.0, 8000.0, 4500.0],
        ]
    )
    f = np.array([0.005, 0.007, 0.01, 0.015])
    observed = inv.calc_norm_compliance(4000, f, model)
    result = inv.invert_compliance_beta(
        observed,
        f,
        4000,
        starting_model=model,
        s=np.full(len(f), 1e-12),
        iteration=100,
        alpha=0,
        seed=0,
        return_profiles=False,
    )
    print(f"ComPy 2.0: {result[0].shape[-1]} saved states; acceptance {result[-1]:.3f}")
    print("This short synthetic run demonstrates the API, not posterior convergence.")


if __name__ == "__main__":
    main()
