"""Install a built wheel into a disposable environment and smoke-test safe CLIs.

Reuses installed scientific dependencies, not the editable project. Runs from
outside the checkout with Python isolated mode. No data/configuration is loaded,
no archive is inspected and no calibration is executed.
"""

import argparse
import os
import subprocess
import tempfile
import venv
import zipfile
from pathlib import Path

CLI_MODULES = (
    "compy_pressure.inventory",
    "compy_pressure.ppsd",
    "compy_pressure.network",
    "compy_pressure.coverage",
    "compy_pressure.quality_review",
    "compy_pressure.noise_model",
    "compy_pressure.rayleigh_calibration",
    "compy_pressure.microseism_diagnostic",
    "compy_pressure.local_p_check",
)
REQUIRED_RESOURCES = (
    "bruitfm/reference_data/pressure_references.json",
    "bruitfm/reference_data/brown_hydroacoustic_fullband.csv",
    "bruitfm/reference_data/rr52_historical_rayleigh_dispersion.csv",
)


def check_wheel(wheel):
    wheel = Path(wheel).resolve(strict=True)
    if wheel.suffix != ".whl":
        raise ValueError("Expected a built .whl file")
    with zipfile.ZipFile(wheel) as archive:
        missing = set(REQUIRED_RESOURCES) - set(archive.namelist())
    if missing:
        raise ValueError(f"Wheel lacks required reference resources: {sorted(missing)}")
    with tempfile.TemporaryDirectory(prefix="bruitfm-wheel-smoke-") as temporary:
        root = Path(temporary)
        environment = root / "venv"
        venv.EnvBuilder(with_pip=True, system_site_packages=True).create(environment)
        python = environment / (
            "Scripts/python.exe" if os.name == "nt" else "bin/python"
        )
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.update(
            MPLBACKEND="Agg",
            MPLCONFIGDIR=str(root / "mpl"),
            NUMBA_CACHE_DIR=str(root / "numba"),
            OPENBLAS_NUM_THREADS="1",
        )
        subprocess.run(
            [
                str(python),
                "-I",
                "-m",
                "pip",
                "install",
                "--no-deps",
                "--force-reinstall",
                str(wheel),
            ],
            check=True,
            cwd=root,
            env=env,
        )
        probe = (
            "import sys; from pathlib import Path; import compy_pressure; "
            "assert Path(compy_pressure.__file__).resolve().is_relative_to("
            "Path(sys.prefix).resolve()), "
            "'Smoke test imported the checkout rather than the installed wheel'; "
            "from importlib.resources import files; "
            "assert files('bruitfm').joinpath("
            "'reference_data/pressure_references.json').is_file(); "
            "from compy_pressure.phase import circular_median_degrees; "
            "assert abs(abs(circular_median_degrees([179,-179])) - 180) < 1e-9"
        )
        subprocess.run([str(python), "-I", "-c", probe], check=True, cwd=root, env=env)
        for module in CLI_MODULES:
            checked = subprocess.run(
                [str(python), "-I", "-m", module, "--help"],
                cwd=root,
                env=env,
                capture_output=True,
                text=True,
            )
            if checked.returncode:
                raise RuntimeError(f"Installed CLI failed: {module}\n{checked.stderr}")
    print(
        f"PASS: installed wheel, packaged references and {len(CLI_MODULES)} safe CLIs"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    check_wheel(parser.parse_args().wheel)


if __name__ == "__main__":
    main()
