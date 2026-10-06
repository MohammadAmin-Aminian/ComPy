# Contributing to ComPy

Contributions that improve correctness, reproducibility, documentation, or portability are welcome.

## Development workflow

1. Create a virtual environment and install the project with `python -m pip install -e ".[dev]"`.
2. Make changes on a branch.
3. Run `pytest` before opening a pull request.
4. For changes to scientific algorithms, include a regression or analytical test and explain the physical/statistical justification.
5. Do not change empirical quality-control thresholds or station-specific assumptions without documenting the reason and expected scientific impact.

## Scientific changes

ComPy accompanies a published methodology. A code change can be software-correct while changing scientific results. Pull requests that affect compliance, calibration, inversion, uncertainty, or forward modelling should therefore report whether numerical outputs change and, where possible, compare against a known station or synthetic case.
