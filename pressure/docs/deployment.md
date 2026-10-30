# Deploy the reviewed BRUIT-FM workflow

BRUIT-FM is a Python research/processing workflow, not a web service. Deployment
means installing a reviewed checkout/package and configuring the owner's local
archive. Software readiness is not absolute calibration acceptance, a complete
all-data campaign, approved station ratings or a finished ocean-noise model.

## Supported setup

Use Python **3.13** on Linux or macOS. The existing locking implementation uses
POSIX `fcntl`; Windows is not a verified deployment target. The scientific and
optional Rayleigh/dispersion dependencies are declared in `pyproject.toml`.

From a clone of this repository, create a dedicated environment (do not replace
an existing environment or user directory):

```bash
python3.13 -m venv .venv
.venv/bin/python -m pip install -e '.[dev,rayleigh,dispersion]'
.venv/bin/python -m bruitfm.network --help
.venv/bin/python -m bruitfm.coverage --help
```

Editable installation is intentional for source-checkout scripts and their
project-relative docs/configuration. A built wheel is also checked in CI for
package imports, bundled reference data and nine module help entry points.
The wheel alone does **not** ship repository scripts, private archives or saved
campaigns. Do not assume those are present after a package-only installation.

## Deployment checks (no real-data campaign)

```bash
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m pytest -q
```

CI additionally builds a wheel and runs
`scripts/check_wheel_install.py WHEEL_PATH` in a disposable environment, from
outside the checkout using Python isolated mode. The check reuses the installed
scientific dependencies, verifies it imports the wheel rather than an editable
checkout, checks reference files, and invokes `--help` only. It does not read
waveforms or start PPSD/calibration jobs. This is an installation smoke test,
not a platform-wide scientific validation or a fully locked dependency build.

## Before using the real archive on a new computer

- Raw data, inventories, masters, cached catalogues and full saved results are
  intentionally not in Git. Copy them separately without overwriting originals;
  preserve filenames, exact channel/epoch identities and source provenance.
- Read `docs/full_data_coverage.md` and `docs/full_data_coverage_handoff.md` before
  planning missing work. Absolute paths and size/mtime provenance may need
  explicit reconciliation; never copy task completion flags blindly.
- Use copied/configured source definitions for the new machine. Keep historical
  result pointers/hashes intact rather than editing saved JSONs to make them fit.
- Set `BRUITFM_PYTHON` to the new environment's Python when using shell wrappers;
  their historical default points at the owner's macOS environment. Set
  `BRUITFM_SOURCE_ROOT` to the configured source checkout as needed.
- A plan/dry run writes only derived plans/cache and executes no data campaign.
  Inspect missing inputs and blocked epochs before deliberately requesting any
  `--execute` stage. Deployment checks do not authorize such execution.

## Preserved scientific limits

No gain is applied to raw data, response files, saved PPSDs or masters by this
readiness work. Pole reflection and weighted fitting stay opt-in experiments;
local-P pre-multiple values stay diagnostic; microseism gains stay conditional.
No new coherence, phase-family, Godin-frequency, nearby-normalization or station
acceptance policy was adopted. Existing scores and all frozen gains are retained.

Fresh phase summaries use the versioned circular L1 median instead of a linear
median of wrapped angles. This corrects reporting at the +/-180-degree boundary;
it does not alter bin-level phase/coherence QC or fit objectives. Historical
tables keep their original values and are identified as legacy summaries when
reused. Native resampling now refuses negative fractional record-boundary tears
rather than bridging them with interpolation; do not rerun completed campaigns
merely because that code fingerprint changed.

See [branch dispositions and verification](reviews/deployment_readiness_20261021.md).
