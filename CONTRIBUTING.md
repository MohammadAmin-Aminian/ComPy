# Contributing

Install the editable development environment with `python -m pip install -e '.[dev]'`.
Before submitting changes, run pytest, Ruff checks/formatting and the package build
as documented in the README. Add a behavioral regression test for a numerical fix.

For scientific changes, state the equation, source, sign convention, physical
units and expected effect on outputs. Prefer analytical limits and independent
reference calculations to assertions that simply duplicate the implementation.
Do not update a numerical fixture merely to make a failing test pass.

Keep preprocessing separate from inversion, preserve input data, use local
random generators, avoid network or filesystem side effects at import, and
make empty quality selections explicit. Do not hide computational failures by
returning apparently valid zero-valued or unprocessed results.

Use small synthetic input for issue reports. Include Python/dependency versions,
parameters and the traceback. For real data, provide the station/response epoch
and describe access conditions without committing restricted waveforms.
