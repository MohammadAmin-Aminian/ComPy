"""Output contracts shared by the experimental local-P command-line tools."""

from pathlib import Path


def local_p_output(path):
    """Resolve an output below LOCAL_P, never inside raw data or exports."""
    path = Path(path).expanduser().resolve()
    forbidden = (
        Path.home() / "Desktop/Data_Bruit_FM",
        Path.home() / "Projects/bruit-fm/exports",
    )
    if any(path == p or p in path.parents for p in forbidden):
        raise ValueError(f"Output lies inside protected inputs: {path}")
    if not any(
        p.name == "LOCAL_P" and p.parent.name == "results" for p in path.parents
    ):
        raise ValueError(f"Output must be below results/LOCAL_P/: {path}")
    return path


def fresh_file(path):
    """Check before expensive work; callers must also open with exclusive 'x'."""
    original = Path(path).expanduser()
    path = local_p_output(original)
    if original.is_symlink() or path.exists():
        raise FileExistsError(f"Refusing to overwrite {path}")
    return path
