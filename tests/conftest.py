import sys
from pathlib import Path


def _ensure_geometry_path() -> None:
    project_root = Path(__file__).resolve().parents[1]
    lib_path = project_root / "src" / "lib"
    if str(lib_path) not in sys.path:
        sys.path.insert(0, str(lib_path))


_ensure_geometry_path()
