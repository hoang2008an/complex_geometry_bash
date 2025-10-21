from pathlib import Path

import pytest

from geometry_engine import GeometryEngine
from geometry_cli import load_script, run_steps


@pytest.mark.usefixtures("capsys")
def test_incenter_arc_perpendicular_demo(capsys):
    engine = GeometryEngine()
    script_path = Path(__file__).resolve().parent.parent / "scripts" / "incenter_arc_demo.json"
    steps = load_script(script_path)
    run_steps(engine, steps)

    output = capsys.readouterr().out
    assert "N = 0" in output
