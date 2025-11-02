import pytest

from geometry_engine import GeometryEngine
from geometry_cli import _demo_steps, run_steps


@pytest.mark.usefixtures("capsys")
def test_demo_midpoint_on_euler_circle(capsys):
    engine = GeometryEngine()
    steps = _demo_steps()

    run_steps(engine, steps)
    captured = capsys.readouterr().out

    assert r"|M_XP - (A+B+C)/2|^2: \frac{1}{4}" in captured
    assert "N = 0" in captured
