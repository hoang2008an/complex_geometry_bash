import os
import shutil
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


pytestmark = pytest.mark.skipif(shutil.which("sage") is None, reason="Sage is not installed")


def _sage_env() -> dict[str, str]:
    env = os.environ.copy()
    env["DOT_SAGE"] = "/tmp/sage-dot-sage"
    return env


def test_symbolic_ring_perp_bisector_oh_closed_form() -> None:
    result = subprocess.run(
        [
            "sage",
            "-python",
            "-c",
            """
from sage.all import SR, factor

zA, zE, zF, zT = SR.var("z_A z_E z_F z_T")
unit_conjugates = {zA: 1 / zA, zE: 1 / zE, zF: 1 / zF, zT: 1 / zT}

def conjugate(expr):
    return expr.subs(unit_conjugates)

H = (-zA * zE - zA * zF + zA * zT + zE * zT + zF * zT) / zT
M_OH = (H + zT) / 2
U_OEF = zE + zF

perpendicular = (
    (M_OH - U_OEF) * (conjugate(zT) - conjugate(H))
    + (conjugate(M_OH) - conjugate(U_OEF)) * (zT - H)
)

assert factor(perpendicular.numerator()) == 0
print("symbolic-ring N = 0")
""",
        ],
        cwd=PROJECT_ROOT,
        env=_sage_env(),
        text=True,
        capture_output=True,
        check=True,
    )
    assert "symbolic-ring N = 0" in result.stdout


def test_sr_engine_perp_bisector_oh_step_by_step() -> None:
    result = subprocess.run(
        [
            "sage",
            "-python",
            "-c",
            """
from geometry_sage_sr_engine import SageSRGeometryEngine

e = SageSRGeometryEngine()
for label in ("A", "T", "E", "F"):
    e.add_point(label)
    e.add_unit_circle(label)

e.line_circle_intersection("A", "F", "T", "A", "B", avoid=["A"])
e.line_circle_intersection("A", "E", "T", "A", "C", avoid=["A"])
e.orthocenter_via_altitudes("A", "B", "C", "H")
e.midpoint("T", "H", "M_OH")
e.reflect_point_over_line("T", "E", "F", "O_ref")
e.circumcenter("O_ref", "E", "F", "U_OEF")

n, d = e.constraint_conjugate_free("perpendicular", ["M_OH", "U_OEF", "T", "H"])[0]
assert n == 0
assert d != 0
print(e.z("H"))
""",
        ],
        cwd=PROJECT_ROOT,
        env={**_sage_env(), "PYTHONPATH": f"{PROJECT_ROOT / 'src' / 'lib'}:{PROJECT_ROOT / 'src'}"},
        text=True,
        capture_output=True,
        check=True,
    )
    assert "z_A" in result.stdout


def test_sr_engine_symbolic_sixty_degree_angle() -> None:
    result = subprocess.run(
        [
            "sage",
            "-python",
            "-c",
            """
from sage.all import pi
from geometry_sage_sr_engine import SageSRGeometryEngine

e = SageSRGeometryEngine()
for label in ("A", "B", "C"):
    e.add_point(label)
e.set_point_value("B", "0", "0")
e.set_point_value("A", "exp(I*pi/3)", "exp(-I*pi/3)")
e.set_point_value("C", "1", "1")

n, d = e.constraint_conjugate_free("angle", ["A", "B", "C"], angle=pi / 3)[0]
assert n == 0
assert d != 0
print("SR 60 degree N = 0")
""",
        ],
        cwd=PROJECT_ROOT,
        env={**_sage_env(), "PYTHONPATH": f"{PROJECT_ROOT / 'src' / 'lib'}:{PROJECT_ROOT / 'src'}"},
        text=True,
        capture_output=True,
        check=True,
    )
    assert "SR 60 degree N = 0" in result.stdout
