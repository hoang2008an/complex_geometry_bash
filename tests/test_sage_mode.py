import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


pytestmark = pytest.mark.skipif(shutil.which("sage") is None, reason="Sage is not installed")


def _sage_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{PROJECT_ROOT / 'src' / 'lib'}:{PROJECT_ROOT / 'src'}"
    env["DOT_SAGE"] = "/tmp/sage-dot-sage"
    return env


def _run_sage_snippet(source: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["sage", "-python", "-c", source],
        cwd=PROJECT_ROOT,
        env=_sage_env(),
        text=True,
        capture_output=True,
        check=True,
    )


def test_sage_midpoint_collinearity() -> None:
    _run_sage_snippet(
        """
from geometry_sage_engine import SageGeometryEngine
e = SageGeometryEngine(point_names=["A", "B", "M"])
e.add_point("A")
e.add_point("B")
e.midpoint("A", "B", "M")
n, d = e.constraint_conjugate_free("collinear", ["A", "B", "M"])[0]
assert n == 0
assert d != 0
"""
    )


def test_sage_unit_circle_learned_rules() -> None:
    result = _run_sage_snippet(
        """
from geometry_sage_engine import SageGeometryEngine
e = SageGeometryEngine(point_names=["A"])
e.add_point("A")
e.add_unit_circle("A")
print(e.learned_rules()["zb_A"])
"""
    )
    assert "1/z_A" in result.stdout


def test_sage_numeric_circumcenter_projection_and_centroid() -> None:
    _run_sage_snippet(
        """
from geometry_sage_engine import SageGeometryEngine
e = SageGeometryEngine(point_names=["A", "B", "C", "U", "P", "H", "G"])
e.set_point_value("A", "0", "0")
e.set_point_value("B", "2", "2")
e.set_point_value("C", "2*1j", "-2*1j")
e.circumcenter("A", "B", "C", "U")
assert e.z("U") == e.expr("1 + 1j")
assert e.zb("U") == e.expr("1 - 1j")
e.set_point_value("P", "2 + 3*1j", "2 - 3*1j")
e.project_point_to_line("P", "A", "B", "H")
assert e.z("H") == 2
assert e.zb("H") == 2
e.centroid("A", "B", "C", "G")
assert e.z("G") == e.expr("(2 + 2*1j) / 3")
"""
    )


def test_sage_fermat_fixed_algebraic_constants() -> None:
    _run_sage_snippet(
        """
from geometry_sage_engine import SageGeometryEngine
e = SageGeometryEngine(point_names=["A", "B", "C", "F1", "F2"])
e.set_point_value("A", "0", "0")
e.set_point_value("B", "1", "1")
e.set_point_value("C", "1j", "-1j")
e.add_fermat_points("A", "B", "C", "F1", "F2")
for label in ("F1", "F2"):
    assert e.z(label).denominator() != 0
    assert e.zb(label).denominator() != 0
"""
    )


def test_sage_line_circle_secant_known_endpoint() -> None:
    result = _run_sage_snippet(
        """
from geometry_sage_engine import SageGeometryEngine
e = SageGeometryEngine(point_names=["A", "F", "T", "B"])
for label in ("A", "F", "T"):
    e.add_point(label)
    e.add_unit_circle(label)
e.line_circle_intersection("A", "F", "T", "A", "B", avoid=["A"])
n_line, d_line = e.constraint_conjugate_free("collinear", ["A", "F", "B"])[0]
assert n_line == 0
assert d_line != 0
assert e.squared_distance("B", "T") - e.squared_distance("A", "T") == 0
print(e.z("B"))
"""
    )
    assert "z_F*z_T" in result.stdout


def test_sage_cli_supported_script_succeeds(tmp_path: Path) -> None:
    script = tmp_path / "sage_supported.json"
    script.write_text(
        json.dumps(
            {
                "steps": [
                    {"op": "add_point", "name": "A", "unit_circle": True},
                    {"op": "add_point", "name": "B", "unit_circle": True},
                    {"op": "midpoint", "args": ["A", "B", "M"]},
                    {"op": "constraint_check", "constraint": "collinear", "args": ["A", "B", "M"]},
                    {"op": "learned_rules"},
                ]
            }
        ),
        encoding="utf-8",
    )
    result = subprocess.run(
        ["sage", "-python", "-m", "geometry_sage_cli", "run", str(script)],
        cwd=PROJECT_ROOT,
        env=_sage_env(),
        text=True,
        capture_output=True,
        check=True,
    )
    assert "N = 0" in result.stdout
    assert "zb_A -> 1/z_A" in result.stdout


def test_sage_cli_qq_supports_p4_json() -> None:
    result = subprocess.run(
        [
            "sage",
            "-python",
            "-m",
            "geometry_sage_cli",
            "run",
            "--field",
            "QQ",
            str(PROJECT_ROOT / "scripts" / "P4_23_10_25.json"),
        ],
        cwd=PROJECT_ROOT,
        env=_sage_env(),
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    assert "N = 0" in result.stdout
    assert "D = 1" in result.stdout


@pytest.mark.parametrize("op", ["circle_intersection"])
def test_sage_cli_branch_ops_fail_clearly(tmp_path: Path, op: str) -> None:
    script = tmp_path / "sage_unsupported.json"
    script.write_text(
        json.dumps(
            {
                "steps": [
                    {"op": "add_point", "name": "A"},
                    {"op": "add_point", "name": "B"},
                    {"op": "add_point", "name": "C"},
                    {"op": "add_point", "name": "D"},
                    {"op": op, "args": ["A", "B", "C", "D", "X"]},
                ]
            }
        ),
        encoding="utf-8",
    )
    result = subprocess.run(
        ["sage", "-python", "-m", "geometry_sage_cli", "run", str(script)],
        cwd=PROJECT_ROOT,
        env=_sage_env(),
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode != 0
    assert "square-root branch handling" in result.stderr or "square-root branch handling" in result.stdout
