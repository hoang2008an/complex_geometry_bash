import sympy as sp

from geometry_engine import Circle, GeometryEngine, Line


def _assign(engine: GeometryEngine, name: str, value: complex) -> None:
    engine.add_point(name)
    expr = sp.nsimplify(value)
    engine._set_point_assignment(name, expr, sp.conjugate(expr))


def test_line_through_points_matches_collinearity():
    engine = GeometryEngine()
    _assign(engine, "P", 0)
    _assign(engine, "Q", 1)
    _assign(engine, "R", 2)
    _assign(engine, "S", 2 + sp.I)

    line = engine.line_through_points("P", "Q")
    assert line.points == ("P", "Q")

    value_R = engine.line_value(line, "R")
    assert sp.simplify(value_R) == 0  # R lies on PQ

    value_S = engine.line_value(line, "S")
    assert sp.simplify(value_S) != 0  # S is off the line


def test_add_point_on_line_records_constraint():
    engine = GeometryEngine()
    engine.add_point("P")
    engine.add_point("Q")
    engine.add_point("X")
    line = engine.line_through_points("P", "Q")
    initial_constraints = len(engine.constraints)
    engine.add_point_on_line(line, "X")
    assert len(engine.constraints) == initial_constraints + 1


def test_circle_from_center_point_and_membership():
    engine = GeometryEngine()
    _assign(engine, "O", 0)
    _assign(engine, "A", 1)
    _assign(engine, "B", sp.exp(sp.I * sp.pi / 3))

    circle = engine.circle_from_center_point("O", "A")
    assert isinstance(circle, Circle)
    assert circle.center == "O"
    assert circle.reference_points == ("A",)

    value_A = engine.circle_value(circle, "A")
    value_B = engine.circle_value(circle, "B")
    assert sp.simplify(value_A) == 0
    assert sp.simplify(value_B) == 0


def test_add_point_on_circle_records_constraint():
    engine = GeometryEngine()
    engine.add_point("O")
    engine.add_point("A")
    circle = engine.circle_from_center_point("O", "A")
    engine.add_point("X")
    initial_constraints = len(engine.constraints)
    engine.add_point_on_circle(circle, "X")
    assert len(engine.constraints) == initial_constraints + 1
