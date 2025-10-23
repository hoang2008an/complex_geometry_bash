import sympy as sp

from geometry_engine import GeometryEngine


def test_line_intersection_places_point_on_both_lines() -> None:
    engine = GeometryEngine()
    for label in ["A", "B", "C", "D"]:
        engine.add_point(label)

    line_ab = engine.line_through_points("A", "B")
    line_cd = engine.line_through_points("C", "D")

    engine.line_intersection(line_ab, line_cd, "X")

    z_x = engine.z("X")
    zb_x = engine.zb("X")

    assert sp.simplify(line_ab.evaluate(z_x, zb_x)) == 0
    assert sp.simplify(line_cd.evaluate(z_x, zb_x)) == 0


def test_symbolic_line_constraints_perpendicular() -> None:
    engine = GeometryEngine()
    engine.add_point("P")
    engine.add_point("I")
    engine.add_point("Q")

    base_line = engine.line_through_points("P", "I")
    symbolic_line = engine.symbolic_line("L")
    engine.add_point_on_line(symbolic_line, "Q")
    engine.add_perpendicular_lines(symbolic_line, base_line)

    alpha, beta, gamma = engine.line_parameters["L"]
    assert alpha in engine.value_subs
    assert beta in engine.value_subs
    assert gamma in engine.value_subs

    z_q = engine.z_symbol("Q")
    zb_q = engine.zb_symbol("Q")
    substituted = symbolic_line.alpha * z_q + symbolic_line.beta * zb_q + symbolic_line.gamma
    substituted = engine._apply_all(substituted)
    assert sp.simplify(substituted) == 0
