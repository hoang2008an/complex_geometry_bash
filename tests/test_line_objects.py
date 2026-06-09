import sympy as sp

from geometry_engine import GeometryEngine


def _assign_value(engine: GeometryEngine, name: str, value: sp.Expr) -> None:
    engine.add_point(name)
    symp_value = sp.sympify(value)
    engine.value_subs[engine.z_symbol(name)] = symp_value
    engine.value_subs[engine.zb_symbol(name)] = sp.conjugate(symp_value)


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


def test_registered_line_name_can_be_used_anywhere_line_object_is_expected() -> None:
    engine = GeometryEngine()
    for name, value in {"A": 0, "B": 2, "C": 1 + sp.I, "D": 1 - sp.I}.items():
        _assign_value(engine, name, value)

    line_ab = engine.register_line("AB", engine.line_through_points("A", "B"))
    line_cd = engine.register_line("CD", engine.line_through_points("C", "D"))

    engine.line_intersection("AB", line_cd, "X")

    assert sp.simplify(engine.line_value("AB", "X")) == 0
    assert sp.simplify(engine.line_value(line_cd, "X")) == 0


def test_circle_objects_support_named_and_direct_object_usage() -> None:
    engine = GeometryEngine()
    for name, value in {"A": 1, "B": -1, "C": sp.I, "D": -sp.I, "P": 2}.items():
        _assign_value(engine, name, value)

    omega = engine.circle_from_three_points("omega", "A", "B", "C")

    assert sp.simplify(engine._apply_all(omega.evaluate(engine, engine.z("D"), engine.zb("D")))) == 0

    engine.add_point_on_circle("omega", "D")
    engine.tangent_lines_from_point_to_circle("P", omega, ["t1", "t2"], tangent_point_names=["T1", "T2"])

    for tangent_name, tangent_point in [("t1", "T1"), ("t2", "T2")]:
        assert sp.simplify(engine.line_value(tangent_name, "P")) == 0
        assert sp.simplify(engine.line_value(tangent_name, tangent_point)) == 0
        assert sp.simplify(engine._apply_all(omega.evaluate(engine, engine.z(tangent_point), engine.zb(tangent_point)))) == 0


def test_circle_and_line_object_intersection_methods_accept_names_and_objects() -> None:
    engine = GeometryEngine()
    for name, value in {"A": 1, "B": -1, "C": sp.I, "E": 2 * sp.I, "P1": -2, "P2": 2}.items():
        _assign_value(engine, name, value)

    omega = engine.circle_from_three_points("omega", "A", "B", "C")
    axis = engine.register_line("axis", engine.line_through_points("P1", "P2"))

    engine.line_circle_object_intersections("axis", omega, ["X1", "X2"])
    assert sp.simplify(engine.line_value(axis, "X1")) == 0
    assert sp.simplify(engine.line_value("axis", "X2")) == 0
    assert sp.simplify(engine._apply_all(omega.evaluate(engine, engine.z("X1"), engine.zb("X1")))) == 0
    assert sp.simplify(engine._apply_all(omega.evaluate(engine, engine.z("X2"), engine.zb("X2")))) == 0

    omega_2 = engine.circle_from_three_points("omega2", "A", "B", "E")
    engine.circle_object_intersections("omega", omega_2, ["Y1", "Y2"])
    for name in ["Y1", "Y2"]:
        assert sp.simplify(engine._apply_all(omega.evaluate(engine, engine.z(name), engine.zb(name)))) == 0
        assert sp.simplify(engine._apply_all(omega_2.evaluate(engine, engine.z(name), engine.zb(name)))) == 0

    radical = engine.radical_line_from_circles("omega", omega_2, "rad")
    assert sp.simplify(engine.line_value(radical, "A")) == 0
    assert sp.simplify(engine.line_value("rad", "B")) == 0
