import sympy as sp

from geometry_engine import GeometryEngine


def test_project_point_to_line_sets_coordinates():
    engine = GeometryEngine()
    engine.add_point("A")
    engine.add_point("B")
    engine.add_point("P")

    engine._set_point_assignment("A", sp.Integer(0), sp.Integer(0))
    engine._set_point_assignment("B", sp.Integer(1), sp.Integer(1))
    engine._set_point_assignment("P", sp.Integer(2) + 3 * sp.I, sp.Integer(2) - 3 * sp.I)

    engine.project_point_to_line("P", "A", "B", "H")

    assert sp.simplify(engine.z("H") - sp.Integer(2)) == 0
    assert sp.simplify(engine.zb("H") - sp.Integer(2)) == 0

    perp_eq, collinear_eq = engine.projection_to_line_polys("P", "A", "B", "H", raw=True)
    substitution_map = {
        engine.z_symbol(label): engine.z(label)
        for label in ("A", "B", "P", "H")
    }
    substitution_map.update({
        engine.zb_symbol(label): engine.zb(label)
        for label in ("A", "B", "P", "H")
    })
    assert sp.simplify(perp_eq.subs(substitution_map)) == 0
    assert sp.simplify(collinear_eq.subs(substitution_map)) == 0


def test_add_projection_to_line_registers_constraints():
    engine = GeometryEngine()
    engine.add_point("A")
    engine.add_point("B")
    engine.add_point("P")
    engine.add_point("H")

    count_before = len(engine.constraints)
    engine.add_projection_to_line("P", "A", "B", "H")
    assert len(engine.constraints) == count_before + 2
