import sympy as sp

from geometry_engine import GeometryEngine


def assign_point(engine: GeometryEngine, name: str, value: complex) -> None:
    engine.add_point(name)
    z_expr = sp.nsimplify(value)
    zb_expr = sp.conjugate(z_expr)
    engine._set_point_assignment(name, z_expr, zb_expr)


def test_triangle_similarity_directed_polys_zero():
    engine = GeometryEngine()
    assign_point(engine, "A", 0)
    assign_point(engine, "B", 1)
    assign_point(engine, "C", 1 + 2 * sp.I)

    assign_point(engine, "D", 5)
    assign_point(engine, "E", 7)
    assign_point(engine, "F", 7 + 4 * sp.I)

    polys = engine.triangle_similarity_polys("A", "B", "C", "D", "E", "F", directed=True)
    for poly in polys:
        assert sp.simplify(engine._apply_all(poly)) == 0


def test_triangle_similarity_undirected_detects_orientation_flip():
    engine = GeometryEngine()
    assign_point(engine, "A", 0)
    assign_point(engine, "B", 1)
    assign_point(engine, "C", 1 + 2 * sp.I)

    assign_point(engine, "D", 5)
    assign_point(engine, "E", 7)
    assign_point(engine, "F", 7 - 4 * sp.I)

    directed_polys = engine.triangle_similarity_polys("A", "B", "C", "D", "E", "F", directed=True)
    evaluated = [sp.simplify(engine._apply_all(poly)) for poly in directed_polys]
    assert any(expr != 0 for expr in evaluated)

    mirror_polys = engine.triangle_similarity_polys("A", "B", "C", "D", "E", "F", directed=False)
    for poly in mirror_polys:
        assert sp.simplify(engine._apply_all(poly)) == 0


def test_add_triangle_similarity_appends_constraints():
    engine = GeometryEngine()
    for label in ["A", "B", "C", "D", "E", "F"]:
        engine.add_point(label)
    initial = len(engine.constraints)
    engine.add_triangle_similarity("A", "B", "C", "D", "E", "F", directed=True)
    assert len(engine.constraints) == initial + 2


def test_triangle_congruence_polys_directed():
    engine = GeometryEngine()
    assign_point(engine, "A", 0)
    assign_point(engine, "B", 1)
    assign_point(engine, "C", 1 + 2 * sp.I)

    assign_point(engine, "D", 5)
    assign_point(engine, "E", 6)
    assign_point(engine, "F", 6 + 2 * sp.I)

    polys = engine.triangle_congruence_polys("A", "B", "C", "D", "E", "F", directed=True)
    for poly in polys:
        assert sp.simplify(engine._apply_all(poly)) == 0


def test_triangle_congruence_polys_undirected():
    engine = GeometryEngine()
    assign_point(engine, "A", 0)
    assign_point(engine, "B", 1)
    assign_point(engine, "C", 1 + 2 * sp.I)

    assign_point(engine, "D", 5)
    assign_point(engine, "E", 6)
    assign_point(engine, "F", 6 - 2 * sp.I)

    directed_polys = engine.triangle_congruence_polys("A", "B", "C", "D", "E", "F", directed=True)
    evaluated = [sp.simplify(engine._apply_all(poly)) for poly in directed_polys]
    assert any(expr != 0 for expr in evaluated)

    undirected_polys = engine.triangle_congruence_polys("A", "B", "C", "D", "E", "F", directed=False)
    for poly in undirected_polys:
        assert sp.simplify(engine._apply_all(poly)) == 0


def test_add_triangle_congruence_appends_constraints():
    engine = GeometryEngine()
    for label in ["A", "B", "C", "D", "E", "F"]:
        engine.add_point(label)
    initial = len(engine.constraints)
    engine.add_triangle_congruence("A", "B", "C", "D", "E", "F", directed=True)
    assert len(engine.constraints) == initial + 3
