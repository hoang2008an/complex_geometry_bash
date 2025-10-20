import sympy as sp

from geometry_engine import GeometryEngine


def _assign_value(engine: GeometryEngine, name: str, value: complex) -> None:
    engine.add_point(name)
    symp_value = sp.nsimplify(value)
    engine.value_subs[engine.z_symbol(name)] = symp_value
    engine.value_subs[engine.zb_symbol(name)] = sp.conjugate(symp_value)


def test_point_reflection_constraint() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "O", 1 + 2 * sp.I)
    _assign_value(engine, "P", 3 - sp.I)

    expected_q = 2 * engine.z("O") - engine.z("P")
    _assign_value(engine, "Q", sp.simplify(expected_q))

    results = engine.constraint_conjugate_free("point_reflection", ["P", "O", "Q"])
    assert len(results) == 2
    for numerator, denominator in results:
        assert sp.simplify(numerator) == 0
        assert denominator != 0


def test_line_reflection_constraint() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 0)
    _assign_value(engine, "B", 4)
    _assign_value(engine, "P", 1 + 3 * sp.I)
    _assign_value(engine, "Q", 1 - 3 * sp.I)
    midpoint_label = "__mid_{}_{}".format(*sorted(("P", "Q")))
    midpoint_value = sp.simplify((engine.z("P") + engine.z("Q")) / 2)
    _assign_value(engine, midpoint_label, midpoint_value)

    results = engine.constraint_conjugate_free("line_reflection", ["P", "A", "B", "Q"])
    assert len(results) == 4
    for numerator, denominator in results:
        assert sp.simplify(numerator) == 0
        assert denominator != 0


def test_reflect_point_over_point_construction() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "O", -1 + sp.I)
    _assign_value(engine, "P", 2 + 3 * sp.I)

    engine.reflect_point_over_point("P", "O", "Q")

    expected = sp.simplify(2 * engine.z("O") - engine.z("P"))
    assert sp.simplify(engine.z("Q") - expected) == 0
    assert sp.simplify(engine.zb("Q") - sp.conjugate(expected)) == 0


def test_reflect_point_over_line_construction() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 0)
    _assign_value(engine, "B", 2 + sp.I)
    _assign_value(engine, "P", 3 + 2 * sp.I)

    engine.reflect_point_over_line("P", "A", "B", "Q")

    results = engine.constraint_conjugate_free("line_reflection", ["P", "A", "B", "Q"])
    assert len(results) == 4
    for numerator, denominator in results:
        assert sp.simplify(numerator) == 0
        assert denominator != 0


def test_midpoint_constraint() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "P", 1 + sp.I)
    _assign_value(engine, "Q", 3 + 5 * sp.I)
    midpoint_value = sp.simplify((engine.z("P") + engine.z("Q")) / 2)
    _assign_value(engine, "M", midpoint_value)

    results = engine.constraint_conjugate_free("midpoint", ["P", "Q", "M"])
    assert len(results) == 2
    for numerator, denominator in results:
        assert sp.simplify(numerator) == 0
        assert denominator != 0


def test_midpoint_construction() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "P", -2 + 4 * sp.I)
    _assign_value(engine, "Q", 6 - 2 * sp.I)

    engine.midpoint("P", "Q", "M")

    expected = sp.simplify((engine.z("P") + engine.z("Q")) / 2)
    assert sp.simplify(engine.z("M") - expected) == 0
    assert sp.simplify(engine.zb("M") - sp.conjugate(expected)) == 0
