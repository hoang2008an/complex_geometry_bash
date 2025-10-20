import sympy as sp
import pytest

from geometry_engine import GeometryEngine, GeometryError


def _assign_value(engine: GeometryEngine, name: str, value: sp.Expr) -> None:
    engine.add_point(name)
    symp_value = sp.sympify(value)
    engine.value_subs[engine.z_symbol(name)] = symp_value
    engine.value_subs[engine.zb_symbol(name)] = sp.conjugate(symp_value)


def test_constraint_conjugate_free_collinear() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 1)
    _assign_value(engine, "B", 2)
    _assign_value(engine, "C", 3)

    results = engine.constraint_conjugate_free("collinear", ["A", "B", "C"])
    assert len(results) == 1
    numerator, denominator = results[0]
    assert sp.simplify(numerator) == 0
    assert denominator != 0


def test_constraint_conjugate_free_perpendicular() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 0)
    _assign_value(engine, "B", 1)
    _assign_value(engine, "C", 0)
    _assign_value(engine, "D", sp.I)

    results = engine.constraint_conjugate_free("perpendicular", ["A", "B", "C", "D"])
    assert len(results) == 1
    numerator, denominator = results[0]
    assert sp.simplify(numerator) == 0
    assert denominator != 0


def test_constraint_conjugate_free_concyclic() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 1)
    _assign_value(engine, "B", -1)
    _assign_value(engine, "C", sp.I)
    _assign_value(engine, "D", -sp.I)

    results = engine.constraint_conjugate_free("concyclic", ["A", "B", "C", "D"])
    assert len(results) == 1
    numerator, denominator = results[0]
    assert sp.simplify(numerator) == 0
    assert denominator != 0


def test_constraint_conjugate_free_angle_value() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "B", 0)
    _assign_value(engine, "A", 1)
    _assign_value(engine, "C", sp.I)

    results = engine.constraint_conjugate_free("angle", ["A", "B", "C"], angle=sp.pi / 2)
    assert len(results) == 1
    numerator, denominator = results[0]
    assert sp.simplify(numerator) == 0
    assert denominator != 0


def test_constraint_conjugate_free_circumcenter() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 1)
    _assign_value(engine, "B", sp.I)
    _assign_value(engine, "C", -1)
    _assign_value(engine, "U", 0)

    results = engine.constraint_conjugate_free("circumcenter", ["A", "B", "C", "U"])
    assert len(results) == 2
    for numerator, denominator in results:
        assert sp.simplify(numerator) == 0
        assert denominator != 0


def test_constraint_conjugate_free_unknown_constraint() -> None:
    engine = GeometryEngine()
    with pytest.raises(GeometryError):
        engine.constraint_conjugate_free("nonexistent", [])
