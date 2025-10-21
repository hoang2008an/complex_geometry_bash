import sympy as sp
import pytest

from geometry_engine import GeometryEngine, GeometryError


def _prepare_unit_triangle_vertices(engine: GeometryEngine) -> None:
    for label in ["A", "B", "C"]:
        engine.add_point(label)
        engine.add_unit_circle(label)


def test_set_main_unit_triangle_assignments():
    engine = GeometryEngine()
    _prepare_unit_triangle_vertices(engine)
    engine.set_main_unit_triangle("A", "B", "C")

    root_a = engine.z_symbol("__unit_root_A")
    root_b = engine.z_symbol("__unit_root_B")
    root_c = engine.z_symbol("__unit_root_C")
    assert sp.simplify(engine.z("A") - root_a ** 2) == 0
    assert sp.simplify(engine.z("B") - root_b ** 2) == 0
    assert sp.simplify(engine.z("C") - root_c ** 2) == 0
    assert sp.simplify(engine.zb("A") - sp.simplify(root_a ** -2)) == 0
    assert sp.simplify(engine.zb("B") - sp.simplify(root_b ** -2)) == 0
    assert sp.simplify(engine.zb("C") - sp.simplify(root_c ** -2)) == 0


def test_main_triangle_incenter_and_excenters():
    engine = GeometryEngine()
    _prepare_unit_triangle_vertices(engine)
    engine.set_main_unit_triangle("A", "B", "C")
    engine.main_triangle_incenter("I")
    engine.main_triangle_excenter("A", "I_A")
    engine.main_triangle_excenter("B", "I_B")
    engine.main_triangle_excenter("C", "I_C")

    root_a = engine.z_symbol("__unit_root_A")
    root_b = engine.z_symbol("__unit_root_B")
    root_c = engine.z_symbol("__unit_root_C")
    expected_incenter = -(root_a * root_b + root_b * root_c + root_c * root_a)
    expected_incenter_conj = -(
        (root_a ** -1) * (root_b ** -1)
        + (root_b ** -1) * (root_c ** -1)
        + (root_c ** -1) * (root_a ** -1)
    )
    assert sp.simplify(engine.z("I") - expected_incenter) == 0
    assert sp.simplify(engine.zb("I") - expected_incenter_conj) == 0

    expected_A = root_a * root_b + root_b * root_c - root_c * root_a
    expected_B = root_a * root_b - root_b * root_c + root_c * root_a
    expected_C = -root_a * root_b + root_b * root_c + root_c * root_a

    expected_A_conj = (root_a ** -1) * (root_b ** -1) + (root_b ** -1) * (root_c ** -1) - (root_c ** -1) * (root_a ** -1)
    expected_B_conj = (root_a ** -1) * (root_b ** -1) - (root_b ** -1) * (root_c ** -1) + (root_c ** -1) * (root_a ** -1)
    expected_C_conj = -(root_a ** -1) * (root_b ** -1) + (root_b ** -1) * (root_c ** -1) + (root_c ** -1) * (root_a ** -1)

    assert sp.simplify(engine.z("I_A") - expected_A) == 0
    assert sp.simplify(engine.z("I_B") - expected_B) == 0
    assert sp.simplify(engine.z("I_C") - expected_C) == 0
    assert sp.simplify(engine.zb("I_A") - expected_A_conj) == 0
    assert sp.simplify(engine.zb("I_B") - expected_B_conj) == 0
    assert sp.simplify(engine.zb("I_C") - expected_C_conj) == 0


def test_main_triangle_arc_midpoints():
    engine = GeometryEngine()
    _prepare_unit_triangle_vertices(engine)
    engine.set_main_unit_triangle("A", "B", "C")
    engine.main_triangle_arc_midpoint("A", "M_A")
    engine.main_triangle_arc_midpoint("A", "M_A_star", containing_vertex=True)
    engine.main_triangle_arc_midpoint("B", "M_B")
    engine.main_triangle_arc_midpoint("C", "M_C")

    root_a = engine.z_symbol("__unit_root_A")
    root_b = engine.z_symbol("__unit_root_B")
    root_c = engine.z_symbol("__unit_root_C")
    assert sp.simplify(engine.z("M_A") + root_b * root_c) == 0
    assert sp.simplify(engine.z("M_A_star") - root_b * root_c) == 0
    assert sp.simplify(engine.z("M_B") + root_c * root_a) == 0
    assert sp.simplify(engine.z("M_C") + root_a * root_b) == 0

    assert sp.simplify(engine.z("M_A") * engine.zb("M_A") - 1) == 0
    assert sp.simplify(engine.z("M_B") * engine.zb("M_B") - 1) == 0
    assert sp.simplify(engine.z("M_C") * engine.zb("M_C") - 1) == 0


def test_main_triangle_operations_require_setup():
    engine = GeometryEngine()
    with pytest.raises(GeometryError):
        engine.main_triangle_incenter("I")


def test_set_main_unit_triangle_requires_unit_circle_vertices():
    engine = GeometryEngine()
    for label in ["A", "B", "C"]:
        engine.add_point(label)
    with pytest.raises(GeometryError):
        engine.set_main_unit_triangle("A", "B", "C")


def test_set_main_unit_triangle_with_custom_roots():
    engine = GeometryEngine()
    _prepare_unit_triangle_vertices(engine)
    engine.set_main_unit_triangle("A", "B", "C", root_names=("x", "y", "z"))

    root_a = engine.z_symbol("x")
    root_b = engine.z_symbol("y")
    root_c = engine.z_symbol("z")

    assert sp.simplify(engine.z("A") - root_a ** 2) == 0
    assert sp.simplify(engine.z("B") - root_b ** 2) == 0
    assert sp.simplify(engine.z("C") - root_c ** 2) == 0
    assert sp.simplify(engine.z("A") * engine.zb("A") - 1) == 0


def test_set_main_unit_triangle_custom_roots_validation():
    engine = GeometryEngine()
    _prepare_unit_triangle_vertices(engine)
    with pytest.raises(GeometryError):
        engine.set_main_unit_triangle("A", "B", "C", root_names=("X", "X", "Z"))
    with pytest.raises(GeometryError):
        engine.set_main_unit_triangle("A", "B", "C", root_names=("A", "Y", "Z"))
