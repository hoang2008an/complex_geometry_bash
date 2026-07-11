import cmath

import pytest
import sympy as sp

from geometry_engine import GeometryEngine, GeometryError


def _assign_value(engine: GeometryEngine, name: str, value: complex) -> None:
    engine.add_point(name)
    symp_value = sp.nsimplify(value)
    engine.value_subs[engine.z_symbol(name)] = symp_value
    engine.value_subs[engine.zb_symbol(name)] = sp.conjugate(symp_value)


def _isogonal_direction(a: complex, b: complex, c: complex, d: complex) -> complex:
    base = (b - a) * (c - a)
    if base == 0:
        raise ValueError("Triangle vertices must be non-collinear.")
    d_vec = d - a
    if d_vec == 0:
        raise ValueError("Point D must be distinct from A.")
    const = base / base.conjugate()
    k = const * d_vec.conjugate() / d_vec
    angle = cmath.phase(k)
    return cmath.exp(0.5j * angle)


def _circumcenter(a: complex, b: complex, c: complex) -> complex:
    ax, ay = a.real, a.imag
    bx, by = b.real, b.imag
    cx, cy = c.real, c.imag
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-12:
        raise ValueError("Triangle is degenerate; circumcenter undefined.")
    ax2 = ax * ax + ay * ay
    bx2 = bx * bx + by * by
    cx2 = cx * cx + cy * cy
    ux = (ax2 * (by - cy) + bx2 * (cy - ay) + cx2 * (ay - by)) / d
    uy = (ax2 * (cx - bx) + bx2 * (ax - cx) + cx2 * (bx - ax)) / d
    return complex(ux, uy)


def _orthocenter(a: complex, b: complex, c: complex, o: complex) -> complex:
    return a + b + c - 2 * o


def _assert_complex_close(expr: sp.Expr, tol: float = 1e-12) -> None:
    numeric = complex(sp.N(expr, 20))
    assert abs(numeric) < tol


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


def test_isogonal_reflection_constraint() -> None:
    engine = GeometryEngine()
    a_val = 0.25 + 0.6j
    b_val = -0.9 + 0.1j
    c_val = 0.4 - 0.8j
    d_val = 0.3 + 0.2j

    direction = _isogonal_direction(a_val, b_val, c_val, d_val)
    e_val = a_val + 1.7 * direction

    for label, value in (("A", a_val), ("B", b_val), ("C", c_val), ("D", d_val), ("E", e_val)):
        _assign_value(engine, label, value)

    results = engine.constraint_conjugate_free("isogonal_reflection", ["A", "B", "C", "D", "E"])
    assert len(results) == 1
    numerator, denominator = results[0]
    ratio = sp.simplify(numerator / denominator)
    if ratio != 0:
        magnitude = float(sp.Abs(ratio.evalf(20)))
        assert magnitude < 1e-12
    else:
        assert ratio == 0
    assert denominator != 0


def test_add_isogonal_reflection_adds_constraint() -> None:
    engine = GeometryEngine()
    a_val = -0.1 + 0.5j
    b_val = 0.9 + 0.3j
    c_val = -0.7 - 0.2j
    d_val = 0.2 + 0.6j
    direction = _isogonal_direction(a_val, b_val, c_val, d_val)
    e_val = a_val + 0.8 * direction

    for label, value in (("A", a_val), ("B", b_val), ("C", c_val), ("D", d_val), ("E", e_val)):
        _assign_value(engine, label, value)

    initial_count = len(engine.constraints)
    engine.add_isogonal_reflection("A", "B", "C", "D", "E")
    assert len(engine.constraints) == initial_count + 1

    stored = engine.constraints[-1]
    expected = engine.isogonal_reflection_poly("A", "B", "C", "D", "E")
    assert sp.simplify(stored - expected) == 0


def test_calculate_point_resolves_collinear_isogonal_point() -> None:
    engine = GeometryEngine()
    for label in ("A", "B", "C", "E", "F"):
        engine.add_point(label)
    for label in ("A", "B", "C", "E"):
        engine.add_unit_circle(label)

    engine.add_collinear("B", "C", "F")
    engine.add_isogonal_reflection("A", "B", "C", "F", "E")
    z_f, zb_f = engine.calculate_point("F")

    assert engine._has_assignment("F")
    assert engine.z_symbol("F") not in z_f.free_symbols
    assert engine.zb_symbol("F") not in z_f.free_symbols
    assert engine.z_symbol("F") not in zb_f.free_symbols
    assert engine.zb_symbol("F") not in zb_f.free_symbols


def test_isogonal_reflection_invalid_vertices() -> None:
    engine = GeometryEngine()
    assignments = {
        "A": 0.1 + 0.2j,
        "B": 0.4 + 0.5j,
        "C": -0.3 + 0.7j,
        "D": 0.6 - 0.1j,
        "E": -0.2 + 0.9j,
    }
    for label, value in assignments.items():
        _assign_value(engine, label, value)

    with pytest.raises(GeometryError):
        engine.isogonal_reflection_poly("A", "A", "C", "D", "E")
    assert sp.simplify(engine.isogonal_reflection_poly("A", "B", "C", "A", "E")) == 0
    assert sp.simplify(engine.isogonal_reflection_poly("A", "B", "C", "D", "A")) == 0


def test_isogonal_conjugate_known_pairs() -> None:
    engine = GeometryEngine()
    assignments = {
        "A": -0.6 + 0.9j,
        "B": 0.8 + 0.3j,
        "C": -0.1 - 0.7j,
    }
    for label, value in assignments.items():
        _assign_value(engine, label, value)

    engine.centroid("A", "B", "C", "G")
    engine.lemoine_point("A", "B", "C", "L")
    engine.isogonal_conjugate_point("A", "B", "C", "G", "K")

    _assert_complex_close(engine.z("K") - engine.z("L"))
    _assert_complex_close(engine.zb("K") - engine.zb("L"))

    o_val = _circumcenter(
        assignments["A"],
        assignments["B"],
        assignments["C"],
    )
    h_val = _orthocenter(
        assignments["A"],
        assignments["B"],
        assignments["C"],
        o_val,
    )
    _assign_value(engine, "O", o_val)
    _assign_value(engine, "H", h_val)
    engine.isogonal_conjugate_point("A", "B", "C", "O", "H_iso")

    _assert_complex_close(engine.z("H_iso") - engine.z("H"))
    _assert_complex_close(engine.zb("H_iso") - engine.zb("H"))

    engine.isogonal_conjugate_point("A", "B", "C", "H", "O_iso")
    _assert_complex_close(engine.z("O_iso") - engine.z("O"))
    _assert_complex_close(engine.zb("O_iso") - engine.zb("O"))


def test_add_isogonal_conjugate_constraints() -> None:
    engine = GeometryEngine()
    assignments = {
        "A": 0.2 + 0.1j,
        "B": -0.7 + 0.4j,
        "C": 0.5 - 0.8j,
        "P": -0.3 + 0.6j,
        "Q": 0.9 + 0.2j,
    }
    for label, value in assignments.items():
        _assign_value(engine, label, value)

    initial = len(engine.constraints)
    engine.add_isogonal_conjugate("A", "B", "C", "P", "Q")
    assert len(engine.constraints) == initial + 3

    expected = engine.isogonal_conjugate_polys("A", "B", "C", "P", "Q")
    stored = engine.constraints[-3:]
    for stored_eq, expected_eq in zip(stored, expected):
        assert sp.simplify(stored_eq - expected_eq) == 0


def test_isogonal_conjugate_p4_problem() -> None:
    engine = GeometryEngine()

    for label in ("A", "B", "C"):
        engine.add_point(label)
        engine.add_unit_circle(label)

    engine.set_main_unit_triangle("A", "B", "C")

    root_assignments = {
        "__unit_root_A": cmath.exp(1j * 0.37),
        "__unit_root_B": cmath.exp(1j * 1.28),
        "__unit_root_C": cmath.exp(1j * 2.49),
    }
    for label, value in root_assignments.items():
        _assign_value(engine, label, value)

    engine.main_triangle_incenter("I")

    engine.intersection_of_lines("A", "I", "B", "C", "A_prime")
    engine.intersection_of_lines("B", "I", "C", "A", "B_prime")
    engine.intersection_of_lines("C", "I", "A", "B", "C_prime")

    engine.isogonal_conjugate_point("I", "B_prime", "C_prime", "A", "A_star")
    engine.isogonal_conjugate_point("I", "C_prime", "A_prime", "B", "B_star")
    engine.isogonal_conjugate_point("I", "A_prime", "B_prime", "C", "C_star")

    engine.isogonal_conjugate_point("A_star", "B_star", "C_star", "I", "I_star")

    results = engine.constraint_conjugate_free("collinear", ["I_star", "O", "I"])
    assert len(results) == 1
    numerator, denominator = results[0]
    _assert_complex_close(engine._apply_all(numerator))
    den_eval = engine._apply_all(denominator)
    assert abs(complex(sp.N(den_eval, 20))) > 1e-8
