import sympy as sp

from geometry_engine import GeometryEngine


def _assign_value(engine: GeometryEngine, name: str, value: complex) -> None:
    engine.add_point(name)
    symp_value = sp.nsimplify(value)
    engine.value_subs[engine.z_symbol(name)] = symp_value
    engine.value_subs[engine.zb_symbol(name)] = sp.conjugate(symp_value)


def _vec(z: sp.Expr) -> sp.Matrix:
    return sp.Matrix([sp.re(z), sp.im(z)])


def _line_intersection(p: sp.Expr, d: sp.Expr, q: sp.Expr, e: sp.Expr) -> sp.Expr:
    matrix = sp.Matrix([[sp.re(d), -sp.re(e)], [sp.im(d), -sp.im(e)]])
    rhs = sp.Matrix([sp.re(q - p), sp.im(q - p)])
    solution = matrix.LUsolve(rhs)
    t = solution[0]
    point_vec = _vec(p) + t * _vec(d)
    return sp.simplify(point_vec[0] + sp.I * point_vec[1])


def _tangent_intersection(zp: sp.Expr, zq: sp.Expr, zu: sp.Expr) -> sp.Expr:
    direction_p = sp.I * (zp - zu)
    direction_q = sp.I * (zq - zu)
    return _line_intersection(zp, direction_p, zq, direction_q)


def test_euler_center_midpoint() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 2 + 0j)
    _assign_value(engine, "B", -1 + 3 * sp.I)
    _assign_value(engine, "C", 1 + 4 * sp.I)

    circ_label = "__circ_expected"
    orth_label = "__orth_expected"
    engine.circumcenter("A", "B", "C", circ_label)
    engine.orthocenter_via_altitudes("A", "B", "C", orth_label)

    engine.euler_center("A", "B", "C", "N")

    expected_z = sp.simplify((engine.z(circ_label) + engine.z(orth_label)) / 2)
    expected_zb = sp.simplify((engine.zb(circ_label) + engine.zb(orth_label)) / 2)

    assert sp.simplify(engine.z("N") - expected_z) == 0
    assert sp.simplify(engine.zb("N") - expected_zb) == 0


def test_lemoine_point_via_tangents() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 2 + 0j)
    _assign_value(engine, "B", -1 + 3 * sp.I)
    _assign_value(engine, "C", 1 + 4 * sp.I)

    circ_label = "__circ_expected"
    engine.circumcenter("A", "B", "C", circ_label)

    engine.lemoine_point("A", "B", "C", "K")

    zA, zB, zC = engine.z("A"), engine.z("B"), engine.z("C")
    zU = engine.z(circ_label)

    tangent_bc = _tangent_intersection(zB, zC, zU)
    tangent_ca = _tangent_intersection(zC, zA, zU)

    expected = _line_intersection(zA, tangent_bc - zA, zB, tangent_ca - zB)
    assert sp.simplify(engine.z("K") - expected) == 0


def test_fermat_point_angles() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 2 + 0j)
    _assign_value(engine, "B", -1 + 3 * sp.I)
    _assign_value(engine, "C", 1 + 4 * sp.I)

    engine.add_fermat_points("A", "B", "C", "F1", "F2")

    sixty = sp.pi / 3
    one_twenty = 2 * sp.pi / 3

    poly_f1_sixty = engine.angle_value_poly("A", "F1", "B", sixty)
    poly_f1_one_twenty = engine.angle_value_poly("C", "F1", "B", one_twenty)
    poly_f2_one_twenty = engine.angle_value_poly("A", "F2", "B", one_twenty)

    assert sp.simplify(engine._apply_all(poly_f1_sixty)) == 0
    assert sp.simplify(engine._apply_all(poly_f1_one_twenty)) == 0
    assert sp.simplify(engine._apply_all(poly_f2_one_twenty)) == 0


def test_angle_bisector_either_poly() -> None:
    engine = GeometryEngine()
    _assign_value(engine, "A", 0)
    _assign_value(engine, "B", 1)
    _assign_value(engine, "C", sp.exp(sp.I * sp.pi / 3))

    _assign_value(engine, "D_int", sp.Rational(3, 2) * sp.exp(sp.I * sp.pi / 6))
    _assign_value(engine, "D_ext", 2 * sp.exp(sp.I * (sp.pi / 6 + sp.pi / 2)))
    _assign_value(engine, "D_off", sp.Rational(3, 2) * sp.exp(sp.I * sp.pi / 3))

    poly_int = engine.angle_bisector_either_poly("A", "B", "C", "D_int")
    poly_ext = engine.angle_bisector_either_poly("A", "B", "C", "D_ext")
    poly_off = engine.angle_bisector_either_poly("A", "B", "C", "D_off")

    assert sp.simplify(engine._apply_all(poly_int)) == 0
    assert sp.simplify(engine._apply_all(poly_ext)) == 0
    assert sp.simplify(engine._apply_all(poly_off)) != 0
