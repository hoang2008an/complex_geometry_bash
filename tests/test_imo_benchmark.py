import cmath
import math

import sympy as sp

from geometry_engine import GeometryEngine


def _line_equation_value(
    p: sp.Expr,
    pb: sp.Expr,
    q: sp.Expr,
    qb: sp.Expr,
    r: sp.Expr,
    rb: sp.Expr,
) -> sp.Expr:
    alpha = qb - pb
    beta = -(q - p)
    gamma = (q - p) * pb - (qb - pb) * p
    return sp.expand(alpha * r + beta * rb + gamma)


def _imo_configuration():
    engine = GeometryEngine()
    for name in ("A", "B", "C"):
        engine.add_point(name)
        engine.add_unit_circle(name)

    engine.set_main_unit_triangle("A", "B", "C", root_names=("x", "y", "z"))
    engine.main_triangle_incenter("I")
    engine.main_triangle_arc_midpoint("A", "P")
    engine.midpoint("A", "C", "K")
    engine.midpoint("A", "B", "L")
    engine.reflect_point_over_point("A", "I", "A2")
    engine.reflect_point_over_point("B", "I", "B2")
    engine.reflect_point_over_point("C", "I", "C2")

    x = engine.z_symbol("x")
    y = engine.z_symbol("y")
    z = engine.z_symbol("z")

    x_point = sp.cancel((2 * x * y**2 + 2 * x * y * z + x * z**2 + 2 * y**2 * z + y * z**2) / (x - y))
    xb_point = sp.cancel(-(x * y + 2 * x * z + y**2 + 2 * y * z + 2 * z**2) / (y * z**2 * (x - y)))
    y_point = sp.cancel((x * y**2 + 2 * x * y * z + 2 * x * z**2 + y**2 * z + 2 * y * z**2) / (x - z))
    yb_point = sp.cancel(-(2 * x * y + x * z + 2 * y**2 + 2 * y * z + z**2) / (y**2 * z * (x - z)))

    return engine, x_point, xb_point, y_point, yb_point


def _unsigned_angle(u: complex, v: complex, w: complex) -> float:
    angle = math.degrees(cmath.phase((u - v) / (w - v)))
    if angle < 0:
        angle = -angle
    if angle > 180:
        angle = 360 - angle
    return angle


def test_imo_incenter_tangent_benchmark_symbolic() -> None:
    engine, x_point, xb_point, y_point, yb_point = _imo_configuration()

    z_b = engine.z("B")
    zb_b = engine.zb("B")
    z_c = engine.z("C")
    zb_c = engine.zb("C")
    z_a2 = engine.z("A2")
    zb_a2 = engine.zb("A2")
    z_b2 = engine.z("B2")
    zb_b2 = engine.zb("B2")
    z_c2 = engine.z("C2")
    zb_c2 = engine.zb("C2")

    assert sp.simplify(_line_equation_value(z_b, zb_b, z_c, zb_c, x_point, xb_point)) == 0
    assert sp.simplify(_line_equation_value(z_a2, zb_a2, z_c2, zb_c2, x_point, xb_point)) == 0
    assert sp.simplify(_line_equation_value(z_b, zb_b, z_c, zb_c, y_point, yb_point)) == 0
    assert sp.simplify(_line_equation_value(z_a2, zb_a2, z_b2, zb_b2, y_point, yb_point)) == 0

    z_k = engine.z("K")
    zb_k = engine.zb("K")
    z_i = engine.z("I")
    zb_i = engine.zb("I")
    z_l = engine.z("L")
    zb_l = engine.zb("L")
    z_p = engine.z("P")
    zb_p = engine.zb("P")

    angle_kil = sp.cancel((z_k - z_i) * (zb_l - zb_i) / ((zb_k - zb_i) * (z_l - z_i)))
    angle_ypx = sp.cancel((y_point - z_p) * (xb_point - zb_p) / ((yb_point - zb_p) * (x_point - z_p)))

    assert sp.simplify(sp.together(angle_kil * angle_ypx - 1)) == 0


def test_imo_incenter_tangent_benchmark_numeric_branch() -> None:
    engine, x_point, _, y_point, _ = _imo_configuration()

    substitutions = {
        engine.z_symbol("x"): sp.N(sp.exp(sp.I * sp.Rational(9, 20))),
        engine.z_symbol("y"): sp.Integer(1),
        engine.z_symbol("z"): sp.N(sp.exp(sp.I * sp.Rational(11, 10))),
    }

    z_a = complex(sp.N(engine.z("A").subs(substitutions)))
    z_b = complex(sp.N(engine.z("B").subs(substitutions)))
    z_c = complex(sp.N(engine.z("C").subs(substitutions)))
    z_i = complex(sp.N(engine.z("I").subs(substitutions)))
    z_k = complex(sp.N(engine.z("K").subs(substitutions)))
    z_l = complex(sp.N(engine.z("L").subs(substitutions)))
    z_p = complex(sp.N(engine.z("P").subs(substitutions)))
    z_x = complex(sp.N(x_point.subs(substitutions)))
    z_y = complex(sp.N(y_point.subs(substitutions)))

    assert abs(z_a - z_b) < abs(z_a - z_c) < abs(z_b - z_c)

    total_angle = _unsigned_angle(z_k, z_i, z_l) + _unsigned_angle(z_y, z_p, z_x)
    assert math.isclose(total_angle, 180.0, rel_tol=0.0, abs_tol=1e-8)
