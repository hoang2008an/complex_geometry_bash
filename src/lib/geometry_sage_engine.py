"""
Isolated Sage backend for rational complex-geometry bashing.

This module intentionally does not import or mutate the SymPy engine.  Sage mode
v1 is limited to rational/linear exact workflows over a fraction field with
QQbar coefficients.  General branch-producing constructions are deliberately
outside this module's supported surface.  A rational secant case is supported
when one provided line point is already known to lie on the circle.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from sage.all import I, QQ, QQbar, FractionField, PolynomialRing


class SageGeometryError(RuntimeError):
    """Raised when a Sage-mode construction or check fails."""


@dataclass(frozen=True)
class SagePointRecord:
    z_name: str
    zb_name: str


@dataclass(frozen=True)
class SageUnitTriangleConfig:
    points: Dict[str, str]
    roots: Dict[str, str]


@dataclass(frozen=True)
class SageLine:
    alpha: Any
    beta: Any
    gamma: Any
    context: Tuple[str, ...] = ()
    point_on: Optional[str] = None

    def evaluate(self, z_expr: Any, zb_expr: Any) -> Any:
        return self.alpha * z_expr + self.beta * zb_expr + self.gamma


class SageGeometryEngine:
    """
    Rational-function complex geometry engine backed by Sage fraction fields.

    Parameters
    ----------
    point_names:
        Optional iterable of point labels to pre-register.  CLI scripts should
        pass all labels up front so the fraction field contains every coordinate
        variable needed by later operations.
    """

    def __init__(self, point_names: Optional[Iterable[str]] = None, *, coefficient_field: str = "QQbar") -> None:
        if coefficient_field not in {"QQbar", "QQ"}:
            raise SageGeometryError("coefficient_field must be 'QQbar' or 'QQ'.")
        self.coefficient_field = coefficient_field
        self._base_field = QQbar if coefficient_field == "QQbar" else QQ
        labels = ["O"]
        if point_names is not None:
            labels.extend(str(name) for name in point_names)
        self._declared_names: List[str] = []
        self.points: Dict[str, SagePointRecord] = {}
        self.constraints: List[Any] = []
        self.value_subs: Dict[str, Any] = {}
        self.point_assignments: Dict[str, Any] = {}
        self.learned_subs: Dict[str, Any] = {}
        self.unit_circle_points: set[str] = set()
        self.lines: Dict[str, SageLine] = {}
        self._main_unit_triangle: Optional[SageUnitTriangleConfig] = None
        self._field = None
        self._poly_ring = None
        self._gens: Dict[str, Any] = {}
        self._apply_cache: Dict[Any, Any] = {}

        for label in labels:
            if label not in self._declared_names:
                self._declared_names.append(label)
        self._rebuild_field()
        for label in self._declared_names:
            self._register_point_record(label)
        self._register_origin()

    # ------------------------------------------------------------------
    # Field and coercion helpers
    # ------------------------------------------------------------------
    def _variable_names(self) -> List[str]:
        names: List[str] = []
        for label in self._declared_names:
            names.extend([f"z_{label}", f"zb_{label}"])
        return names

    def _rebuild_field(self) -> None:
        names = self._variable_names() or ["z_O", "zb_O"]
        old_values = {
            "constraints": [str(expr) for expr in self.constraints],
            "value_subs": {name: str(expr) for name, expr in self.value_subs.items()},
            "point_assignments": {name: str(expr) for name, expr in self.point_assignments.items()},
            "learned_subs": {name: str(expr) for name, expr in self.learned_subs.items()},
        }

        self._poly_ring = PolynomialRing(self._base_field, names=names)
        self._field = FractionField(self._poly_ring)
        self._gens = {name: self._field(gen) for name, gen in zip(names, self._poly_ring.gens())}

        self.constraints = [self.expr(text) for text in old_values["constraints"]]
        self.value_subs = {name: self.expr(text) for name, text in old_values["value_subs"].items()}
        self.point_assignments = {name: self.expr(text) for name, text in old_values["point_assignments"].items()}
        self.learned_subs = {name: self.expr(text) for name, text in old_values["learned_subs"].items()}
        self._apply_cache.clear()

    def expr(self, value: Any) -> Any:
        """Coerce a Python/Sage value into the active fraction field."""
        if isinstance(value, str):
            text = value.replace("1j", "I").replace("j", "I")
            if not any(name in text for name in self._gens):
                try:
                    constant = eval(text, {"__builtins__": {}}, {"I": I, "QQbar": QQbar, "QQ": QQ})
                    return self._field(self._base_field(constant))
                except (NameError, SyntaxError, TypeError, ValueError):
                    pass
            return self._field(text)
        return self._field(value)

    def root_of_unity_3(self) -> Any:
        """Return a primitive third root of unity as a QQbar coefficient."""
        return self.expr(QQbar.zeta(3))

    def _substitution_items(self) -> List[Tuple[Any, Any]]:
        items: List[Tuple[Any, Any]] = []
        for mapping in (self.value_subs, self.point_assignments, self.learned_subs):
            for symbol_name, expr in mapping.items():
                items.append((self._gens[symbol_name], expr))
        return items

    def _apply_all(self, expr: Any) -> Any:
        current = self.expr(expr)
        cached = self._apply_cache.get(current)
        if cached is not None:
            return cached
        substitutions = dict(self._substitution_items())
        if not substitutions:
            return current
        for _ in range(8):
            updated = current.subs(substitutions)
            updated = self.expr(updated)
            if updated == current:
                self._apply_cache[expr] = updated
                self._apply_cache[updated] = updated
                return updated
            current = updated
        self._apply_cache[expr] = current
        self._apply_cache[current] = current
        return current

    def _numerator(self, expr: Any) -> Any:
        return self.expr(expr).numerator()

    def _is_zero(self, expr: Any) -> bool:
        return self.expr(expr) == 0

    # ------------------------------------------------------------------
    # Point management
    # ------------------------------------------------------------------
    def add_point(self, name: str) -> None:
        if name in self.points:
            return
        if name not in self._declared_names:
            self._declared_names.append(name)
            self._rebuild_field()
        self._register_point_record(name)

    def _register_point_record(self, name: str) -> None:
        self.points[name] = SagePointRecord(z_name=f"z_{name}", zb_name=f"zb_{name}")

    def ensure_point(self, name: str) -> None:
        if name not in self.points:
            raise SageGeometryError(f"Point '{name}' is not registered.")

    def z_symbol(self, name: str) -> Any:
        self.ensure_point(name)
        return self._gens[self.points[name].z_name]

    def zb_symbol(self, name: str) -> Any:
        self.ensure_point(name)
        return self._gens[self.points[name].zb_name]

    def z(self, name: str) -> Any:
        record = self._record(name)
        return self._apply_all(self.point_assignments.get(record.z_name, self.z_symbol(name)))

    def zb(self, name: str) -> Any:
        record = self._record(name)
        return self._apply_all(self.point_assignments.get(record.zb_name, self.zb_symbol(name)))

    def _record(self, name: str) -> SagePointRecord:
        self.ensure_point(name)
        return self.points[name]

    def _register_origin(self) -> None:
        self.add_point("O")
        record = self._record("O")
        self.value_subs[record.z_name] = self.expr(0)
        self.value_subs[record.zb_name] = self.expr(0)

    def set_point_value(self, name: str, z_value: Any, zb_value: Optional[Any] = None) -> None:
        """Fix a point to an exact algebraic value."""
        self.add_point(name)
        z_expr = self.expr(z_value)
        if zb_value is None:
            zb_expr = self.expr(QQbar(z_expr).conjugate()) if not z_expr.denominator().variables() else z_expr
        else:
            zb_expr = self.expr(zb_value)
        record = self._record(name)
        self.value_subs[record.z_name] = z_expr
        self.value_subs[record.zb_name] = zb_expr
        self._apply_cache.clear()

    def _set_point_assignment(self, name: str, z_expr: Any, zb_expr: Any, *, apply: bool = True) -> None:
        self.add_point(name)
        record = self._record(name)
        if apply:
            self.point_assignments[record.z_name] = self._apply_all(z_expr)
            self.point_assignments[record.zb_name] = self._apply_all(zb_expr)
        else:
            self.point_assignments[record.z_name] = self.expr(z_expr)
            self.point_assignments[record.zb_name] = self.expr(zb_expr)
        self._apply_cache.clear()

    def add_unit_circle(self, name: str) -> None:
        self.ensure_point(name)
        if name in self.unit_circle_points:
            return
        z = self.z_symbol(name)
        zb = self.zb_symbol(name)
        self.constraints.append(z * zb - 1)
        self.learned_subs[self._record(name).zb_name] = 1 / z
        self._apply_cache.clear()
        self.unit_circle_points.add(name)

    # ------------------------------------------------------------------
    # Constraint handling
    # ------------------------------------------------------------------
    def add_constraint(self, constraint: Any) -> None:
        expr = self.expr(constraint)
        self.constraints.append(expr)
        self._auto_learn_from_constraint(self._apply_all(expr))
        self._apply_cache.clear()

    def _auto_learn_from_constraint(self, constraint: Any) -> None:
        if self._is_zero(constraint):
            return
        numerator = self._numerator(constraint)
        present = []
        for label, record in self.points.items():
            symbol = self._poly_ring(record.zb_name)
            if numerator.degree(symbol) > 0:
                present.append(record.zb_name)
        if len(present) != 1:
            return

        target_name = present[0]
        target_poly_gen = self._poly_ring(target_name)
        if numerator.degree(target_poly_gen) != 1:
            return
        a = numerator.monomial_coefficient(target_poly_gen)
        b = numerator.monomial_coefficient(self._poly_ring(1))
        if a == 0:
            return
        rhs = self.expr(-b / a)
        if any(self._poly_ring(name) in rhs.numerator().variables() for name in self._zb_names()):
            return
        self.learned_subs[target_name] = rhs
        self._apply_cache.clear()

    def _zb_names(self) -> List[str]:
        return [record.zb_name for record in self.points.values()]

    # ------------------------------------------------------------------
    # Predicate polynomials
    # ------------------------------------------------------------------
    def collinear_poly(self, A: str, B: str, C: str) -> Any:
        zA, zB, zC = self.z_symbol(A), self.z_symbol(B), self.z_symbol(C)
        zbA, zbB, zbC = self.zb_symbol(A), self.zb_symbol(B), self.zb_symbol(C)
        return (zA - zC) * (zbB - zbC) - (zbA - zbC) * (zB - zC)

    def perpendicular_poly(self, A: str, B: str, C: str, D: str) -> Any:
        zA, zB = self.z_symbol(A), self.z_symbol(B)
        zC, zD = self.z_symbol(C), self.z_symbol(D)
        zbA, zbB = self.zb_symbol(A), self.zb_symbol(B)
        zbC, zbD = self.zb_symbol(C), self.zb_symbol(D)
        return (zA - zB) * (zbC - zbD) + (zbA - zbB) * (zC - zD)

    def angle_value_poly(self, A: str, B: str, C: str, angle: Any) -> Any:
        w = self._angle_unit(angle)
        wc = 1 / w
        zA, zB, zC = self.z_symbol(A), self.z_symbol(B), self.z_symbol(C)
        zbA, zbB, zbC = self.zb_symbol(A), self.zb_symbol(B), self.zb_symbol(C)
        return wc * (zA - zB) * (zbC - zbB) - w * (zbA - zbB) * (zC - zB)

    def angle_bisector_either_poly(self, A: str, B: str, C: str, D: str) -> Any:
        zA, zB, zC, zD = self.z_symbol(A), self.z_symbol(B), self.z_symbol(C), self.z_symbol(D)
        zbA, zbB, zbC, zbD = self.zb_symbol(A), self.zb_symbol(B), self.zb_symbol(C), self.zb_symbol(D)
        return (zD - zA) ** 2 * (zbB - zbA) * (zbC - zbA) - (zbD - zbA) ** 2 * (zB - zA) * (zC - zA)

    def midpoint_polys(self, P: str, Q: str, M: str) -> Tuple[Any, Any]:
        return 2 * self.z_symbol(M) - self.z_symbol(P) - self.z_symbol(Q), 2 * self.zb_symbol(M) - self.zb_symbol(P) - self.zb_symbol(Q)

    def centroid_polys(self, A: str, B: str, C: str, G: str) -> Tuple[Any, Any]:
        return self.z_symbol(A) + self.z_symbol(B) + self.z_symbol(C) - 3 * self.z_symbol(G), self.zb_symbol(A) + self.zb_symbol(B) + self.zb_symbol(C) - 3 * self.zb_symbol(G)

    def circumcenter_polys(self, A: str, B: str, C: str, U: str) -> Tuple[Any, Any]:
        zU, zbU = self.z_symbol(U), self.zb_symbol(U)

        def dist(P: str) -> Any:
            return (self.z_symbol(P) - zU) * (self.zb_symbol(P) - zbU)

        return dist(A) - dist(B), dist(A) - dist(C)

    def projection_to_line_polys(self, P: str, A: str, B: str, H: str) -> Tuple[Any, Any]:
        if A == B:
            raise SageGeometryError("Projection requires two distinct points to define line AB.")
        return self.perpendicular_poly(A, B, P, H), self.collinear_poly(H, A, B)

    def line_reflection_polys(self, P: str, A: str, B: str, Q: str) -> Tuple[Any, Any]:
        midpoint_z = (self.z_symbol(P) + self.z_symbol(Q)) / 2
        midpoint_zb = (self.zb_symbol(P) + self.zb_symbol(Q)) / 2
        zA, zB = self.z_symbol(A), self.z_symbol(B)
        zbA, zbB = self.zb_symbol(A), self.zb_symbol(B)
        eq_mid_on_line = (midpoint_z - zB) * (zbA - zbB) - (midpoint_zb - zbB) * (zA - zB)
        return self.perpendicular_poly(A, B, P, Q), eq_mid_on_line

    def isogonal_reflection_poly(self, A: str, B: str, C: str, D: str, E: str) -> Any:
        for label in (A, B, C, D, E):
            self.ensure_point(label)
        zA, zB, zC, zD, zE = (
            self.z_symbol(A),
            self.z_symbol(B),
            self.z_symbol(C),
            self.z_symbol(D),
            self.z_symbol(E),
        )
        zbA, zbB, zbC, zbD, zbE = (
            self.zb_symbol(A),
            self.zb_symbol(B),
            self.zb_symbol(C),
            self.zb_symbol(D),
            self.zb_symbol(E),
        )
        left = (zD - zA) * (zE - zA) * (zbB - zbA) * (zbC - zbA)
        right = (zB - zA) * (zC - zA) * (zbD - zbA) * (zbE - zbA)
        return left - right

    def isogonal_conjugate_polys(self, A: str, B: str, C: str, P: str, Q: str) -> Tuple[Any, Any, Any]:
        return (
            self.isogonal_reflection_poly(A, B, C, P, Q),
            self.isogonal_reflection_poly(B, C, A, P, Q),
            self.isogonal_reflection_poly(C, A, B, P, Q),
        )

    # ------------------------------------------------------------------
    # Constraint guards
    # ------------------------------------------------------------------
    def add_collinear(self, A: str, B: str, C: str) -> None:
        self.add_constraint(self.collinear_poly(A, B, C))

    def add_perpendicular(self, A: str, B: str, C: str, D: str) -> None:
        self.add_constraint(self.perpendicular_poly(A, B, C, D))

    def add_angle_value(self, A: str, B: str, C: str, angle: Any) -> None:
        self.add_constraint(self.angle_value_poly(A, B, C, angle))

    def add_angle_bisector_either(self, A: str, B: str, C: str, D: str) -> None:
        self.add_constraint(self.angle_bisector_either_poly(A, B, C, D))

    def add_midpoint(self, P: str, Q: str, M: str) -> None:
        self.add_point(M)
        for eq in self.midpoint_polys(P, Q, M):
            self.add_constraint(eq)

    def add_centroid_constraint(self, A: str, B: str, C: str, G: str) -> None:
        self.add_point(G)
        for eq in self.centroid_polys(A, B, C, G):
            self.add_constraint(eq)

    def add_circumcenter(self, A: str, B: str, C: str, U: str) -> None:
        self.add_point(U)
        for eq in self.circumcenter_polys(A, B, C, U):
            self.add_constraint(eq)

    def add_projection_to_line(self, P: str, A: str, B: str, H: str) -> None:
        self.add_point(H)
        for eq in self.projection_to_line_polys(P, A, B, H):
            self.add_constraint(eq)

    def add_line_reflection(self, P: str, A: str, B: str, Q: str) -> None:
        self.add_point(Q)
        for eq in self.line_reflection_polys(P, A, B, Q):
            self.add_constraint(eq)

    def add_isogonal_conjugate(self, A: str, B: str, C: str, P: str, Q: str) -> None:
        for eq in self.isogonal_conjugate_polys(A, B, C, P, Q):
            self.add_constraint(eq)

    # ------------------------------------------------------------------
    # Line objects
    # ------------------------------------------------------------------
    def line_from_coefficients(
        self,
        alpha: Any,
        beta: Any,
        gamma: Any,
        *,
        context: Optional[Sequence[str]] = None,
        point_on: Optional[str] = None,
    ) -> SageLine:
        return SageLine(
            alpha=self.expr(alpha),
            beta=self.expr(beta),
            gamma=self.expr(gamma),
            context=tuple(context or ()),
            point_on=point_on,
        )

    def symbolic_line(self, name: str) -> SageLine:
        return SageLine(alpha=None, beta=None, gamma=None, context=(name,))

    def register_line(self, name: str, line: SageLine) -> SageLine:
        if name in self.lines:
            raise SageGeometryError(f"Line '{name}' is already defined.")
        self.lines[name] = line
        return line

    def get_line(self, name: str) -> SageLine:
        try:
            return self.lines[name]
        except KeyError as exc:
            raise SageGeometryError(f"Line '{name}' is not defined.") from exc

    def _coerce_line(self, line: SageLine | str) -> SageLine:
        if isinstance(line, SageLine):
            return line
        if isinstance(line, str):
            return self.get_line(line)
        raise SageGeometryError("Expected a line label or line object.")

    def line_through_points(self, P: str, Q: str) -> SageLine:
        self.ensure_point(P)
        self.ensure_point(Q)
        if P == Q:
            raise SageGeometryError("Line requires two distinct points.")
        zP, zQ = self.z_symbol(P), self.z_symbol(Q)
        zbP, zbQ = self.zb_symbol(P), self.zb_symbol(Q)
        alpha = zbQ - zbP
        beta = -(zQ - zP)
        gamma = (zQ - zP) * zbP - (zbQ - zbP) * zP
        return self.line_from_coefficients(alpha, beta, gamma, context=(P, Q))

    def add_point_on_line(self, line: SageLine | str, point: str) -> None:
        self.ensure_point(point)
        line_obj = self._coerce_line(line)
        if line_obj.alpha is None:
            label = line_obj.context[0] if line_obj.context else None
            updated = SageLine(None, None, None, line_obj.context, point_on=point)
            if label and label in self.lines:
                self.lines[label] = updated
            return
        self.add_constraint(line_obj.evaluate(self.z_symbol(point), self.zb_symbol(point)))

    def add_perpendicular_lines(self, line1: SageLine | str, line2: SageLine | str) -> None:
        line1_obj = self._coerce_line(line1)
        line2_obj = self._coerce_line(line2)

        materialized = self._materialize_perpendicular_line(line1_obj, line2_obj)
        if materialized is not None:
            return
        materialized = self._materialize_perpendicular_line(line2_obj, line1_obj)
        if materialized is not None:
            return

        if line1_obj.alpha is None or line2_obj.alpha is None:
            raise SageGeometryError("Perpendicular symbolic lines need a known point and a concrete reference line.")
        a1, b1 = self._line_normal_components(line1_obj)
        a2, b2 = self._line_normal_components(line2_obj)
        constraint = a1 * a2 + b1 * b2
        if constraint != 0:
            self.add_constraint(constraint)

    def _materialize_perpendicular_line(self, free_line: SageLine, reference_line: SageLine) -> Optional[SageLine]:
        if free_line.alpha is not None:
            return None
        if reference_line.alpha is None or free_line.point_on is None:
            return None
        label = free_line.context[0] if free_line.context else None
        if label is None:
            return None
        point = free_line.point_on
        alpha = -reference_line.alpha
        beta = reference_line.beta
        gamma = -alpha * self.z_symbol(point) - beta * self.zb_symbol(point)
        concrete = self.line_from_coefficients(alpha, beta, gamma, context=free_line.context, point_on=point)
        self.lines[label] = concrete
        return concrete

    def _line_normal_components(self, line: SageLine) -> Tuple[Any, Any]:
        if line.alpha is None:
            raise SageGeometryError("Symbolic line has not been determined yet.")
        return line.alpha + line.beta, I * (line.beta - line.alpha)

    def line_intersection(self, line1: SageLine | str, line2: SageLine | str, name: str) -> Tuple[Any, Any]:
        self.add_point(name)
        line1_obj = self._coerce_line(line1)
        line2_obj = self._coerce_line(line2)
        if line1_obj.alpha is None or line2_obj.alpha is None:
            raise SageGeometryError("Line intersection requires concrete lines.")
        z_var, zb_var = self.z_symbol(name), self.zb_symbol(name)
        eq1 = self._apply_all(line1_obj.evaluate(z_var, zb_var))
        eq2 = self._apply_all(line2_obj.evaluate(z_var, zb_var))
        solution = self._solve_linear_point_system(eq1, eq2, z_var, zb_var, "Line intersection")
        self._set_point_assignment(name, *solution)
        return self.z(name), self.zb(name)

    # ------------------------------------------------------------------
    # Constructions
    # ------------------------------------------------------------------
    def midpoint(self, P: str, Q: str, M: str) -> Any:
        self._set_point_assignment(M, (self.z(P) + self.z(Q)) / 2, (self.zb(P) + self.zb(Q)) / 2)
        return self.z(M)

    def centroid(self, A: str, B: str, C: str, G: str) -> Any:
        self._set_point_assignment(G, (self.z(A) + self.z(B) + self.z(C)) / 3, (self.zb(A) + self.zb(B) + self.zb(C)) / 3)
        return self.z(G)

    def circumcenter(self, A: str, B: str, C: str, U: str) -> Any:
        self.add_point(U)
        zU, zbU = self.z_symbol(U), self.zb_symbol(U)
        eq1, eq2 = self.circumcenter_polys(A, B, C, U)
        self._set_point_assignment(U, *self._solve_linear_point_system(self._apply_all(eq1), self._apply_all(eq2), zU, zbU, "Circumcenter"))
        return self.z(U)

    def orthocenter_via_altitudes(self, A: str, B: str, C: str, H: str) -> Tuple[Any, Any]:
        self.add_point(H)
        zH, zbH = self.z_symbol(H), self.zb_symbol(H)
        eq1 = self._apply_all(self.perpendicular_poly(A, H, B, C))
        eq2 = self._apply_all(self.perpendicular_poly(B, H, A, C))
        solution = self._solve_linear_point_system(eq1, eq2, zH, zbH, "Orthocenter")
        self._set_point_assignment(H, *solution)
        return self.z(H), self.zb(H)

    def intersection_of_lines(self, A: str, B: str, C: str, D: str, X: str) -> Tuple[Any, Any]:
        self.add_point(X)
        zX, zbX = self.z_symbol(X), self.zb_symbol(X)
        eq1 = self._apply_all(self.collinear_poly(X, A, B))
        eq2 = self._apply_all(self.collinear_poly(X, C, D))
        solution = self._solve_linear_point_system(eq1, eq2, zX, zbX, "Line intersection")
        self._set_point_assignment(X, *solution)
        return self.z(X), self.zb(X)

    def project_point_to_line(self, P: str, A: str, B: str, H: str) -> Tuple[Any, Any]:
        self.add_point(H)
        zH, zbH = self.z_symbol(H), self.zb_symbol(H)
        eq1, eq2 = self.projection_to_line_polys(P, A, B, H)
        solution = self._solve_linear_point_system(self._apply_all(eq1), self._apply_all(eq2), zH, zbH, "Projection")
        self._set_point_assignment(H, *solution)
        return self.z(H), self.zb(H)

    def reflect_point_over_line(self, P: str, A: str, B: str, Q: str) -> Any:
        self.add_point(Q)
        zQ, zbQ = self.z_symbol(Q), self.zb_symbol(Q)
        eq1, eq2 = self.line_reflection_polys(P, A, B, Q)
        solution = self._solve_linear_point_system(self._apply_all(eq1), self._apply_all(eq2), zQ, zbQ, "Line reflection")
        self._set_point_assignment(Q, *solution)
        return self.z(Q)

    def line_circle_intersection(
        self,
        line_point1: str,
        line_point2: str,
        center: str,
        radius_point: str,
        name: str,
        *,
        avoid: Optional[Sequence[str]] = None,
    ) -> Tuple[Any, Any]:
        """
        Construct a line-circle intersection in the rational secant case.

        This intentionally does not solve a general quadratic.  It requires one
        provided line endpoint to already lie on the circle, then gets the other
        intersection by Vieta.
        """
        self.add_point(name)
        for label in (line_point1, line_point2, center, radius_point):
            self.ensure_point(label)

        radius_sq = self.squared_distance(center, radius_point)
        candidates: List[Tuple[Any, Any]] = []

        ordered_pairs: List[Tuple[str, str, bool]] = []
        if line_point1 == radius_point:
            ordered_pairs.append((line_point1, line_point2, True))
        if line_point2 == radius_point:
            ordered_pairs.append((line_point2, line_point1, True))
        ordered_pairs.extend([(line_point1, line_point2, False), (line_point2, line_point1, False)])

        for base, direction_point, known_on_circle in ordered_pairs:
            if known_on_circle or self._is_zero(self.squared_distance(base, center) - radius_sq):
                base_z, base_zb = self.z(base), self.zb(base)
                dir_z = self.z(direction_point) - base_z
                dir_zb = self.zb(direction_point) - base_zb
                denominator = self.expr(dir_z * dir_zb)
                if denominator == 0:
                    raise SageGeometryError("Line-circle intersection requires distinct line points.")
                numerator = (base_z - self.z(center)) * dir_zb + (base_zb - self.zb(center)) * dir_z
                t_other = self.expr(-numerator / denominator)
                other = (
                    self.expr(base_z + t_other * dir_z),
                    self.expr(base_zb + t_other * dir_zb),
                )
                candidates.extend([(base_z, base_zb), other])
                break

        if not candidates:
            raise SageGeometryError(
                "Sage mode supports line_circle_intersection only when one line endpoint is already on the circle."
            )

        if avoid and len(candidates) == 2:
            avoid_set = set(avoid)
            for endpoint in (line_point1, line_point2, radius_point):
                if endpoint in avoid_set:
                    z_endpoint, zb_endpoint = self.z(endpoint), self.zb(endpoint)
                    if candidates[0] == (z_endpoint, zb_endpoint):
                        self._set_point_assignment(name, candidates[1][0], candidates[1][1], apply=False)
                        return self.z(name), self.zb(name)

        avoid_labels = set(avoid or [])
        avoid_pairs = self._prepare_avoid_pairs(avoid)
        for z_candidate, zb_candidate in candidates:
            if z_candidate == self.z(radius_point) and zb_candidate == self.zb(radius_point) and radius_point in avoid_labels:
                continue
            if not self._matches_any(z_candidate, zb_candidate, avoid_pairs):
                self._set_point_assignment(name, z_candidate, zb_candidate, apply=False)
                return self.z(name), self.zb(name)

        raise SageGeometryError("Line-circle intersection has no candidate left after applying avoid.")

    def circle_intersection(
        self,
        center1: str,
        radius_point1: str,
        center2: str,
        radius_point2: str,
        name: str,
        *,
        avoid: Optional[Sequence[str]] = None,
    ) -> Tuple[Any, Any]:
        """
        Construct the second intersection of two circles when a known common
        point is provided through ``avoid``.
        """
        self.add_point(name)
        for label in (center1, radius_point1, center2, radius_point2):
            self.ensure_point(label)

        radius_sq1 = self.squared_distance(center1, radius_point1)
        radius_sq2 = self.squared_distance(center2, radius_point2)
        for known in avoid or ():
            self.ensure_point(known)
            on_first = self._is_zero(self.squared_distance(known, center1) - radius_sq1)
            on_second = self._is_zero(self.squared_distance(known, center2) - radius_sq2)
            if on_first and on_second:
                return self._circle_intersection_from_known_point(
                    center1,
                    radius_sq1,
                    center2,
                    radius_sq2,
                    known,
                    name,
                    avoid=avoid,
                )
        raise SageGeometryError(
            "Sage mode supports circle_intersection only when avoid names a known common intersection."
        )

    def _circle_intersection_from_known_point(
        self,
        center1: str,
        radius_sq1: Any,
        center2: str,
        radius_sq2: Any,
        known: str,
        name: str,
        *,
        avoid: Optional[Sequence[str]] = None,
    ) -> Tuple[Any, Any]:
        p_z, p_zb = self.z(known), self.zb(known)
        z1, zb1 = self.z(center1), self.zb(center1)
        z2, zb2 = self.z(center2), self.zb(center2)

        alpha = zb2 - zb1
        beta = z2 - z1
        gamma = z1 * zb1 - z2 * zb2 - radius_sq1 + radius_sq2
        dir_z = beta
        dir_zb = -alpha
        denominator = self.expr(dir_z * dir_zb)
        if denominator == 0:
            raise SageGeometryError("Circle intersection radical axis is degenerate.")
        numerator = (p_z - z1) * dir_zb + (p_zb - zb1) * dir_z
        t_other = self.expr(-numerator / denominator)
        candidates = [(p_z, p_zb), (self.expr(p_z + t_other * dir_z), self.expr(p_zb + t_other * dir_zb))]

        avoid_pairs = self._prepare_avoid_pairs(avoid)
        for z_candidate, zb_candidate in candidates:
            if not self._matches_any(z_candidate, zb_candidate, avoid_pairs):
                self._set_point_assignment(name, z_candidate, zb_candidate, apply=False)
                return self.z(name), self.zb(name)
        raise SageGeometryError("Circle intersection has no candidate left after applying avoid.")

    def set_main_unit_triangle(
        self,
        A: str,
        B: str,
        C: str,
        *,
        root_names: Optional[Sequence[str]] = None,
    ) -> None:
        if len({A, B, C}) != 3:
            raise SageGeometryError("Main unit triangle requires three distinct vertex labels.")
        for vertex in (A, B, C):
            self.ensure_point(vertex)
            if vertex not in self.unit_circle_points:
                raise SageGeometryError(f"Point '{vertex}' must be declared on the unit circle first.")

        roots_sequence = tuple(root_names) if root_names is not None else (
            self._unit_root_name("A"),
            self._unit_root_name("B"),
            self._unit_root_name("C"),
        )
        if len(roots_sequence) != 3 or len(set(roots_sequence)) != 3:
            raise SageGeometryError("root_names must provide three distinct labels.")

        roots: Dict[str, str] = {}
        for canonical, vertex, root in zip(("A", "B", "C"), (A, B, C), roots_sequence):
            self.add_point(root)
            self.add_unit_circle(root)
            self._set_point_assignment(vertex, self.z_symbol(root) ** 2, self.zb_symbol(root) ** 2, apply=False)
            roots[canonical] = root
        self._main_unit_triangle = SageUnitTriangleConfig(
            points={"A": A, "B": B, "C": C},
            roots=roots,
        )

    def main_triangle_incenter(self, name: str) -> Any:
        config = self._require_main_unit_triangle()
        x = self.z_symbol(config.roots["A"])
        y = self.z_symbol(config.roots["B"])
        z = self.z_symbol(config.roots["C"])
        xb = self.zb_symbol(config.roots["A"])
        yb = self.zb_symbol(config.roots["B"])
        zb = self.zb_symbol(config.roots["C"])
        self._set_point_assignment(name, -(x * y + x * z + y * z), -(xb * yb + xb * zb + yb * zb), apply=False)
        return self.z(name)

    def main_triangle_arc_midpoint(self, which: str, name: str, *, containing_vertex: bool = False) -> Any:
        config = self._require_main_unit_triangle()
        normalized = which.upper()
        if normalized not in config.roots:
            raise SageGeometryError(f"Arc midpoint '{which}' is not valid for the configured main triangle.")
        if normalized == "A":
            root1, root2 = config.roots["B"], config.roots["C"]
        elif normalized == "B":
            root1, root2 = config.roots["C"], config.roots["A"]
        else:
            root1, root2 = config.roots["A"], config.roots["B"]
        sign = self.expr(1 if containing_vertex else -1)
        self._set_point_assignment(
            name,
            sign * self.z_symbol(root1) * self.z_symbol(root2),
            sign * self.zb_symbol(root1) * self.zb_symbol(root2),
            apply=False,
        )
        self.add_unit_circle(name)
        return self.z(name)

    def isogonal_conjugate_point(self, A: str, B: str, C: str, P: str, Q: str) -> Any:
        if len({A, B, C}) != 3:
            raise SageGeometryError("Isogonal conjugate requires three distinct triangle vertices.")
        self.add_point(Q)
        if P in {A, B, C}:
            self._set_point_assignment(Q, self.z(P), self.zb(P), apply=False)
            return self.z(Q)
        zQ, zbQ = self.z_symbol(Q), self.zb_symbol(Q)
        eq1, eq2, eq3 = self.isogonal_conjugate_polys(A, B, C, P, Q)
        solution = self._solve_linear_point_system(self._apply_all(eq1), self._apply_all(eq2), zQ, zbQ, "Isogonal conjugate")
        check = self._apply_all(eq3).subs({zQ: solution[0], zbQ: solution[1]})
        if not self._is_zero(check):
            raise SageGeometryError("Isogonal conjugate construction failed consistency check.")
        self._set_point_assignment(Q, *solution, apply=False)
        return self.z(Q)

    def lemoine_point(self, A: str, B: str, C: str, L: str, circumcenter_name: Optional[str] = None) -> Any:
        circ = circumcenter_name or self._unique_internal_name("LemoineCirc")
        self.circumcenter(A, B, C, circ)
        t_bc = self._unique_internal_name(f"T_{B}{C}")
        t_ca = self._unique_internal_name(f"T_{C}{A}")
        self._tangent_intersection(B, C, circ, t_bc)
        self._tangent_intersection(C, A, circ, t_ca)
        self.intersection_of_lines(A, t_bc, B, t_ca, L)
        return self.z(L)

    def add_fermat_points(self, A: str, B: str, C: str, F1: str, F2: str) -> Tuple[Any, Any]:
        if len({A, B, C}) != 3:
            raise SageGeometryError("Fermat point construction expects three distinct vertices.")
        w = self.root_of_unity_3()
        self._fermat_branch_point(A, B, C, F1, w)
        self._fermat_branch_point(A, B, C, F2, w ** 2)
        return self.z(F1), self.z(F2)

    def _fermat_branch_point(self, A: str, B: str, C: str, label: str, w: Any) -> None:
        zA, zB, zC = self.z(A), self.z(B), self.z(C)
        zbA, zbB, zbC = self.zb(A), self.zb(B), self.zb(C)
        w = self.expr(w)
        w2 = w ** 2
        z_num = zA * (zbB - zbC) + w * zB * (zbC - zbA) + w2 * zC * (zbA - zbB)
        z_den = (zbB - zbC) + w * (zbC - zbA) + w2 * (zbA - zbB)
        wc = 1 / w
        wc2 = 1 / w2
        zb_num = zbA * (zB - zC) + wc * zbB * (zC - zA) + wc2 * zbC * (zA - zB)
        zb_den = (zB - zC) + wc * (zC - zA) + wc2 * (zA - zB)
        if z_den == 0 or zb_den == 0:
            raise SageGeometryError("Fermat point construction is undefined: denominator vanished.")
        self._set_point_assignment(label, z_num / z_den, zb_num / zb_den)

    def _tangent_intersection(self, P: str, Q: str, center: str, name: str) -> Tuple[Any, Any]:
        self.add_point(name)
        zX, zbX = self.z_symbol(name), self.zb_symbol(name)
        eq1 = self._apply_all(self.perpendicular_poly(P, name, P, center))
        eq2 = self._apply_all(self.perpendicular_poly(Q, name, Q, center))
        solution = self._solve_linear_point_system(eq1, eq2, zX, zbX, "Tangent intersection")
        self._set_point_assignment(name, *solution)
        return self.z(name), self.zb(name)

    def _solve_linear_point_system(self, eq1: Any, eq2: Any, z_var: Any, zb_var: Any, label: str) -> Tuple[Any, Any]:
        p1 = self._numerator(eq1)
        p2 = self._numerator(eq2)
        z_poly = self._poly_ring(str(z_var))
        zb_poly = self._poly_ring(str(zb_var))
        if not self._is_linear_in_pair(p1, z_poly, zb_poly) or not self._is_linear_in_pair(p2, z_poly, zb_poly):
            raise SageGeometryError(f"{label} construction only supports linear rational systems in Sage mode v1.")
        a1, b1, c1 = self._linear_coefficients(p1, z_poly, zb_poly)
        a2, b2, c2 = self._linear_coefficients(p2, z_poly, zb_poly)
        det = self.expr(a1 * b2 - a2 * b1)
        if det == 0:
            raise SageGeometryError(f"{label} construction failed: equations do not determine a unique point.")
        return self.expr((b1 * c2 - b2 * c1) / det), self.expr((c1 * a2 - c2 * a1) / det)

    def _is_linear_in_pair(self, poly: Any, z_poly: Any, zb_poly: Any) -> bool:
        return (
            poly.derivative(z_poly).derivative(z_poly) == 0
            and poly.derivative(zb_poly).derivative(zb_poly) == 0
            and poly.derivative(z_poly).derivative(zb_poly) == 0
        )

    def _linear_coefficients(self, poly: Any, z_poly: Any, zb_poly: Any) -> Tuple[Any, Any, Any]:
        zero_subs = {z_poly: 0, zb_poly: 0}
        a = poly.derivative(z_poly).subs(zero_subs)
        b = poly.derivative(zb_poly).subs(zero_subs)
        c = poly.subs(zero_subs)
        return a, b, c

    def _prepare_avoid_pairs(self, avoid: Optional[Sequence[str]]) -> List[Tuple[Any, Any]]:
        pairs: List[Tuple[Any, Any]] = []
        if avoid is None:
            return pairs
        for label in avoid:
            self.ensure_point(label)
            pairs.append((self.z(label), self.zb(label)))
        return pairs

    def _matches_any(self, z_candidate: Any, zb_candidate: Any, pairs: Sequence[Tuple[Any, Any]]) -> bool:
        for z_avoid, zb_avoid in pairs:
            if self._is_zero(z_candidate - z_avoid) and self._is_zero(zb_candidate - zb_avoid):
                return True
        return False

    def _require_main_unit_triangle(self) -> SageUnitTriangleConfig:
        if self._main_unit_triangle is None:
            raise SageGeometryError("Main unit triangle is not configured. Use set_main_unit_triangle first.")
        return self._main_unit_triangle

    def _unit_root_name(self, vertex: str) -> str:
        return f"__unit_root_{vertex}"

    def _unique_internal_name(self, base: str) -> str:
        candidate = f"__{base.replace(' ', '_')}"
        idx = 1
        while candidate in self.points:
            candidate = f"__{base.replace(' ', '_')}_{idx}"
            idx += 1
        self.add_point(candidate)
        return candidate

    # ------------------------------------------------------------------
    # Checks and display
    # ------------------------------------------------------------------
    def constraint_polynomials(self, constraint: str, args: Sequence[str], *, angle: Optional[Any] = None) -> List[Any]:
        normalized = constraint.lower()
        if normalized in {"collinear", "line"}:
            return [self.collinear_poly(*args)]
        if normalized in {"perpendicular", "perp"}:
            return [self.perpendicular_poly(*args)]
        if normalized in {"angle", "angle_value"}:
            if angle is None:
                raise SageGeometryError("Angle constraint requires an angle.")
            return [self.angle_value_poly(*args, angle)]
        if normalized in {"angle_bisector_either", "angle_bisector_any"}:
            return [self.angle_bisector_either_poly(*args)]
        if normalized == "circumcenter":
            return list(self.circumcenter_polys(*args))
        if normalized == "midpoint":
            return list(self.midpoint_polys(*args))
        if normalized == "centroid":
            return list(self.centroid_polys(*args))
        if normalized == "projection":
            return list(self.projection_to_line_polys(*args))
        if normalized in {"isogonal_conjugate", "isogonal_conj", "add_isogonal_conjugate"}:
            return list(self.isogonal_conjugate_polys(*args))
        raise SageGeometryError(f"Unsupported Sage constraint '{constraint}'.")

    def constraint_conjugate_free(self, constraint: str, args: Sequence[str], *, angle: Optional[Any] = None) -> List[Tuple[Any, Any]]:
        results: List[Tuple[Any, Any]] = []
        for poly in self.constraint_polynomials(constraint, args, angle=angle):
            substituted = self._apply_all(poly)
            results.append((substituted.numerator(), substituted.denominator()))
        return results

    def squared_distance(self, P: str, Q: str) -> Any:
        return self._apply_all((self.z(P) - self.z(Q)) * (self.zb(P) - self.zb(Q)))

    def learned_rules(self) -> Dict[str, Any]:
        return dict(self.learned_subs)

    def point_summary(self, names: Optional[Sequence[str]] = None) -> Dict[str, str]:
        selected = sorted(self.points.keys()) if names is None else names
        return {name: str(self.z(name)) for name in selected}

    def constraint_strings(self) -> List[str]:
        return [str(expr) for expr in self.constraints]

    def _angle_unit(self, angle: Any) -> Any:
        if angle in ("pi/3", "60", 60):
            return self.expr(-QQbar.zeta(3) ** 2)
        if angle in ("2*pi/3", "120", 120):
            return self.expr(QQbar.zeta(3))
        if angle in ("-pi/3", "-60", -60):
            return self.expr(-QQbar.zeta(3))
        raise SageGeometryError("Sage mode v1 supports fixed algebraic angles pi/3 and 2*pi/3 only.")


__all__ = ["SageGeometryEngine", "SageGeometryError"]
