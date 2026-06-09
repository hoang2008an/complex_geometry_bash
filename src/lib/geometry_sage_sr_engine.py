"""
Sage symbolic-ring geometry engine.

This backend mirrors the SymPy engine's substitution-oriented style while using
Sage's symbolic ring (SR).  Point creation is cheap: each point gets two SR
variables and no global polynomial/fraction field is rebuilt.

The initial surface intentionally covers the construction/check workflow needed
by the current unit-circle scripts:

* point registration and fixed values
* unit-circle conjugate rules
* collinear/perpendicular/angle predicates
* midpoint, centroid, circumcenter, orthocenter, projection
* reflection over a line
* rational secant line-circle intersections with avoid
* constraint checks and summaries
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

from sage.all import I, SR, factor


class SageSRGeometryError(RuntimeError):
    """Raised when the SR geometry backend cannot complete an operation."""


@dataclass(frozen=True)
class SageSRPointRecord:
    z: Any
    zb: Any


@dataclass(frozen=True)
class SageSRUnitTriangleConfig:
    points: Dict[str, str]
    roots: Dict[str, str]


class SageSRGeometryEngine:
    """Complex geometry engine backed by Sage's symbolic ring."""

    def __init__(self) -> None:
        self.points: Dict[str, SageSRPointRecord] = {}
        self.constraints: List[Any] = []
        self.value_subs: Dict[Any, Any] = {}
        self.point_assignments: Dict[Any, Any] = {}
        self.learned_subs: Dict[Any, Any] = {}
        self.unit_circle_points: set[str] = set()
        self._main_unit_triangle: Optional[SageSRUnitTriangleConfig] = None
        self._apply_cache: Dict[Any, Any] = {}

        self.add_point("O")
        self._register_origin()

    # ------------------------------------------------------------------
    # Coercion and substitution
    # ------------------------------------------------------------------
    def expr(self, value: Any) -> Any:
        if isinstance(value, str):
            return SR(value.replace("1j", "I").replace("j", "I"))
        return SR(value)

    def _substitution_map(self) -> Dict[Any, Any]:
        subs: Dict[Any, Any] = {}
        subs.update(self.value_subs)
        subs.update(self.point_assignments)
        subs.update(self.learned_subs)
        return subs

    def _apply_all(self, expr: Any) -> Any:
        current = self.expr(expr)
        cached = self._apply_cache.get(current)
        if cached is not None:
            return cached
        substitutions = self._substitution_map()
        if not substitutions:
            return current
        for _ in range(8):
            updated = current.subs(substitutions)
            if updated == current:
                self._apply_cache[current] = updated
                return updated
            current = updated
        self._apply_cache[current] = current
        return current

    def _clear_cache(self) -> None:
        self._apply_cache.clear()

    def _is_zero(self, expr: Any) -> bool:
        return factor(self.expr(expr)) == 0

    def _reduce_rational(self, expr: Any) -> Any:
        value = self.expr(expr)
        numerator = factor(value.numerator())
        denominator = factor(value.denominator())
        if denominator == 1:
            return numerator
        return numerator / denominator

    # ------------------------------------------------------------------
    # Point management
    # ------------------------------------------------------------------
    def add_point(self, name: str) -> None:
        if name in self.points:
            return
        z = SR.var(f"z_{name}")
        zb = SR.var(f"zb_{name}")
        self.points[name] = SageSRPointRecord(z=z, zb=zb)

    def ensure_point(self, name: str) -> None:
        if name not in self.points:
            raise SageSRGeometryError(f"Point '{name}' is not registered.")

    def z_symbol(self, name: str) -> Any:
        self.ensure_point(name)
        return self.points[name].z

    def zb_symbol(self, name: str) -> Any:
        self.ensure_point(name)
        return self.points[name].zb

    def z(self, name: str) -> Any:
        symbol = self.z_symbol(name)
        return self._apply_all(self.point_assignments.get(symbol, symbol))

    def zb(self, name: str) -> Any:
        symbol = self.zb_symbol(name)
        return self._apply_all(self.point_assignments.get(symbol, symbol))

    def _register_origin(self) -> None:
        self.value_subs[self.z_symbol("O")] = SR(0)
        self.value_subs[self.zb_symbol("O")] = SR(0)
        self._clear_cache()

    def set_point_value(self, name: str, z_value: Any, zb_value: Optional[Any] = None) -> None:
        self.add_point(name)
        z_expr = self.expr(z_value)
        zb_expr = self.expr(zb_value) if zb_value is not None else z_expr.conjugate()
        self.value_subs[self.z_symbol(name)] = z_expr
        self.value_subs[self.zb_symbol(name)] = zb_expr
        self._clear_cache()

    def _set_point_assignment(self, name: str, z_expr: Any, zb_expr: Any, *, apply: bool = True) -> None:
        self.add_point(name)
        z_value = self._apply_all(z_expr) if apply else self.expr(z_expr)
        zb_value = self._apply_all(zb_expr) if apply else self.expr(zb_expr)
        self.point_assignments[self.z_symbol(name)] = self._reduce_rational(z_value)
        self.point_assignments[self.zb_symbol(name)] = self._reduce_rational(zb_value)
        self._clear_cache()

    def add_unit_circle(self, name: str) -> None:
        self.ensure_point(name)
        if name in self.unit_circle_points:
            return
        z = self.z_symbol(name)
        zb = self.zb_symbol(name)
        self.constraints.append(z * zb - 1)
        self.learned_subs[zb] = 1 / z
        self.unit_circle_points.add(name)
        self._clear_cache()

    # ------------------------------------------------------------------
    # Predicates and constraints
    # ------------------------------------------------------------------
    def add_constraint(self, constraint: Any) -> None:
        self.constraints.append(self.expr(constraint))
        self._clear_cache()

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
        theta = self.expr(angle)
        w = (I * theta).exp()
        wc = (-I * theta).exp()
        zA, zB, zC = self.z_symbol(A), self.z_symbol(B), self.z_symbol(C)
        zbA, zbB, zbC = self.zb_symbol(A), self.zb_symbol(B), self.zb_symbol(C)
        return wc * (zA - zB) * (zbC - zbB) - w * (zbA - zbB) * (zC - zB)

    def midpoint_polys(self, P: str, Q: str, M: str) -> Tuple[Any, Any]:
        return (
            2 * self.z_symbol(M) - self.z_symbol(P) - self.z_symbol(Q),
            2 * self.zb_symbol(M) - self.zb_symbol(P) - self.zb_symbol(Q),
        )

    def centroid_polys(self, A: str, B: str, C: str, G: str) -> Tuple[Any, Any]:
        return (
            self.z_symbol(A) + self.z_symbol(B) + self.z_symbol(C) - 3 * self.z_symbol(G),
            self.zb_symbol(A) + self.zb_symbol(B) + self.zb_symbol(C) - 3 * self.zb_symbol(G),
        )

    def circumcenter_polys(self, A: str, B: str, C: str, U: str) -> Tuple[Any, Any]:
        zU, zbU = self.z_symbol(U), self.zb_symbol(U)

        def dist(P: str) -> Any:
            return (self.z_symbol(P) - zU) * (self.zb_symbol(P) - zbU)

        return dist(A) - dist(B), dist(A) - dist(C)

    def projection_to_line_polys(self, P: str, A: str, B: str, H: str) -> Tuple[Any, Any]:
        if A == B:
            raise SageSRGeometryError("Projection requires two distinct points to define line AB.")
        return self.perpendicular_poly(A, B, P, H), self.collinear_poly(H, A, B)

    def line_reflection_polys(self, P: str, A: str, B: str, Q: str) -> Tuple[Any, Any]:
        midpoint_z = (self.z_symbol(P) + self.z_symbol(Q)) / 2
        midpoint_zb = (self.zb_symbol(P) + self.zb_symbol(Q)) / 2
        zA, zB = self.z_symbol(A), self.z_symbol(B)
        zbA, zbB = self.zb_symbol(A), self.zb_symbol(B)
        midpoint_on_line = (midpoint_z - zB) * (zbA - zbB) - (midpoint_zb - zbB) * (zA - zB)
        return self.perpendicular_poly(A, B, P, Q), midpoint_on_line

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

    def add_collinear(self, A: str, B: str, C: str) -> None:
        self.add_constraint(self.collinear_poly(A, B, C))

    def add_perpendicular(self, A: str, B: str, C: str, D: str) -> None:
        self.add_constraint(self.perpendicular_poly(A, B, C, D))

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

    def add_isogonal_conjugate(self, A: str, B: str, C: str, P: str, Q: str) -> None:
        for eq in self.isogonal_conjugate_polys(A, B, C, P, Q):
            self.add_constraint(eq)

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
                denominator = dir_z * dir_zb
                if self._is_zero(denominator):
                    raise SageSRGeometryError("Line-circle intersection requires distinct line points.")
                numerator = (base_z - self.z(center)) * dir_zb + (base_zb - self.zb(center)) * dir_z
                t_other = -numerator / denominator
                other = (base_z + t_other * dir_z, base_zb + t_other * dir_zb)
                candidates.extend([(base_z, base_zb), other])
                break

        if not candidates:
            raise SageSRGeometryError(
                "SR mode supports line_circle_intersection only when one line endpoint is already on the circle."
            )

        avoid_set = set(avoid or [])
        if len(candidates) == 2:
            for endpoint in (line_point1, line_point2, radius_point):
                if endpoint in avoid_set and candidates[0] == (self.z(endpoint), self.zb(endpoint)):
                    self._set_point_assignment(name, candidates[1][0], candidates[1][1], apply=False)
                    return self.z(name), self.zb(name)

        avoid_pairs = self._prepare_avoid_pairs(avoid)
        for z_candidate, zb_candidate in candidates:
            if not self._matches_any(z_candidate, zb_candidate, avoid_pairs):
                self._set_point_assignment(name, z_candidate, zb_candidate, apply=False)
                return self.z(name), self.zb(name)
        raise SageSRGeometryError("Line-circle intersection has no candidate left after applying avoid.")

    def set_main_unit_triangle(
        self,
        A: str,
        B: str,
        C: str,
        *,
        root_names: Optional[Sequence[str]] = None,
    ) -> None:
        if len({A, B, C}) != 3:
            raise SageSRGeometryError("Main unit triangle requires three distinct vertex labels.")
        for vertex in (A, B, C):
            self.ensure_point(vertex)
            if vertex not in self.unit_circle_points:
                raise SageSRGeometryError(f"Point '{vertex}' must be declared on the unit circle first.")

        roots_sequence = tuple(root_names) if root_names is not None else (
            self._unit_root_name("A"),
            self._unit_root_name("B"),
            self._unit_root_name("C"),
        )
        if len(roots_sequence) != 3 or len(set(roots_sequence)) != 3:
            raise SageSRGeometryError("root_names must provide three distinct labels.")

        roots: Dict[str, str] = {}
        for canonical, vertex, root in zip(("A", "B", "C"), (A, B, C), roots_sequence):
            self.add_point(root)
            self.add_unit_circle(root)
            self._set_point_assignment(vertex, self.z_symbol(root) ** 2, self.zb_symbol(root) ** 2, apply=False)
            roots[canonical] = root
        self._main_unit_triangle = SageSRUnitTriangleConfig(
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

    def isogonal_conjugate_point(self, A: str, B: str, C: str, P: str, Q: str) -> Any:
        if len({A, B, C}) != 3:
            raise SageSRGeometryError("Isogonal conjugate requires three distinct triangle vertices.")
        self.add_point(Q)
        if P in {A, B, C}:
            self._set_point_assignment(Q, self.z(P), self.zb(P), apply=False)
            return self.z(Q)
        zQ, zbQ = self.z_symbol(Q), self.zb_symbol(Q)
        eq1, eq2, eq3 = self.isogonal_conjugate_polys(A, B, C, P, Q)
        solution = self._solve_linear_point_system(self._apply_all(eq1), self._apply_all(eq2), zQ, zbQ, "Isogonal conjugate")
        check = self._apply_all(eq3).subs({zQ: solution[0], zbQ: solution[1]})
        if not self._is_zero(check):
            raise SageSRGeometryError("Isogonal conjugate construction failed consistency check.")
        self._set_point_assignment(Q, *solution, apply=False)
        return self.z(Q)

    # ------------------------------------------------------------------
    # Linear solving and checks
    # ------------------------------------------------------------------
    def _solve_linear_point_system(self, eq1: Any, eq2: Any, z_var: Any, zb_var: Any, label: str) -> Tuple[Any, Any]:
        p1 = self.expr(eq1).numerator()
        p2 = self.expr(eq2).numerator()
        if not self._is_linear_in_pair(p1, z_var, zb_var) or not self._is_linear_in_pair(p2, z_var, zb_var):
            raise SageSRGeometryError(f"{label} construction only supports linear systems.")
        a1, b1, c1 = self._linear_coefficients(p1, z_var, zb_var)
        a2, b2, c2 = self._linear_coefficients(p2, z_var, zb_var)
        det = a1 * b2 - a2 * b1
        if self._is_zero(det):
            raise SageSRGeometryError(f"{label} construction failed: equations do not determine a unique point.")
        return (b1 * c2 - b2 * c1) / det, (c1 * a2 - c2 * a1) / det

    def _is_linear_in_pair(self, expr: Any, z_var: Any, zb_var: Any) -> bool:
        return (
            self.expr(expr).derivative(z_var, 2) == 0
            and self.expr(expr).derivative(zb_var, 2) == 0
            and self.expr(expr).derivative(z_var).derivative(zb_var) == 0
        )

    def _linear_coefficients(self, expr: Any, z_var: Any, zb_var: Any) -> Tuple[Any, Any, Any]:
        zero_subs = {z_var: 0, zb_var: 0}
        a = self.expr(expr).derivative(z_var).subs(zero_subs)
        b = self.expr(expr).derivative(zb_var).subs(zero_subs)
        c = self.expr(expr).subs(zero_subs)
        return a, b, c

    def _prepare_avoid_pairs(self, avoid: Optional[Sequence[str]]) -> List[Tuple[Any, Any]]:
        if avoid is None:
            return []
        return [(self.z(label), self.zb(label)) for label in avoid]

    def _matches_any(self, z_candidate: Any, zb_candidate: Any, pairs: Sequence[Tuple[Any, Any]]) -> bool:
        for z_avoid, zb_avoid in pairs:
            if self._is_zero(z_candidate - z_avoid) and self._is_zero(zb_candidate - zb_avoid):
                return True
        return False

    def _require_main_unit_triangle(self) -> SageSRUnitTriangleConfig:
        if self._main_unit_triangle is None:
            raise SageSRGeometryError("Main unit triangle is not configured. Use set_main_unit_triangle first.")
        return self._main_unit_triangle

    def _unit_root_name(self, vertex: str) -> str:
        return f"__unit_root_{vertex}"

    def constraint_polynomials(self, constraint: str, args: Sequence[str], *, angle: Optional[Any] = None) -> List[Any]:
        normalized = constraint.lower()
        if normalized in {"collinear", "line"}:
            return [self.collinear_poly(*args)]
        if normalized in {"perpendicular", "perp"}:
            return [self.perpendicular_poly(*args)]
        if normalized in {"angle", "angle_value"}:
            if angle is None:
                raise SageSRGeometryError("Angle constraint requires an angle.")
            return [self.angle_value_poly(*args, angle)]
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
        raise SageSRGeometryError(f"Unsupported SR constraint '{constraint}'.")

    def constraint_conjugate_free(self, constraint: str, args: Sequence[str], *, angle: Optional[Any] = None) -> List[Tuple[Any, Any]]:
        results: List[Tuple[Any, Any]] = []
        for poly in self.constraint_polynomials(constraint, args, angle=angle):
            substituted = self._apply_all(poly)
            numerator = factor(substituted.numerator())
            denominator = substituted.denominator()
            results.append((numerator, denominator))
        return results

    def squared_distance(self, P: str, Q: str) -> Any:
        return self._apply_all((self.z(P) - self.z(Q)) * (self.zb(P) - self.zb(Q)))

    def point_summary(self, names: Optional[Sequence[str]] = None) -> Dict[str, str]:
        selected = sorted(self.points.keys()) if names is None else names
        return {name: str(self.z(name)) for name in selected}

    def constraint_strings(self) -> List[str]:
        return [str(expr) for expr in self.constraints]

    def learned_rules(self) -> Dict[str, Any]:
        return {str(symbol): expr for symbol, expr in self.learned_subs.items()}


__all__ = ["SageSRGeometryEngine", "SageSRGeometryError"]
