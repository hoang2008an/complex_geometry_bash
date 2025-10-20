"""
Symbolic geometry engine for complex-number bashing in the circumcenter model.

The implementation follows the behavior spec outlined in AGENTS.md and the
user-provided requirements.  Points are represented by independent symbols
(z_X, zb_X).  Constraints are maintained as polynomials, and the engine learns
conjugate substitution rules whenever a single-conjugate equation is detected.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import sympy as sp


class GeometryError(RuntimeError):
    """Raised when a construction or solve step fails."""


@dataclass
class PointRecord:
    z: sp.Symbol
    zb: sp.Symbol


class GeometryEngine:
    """
    Geometry Engine working in the complex plane with symbolic coordinates.

    Each point has two independent symbols (z_X, zb_X).  The engine stores
    constraints as polynomials and maintains a substitution system composed of:

    * Value substitutions (e.g. origin at 0)
    * Assignments produced by constructions (rational expressions)
    * Unit-circle substitutions (zb_U -> 1/z_U)
    * Auto-learned conjugate collapses from single-conjugate constraints
    """

    def __init__(self) -> None:
        self.points: Dict[str, PointRecord] = {}
        self.constraints: List[sp.Expr] = []
        self.value_subs: Dict[sp.Symbol, sp.Expr] = {}
        self.point_assignments: Dict[sp.Symbol, sp.Expr] = {}
        self.learned_subs: Dict[sp.Symbol, sp.Expr] = {}

        # Declare origin O with fixed coordinates at 0
        self.add_point("O")
        self._register_origin()

    # ------------------------------------------------------------------
    # Point management helpers
    # ------------------------------------------------------------------
    def add_point(self, name: str) -> None:
        """
        Register a point with independent (z, zb) symbols.

        Parameters
        ----------
        name:
            Label of the point (must be unique).
        """
        if name in self.points:
            return

        z_symbol = sp.Symbol(f"z_{name}")
        zb_symbol = sp.Symbol(f"zb_{name}")
        self.points[name] = PointRecord(z=z_symbol, zb=zb_symbol)


    def ensure_point(self, name: str) -> None:
        """Ensure a point is registered; raise an informative error otherwise."""
        if name not in self.points:
            raise GeometryError(f"Point '{name}' is not registered.")

    def _register_origin(self) -> None:
        """Apply the origin substitution rules (O fixed at 0)."""
        self.ensure_point("O")
        origin = self.points["O"]
        self.value_subs[origin.z] = sp.Integer(0)
        self.value_subs[origin.zb] = sp.Integer(0)

    def add_unit_circle(self, name: str) -> None:
        """Declare that point name lies on the unit circle."""
        self.ensure_point(name)
        record = self.points[name]
        if record.zb in self.learned_subs:
            return
        constraint = record.z * record.zb - 1
        self.add_constraint(constraint)

    # ------------------------------------------------------------------
    # Symbol and expression helpers
    # ------------------------------------------------------------------
    def z_symbol(self, name: str) -> sp.Symbol:
        self.ensure_point(name)
        return self.points[name].z

    def zb_symbol(self, name: str) -> sp.Symbol:
        self.ensure_point(name)
        return self.points[name].zb

    def z(self, name: str) -> sp.Expr:
        """Return the (possibly substituted) expression for z_name."""
        symbol = self.z_symbol(name)
        expr = self.point_assignments.get(symbol, symbol)
        return self._apply_all(expr)

    def zb(self, name: str) -> sp.Expr:
        """Return the (possibly substituted) expression for zb_name."""
        symbol = self.zb_symbol(name)
        expr = self.point_assignments.get(symbol, symbol)
        return self._apply_all(expr)

    def _substitution_map(self) -> Dict[sp.Symbol, sp.Expr]:
        """Compose the global substitution map in precedence order."""
        subs: Dict[sp.Symbol, sp.Expr] = {}
        subs.update(self.value_subs)
        subs.update(self.point_assignments)
        subs.update(self.learned_subs)
        return subs

    def _apply_all(self, expr: sp.Expr) -> sp.Expr:
        """
        Apply all known substitutions (value, assignments, unit-circle, learned).

        Multiple passes are used to ensure transitive substitutions settle.
        """
        if expr is None:
            return expr

        current = sp.simplify(expr)
        if current == 0:
            return sp.Integer(0)

        subs = self._substitution_map()
        for _ in range(4):
            updated = current.subs(subs, simultaneous=False)
            updated = sp.simplify(updated)
            if updated == current:
                break
            current = updated
        return sp.simplify(current)

    # ------------------------------------------------------------------
    # Constraint handling and auto-learning
    # ------------------------------------------------------------------
    def add_constraint(self, constraint: sp.Expr) -> None:
        """
        Record a polynomial constraint (assumed == 0) and attempt auto-learning.
        """
        raw_constraint = sp.expand(constraint)
        self.constraints.append(raw_constraint)
        processed = self._apply_all(raw_constraint)
        self._auto_learn_from_constraint(processed)

    def _auto_learn_from_constraint(self, constraint: sp.Expr) -> None:
        """
        Implement the single-equation conjugate collapse:

        1. Apply substitutions (already handled before getting here).
        2. Convert to a numerator-only polynomial (clear denominators).
        3. If the numerator contains exactly one conjugate symbol zb_X,
           solve for it and record the substitution if conjugate-free.
        """
        if constraint == 0:
            return

        together_expr = sp.together(constraint)
        numerator, _ = sp.fraction(together_expr)
        numerator = sp.expand(numerator)

        zb_symbols = [symbol for symbol in numerator.free_symbols if isinstance(symbol, sp.Symbol) and str(symbol).startswith("zb_")]
        if len(zb_symbols) != 1:
            return

        target = zb_symbols[0]
        solution = sp.solve(sp.Eq(numerator, 0), target, dict=True)
        if not solution:
            return

        rhs = sp.simplify(solution[0][target])
        if any(str(sym).startswith("zb_") for sym in rhs.free_symbols):
            return

        # Record the learned substitution and shrink existing constraints.
        self.learned_subs[target] = rhs

    # ------------------------------------------------------------------
    # Predicate polynomials
    # ------------------------------------------------------------------
    def collinear_poly(self, A: str, B: str, C: str, *, raw: bool = False) -> sp.Expr:
        """
        Return the polynomial witnessing collinearity of points A, B, C.
        """
        zA, zB, zC = self.z_symbol(A), self.z_symbol(B), self.z_symbol(C)
        zbA, zbB, zbC = self.zb_symbol(A), self.zb_symbol(B), self.zb_symbol(C)
        poly = (zA - zC) * (zbB - zbC) - (zbA - zbC) * (zB - zC)
        return poly if raw else sp.expand(poly)

    def perpendicular_poly(self, A: str, B: str, C: str, D: str, *, raw: bool = False) -> sp.Expr:
        """
        Return the polynomial encoding AB ⟂ CD.
        """
        zA, zB = self.z_symbol(A), self.z_symbol(B)
        zC, zD = self.z_symbol(C), self.z_symbol(D)
        zbA, zbB = self.zb_symbol(A), self.zb_symbol(B)
        zbC, zbD = self.zb_symbol(C), self.zb_symbol(D)
        poly = (zA - zB) * (zbC - zbD) + (zbA - zbB) * (zC - zD)
        return poly if raw else sp.expand(poly)

    def concyclic_poly(self, A: str, B: str, C: str, D: str, *, raw: bool = False) -> sp.Expr:
        """
        Return a polynomial enforcing that A, B, C, D lie on a common circle.

        The polynomial is obtained by eliminating the circle center U from the
        system |A-U|^2 = |B-U|^2, |A-U|^2 = |C-U|^2, |A-U|^2 = |D-U|^2. The first
        two equations are solved for U; the solution is substituted into the
        third, and denominators are cleared to yield a single polynomial in the
        point coordinates.
        """
        zA, zB, zC, zD = (self.z_symbol(A), self.z_symbol(B),
                          self.z_symbol(C), self.z_symbol(D))
        zbA, zbB, zbC, zbD = (self.zb_symbol(A), self.zb_symbol(B),
                              self.zb_symbol(C), self.zb_symbol(D))

        zU = sp.Symbol(f"z_circ_{A}{B}{C}", complex=True)
        zbU = sp.Symbol(f"zb_circ_{A}{B}{C}", complex=True)

        eq1 = sp.expand((zA - zU) * (zbA - zbU) - (zB - zU) * (zbB - zbU))
        eq2 = sp.expand((zA - zU) * (zbA - zbU) - (zC - zU) * (zbC - zbU))

        solutions = sp.solve((sp.Eq(eq1, 0), sp.Eq(eq2, 0)), (zU, zbU), dict=True)
        if not solutions:
            raise GeometryError("Failed to derive circumcenter for concyclicity test.")

        solution = solutions[0]
        zU_expr = sp.simplify(solution[zU])
        zbU_expr = sp.simplify(solution[zbU])

        eq3 = sp.simplify((zA - zU_expr) * (zbA - zbU_expr) - (zD - zU_expr) * (zbD - zbU_expr))
        numerator, denominator = sp.fraction(sp.simplify(eq3))
        poly = sp.expand(numerator)
        return poly if raw else sp.expand(poly)

    def angle_value_poly(self, A: str, B: str, C: str, angle_radians: sp.Expr, *, raw: bool = False) -> sp.Expr:
        """
        Return polynomial enforcing that angle ABC equals the provided angle (in radians).

        Uses the identity that (z_A - z_B)/(z_C - z_B) / exp(i*theta) is real.
        """
        zA, zB, zC = self.z_symbol(A), self.z_symbol(B), self.z_symbol(C)
        zbA, zbB, zbC = self.zb_symbol(A), self.zb_symbol(B), self.zb_symbol(C)
        w = sp.exp(sp.I * angle_radians)
        w_conj = sp.exp(-sp.I * angle_radians)
        poly = w_conj * (zA - zB) * (zbC - zbB) - w * (zbA - zbB) * (zC - zB)
        return poly if raw else sp.expand(poly)

    def circumcenter_polys(self, A: str, B: str, C: str, U: str, *, raw: bool = False) -> Tuple[sp.Expr, sp.Expr]:
        """
        Return the two polynomials that enforce U as the circumcenter of triangle ABC.

        |A-U|^2 = |B-U|^2 and |A-U|^2 = |C-U|^2, expressed as polynomial equalities.
        """
        zU = self.z_symbol(U)
        zbU = self.zb_symbol(U)

        def squared_distance(zX: sp.Expr, zbX: sp.Expr) -> sp.Expr:
            return (zX - zU) * (zbX - zbU)

        zA, zB, zC = self.z_symbol(A), self.z_symbol(B), self.z_symbol(C)
        zbA, zbB, zbC = self.zb_symbol(A), self.zb_symbol(B), self.zb_symbol(C)

        dist_A = squared_distance(zA, zbA)
        dist_B = squared_distance(zB, zbB)
        dist_C = squared_distance(zC, zbC)

        eq1 = dist_A - dist_B
        eq2 = dist_A - dist_C
        if not raw:
            eq1 = sp.expand(eq1)
            eq2 = sp.expand(eq2)
        return eq1, eq2

    def midpoint_polys(self, P: str, Q: str, M: str, *, raw: bool = False) -> Tuple[sp.Expr, sp.Expr]:
        """
        Polynomials enforcing that M is the midpoint of segment PQ.
        """
        zP, zQ, zM = self.z_symbol(P), self.z_symbol(Q), self.z_symbol(M)
        zbP, zbQ, zbM = self.zb_symbol(P), self.zb_symbol(Q), self.zb_symbol(M)
        eq_z = 2 * zM - zP - zQ
        eq_zb = 2 * zbM - zbP - zbQ
        if not raw:
            eq_z = sp.expand(eq_z)
            eq_zb = sp.expand(eq_zb)
        return eq_z, eq_zb

    def point_reflection_polys(self, P: str, O: str, Q: str, *, raw: bool = False) -> Tuple[sp.Expr, sp.Expr]:
        """
        Polynomials enforcing that Q is the reflection of P across point O.
        """
        zP, zO, zQ = self.z_symbol(P), self.z_symbol(O), self.z_symbol(Q)
        zbP, zbO, zbQ = self.zb_symbol(P), self.zb_symbol(O), self.zb_symbol(Q)
        eq_z = zQ + zP - 2 * zO
        eq_zb = zbQ + zbP - 2 * zbO
        if not raw:
            eq_z = sp.expand(eq_z)
            eq_zb = sp.expand(eq_zb)
        return eq_z, eq_zb

    def line_reflection_polys(self, P: str, A: str, B: str, Q: str, *, raw: bool = False) -> Tuple[sp.Expr, ...]:
        """
        Polynomials enforcing that Q is the reflection of P across line AB.
        """
        eq_perp = self.perpendicular_poly(A, B, P, Q, raw=True)

        midpoint_label = self._midpoint_label(P, Q)
        self.add_point(midpoint_label)
        eq_mid_z, eq_mid_zb = self.midpoint_polys(P, Q, midpoint_label, raw=True)
        collinear_midpoint = self.collinear_poly(midpoint_label, A, B, raw=True)

        if not raw:
            eq_perp = sp.expand(eq_perp)
            eq_mid_z = sp.expand(eq_mid_z)
            eq_mid_zb = sp.expand(eq_mid_zb)
            collinear_midpoint = sp.expand(collinear_midpoint)

        return eq_perp, eq_mid_z, eq_mid_zb, collinear_midpoint

    # ------------------------------------------------------------------
    # Constraint guards
    # ------------------------------------------------------------------
    def add_collinear(self, A: str, B: str, C: str) -> None:
        """Store the collinearity constraint for (A, B, C)."""
        poly = self.collinear_poly(A, B, C)
        self.add_constraint(poly)

    def add_perpendicular(self, A: str, B: str, C: str, D: str) -> None:
        """Store the perpendicularity constraint AB ⟂ CD."""
        poly = self.perpendicular_poly(A, B, C, D)
        self.add_constraint(poly)

    def add_concyclic(self, A: str, B: str, C: str, D: str) -> None:
        """Store the concyclicity constraint for points A, B, C, D."""
        poly = self.concyclic_poly(A, B, C, D)
        self.add_constraint(poly)

    def add_angle_value(self, A: str, B: str, C: str, angle_radians: sp.Expr) -> None:
        """Store the constraint that angle ABC equals the specified angle (in radians)."""
        poly = self.angle_value_poly(A, B, C, angle_radians)
        self.add_constraint(poly)

    def add_circumcenter(self, A: str, B: str, C: str, U: str) -> None:
        """
        Store the circumcenter constraints ensuring U is the circumcenter of triangle ABC.
        """
        self.add_point(U)
        eq1, eq2 = self.circumcenter_polys(A, B, C, U)
        self.add_constraint(eq1)
        self.add_constraint(eq2)

    def add_midpoint(self, P: str, Q: str, M: str) -> None:
        """Store the constraint that M is the midpoint of segment PQ."""
        self.add_point(M)
        eq_z, eq_zb = self.midpoint_polys(P, Q, M)
        self.add_constraint(eq_z)
        self.add_constraint(eq_zb)

    def add_point_reflection(self, P: str, O: str, Q: str) -> None:
        """Store the constraint that Q is the reflection of P across point O."""
        eq_z, eq_zb = self.point_reflection_polys(P, O, Q)
        self.add_constraint(eq_z)
        self.add_constraint(eq_zb)

    def add_line_reflection(self, P: str, A: str, B: str, Q: str) -> None:
        """Store the constraint that Q is the reflection of P across line AB."""
        for eq in self.line_reflection_polys(P, A, B, Q):
            self.add_constraint(eq)

    def circumcenter(self, A: str, B: str, C: str, U: str) -> sp.Expr:
        """
        Construct the circumcenter U of triangle ABC by solving the circumcenter equations.
        """
        self.add_point(U)
        zU, zbU = self.z_symbol(U), self.zb_symbol(U)
        if zU in self.point_assignments and zbU in self.point_assignments:
            return self.z(U)
        eq1, eq2 = self.circumcenter_polys(A, B, C, U, raw=True)
        eq1 = self._apply_all(eq1)
        eq2 = self._apply_all(eq2)
        solutions = sp.solve([eq1, eq2], (zU, zbU), dict=True)
        if not solutions:
            raise GeometryError("Circumcenter construction failed: system is underdetermined.")

        solution = solutions[0]
        self._set_point_assignment(U, solution[zU], solution[zbU])
        return self.z(U)

    def midpoint(self, P: str, Q: str, M: str) -> sp.Expr:
        """
        Construct the midpoint M of segment PQ.
        """
        self.add_point(M)
        z_expr = sp.simplify((self.z(P) + self.z(Q)) / 2)
        zb_expr = sp.simplify((self.zb(P) + self.zb(Q)) / 2)
        self._set_point_assignment(M, z_expr, zb_expr)
        return self.z(M)

    def euler_center(
        self,
        A: str,
        B: str,
        C: str,
        N: str,
        circumcenter_name: Optional[str] = None,
        orthocenter_name: Optional[str] = None,
    ) -> sp.Expr:
        """
        Construct the Euler (nine-point) center N as midpoint of circumcenter and orthocenter.

        Optional circumcenter/orthocenter labels allow callers to reuse existing points; otherwise
        temporary internal names are introduced for the construction.
        """
        circ_label = circumcenter_name or self._unique_internal_name("EulerCirc")
        orth_label = orthocenter_name or self._unique_internal_name("EulerOrth")

        self.circumcenter(A, B, C, circ_label)
        self.orthocenter_via_altitudes(A, B, C, orth_label)
        self.add_point(N)

        z_expr = sp.simplify((self.z(circ_label) + self.z(orth_label)) / 2)
        zb_expr = sp.simplify((self.zb(circ_label) + self.zb(orth_label)) / 2)
        self._set_point_assignment(N, z_expr, zb_expr)
        return self.z(N)

    def lemoine_point(
        self,
        A: str,
        B: str,
        C: str,
        L: str,
        circumcenter_name: Optional[str] = None,
    ) -> sp.Expr:
        """
        Construct the Lemoine (symmedian) point via intersections of circumcircle tangents.
        """
        circ_label = circumcenter_name or self._unique_internal_name("LemoineCirc")

        self.circumcenter(A, B, C, circ_label)

        tangent_BC = self._unique_internal_name(f"T_{B}{C}")
        self._tangent_intersection(B, C, circ_label, tangent_BC)

        tangent_CA = self._unique_internal_name(f"T_{C}{A}")
        self._tangent_intersection(C, A, circ_label, tangent_CA)

        self.intersection_of_lines(A, tangent_BC, B, tangent_CA, L)
        return self.z(L)

    def reflect_point_over_point(self, P: str, O: str, Q: str) -> sp.Expr:
        """
        Construct the reflection of point P across point O and assign it to Q.
        """
        self.add_point(Q)
        z_expr = sp.simplify(2 * self.z(O) - self.z(P))
        zb_expr = sp.simplify(2 * self.zb(O) - self.zb(P))
        self._set_point_assignment(Q, z_expr, zb_expr)
        return self.z(Q)

    def reflect_point_over_line(self, P: str, A: str, B: str, Q: str) -> sp.Expr:
        """
        Construct the reflection of point P across line AB and assign it to Q.
        """
        self.add_point(Q)
        midpoint_label = self._midpoint_label(P, Q)
        self.add_point(midpoint_label)

        eq_perp, eq_mid_z, eq_mid_zb, eq_collinear = self.line_reflection_polys(P, A, B, Q, raw=True)
        equations = [eq_perp, eq_mid_z, eq_mid_zb, eq_collinear]

        zQ, zbQ = self.z_symbol(Q), self.zb_symbol(Q)
        zM, zbM = self.z_symbol(midpoint_label), self.zb_symbol(midpoint_label)
        processed = [self._apply_all(eq) for eq in equations]
        solutions = sp.solve(processed, (zQ, zbQ, zM, zbM), dict=True)
        if not solutions:
            raise GeometryError("Line reflection construction failed: system has no solution.")

        solution = solutions[0]
        self._set_point_assignment(Q, solution[zQ], solution[zbQ])
        self._set_point_assignment(midpoint_label, solution[zM], solution[zbM])
        return self.z(Q)

    # ------------------------------------------------------------------
    # Constructions
    # ------------------------------------------------------------------
    def orthocenter_via_altitudes(self, T1: str, T2: str, T3: str, H: str) -> Tuple[sp.Expr, sp.Expr]:
        """
        Construct the orthocenter H of triangle T1T2T3 by solving the two altitude constraints.
        """
        self.add_point(H)
        zH, zbH = self.z_symbol(H), self.zb_symbol(H)

        eq1 = self._apply_all(self.perpendicular_poly(T1, H, T2, T3))
        eq2 = self._apply_all(self.perpendicular_poly(T2, H, T1, T3))
        solutions = sp.solve([eq1, eq2], (zH, zbH), dict=True)
        if not solutions:
            raise GeometryError("Orthocenter construction failed: system is underdetermined.")

        solution = solutions[0]
        self._set_point_assignment(H, solution[zH], solution[zbH])
        return self.z(H), self.zb(H)

    def intersection_of_lines(self, A: str, B: str, C: str, D: str, X: str) -> Tuple[sp.Expr, sp.Expr]:
        """
        Intersection point X of lines AB and CD.
        """
        self.add_point(X)
        zX, zbX = self.z_symbol(X), self.zb_symbol(X)

        eq1 = self._apply_all(self.collinear_poly(X, A, B))
        eq2 = self._apply_all(self.collinear_poly(X, C, D))
        solutions = sp.solve([eq1, eq2], (zX, zbX), dict=True)
        if not solutions:
            raise GeometryError("Intersection construction failed: lines do not determine a unique point.")

        solution = solutions[0]
        self._set_point_assignment(X, solution[zX], solution[zbX])
        return self.z(X), self.zb(X)

    def line_circle_intersection(
        self,
        line_point1: str,
        line_point2: str,
        center: str,
        radius_point: str,
        name: str,
        *,
        avoid: Optional[Sequence[str]] = None,
    ) -> Tuple[sp.Expr, sp.Expr]:
        """
        Construct the intersection of line (line_point1 line_point2) with the circle centered at `center`
        passing through `radius_point`.

        Parameters
        ----------
        line_point1, line_point2:
            Points defining the target line.
        center:
            Circle center.
        radius_point:
            Point fixing the circle radius (must be distinct from center).
        name:
            Label to assign the constructed intersection.
        avoid:
            Optional list of point labels to skip if multiple intersections exist.
        """
        self.add_point(name)
        self.ensure_point(line_point1)
        self.ensure_point(line_point2)
        self.ensure_point(center)
        self.ensure_point(radius_point)

        z_center, zb_center = self.z_symbol(center), self.zb_symbol(center)
        z_radius, zb_radius = self.z_symbol(radius_point), self.zb_symbol(radius_point)
        z_var, zb_var = self.z_symbol(name), self.zb_symbol(name)

        eq_line = self.collinear_poly(name, line_point1, line_point2, raw=True)
        eq_circle = (z_var - z_center) * (zb_var - zb_center) - (z_radius - z_center) * (zb_radius - zb_center)

        return self._solve_point_from_equations(name, [eq_line, eq_circle], avoid=avoid)

    def centroid(self, A: str, B: str, C: str, G: str) -> sp.Expr:
        """
        Assign centroid G of triangle ABC.
        """
        self.add_point(G)
        z_expr = sp.simplify((self.z(A) + self.z(B) + self.z(C)) / 3)
        zb_expr = sp.simplify((self.zb(A) + self.zb(B) + self.zb(C)) / 3)
        self._set_point_assignment(G, z_expr, zb_expr)
        return self.z(G)

    def squared_distance(self, P: str, Q: str) -> sp.Expr:
        """
        Return the squared distance |P - Q|^2 using current substitutions.
        """
        self.ensure_point(P)
        self.ensure_point(Q)
        diff_z = sp.simplify(self.z(P) - self.z(Q))
        diff_zb = sp.simplify(self.zb(P) - self.zb(Q))
        return sp.simplify(diff_z * diff_zb)

    def _set_point_assignment(self, name: str, z_expr: sp.Expr, zb_expr: sp.Expr) -> None:
        """Register rational expressions for a constructed point."""
        z_symbol = self.z_symbol(name)
        zb_symbol = self.zb_symbol(name)
        self.point_assignments[z_symbol] = sp.simplify(sp.together(z_expr))
        self.point_assignments[zb_symbol] = sp.simplify(sp.together(zb_expr))

    def _unique_internal_name(self, base: str) -> str:
        """Return a point name unlikely to clash with user-provided labels."""
        sanitized = base.replace(" ", "_")
        candidate = f"__{sanitized}"
        index = 1
        while candidate in self.points:
            candidate = f"__{sanitized}_{index}"
            index += 1
        return candidate

    def _midpoint_label(self, P: str, Q: str) -> str:
        """Deterministic internal name for the midpoint of segment PQ."""
        ordered = tuple(sorted((P, Q)))
        return f"__mid_{ordered[0]}_{ordered[1]}"

    def _solve_point_from_equations(
        self,
        name: str,
        equations: Sequence[sp.Expr],
        *,
        avoid: Optional[Sequence[str]] = None,
    ) -> Tuple[sp.Expr, sp.Expr]:
        """Solve a system for point `name` using the provided equations."""
        self.add_point(name)
        z_var, zb_var = self.z_symbol(name), self.zb_symbol(name)
        processed = [self._apply_all(eq) for eq in equations]
        solutions = sp.solve(processed, (z_var, zb_var), dict=True)
        if not solutions:
            raise GeometryError(f"Construction for point '{name}' failed: system has no solution.")
        chosen = None
        avoid_pairs: List[Tuple[sp.Expr, sp.Expr]] = []
        if avoid:
            for label in avoid:
                self.ensure_point(label)
                avoid_pairs.append((self.z(label), self.zb(label)))
        for candidate in solutions:
            z_candidate = sp.simplify(candidate[z_var])
            zb_candidate = sp.simplify(candidate[zb_var])
            if avoid_pairs:
                skip = False
                for z_avoid, zb_avoid in avoid_pairs:
                    if sp.simplify(z_candidate - z_avoid) == 0 and sp.simplify(zb_candidate - zb_avoid) == 0:
                        skip = True
                        break
                if skip:
                    continue
            chosen = candidate
            break
        if chosen is None:
            chosen = solutions[0]
        self._set_point_assignment(name, chosen[z_var], chosen[zb_var])
        return self.z(name), self.zb(name)

    def _tangent_intersection(
        self,
        P: str,
        Q: str,
        center: str,
        name: str,
    ) -> Tuple[sp.Expr, sp.Expr]:
        """
        Solve for the intersection of tangents at P and Q to the circumcircle centered at `center`.
        """
        self.add_point(name)
        eq1 = self.perpendicular_poly(P, name, P, center)
        eq2 = self.perpendicular_poly(Q, name, Q, center)
        return self._solve_point_from_equations(name, [eq1, eq2])

    # ------------------------------------------------------------------
    # Constraint inspection helpers
    # ------------------------------------------------------------------
    def _conjugate_free_expr(self, expr: sp.Expr) -> Tuple[sp.Expr, sp.Expr]:
        """
        Reduce a constraint polynomial and return the simplified numerator/denominator.

        The result is simplified as much as possible; residual conjugate symbols may
        remain if substitutions cannot eliminate them completely.
        """
        substituted = self._apply_all(expr)
        together_expr = sp.together(substituted)
        numerator, denominator = sp.fraction(together_expr)
        numerator = sp.expand(numerator)
        denominator = sp.expand(denominator)

        return numerator, denominator

    def _constraint_polynomials(
        self,
        constraint: str,
        args: Sequence[str],
        *,
        angle: Optional[sp.Expr] = None,
    ) -> List[sp.Expr]:
        """
        Resolve the requested constraint into the underlying polynomial(s).
        """
        normalized = constraint.lower()
        if normalized in {"collinear", "line"}:
            if len(args) != 3:
                raise GeometryError("Collinear constraint expects exactly three point labels.")
            return [self.collinear_poly(*args)]
        if normalized in {"perpendicular", "perp"}:
            if len(args) != 4:
                raise GeometryError("Perpendicular constraint expects exactly four point labels.")
            return [self.perpendicular_poly(*args)]
        if normalized in {"concyclic", "cyclic"}:
            if len(args) != 4:
                raise GeometryError("Concyclic constraint expects exactly four point labels.")
            return [self.concyclic_poly(*args)]
        if normalized in {"angle", "angle_value"}:
            if len(args) != 3:
                raise GeometryError("Angle constraint expects three point labels (A, B, C).")
            if angle is None:
                raise GeometryError("Angle constraint requires an 'angle' expression.")
            return [self.angle_value_poly(*args, angle, raw=False)]
        if normalized == "circumcenter":
            if len(args) != 4:
                raise GeometryError("Circumcenter constraint expects four point labels (A, B, C, U).")
            eq1, eq2 = self.circumcenter_polys(*args)
            return [eq1, eq2]
        if normalized in {"midpoint"}:
            if len(args) != 3:
                raise GeometryError("Midpoint constraint expects three point labels (P, Q, M).")
            eq_z, eq_zb = self.midpoint_polys(*args)
            return [eq_z, eq_zb]
        if normalized in {"point_reflection", "reflection_point", "reflect_point"}:
            if len(args) != 3:
                raise GeometryError("Point reflection expects three point labels (P, O, Q).")
            eq_z, eq_zb = self.point_reflection_polys(*args)
            return [eq_z, eq_zb]
        if normalized in {"line_reflection", "reflection_line", "reflect_line"}:
            if len(args) != 4:
                raise GeometryError("Line reflection expects four point labels (P, A, B, Q).")
            return list(self.line_reflection_polys(*args))

        raise GeometryError(f"Unsupported constraint '{constraint}'.")

    def constraint_conjugate_free(
        self,
        constraint: str,
        args: Sequence[str],
        *,
        angle: Optional[sp.Expr] = None,
    ) -> List[Tuple[sp.Expr, sp.Expr]]:
        """
        Return conjugate-free numerator/denominator pairs for the chosen constraint.
        """
        polys = self._constraint_polynomials(constraint, args, angle=angle)
        return [self._conjugate_free_expr(poly) for poly in polys]

    def perpendicular_conjugate_free(self, A: str, B: str, C: str, D: str) -> Tuple[sp.Expr, sp.Expr]:
        """
        Backwards-compatible helper for AB ⟂ CD checks.
        """
        return self.constraint_conjugate_free("perpendicular", [A, B, C, D])[0]

    # ------------------------------------------------------------------
    # Introspection utilities
    # ------------------------------------------------------------------
    def learned_rules(self) -> Dict[str, sp.Expr]:
        """Return learned conjugate substitutions as a plain dictionary."""
        return {str(symbol): sp.simplify(expr) for symbol, expr in self.learned_subs.items()}

    def point_summary(self, names: Optional[Sequence[str]] = None) -> Dict[str, sp.Expr]:
        """Return simplified z-coordinate expressions for selected points."""
        if names is None:
            names = sorted(self.points.keys())
        summary: Dict[str, sp.Expr] = {}
        for name in names:
            summary[name] = sp.simplify(self.z(name))
        return summary

    def constraint_strings(self) -> List[str]:
        """Return stored constraints in their recorded polynomial form."""
        return [str(sp.expand(constraint)) for constraint in self.constraints]


__all__ = [
    "GeometryEngine",
    "GeometryError",
]
