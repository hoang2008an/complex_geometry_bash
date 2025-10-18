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
        Return the polynomial enforcing that A, B, C, D lie on a common circle.

        Derived from the cross-ratio being real:
        (z_A - z_C)(zb_A - zb_C)(z_B - z_D)(zb_B - zb_D)
        -
        (z_A - z_D)(zb_A - zb_D)(z_B - z_C)(zb_B - zb_C) = 0
        """
        zA, zB, zC, zD = (self.z_symbol(A), self.z_symbol(B),
                          self.z_symbol(C), self.z_symbol(D))
        zbA, zbB, zbC, zbD = (self.zb_symbol(A), self.zb_symbol(B),
                              self.zb_symbol(C), self.zb_symbol(D))
        term1 = (zA - zC) * (zbA - zbC) * (zB - zD) * (zbB - zbD)
        term2 = (zA - zD) * (zbA - zbD) * (zB - zC) * (zbB - zbC)
        poly = term1 - term2
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

    def centroid(self, A: str, B: str, C: str, G: str) -> sp.Expr:
        """
        Assign centroid G of triangle ABC.
        """
        self.add_point(G)
        z_expr = sp.simplify((self.z(A) + self.z(B) + self.z(C)) / 3)
        zb_expr = sp.simplify((self.zb(A) + self.zb(B) + self.zb(C)) / 3)
        self._set_point_assignment(G, z_expr, zb_expr)
        return self.z(G)

    def _set_point_assignment(self, name: str, z_expr: sp.Expr, zb_expr: sp.Expr) -> None:
        """Register rational expressions for a constructed point."""
        z_symbol = self.z_symbol(name)
        zb_symbol = self.zb_symbol(name)
        self.point_assignments[z_symbol] = sp.simplify(sp.together(z_expr))
        self.point_assignments[zb_symbol] = sp.simplify(sp.together(zb_expr))

    # ------------------------------------------------------------------
    # Conjugate-free perpendicular check
    # ------------------------------------------------------------------
    def perpendicular_conjugate_free(self, A: str, B: str, C: str, D: str) -> Tuple[sp.Expr, sp.Expr]:
        """
        Return the conjugate-free numerator and denominator witnessing AB ⟂ CD.
        """
        raw_poly = self.perpendicular_poly(A, B, C, D)
        substituted = self._apply_all(raw_poly)
        together_expr = sp.together(substituted)
        numerator, denominator = sp.fraction(together_expr)
        numerator = sp.expand(numerator)
        denominator = sp.expand(denominator)

        zb_symbols = [
            symbol
            for symbol in numerator.free_symbols | denominator.free_symbols
            if isinstance(symbol, sp.Symbol) and str(symbol).startswith("zb_")
        ]
        if zb_symbols:
            raise GeometryError(f"Perpendicular check still depends on conjugates: {zb_symbols}")

        return numerator, denominator

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
