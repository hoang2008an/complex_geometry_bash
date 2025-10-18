"""
Command-line interface for the complex geometry engine.

The CLI operates on JSON scripts describing a sequence of steps.  Each step is
an object with an `op` field and optional arguments required by the operation.

Example script (abbreviated):

{
  "steps": [
    {"op": "add_point", "name": "A", "unit_circle": true},
    {"op": "add_collinear", "args": ["O", "A", "P"]},
    {"op": "orthocenter", "args": ["A", "B", "C", "H"]},
    {"op": "print_points", "names": ["A", "B", "C", "H"]},
    {"op": "perp_check", "args": ["A", "Q", "P", "H"]}
  ]
}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

import sympy as sp

from .geometry_engine import GeometryEngine, GeometryError


def _print_point_summary(summary: Dict[str, sp.Expr]) -> None:
    for name in sorted(summary.keys()):
        print(f"{name}: {sp.simplify(summary[name])}")


def _print_learned_rules(rules: Dict[str, sp.Expr]) -> None:
    if not rules:
        print("No learned conjugate rules yet.")
        return
    for symbol, expr in sorted(rules.items()):
        print(f"{symbol} -> {sp.simplify(expr)}")


def _print_constraints(constraints: Iterable[str]) -> None:
    for idx, constraint in enumerate(constraints, start=1):
        print(f"[{idx}] {constraint} = 0")


def run_steps(engine: GeometryEngine, steps: Sequence[Dict[str, Any]]) -> None:
    for index, step in enumerate(steps, start=1):
        op = step.get("op")
        if not op:
            raise GeometryError(f"Step {index} missing 'op' field.")

        try:
            if op == "add_point":
                name = step["name"]
                engine.add_point(name)
                if step.get("unit_circle"):
                    engine.add_unit_circle(name)
            elif op == "add_unit_circle":
                name = step["name"]
                engine.add_unit_circle(name)
            elif op == "add_collinear":
                engine.add_collinear(*step["args"])
            elif op == "add_perpendicular":
                engine.add_perpendicular(*step["args"])
            elif op == "add_concyclic":
                engine.add_concyclic(*step["args"])
            elif op == "add_angle_value":
                A, B, C = step["args"]
                angle_expr = sp.sympify(step["angle"])
                engine.add_angle_value(A, B, C, angle_expr)
            elif op == "add_circumcenter":
                engine.add_circumcenter(*step["args"])
            elif op == "orthocenter":
                engine.orthocenter_via_altitudes(*step["args"])
            elif op == "intersection":
                engine.intersection_of_lines(*step["args"])
            elif op == "centroid":
                engine.centroid(*step["args"])
            elif op == "print_points":
                names = step.get("names")
                summary = engine.point_summary(names)
                _print_point_summary(summary)
            elif op == "perp_check":
                numerator, denominator = engine.perpendicular_conjugate_free(*step["args"])
                print(f"N = {sp.expand(numerator)}")
                print(f"D = {sp.expand(denominator)}")
            elif op == "learned_rules":
                _print_learned_rules(engine.learned_rules())
            elif op == "print_constraints":
                _print_constraints(engine.constraint_strings())
            else:
                raise GeometryError(f"Unsupported op '{op}' on step {index}.")
        except KeyError as exc:
            raise GeometryError(f"Step {index} missing required key: {exc}") from exc


def load_script(path: Path) -> List[Dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    steps = payload.get("steps")
    if not isinstance(steps, list):
        raise GeometryError("Script must contain a 'steps' list.")
    return steps


def _demo_steps() -> List[Dict[str, Any]]:
    """Workflow mirroring the spec's sample scenario."""
    return [
        {"op": "add_point", "name": "A", "unit_circle": True},
        {"op": "add_point", "name": "B", "unit_circle": True},
        {"op": "add_point", "name": "C", "unit_circle": True},
        {"op": "add_circumcenter", "args": ["A", "B", "C", "O"]},
        {"op": "orthocenter", "args": ["A", "B", "C", "H"]},
        {"op": "add_point", "name": "P"},
        {"op": "add_collinear", "args": ["O", "A", "P"]},
        {"op": "orthocenter", "args": ["P", "A", "B", "K"]},
        {"op": "orthocenter", "args": ["P", "C", "A", "L"]},
        {"op": "intersection", "args": ["B", "L", "C", "K", "Q"]},
        {"op": "centroid", "args": ["A", "B", "C", "G"]},
        {"op": "add_angle_value", "args": ["A", "B", "P"], "angle": "pi/2"},
        {"op": "print_points", "names": ["A", "B", "C", "H", "P", "K", "L", "Q", "G"]},
        {"op": "perp_check", "args": ["A", "Q", "P", "H"]},
        {"op": "learned_rules"},
        {"op": "print_constraints"},
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Complex geometry engine CLI.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Execute operations from a JSON script.")
    run_parser.add_argument("script", type=Path, help="Path to JSON script.")

    demo_parser = subparsers.add_parser("demo", help="Run the bundled demonstration workflow.")
    demo_parser.add_argument("--no-learned", action="store_true", help="Skip printing learned rules.")

    poly_parser = subparsers.add_parser("poly", help="Inspect predicate polynomials.")
    poly_subparsers = poly_parser.add_subparsers(dest="poly_command", required=True)

    concyc_parser = poly_subparsers.add_parser("concyclic", help="Polynomial enforcing four concyclic points.")
    concyc_parser.add_argument("A")
    concyc_parser.add_argument("B")
    concyc_parser.add_argument("C")
    concyc_parser.add_argument("D")
    concyc_parser.add_argument("--raw", action="store_true", help="Return the unexpanded polynomial.")

    circum_parser = poly_subparsers.add_parser("circumcenter", help="Polynomials enforcing U as circumcenter of triangle ABC.")
    circum_parser.add_argument("A")
    circum_parser.add_argument("B")
    circum_parser.add_argument("C")
    circum_parser.add_argument("U")
    circum_parser.add_argument("--raw", action="store_true", help="Return unexpanded polynomials.")

    angle_parser = poly_subparsers.add_parser("angle", help="Polynomial enforcing ∠ABC equals a specified angle (radians).")
    angle_parser.add_argument("A")
    angle_parser.add_argument("B")
    angle_parser.add_argument("C")
    angle_parser.add_argument("angle")
    angle_parser.add_argument("--raw", action="store_true", help="Return the unexpanded polynomial.")

    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    engine = GeometryEngine()

    if args.command == "run":
        steps = load_script(args.script)
        run_steps(engine, steps)
    elif args.command == "demo":
        steps = _demo_steps()
        if args.no_learned:
            steps = [step for step in steps if step.get("op") != "learned_rules"]
        run_steps(engine, steps)
    elif args.command == "poly":
        if args.poly_command == "concyclic":
            for name in {args.A, args.B, args.C, args.D}:
                engine.add_point(name)
            poly = engine.concyclic_poly(args.A, args.B, args.C, args.D, raw=args.raw)
            print(sp.expand(poly) if not args.raw else poly)
        elif args.poly_command == "circumcenter":
            for name in {args.A, args.B, args.C, args.U}:
                engine.add_point(name)
            eq1, eq2 = engine.circumcenter_polys(args.A, args.B, args.C, args.U, raw=args.raw)
            if args.raw:
                print(f"eq1: {eq1}")
                print(f"eq2: {eq2}")
            else:
                print(f"eq1: {sp.expand(eq1)}")
                print(f"eq2: {sp.expand(eq2)}")
        elif args.poly_command == "angle":
            for name in {args.A, args.B, args.C}:
                engine.add_point(name)
            angle_expr = sp.sympify(args.angle)
            poly = engine.angle_value_poly(args.A, args.B, args.C, angle_expr, raw=args.raw)
            print(poly if args.raw else sp.expand(poly))
        else:
            parser.error("Unsupported poly command.")
    else:
        parser.error("Unknown command.")


if __name__ == "__main__":
    main()
