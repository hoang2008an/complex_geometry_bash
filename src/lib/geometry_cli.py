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
    {"op": "constraint_check", "constraint": "perpendicular", "args": ["A", "Q", "P", "H"]}
  ]
}
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import sympy as sp

from geometry_engine import GeometryEngine, GeometryError


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


def _print_constraint_results(results: Sequence[Tuple[sp.Expr, sp.Expr]]) -> None:
    multi = len(results) > 1
    for idx, (numerator, denominator) in enumerate(results, start=1):
        prefix = f"[{idx}] " if multi else ""
        print(f"{prefix}N = {sp.expand(numerator)}")
        print(f"{prefix}D = {sp.expand(denominator)}")


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
            elif op == "add_midpoint":
                engine.add_midpoint(*step["args"])
            elif op == "add_point_reflection":
                engine.add_point_reflection(*step["args"])
            elif op == "add_line_reflection":
                engine.add_line_reflection(*step["args"])
            elif op in {"add_triangle_similarity", "triangle_similarity"}:
                args = step.get("args", [])
                if len(args) != 6:
                    raise GeometryError("add_triangle_similarity expects six point labels (A, B, C, D, E, F).")
                directed_flag = step.get("directed")
                engine.add_triangle_similarity(*args, directed=True if directed_flag is None else bool(directed_flag))
            elif op in {
                "add_triangle_similarity_undirected",
                "triangle_similarity_undirected",
                "triangle_similarity_reflected",
            }:
                args = step.get("args", [])
                if len(args) != 6:
                    raise GeometryError("triangle_similarity_undirected expects six point labels (A, B, C, D, E, F).")
                engine.add_triangle_similarity(*args, directed=False)
            elif op in {"add_triangle_congruence", "triangle_congruence", "triangle_equal"}:
                args = step.get("args", [])
                if len(args) != 6:
                    raise GeometryError("add_triangle_congruence expects six point labels (A, B, C, D, E, F).")
                directed_flag = step.get("directed")
                engine.add_triangle_congruence(*args, directed=True if directed_flag is None else bool(directed_flag))
            elif op in {
                "add_triangle_congruence_undirected",
                "triangle_congruence_undirected",
                "triangle_equal_undirected",
            }:
                args = step.get("args", [])
                if len(args) != 6:
                    raise GeometryError("triangle_congruence_undirected expects six point labels (A, B, C, D, E, F).")
                engine.add_triangle_congruence(*args, directed=False)
            elif op == "orthocenter":
                engine.orthocenter_via_altitudes(*step["args"])
            elif op == "intersection":
                engine.intersection_of_lines(*step["args"])
            elif op == "circumcenter":
                args = step.get("args", [])
                if len(args) != 4:
                    raise GeometryError("circumcenter expects arguments [A, B, C, U].")
                engine.circumcenter(*args)
            elif op in {"set_unit_triangle", "set_main_unit_triangle"}:
                raw_args = step.get("args")
                triangle = step.get("triangle")
                roots = step.get("roots")

                if raw_args and triangle:
                    raise GeometryError("Provide triangle vertices via either 'args' or 'triangle', not both.")

                if raw_args is not None:
                    if not isinstance(raw_args, (list, tuple)):
                        raise GeometryError("set_unit_triangle 'args' must be a list.")
                    if len(raw_args) == 3:
                        triangle = list(raw_args)
                    elif len(raw_args) == 6:
                        triangle = list(raw_args[:3])
                        roots = list(raw_args[3:])
                    else:
                        raise GeometryError("set_unit_triangle expects 3 or 6 entries in 'args'.")

                if triangle is None:
                    raise GeometryError("set_unit_triangle requires triangle vertices.")
                if not isinstance(triangle, (list, tuple)) or len(triangle) != 3:
                    raise GeometryError("Triangle vertices must be a list of three labels.")

                root_names = None
                if roots is not None:
                    if not isinstance(roots, (list, tuple)) or len(roots) != 3:
                        raise GeometryError("roots must provide exactly three labels.")
                    root_names = list(roots)

                engine.set_main_unit_triangle(*triangle, root_names=root_names)
            elif op in {"unit_triangle_incenter", "main_triangle_incenter"}:
                name = step.get("name")
                if not isinstance(name, str):
                    raise GeometryError("unit_triangle_incenter requires a 'name' field.")
                engine.main_triangle_incenter(name)
            elif op in {"unit_triangle_excenter", "main_triangle_excenter"}:
                name = step.get("name")
                which = step.get("which")
                if not isinstance(name, str) or not isinstance(which, str):
                    raise GeometryError("unit_triangle_excenter requires 'name' and 'which' fields.")
                engine.main_triangle_excenter(which, name)
            elif op in {"unit_triangle_arc_midpoint", "main_triangle_arc_midpoint"}:
                name = step.get("name")
                which = step.get("which")
                if not isinstance(name, str) or not isinstance(which, str):
                    raise GeometryError("unit_triangle_arc_midpoint requires 'name' and 'which' fields.")
                containing = bool(step.get("containing_vertex") or step.get("containing"))
                engine.main_triangle_arc_midpoint(which, name, containing_vertex=containing)
            elif op == "line_circle_intersection":
                args = step.get("args", [])
                if len(args) != 5:
                    raise GeometryError(
                        "line_circle_intersection expects args [line_point1, line_point2, center, radius_point, name]."
                    )
                avoid = step.get("avoid")
                if avoid is None:
                    engine.line_circle_intersection(*args)
                else:
                    if isinstance(avoid, str):
                        avoid_list = [avoid]
                    elif isinstance(avoid, (list, tuple)):
                        avoid_list = list(avoid)
                    else:
                        raise GeometryError("avoid must be a string or list of point labels.")
                    engine.line_circle_intersection(*args, avoid=avoid_list)
            elif op == "centroid":
                engine.centroid(*step["args"])
            elif op == "euler_center":
                args = step.get("args", [])
                if len(args) == 4:
                    engine.euler_center(*args)
                elif len(args) == 6:
                    engine.euler_center(*args)
                else:
                    raise GeometryError("Euler center expects 4 or 6 arguments.")
            elif op in {"lemoine_point", "lemoine", "symmedian"}:
                args = step.get("args", [])
                if len(args) == 4:
                    engine.lemoine_point(*args)
                elif len(args) == 5:
                    engine.lemoine_point(*args)
                else:
                    raise GeometryError("Lemoine point expects 4 or 5 arguments.")
            elif op in {"reflect_point_over_point", "point_reflection"}:
                engine.reflect_point_over_point(*step["args"])
            elif op in {"reflect_point_over_line", "line_reflection"}:
                engine.reflect_point_over_line(*step["args"])
            elif op == "midpoint":
                engine.midpoint(*step["args"])
            elif op == "squared_distance":
                args = step.get("args", [])
                if len(args) != 2:
                    raise GeometryError("squared_distance expects exactly two point labels.")
                label = step.get("label") or f"|{args[0]}{args[1]}|^2"
                value = engine.squared_distance(*args)
                print(f"{label}: {sp.simplify(value)}")
            elif op == "print_points":
                names = step.get("names")
                summary = engine.point_summary(names)
                _print_point_summary(summary)
            elif op == "constraint_check":
                constraint = step["constraint"]
                args = step.get("args", [])
                if not isinstance(args, (list, tuple)):
                    raise GeometryError("Constraint check requires an 'args' list.")
                angle_expr = sp.sympify(step["angle"]) if "angle" in step else None
                directed_flag = step.get("directed")
                results = engine.constraint_conjugate_free(
                    constraint,
                    list(args),
                    angle=angle_expr,
                    directed=directed_flag,
                )
                _print_constraint_results(results)
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
    """Demonstration showing the Euler circle passing through the midpoint of XP."""
    return [
        {"op": "add_point", "name": "A", "unit_circle": True},
        {"op": "add_point", "name": "B", "unit_circle": True},
        {"op": "add_point", "name": "C", "unit_circle": True},
        {"op": "add_circumcenter", "args": ["A", "B", "C", "O"]},
        {"op": "orthocenter", "args": ["A", "B", "C", "H"]},
        {"op": "add_point", "name": "X"},
        {"op": "add_collinear", "args": ["B", "C", "X"]},
        {"op": "add_point", "name": "X_perp"},
        {"op": "add_perpendicular", "args": ["X", "X_perp", "B", "C"]},
        {"op": "intersection", "args": ["X", "X_perp", "A", "C", "Y"]},
        {"op": "midpoint", "args": ["B", "C", "M_BC"]},
        {"op": "midpoint", "args": ["A", "C", "M_AC"]},
        {"op": "reflect_point_over_point", "args": ["X", "M_BC", "Z"]},
        {"op": "reflect_point_over_point", "args": ["Y", "M_AC", "T"]},
        {"op": "line_circle_intersection", "args": ["Z", "T", "O", "X", "P"], "avoid": ["Z"]},
        {"op": "add_circumcenter", "args": ["X", "Z", "P", "O"]},
        {"op": "midpoint", "args": ["A", "B", "M_AB"]},
        {"op": "midpoint", "args": ["X", "P", "M_XP"]},
        {"op": "midpoint", "args": ["O", "H", "N_euler"]},
        {"op": "squared_distance", "args": ["M_XP", "N_euler"], "label": "|M_XP - (A+B+C)/2|^2"},
        {"op": "print_points", "names": ["A", "B", "C", "H", "X", "Y", "Z", "T", "P", "M_XP"]},
        {"op": "constraint_check", "constraint": "concyclic", "args": ["M_AB", "M_BC", "M_AC", "M_XP"]},
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
