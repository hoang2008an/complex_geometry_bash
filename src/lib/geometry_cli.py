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

DISPLAY_STYLE = "latex"


def _print_point_summary(summary: Dict[str, str]) -> None:
    for name in sorted(summary.keys()):
        value = summary[name]
        print(f"{name}: {value}")


def _print_learned_rules(rules: Dict[str, str]) -> None:
    if not rules:
        print("No learned conjugate rules yet.")
        return
    for symbol, expr in sorted(rules.items()):
        print(f"{symbol} -> {expr}")


def _print_constraints(constraints: Iterable[str]) -> None:
    for idx, constraint in enumerate(constraints, start=1):
        print(f"[{idx}] {constraint} = 0")


def _print_constraint_results(results: Sequence[Tuple[str, str]]) -> None:
    multi = len(results) > 1
    for idx, (numerator, denominator) in enumerate(results, start=1):
        prefix = f"[{idx}] " if multi else ""
        print(f"{prefix}N = {numerator}")
        print(f"{prefix}D = {denominator}")


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
            elif op in {"add_centroid_constraint", "add_centroid", "centroid_constraint"}:
                args = step.get("args", [])
                if len(args) != 4:
                    raise GeometryError("add_centroid_constraint expects args [A, B, C, G].")
                engine.add_centroid_constraint(*args)
            elif op == "add_point_reflection":
                engine.add_point_reflection(*step["args"])
            elif op == "add_line_reflection":
                engine.add_line_reflection(*step["args"])
            elif op in {"add_projection_to_line", "projection_constraint"}:
                args = step.get("args", [])
                if len(args) != 4:
                    raise GeometryError("add_projection_to_line expects args [P, A, B, H].")
                engine.add_projection_to_line(*args)
            elif op in {"add_isogonal_reflection", "isogonal_reflection", "angle_bisector_reflection"}:
                args = step.get("args", [])
                if len(args) != 5:
                    raise GeometryError("isogonal_reflection expects args [A, B, C, D, E].")
                engine.add_isogonal_reflection(*args)
            elif op == "add_isogonal_conjugate":
                args = step.get("args", [])
                if len(args) != 5:
                    raise GeometryError("add_isogonal_conjugate expects args [A, B, C, P, Q].")
                engine.add_isogonal_conjugate(*args)
            elif op in {"isogonal_conjugate", "isogonal_conj"}:
                args = step.get("args", [])
                if len(args) != 5:
                    raise GeometryError("isogonal_conjugate expects args [A, B, C, P, Q].")
                engine.isogonal_conjugate_point(*args)
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
                "add_angle_bisector_either",
                "angle_bisector_either",
                "angle_bisector_any",
            }:
                args = step.get("args", [])
                if len(args) != 4:
                    raise GeometryError("add_angle_bisector_either expects args [A, B, C, D].")
                engine.add_angle_bisector_either(*args)
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
            elif op in {"project_point_to_line", "point_projection"}:
                args = step.get("args", [])
                if len(args) != 4:
                    raise GeometryError("project_point_to_line expects args [P, A, B, H].")
                engine.project_point_to_line(*args)
            elif op == "circumcenter":
                args = step.get("args", [])
                if len(args) != 4:
                    raise GeometryError("circumcenter expects arguments [A, B, C, U].")
                engine.circumcenter(*args)
            elif op in {"add_fermat_points", "fermat_points"}:
                args = step.get("args", [])
                if len(args) != 5:
                    raise GeometryError("add_fermat_points expects args [A, B, C, F1, F2].")
                engine.add_fermat_points(*args)
            elif op in {"humpty_point", "humpy_point", "humpty"}:
                args = step.get("args", [])
                if len(args) == 4:
                    engine.humpty_point(*args)
                elif len(args) == 5:
                    A, B, C, which, H = args
                    engine.humpty_point(A, B, C, H, which=which)
                else:
                    raise GeometryError("humpty_point expects args [A, B, C, H] or [A, B, C, which, H].")
            elif op in {"dumpty_point", "dumpty"}:
                args = step.get("args", [])
                if len(args) == 4:
                    engine.dumpty_point(*args)
                elif len(args) == 5:
                    A, B, C, which, D = args
                    engine.dumpty_point(A, B, C, D, which=which)
                else:
                    raise GeometryError("dumpty_point expects args [A, B, C, D] or [A, B, C, which, D].")
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
            elif op in {"line_from_points", "line_through"}:
                name = step.get("name")
                args = step.get("args", [])
                if not isinstance(name, str) or len(args) != 2:
                    raise GeometryError("line_from_points expects 'name' and args [P, Q].")
                if name in engine.lines:
                    raise GeometryError(f"Line '{name}' is already defined.")
                engine.register_line(name, engine.line_through_points(*args))
            elif op in {"line_init", "line_new"}:
                name = step.get("name")
                if not isinstance(name, str):
                    raise GeometryError("line_init expects a 'name'.")
                if name in engine.lines:
                    raise GeometryError(f"Line '{name}' is already defined.")
                engine.register_line(name, engine.symbolic_line(name))
            elif op in {"line_from_coefficients", "line_define"}:
                name = step.get("name")
                alpha = step.get("alpha")
                beta = step.get("beta")
                gamma = step.get("gamma")
                if not isinstance(name, str):
                    raise GeometryError("line_from_coefficients requires a 'name'.")
                if alpha is None or beta is None or gamma is None:
                    raise GeometryError("line_from_coefficients requires 'alpha', 'beta', and 'gamma'.")
                if name in engine.lines:
                    raise GeometryError(f"Line '{name}' is already defined.")
                line = engine.line_from_coefficients(sp.sympify(alpha), sp.sympify(beta), sp.sympify(gamma))
                engine.register_line(name, line)
            elif op == "point_on_line":
                line_name = step.get("line")
                point = step.get("point")
                if not isinstance(line_name, str) or not isinstance(point, str):
                    raise GeometryError("point_on_line expects 'line' and 'point'.")
                if line_name not in engine.lines:
                    raise GeometryError(f"Line '{line_name}' is not defined.")
                engine.add_point_on_line(engine.lines[line_name], point)
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
            elif op in {"circle_from_three_points", "circle_three_points"}:
                name = step.get("name")
                args = step.get("args", [])
                if not isinstance(name, str):
                    raise GeometryError("circle_from_three_points requires a 'name'.")
                if len(args) != 3:
                    raise GeometryError("circle_from_three_points expects args [A, B, C].")
                engine.circle_from_three_points(name, *args)
            elif op in {
                "line_circle_object_intersection",
                "line_circle_object_intersections",
                "line_circle_object",
            }:
                line_name = step.get("line")
                circle_name = step.get("circle")
                names = step.get("names")
                if not isinstance(line_name, str):
                    raise GeometryError("line_circle_object_intersection requires 'line'.")
                if not isinstance(circle_name, str):
                    raise GeometryError("line_circle_object_intersection requires 'circle'.")
                if not isinstance(names, (list, tuple)) or not names:
                    raise GeometryError("line_circle_object_intersection requires a non-empty 'names' list.")
                point_names = []
                for entry in names:
                    if not isinstance(entry, str):
                        raise GeometryError("Intersection point labels must be strings.")
                    point_names.append(entry)
                avoid = step.get("avoid")
                if avoid is None:
                    engine.line_circle_object_intersections(line_name, circle_name, point_names)
                else:
                    if isinstance(avoid, str):
                        avoid_list = [avoid]
                    elif isinstance(avoid, (list, tuple)):
                        avoid_list = list(avoid)
                    else:
                        raise GeometryError("avoid must be a string or list of point labels.")
                    engine.line_circle_object_intersections(line_name, circle_name, point_names, avoid=avoid_list)
            elif op in {
                "circle_object_intersection",
                "circle_object_intersections",
                "circle_circle_object",
            }:
                circle1_name = step.get("circle1")
                circle2_name = step.get("circle2")
                names = step.get("names")
                if not isinstance(circle1_name, str):
                    raise GeometryError("circle_object_intersection requires 'circle1'.")
                if not isinstance(circle2_name, str):
                    raise GeometryError("circle_object_intersection requires 'circle2'.")
                if not isinstance(names, (list, tuple)) or not names:
                    raise GeometryError("circle_object_intersection requires a non-empty 'names' list.")
                point_names = []
                for entry in names:
                    if not isinstance(entry, str):
                        raise GeometryError("Intersection point labels must be strings.")
                    point_names.append(entry)
                avoid = step.get("avoid")
                if avoid is None:
                    engine.circle_object_intersections(circle1_name, circle2_name, point_names)
                else:
                    if isinstance(avoid, str):
                        avoid_list = [avoid]
                    elif isinstance(avoid, (list, tuple)):
                        avoid_list = list(avoid)
                    else:
                        raise GeometryError("avoid must be a string or list of point labels.")
                    engine.circle_object_intersections(circle1_name, circle2_name, point_names, avoid=avoid_list)
            elif op in {"radical_line", "circle_radical_line"}:
                line_name = step.get("name")
                circle1_name = step.get("circle1")
                circle2_name = step.get("circle2")
                if not isinstance(line_name, str):
                    raise GeometryError("radical_line requires a 'name' for the resulting line.")
                if not isinstance(circle1_name, str) or not isinstance(circle2_name, str):
                    raise GeometryError("radical_line requires 'circle1' and 'circle2' labels.")
                engine.radical_line_from_circles(circle1_name, circle2_name, line_name)
            elif op in {
                "tangent_lines_point_circle",
                "point_circle_tangents",
                "tangent_line_point_circle",
            }:
                point = step.get("point")
                circle_name = step.get("circle")
                line_entries = step.get("lines")
                if not isinstance(point, str):
                    raise GeometryError("tangent_lines_point_circle requires a 'point' label.")
                if not isinstance(circle_name, str):
                    raise GeometryError("tangent_lines_point_circle requires a 'circle' name.")
                if not isinstance(line_entries, (list, tuple)) or not line_entries:
                    raise GeometryError("tangent_lines_point_circle expects a non-empty 'lines' list.")
                line_names = []
                for entry in line_entries:
                    if not isinstance(entry, str):
                        raise GeometryError("Line labels must be strings.")
                    line_names.append(entry)
                tangent_points = step.get("points") or step.get("tangent_points")
                if tangent_points is not None:
                    if not isinstance(tangent_points, (list, tuple)):
                        raise GeometryError("tangent point labels must be provided as a list.")
                    tangent_point_names = []
                    for label in tangent_points:
                        if not isinstance(label, str):
                            raise GeometryError("Tangent point labels must be strings.")
                        tangent_point_names.append(label)
                else:
                    tangent_point_names = None
                engine.tangent_lines_from_point_to_circle(
                    point,
                    circle_name,
                    line_names,
                    tangent_point_names=tangent_point_names,
                )
            elif op == "circle_intersection":
                args = step.get("args", [])
                if len(args) != 5:
                    raise GeometryError(
                        "circle_intersection expects args [center1, radius_point1, center2, radius_point2, name]."
                    )
                avoid = step.get("avoid")
                if avoid is None:
                    engine.circle_intersection(*args)
                else:
                    if isinstance(avoid, str):
                        avoid_list = [avoid]
                    elif isinstance(avoid, (list, tuple)):
                        avoid_list = list(avoid)
                    else:
                        raise GeometryError("avoid must be a string or list of point labels.")
                    engine.circle_intersection(*args, avoid=avoid_list)
            elif op in {"avoid", "add_avoid", "distinct_points"}:
                args = step.get("args", [])
                if len(args) != 2:
                    raise GeometryError("avoid expects args [P, Q].")
                engine.add_distinct_points(*args)
            elif op in {"perpendicular_lines", "line_perpendicular"}:
                args = step.get("args", [])
                if len(args) != 2:
                    raise GeometryError("perpendicular_lines expects args [line1_name, line2_name].")
                line1_name, line2_name = args
                if not isinstance(line1_name, str) or not isinstance(line2_name, str):
                    raise GeometryError("perpendicular_lines arguments must be line labels.")
                if line1_name not in engine.lines or line2_name not in engine.lines:
                    missing = [name for name in (line1_name, line2_name) if name not in engine.lines]
                    raise GeometryError(f"Lines not defined: {', '.join(missing)}.")
                engine.add_perpendicular_lines(engine.lines[line1_name], engine.lines[line2_name])
            elif op in {"line_intersection", "intersection_lines"}:
                args = step.get("args", [])
                if len(args) != 3:
                    raise GeometryError("line_intersection expects args [line1_name, line2_name, point_label].")
                line1_name, line2_name, point_label = args
                if line1_name not in engine.lines:
                    raise GeometryError(f"Line '{line1_name}' is not defined.")
                if line2_name not in engine.lines:
                    raise GeometryError(f"Line '{line2_name}' is not defined.")
                engine.line_intersection(engine.lines[line1_name], engine.lines[line2_name], point_label)
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
                print(f"{label}: {engine.format_expr(value, style=DISPLAY_STYLE)}")
            elif op == "print_points":
                names = step.get("names")
                summary = engine.point_summary(names, style=DISPLAY_STYLE)
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
                formatted_results = [
                    (
                        engine.format_expr(numerator, style=DISPLAY_STYLE),
                        engine.format_expr(denominator, style=DISPLAY_STYLE),
                    )
                    for numerator, denominator in results
                ]
                _print_constraint_results(formatted_results)
            elif op == "learned_rules":
                _print_learned_rules(engine.display_learned_rules(style=DISPLAY_STYLE))
            elif op == "print_constraints":
                _print_constraints(engine.constraint_strings(style=DISPLAY_STYLE))
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
            expr = poly if args.raw else sp.expand(poly)
            print(engine.format_expr(expr, style=DISPLAY_STYLE))
        elif args.poly_command == "circumcenter":
            for name in {args.A, args.B, args.C, args.U}:
                engine.add_point(name)
            eq1, eq2 = engine.circumcenter_polys(args.A, args.B, args.C, args.U, raw=args.raw)
            if args.raw:
                print(f"eq1: {engine.format_expr(eq1, style=DISPLAY_STYLE)}")
                print(f"eq2: {engine.format_expr(eq2, style=DISPLAY_STYLE)}")
            else:
                print(f"eq1: {engine.format_expr(sp.expand(eq1), style=DISPLAY_STYLE)}")
                print(f"eq2: {engine.format_expr(sp.expand(eq2), style=DISPLAY_STYLE)}")
        elif args.poly_command == "angle":
            for name in {args.A, args.B, args.C}:
                engine.add_point(name)
            angle_expr = sp.sympify(args.angle)
            poly = engine.angle_value_poly(args.A, args.B, args.C, angle_expr, raw=args.raw)
            expr = poly if args.raw else sp.expand(poly)
            print(engine.format_expr(expr, style=DISPLAY_STYLE))
        else:
            parser.error("Unsupported poly command.")
    else:
        parser.error("Unknown command.")


if __name__ == "__main__":
    main()
