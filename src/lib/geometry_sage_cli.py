"""
Command-line interface for the isolated Sage geometry backend.

Run through ``bin/geometry-sage`` so Sage imports are handled by ``sage -python``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

from geometry_sage_engine import SageGeometryEngine, SageGeometryError
from geometry_sage_sr_engine import SageSRGeometryEngine, SageSRGeometryError


UNSUPPORTED_BRANCH_OPS = {
    "line_circle_object_intersection",
    "line_circle_object_intersections",
    "circle_object_intersection",
    "circle_object_intersections",
    "tangent_lines_point_circle",
    "point_circle_tangents",
    "radical_line",
    "circle_radical_line",
}


def _print_point_summary(summary: Dict[str, str]) -> None:
    for name in sorted(summary.keys()):
        print(f"{name}: {summary[name]}")


def _print_learned_rules(rules: Dict[str, Any]) -> None:
    if not rules:
        print("No learned conjugate rules yet.")
        return
    for symbol, expr in sorted(rules.items()):
        print(f"{symbol} -> {expr}")


def _print_constraints(constraints: Iterable[str]) -> None:
    for idx, constraint in enumerate(constraints, start=1):
        print(f"[{idx}] {constraint} = 0")


def _print_constraint_results(results: Sequence[Tuple[Any, Any]]) -> None:
    multi = len(results) > 1
    for idx, (numerator, denominator) in enumerate(results, start=1):
        prefix = f"[{idx}] " if multi else ""
        print(f"{prefix}N = {numerator}")
        print(f"{prefix}D = {denominator}")


def _collect_labels_from_value(value: Any, labels: set[str]) -> None:
    if isinstance(value, str):
        labels.add(value)
    elif isinstance(value, list):
        for item in value:
            _collect_labels_from_value(item, labels)


def collect_point_labels(steps: Sequence[Dict[str, Any]]) -> List[str]:
    labels = {"O"}
    for step in steps:
        op = step.get("op")
        for key in ("name", "point"):
            value = step.get(key)
            if isinstance(value, str):
                labels.add(value)
        args = step.get("args")
        if isinstance(args, list):
            for value in args:
                if isinstance(value, str):
                    labels.add(value)
        names = step.get("names")
        if isinstance(names, list):
            for value in names:
                if isinstance(value, str):
                    labels.add(value)
        if op in {"set_point_value", "assign_value", "fix_point"}:
            name = step.get("name")
            if isinstance(name, str):
                labels.add(name)
    return sorted(labels)


GeometryErrorTypes = (SageGeometryError, SageSRGeometryError)


def _parse_value(engine: Any, raw: Any) -> Any:
    if isinstance(raw, (int, float)):
        return engine.expr(raw)
    if isinstance(raw, str):
        return engine.expr(raw)
    raise SageGeometryError(f"Unsupported Sage value {raw!r}. Use an exact number string.")


def run_steps(engine: Any, steps: Sequence[Dict[str, Any]]) -> None:
    for index, step in enumerate(steps, start=1):
        op = step.get("op")
        if not op:
            raise SageGeometryError(f"Step {index} missing 'op' field.")

        try:
            if op in UNSUPPORTED_BRANCH_OPS:
                raise SageGeometryError(
                    f"Unsupported op '{op}' in Sage mode v1: circle/radical intersections need square-root branch handling."
                )
            if op == "add_point":
                engine.add_point(step["name"])
                if step.get("unit_circle"):
                    engine.add_unit_circle(step["name"])
            elif op == "add_unit_circle":
                engine.add_unit_circle(step["name"])
            elif op in {"set_point_value", "assign_value", "fix_point"}:
                name = step["name"]
                z_value = _parse_value(engine, step["z"])
                zb_value = _parse_value(engine, step["zb"]) if "zb" in step else None
                engine.set_point_value(name, z_value, zb_value)
            elif op == "add_collinear":
                engine.add_collinear(*step["args"])
            elif op == "add_perpendicular":
                engine.add_perpendicular(*step["args"])
            elif op == "add_angle_value":
                A, B, C = step["args"]
                engine.add_angle_value(A, B, C, step["angle"])
            elif op in {"add_angle_bisector_either", "angle_bisector_either", "angle_bisector_any"}:
                engine.add_angle_bisector_either(*step["args"])
            elif op in {"set_unit_triangle", "set_main_unit_triangle"}:
                raw_args = step.get("args")
                triangle = step.get("triangle")
                roots = step.get("roots")
                if raw_args is not None:
                    if len(raw_args) == 3:
                        triangle = raw_args
                    elif len(raw_args) == 6:
                        triangle = raw_args[:3]
                        roots = raw_args[3:]
                    else:
                        raise SageGeometryError("set_unit_triangle expects 3 or 6 entries in args.")
                if not isinstance(triangle, (list, tuple)) or len(triangle) != 3:
                    raise SageGeometryError("set_unit_triangle requires three triangle vertices.")
                engine.set_main_unit_triangle(*triangle, root_names=roots)
            elif op in {"unit_triangle_incenter", "main_triangle_incenter"}:
                name = step.get("name")
                if not isinstance(name, str):
                    raise SageGeometryError("main_triangle_incenter requires a name.")
                engine.main_triangle_incenter(name)
            elif op in {"unit_triangle_arc_midpoint", "main_triangle_arc_midpoint"}:
                name = step.get("name")
                which = step.get("which")
                if not isinstance(name, str) or not isinstance(which, str):
                    raise SageGeometryError("unit_triangle_arc_midpoint requires 'name' and 'which' fields.")
                containing = bool(step.get("containing_vertex") or step.get("containing"))
                engine.main_triangle_arc_midpoint(which, name, containing_vertex=containing)
            elif op in {"line_from_points", "line_through"}:
                name = step.get("name")
                args = step.get("args", [])
                if not isinstance(name, str) or len(args) != 2:
                    raise SageGeometryError("line_from_points expects 'name' and args [P, Q].")
                engine.register_line(name, engine.line_through_points(*args))
            elif op in {"line_init", "line_new"}:
                name = step.get("name")
                if not isinstance(name, str):
                    raise SageGeometryError("line_init expects a 'name'.")
                engine.register_line(name, engine.symbolic_line(name))
            elif op == "point_on_line":
                line_name = step.get("line")
                point = step.get("point")
                if not isinstance(line_name, str) or not isinstance(point, str):
                    raise SageGeometryError("point_on_line expects 'line' and 'point'.")
                engine.add_point_on_line(engine.lines[line_name], point)
            elif op in {"perpendicular_lines", "line_perpendicular"}:
                args = step.get("args", [])
                if len(args) != 2:
                    raise SageGeometryError("perpendicular_lines expects args [line1_name, line2_name].")
                line1_name, line2_name = args
                engine.add_perpendicular_lines(engine.lines[line1_name], engine.lines[line2_name])
            elif op in {"line_intersection", "intersection_lines"}:
                args = step.get("args", [])
                if len(args) != 3:
                    raise SageGeometryError("line_intersection expects args [line1_name, line2_name, point_label].")
                line1_name, line2_name, point_label = args
                engine.line_intersection(engine.lines[line1_name], engine.lines[line2_name], point_label)
            elif op == "add_circumcenter":
                engine.add_circumcenter(*step["args"])
            elif op == "add_midpoint":
                engine.add_midpoint(*step["args"])
            elif op in {"add_centroid_constraint", "add_centroid", "centroid_constraint"}:
                engine.add_centroid_constraint(*step["args"])
            elif op in {"add_projection_to_line", "projection_constraint"}:
                engine.add_projection_to_line(*step["args"])
            elif op == "midpoint":
                engine.midpoint(*step["args"])
            elif op == "centroid":
                engine.centroid(*step["args"])
            elif op == "circumcenter":
                engine.circumcenter(*step["args"])
            elif op == "orthocenter":
                engine.orthocenter_via_altitudes(*step["args"])
            elif op in {"intersection", "line_intersection"}:
                engine.intersection_of_lines(*step["args"])
            elif op in {"project_point_to_line", "point_projection"}:
                engine.project_point_to_line(*step["args"])
            elif op in {"reflect_point_over_line", "line_reflection"}:
                engine.reflect_point_over_line(*step["args"])
            elif op == "line_circle_intersection":
                args = step.get("args", [])
                if len(args) != 5:
                    raise SageGeometryError(
                        "line_circle_intersection expects args [line_point1, line_point2, center, radius_point, name]."
                    )
                avoid = step.get("avoid")
                if isinstance(avoid, str):
                    avoid = [avoid]
                engine.line_circle_intersection(*args, avoid=avoid)
            elif op == "circle_intersection":
                args = step.get("args", [])
                if len(args) != 5:
                    raise SageGeometryError(
                        "circle_intersection expects args [center1, radius_point1, center2, radius_point2, name]."
                    )
                avoid = step.get("avoid")
                if isinstance(avoid, str):
                    avoid = [avoid]
                engine.circle_intersection(*args, avoid=avoid)
            elif op in {"lemoine_point", "lemoine", "symmedian"}:
                engine.lemoine_point(*step["args"])
            elif op in {"add_isogonal_conjugate"}:
                engine.add_isogonal_conjugate(*step["args"])
            elif op in {"isogonal_conjugate", "isogonal_conj"}:
                engine.isogonal_conjugate_point(*step["args"])
            elif op in {"add_fermat_points", "fermat_points"}:
                engine.add_fermat_points(*step["args"])
            elif op == "squared_distance":
                args = step.get("args", [])
                if len(args) != 2:
                    raise SageGeometryError("squared_distance expects exactly two point labels.")
                label = step.get("label") or f"|{args[0]}{args[1]}|^2"
                print(f"{label}: {engine.squared_distance(*args)}")
            elif op == "constraint_check":
                angle = step.get("angle")
                results = engine.constraint_conjugate_free(step["constraint"], step.get("args", []), angle=angle)
                _print_constraint_results(results)
            elif op == "print_points":
                _print_point_summary(engine.point_summary(step.get("names")))
            elif op == "learned_rules":
                _print_learned_rules(engine.learned_rules())
            elif op == "print_constraints":
                _print_constraints(engine.constraint_strings())
            else:
                raise SageGeometryError(f"Unsupported op '{op}' on step {index} in Sage mode v1.")
        except KeyError as exc:
            raise SageGeometryError(f"Step {index} missing required key: {exc}") from exc


def load_script(path: Path) -> List[Dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    steps = payload.get("steps")
    if not isinstance(steps, list):
        raise SageGeometryError("Script must contain a 'steps' list.")
    return steps


def _demo_steps() -> List[Dict[str, Any]]:
    return [
        {"op": "add_point", "name": "A", "unit_circle": True},
        {"op": "add_point", "name": "B", "unit_circle": True},
        {"op": "midpoint", "args": ["A", "B", "M"]},
        {"op": "constraint_check", "constraint": "collinear", "args": ["A", "B", "M"]},
        {"op": "print_points", "names": ["A", "B", "M"]},
        {"op": "learned_rules"},
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Experimental isolated Sage geometry CLI.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Execute Sage-mode operations from a JSON script.")
    run_parser.add_argument("script", type=Path, help="Path to JSON script.")
    run_parser.add_argument("--backend", choices=["fraction", "sr"], default="fraction", help="Sage backend to use.")
    run_parser.add_argument("--field", choices=["QQbar", "QQ"], default="QQbar", help="Coefficient field for the Sage fraction field.")

    demo_parser = subparsers.add_parser("demo", help="Run a small Sage-mode rational workflow.")
    demo_parser.add_argument("--backend", choices=["fraction", "sr"], default="fraction", help="Sage backend to use.")
    demo_parser.add_argument("--field", choices=["QQbar", "QQ"], default="QQbar", help="Coefficient field for the Sage fraction field.")
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run":
        steps = load_script(args.script)
    else:
        steps = _demo_steps()

    if args.backend == "sr":
        engine = SageSRGeometryEngine()
        for name in collect_point_labels(steps):
            engine.add_point(name)
    else:
        engine = SageGeometryEngine(point_names=collect_point_labels(steps), coefficient_field=args.field)
    run_steps(engine, steps)


if __name__ == "__main__":
    try:
        main()
    except GeometryErrorTypes as exc:
        raise SystemExit(str(exc))
