# Complex Geometry Bash

Symbolic complex-geometry playground powered by SymPy. The engine lives under `src/lib/geometry_engine.py`, while a JSON-driven CLI orchestrator sits in `src/lib/geometry_cli.py`. Launchers and tests wire everything together so you can explore triangle centers, constraints, and algebraic identities.

## Quick Start

```bash
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
bin/geometry demo          # run the bundled workflow
bin/geometry run scripts/example.json
bin/geometry poly concyclic A B C D
```

The `bin/geometry` wrapper exports `PYTHONPATH=src/lib:src` and calls `python -m geometry_cli`, so all absolute imports resolve without extra setup. If you invoke the CLI directly, mimic that environment.

## Available Operations

The CLI consumes JSON scripts whose `steps` array contains operations with optional arguments. The table below lists the supported `op` keys and their expected arguments.

| Operation | Arguments | Description |
|-----------|-----------|-------------|
| `add_point` | `name`, optional `unit_circle` flag | Registers a point; optionally constrains it to the unit circle. |
| `add_unit_circle` | `name` | Declares that an existing point lies on the unit circle. |
| `add_collinear` | `[A, B, C]` | Adds the collinearity constraint for the three points. |
| `add_perpendicular` | `[A, B, C, D]` | Constrains line AB perpendicular to CD. |
| `add_concyclic` | `[A, B, C, D]` | Forces four points onto a common circle. |
| `add_angle_value` | `[A, B, C]`, `angle` | Fixes ∠ABC to the provided SymPy expression. |
| `add_circumcenter` | `[A, B, C, U]` | Enforces U as the circumcenter of triangle ABC. |
| `add_midpoint` | `[P, Q, M]` | Constrains M to be the midpoint of segment PQ. |
| `add_point_reflection` | `[P, O, Q]` | Adds Q as the reflection of P across point O. |
| `add_line_reflection` | `[P, A, B, Q]` | Adds Q as the reflection of P across line AB. |
| `orthocenter` | `[A, B, C, H]` | Constructs the altitudes-based orthocenter. |
| `intersection` | `[A, B, C, D, X]` | Intersects lines AB and CD. |
| `centroid` | `[A, B, C, G]` | Assigns the centroid of triangle ABC. |
| `euler_center` | `[A, B, C, N]` or `[A, B, C, N, U, H]` | Midpoint of circumcenter/orthocenter (optional helper labels). |
| `lemoine_point` | `[A, B, C, L]` or `[A, B, C, L, U]` | Symmedian point via tangents (optional circumcenter label). |
| `reflect_point_over_point` | `[P, O, Q]` | Constructs Q as reflection of P about point O. |
| `reflect_point_over_line` | `[P, A, B, Q]` | Constructs Q as reflection of P across line AB. |
| `midpoint` | `[P, Q, M]` | Constructs the midpoint M directly. |
| `print_points` | optional `names` list | Prints simplified z-coordinates for selected points. |
| `constraint_check` | `constraint`, `args`, optional `angle` | Evaluates any supported constraint (`collinear`, `perpendicular`, `concyclic`, `angle`, `circumcenter`, `midpoint`, `point_reflection`, `line_reflection`, …). |
| `learned_rules` | — | Prints auto-learned conjugate substitutions. |
| `print_constraints` | — | Dumps stored constraint polynomials. |

Shortcut aliases: `lemoine`, `symmedian`, `point_reflection`, and `line_reflection` map to their corresponding constructors.

## Constraints & Helpers

- `constraint_check` returns numerator/denominator pairs for supported predicates. Pair it with the `poly` subcommands (`concyclic`, `circumcenter`, `angle`) to inspect raw algebra.
- The engine auto-learns conjugate substitutions from single-conjugate equations; retrieve them with `learned_rules`.

## Testing

```
. .venv/bin/activate
pytest
```

Tests live under `tests/` and cover constraint reduction, reflections, and special-point constructions. Add new suites as you extend the engine or CLI.
