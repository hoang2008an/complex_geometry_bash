# Complex Geometry Bash — Library Usage (Python API)

## What is it?

Complex Geometry Bash is a SymPy-powered symbolic geometry engine for the complex plane. It provides a compact API for
describing geometric configurations (points + constraints), constructing derived objects (centers, reflections,
projections, intersections), and verifying relations via exact symbolic expressions.

## Table of Contents

- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Concepts (brief)](#concepts-brief)
- [API Reference (concise)](#api-reference-concise)

---

## Introduction

Complex Geometry Bash is a **symbolic geometry engine** for the complex plane built on SymPy. You describe a geometric
configuration by:

1. registering points by label,
2. adding geometric constraints (stored internally as polynomial equations),
3. constructing derived points (centers, projections, reflections, intersections),
4. querying symbolic expressions and constraint checks.

The engine also maintains a lightweight substitution system (including auto-learned conjugate substitutions) to simplify
expressions across the workflow.

---

## Quick Start

### 1) Install dependencies

```bash
python -m venv .venv
```

- Windows (PowerShell):

```bash
.venv\Scripts\activate
pip install -r requirements.txt
```

- macOS/Linux:

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

### 2) Import

When running from the repository, ensure `src/lib` is on `PYTHONPATH` (the `bin/geometry` wrapper handles this
automatically).

```python
from geometry_engine import GeometryEngine, GeometryError
```

### 3) Create an engine, add points, add constraints

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()

# Note: O is pre-registered at the origin. Other points must be registered explicitly.
for name in ["A", "B", "C"]:
    eng.add_point(name)

# Example constraints (optional; depends on your problem)
eng.add_unit_circle("A")
eng.add_unit_circle("B")
eng.add_unit_circle("C")
```

### 4) Construct derived points and inspect results

```python
eng.circumcenter("A", "B", "C", "U")
eng.orthocenter_via_altitudes("A", "B", "C", "H")
eng.centroid("A", "B", "C", "G")

print(eng.point_summary(["A", "B", "C", "U", "H", "G"], style="text"))
```

### 5) Verify a relationship

```python
num, den = eng.constraint_conjugate_free("collinear", ["U", "G", "H"])[0]
print("Collinear check numerator =", eng.format_expr(num, style="text"))
```

---

## Concepts (brief)

### Points and symbols

Each point label `X` corresponds to two independent symbols:
- `z_X` (complex coordinate)
- `zb_X` (formal conjugate symbol)

These are not assumed to be complex conjugates unless constraints force that relationship; the engine can learn such
substitutions during solving.

### Constraints

Constraint methods (`add_*`) store polynomial equations (assumed equal to 0) in `engine.constraints` and may trigger
auto-learning of substitutions.

### Lines and circles

Most workflows use point-based constraints directly. The engine also provides lightweight `Line` and `Circle` objects
primarily to support object-based intersections and tangents via named registries.

Example: intersect a stored line with a stored circle.

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()
for name in ["A", "B", "C", "D"]:
    eng.add_point(name)

eng.register_line("lAB", eng.line_through_points("A", "B"))
eng.circle_from_three_points("omega", "A", "C", "D")
eng.line_circle_object_intersections("lAB", "omega", ["X1", "X2"])
```

Object-oriented helpers accept either the registered name or the `Line` / `Circle` instance itself. That applies to
methods such as `line_value(...)`, `line_intersection(...)`, `add_point_on_line(...)`, `add_point_on_circle(...)`,
`line_circle_object_intersections(...)`, `circle_object_intersections(...)`, `tangent_lines_from_point_to_circle(...)`,
and `radical_line_from_circles(...)`.

Current Python signatures use generic object parameters on those methods. If you call them with keywords, use:

- `line_intersection(line1=..., line2=..., name=...)`
- `add_point_on_line(line=..., point=...)`
- `add_point_on_circle(circle=..., point=...)`
- `line_circle_object_intersections(line=..., circle=..., point_names=...)`
- `circle_object_intersections(circle1=..., circle2=..., point_names=...)`
- `tangent_lines_from_point_to_circle(point=..., circle=..., line_names=...)`
- `radical_line_from_circles(circle1=..., circle2=..., line_name=...)`

Examples:

```python
axis = eng.register_line("axis", eng.line_through_points("P1", "P2"))
omega = eng.circle_from_three_points("omega", "A", "B", "C")

eng.line_circle_object_intersections("axis", "omega", ["X1", "X2"])
eng.line_circle_object_intersections(axis, omega, ["Y1", "Y2"])
eng.line_circle_object_intersections(line=axis, circle=omega, point_names=["Z1", "Z2"])
```

---

## API Reference (concise)

### Core engine

- `GeometryEngine()`
- `GeometryError`

### Point management

- `add_point(name)`
- `add_unit_circle(name)`
- `add_distinct_points(P, Q)`
- `z(name)`, `zb(name)`

### Constraints (store relations)

- `add_collinear(A, B, C)`
- `add_perpendicular(A, B, C, D)`
- `add_concyclic(A, B, C, D)`
- `add_angle_value(A, B, C, angle_radians)`
- `add_triangle_similarity(A, B, C, D, E, F, directed=True)`
- `add_triangle_congruence(A, B, C, D, E, F, directed=True)`
- `add_point_reflection(P, O, Q)`
- `add_line_reflection(P, A, B, Q)`
- `add_projection_to_line(P, A, B, H)`
- `add_isogonal_reflection(A, B, C, D, E)`
- `add_isogonal_conjugate(A, B, C, P, Q)`
- `add_angle_bisector_either(A, B, C, D)`

### Constructions (assign points)

- `circumcenter(A, B, C, U)`
- `orthocenter_via_altitudes(A, B, C, H)`
- `centroid(A, B, C, G)`
- `euler_center(A, B, C, N, circumcenter_name=None, orthocenter_name=None)`
- `lemoine_point(A, B, C, L, circumcenter_name=None)`
- `midpoint(P, Q, M)`
- `reflect_point_over_point(P, O, Q)`
- `reflect_point_over_line(P, A, B, Q)`
- `project_point_to_line(P, A, B, H)`
- `isogonal_conjugate_point(A, B, C, P, Q)`
- `add_fermat_points(A, B, C, F1, F2)`

### Lines and circles

- **Line objects**: `line_through_points(P, Q)`, `symbolic_line(name)`, `line_intersection(line1, line2, name)`
- **Point incidence**: `add_point_on_line(line_or_name, point)`, `add_point_on_circle(circle_or_name, point)`
- **Line registry**: `register_line(name, line)`, `get_line(name)`
- **Circle registry**: `circle_from_three_points(name, A, B, C)`, `register_circle(name, circle)`, `get_circle(name)`
- **Intersections**: `line_circle_intersection(...)`, `circle_intersection(...)`,
  `line_circle_object_intersections(line, circle, point_names, avoid=None)`,
  `circle_object_intersections(circle1, circle2, point_names, avoid=None)`
- **Tangents / advanced**: `tangent_lines_from_point_to_circle(point, circle, line_names, tangent_point_names=None)`,
  `radical_line_from_circles(circle1, circle2, line_name)`

### Query and debugging

- `squared_distance(P, Q)`
- `format_expr(expr, style="latex"|"text")`
- `point_summary(names=None, style="latex"|"text")`
- `constraint_strings(style="text")`
- `learned_rules()`, `display_learned_rules(style="latex"|"text")`
- `constraint_conjugate_free(constraint_type, args)`

### Unit triangle helpers

These helpers are a core part of the unit-circle workflow. In particular, they are the intended entry points for
constructing the incenter/excenters and several related points when working with a triangle on the unit circle.

- `set_main_unit_triangle(A, B, C, root_names=None)`
- `main_triangle_incenter(name)`
- `main_triangle_excenter(which, name)`
- `main_triangle_arc_midpoint(which, name, containing_vertex=False)`
