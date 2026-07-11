> **⚠️ Unofficial Guide**
>
> This guide was mostly AI-generated from an older commit and is **not** official documentation. It may contain outdated or inaccurate information. Please use it only as an informal reference.


# Complex Geometry Bash — Library Usage Guide

This comprehensive guide explains how to use the Complex Geometry Bash library as a Python API. It covers all public methods and provides practical examples for common geometric constructions.

## Overview

The library provides a symbolic geometry engine for the complex plane, designed for olympiad-style geometry problems and algebraic verification of geometric relations.

### Key Features

- **Symbolic Computation**: All coordinates are complex numbers with automatic conjugate handling
- **Constraint System**: Geometric relations expressed as polynomial equations
- **Auto-Learning**: Engine automatically discovers conjugate substitutions from constraints
- **Black-Box API**: Focus on geometric description, let the engine handle algebra
- **Extensible**: Easy to add new constructions while maintaining API stability

### Quick Import

```python
from geometry_engine import GeometryEngine, GeometryError
```

### Basic Usage Pattern

```python
# Create engine
eng = GeometryEngine()

# Define points
eng.add_point("A")
eng.add_point("B")
eng.add_point("C")

# Add constraints
eng.add_unit_circle("A")  # A on unit circle
eng.add_collinear("A", "B", "C")  # A, B, C collinear

# Query results
z_A = eng.z("A")  # Get coordinate expression
```

## Table of Contents

- [Concepts and Mental Model](#concepts-and-mental-model)
- [Quick Start](#quick-start)
- [Core API Reference](#core-api-overview)
- [Common Usage Recipes](#common-usage-recipes)
- [Error Handling](#error-handling)
- [Import and Configuration](#how-to-import-and-configure-in-your-own-project)
- [Method Reference](#minimal-reference-of-key-methods)


## Concepts and Mental Model

The engine operates in the complex plane using a symbolic approach optimized for olympiad geometry.

### Coordinate System

Each geometric point `X` is represented by **two independent complex symbols**:
- `z_X` — the point's complex coordinate
- `zb_X` — the formal conjugate symbol (not necessarily `conjugate(z_X)`)

This dual-symbol approach allows the engine to handle:
- Points not constrained to the real axis
- Complex conjugate relationships discovered through constraints
- Algebraic simplifications via conjugate substitutions

### Constraint System

Geometric relations are encoded as **polynomial equations** in the point symbols:

- **Collinearity**: `(z_A - z_C)(zb_B - zb_C) - (zb_A - zb_C)(z_B - z_C) = 0`
- **Perpendicularity**: `(z_A - z_B)(zb_C - zb_D) + (zb_A - zb_B)(z_C - z_D) = 0`
- **Concyclicity**: Complex polynomial eliminating the circumcenter
- **Unit Circle**: `z_X * zb_X - 1 = 0`

## Conjugate Substitution

The engine automatically discovers **conjugate substitutions** when the given constraints uniquely determine a conjugate symbol in terms of non-conjugate symbols. For example, if a constraint implies

```text
zb_A = 1 / z_A
```

the engine records this substitution and automatically applies it to all subsequent expressions.

> **Note**
>
> This behavior was recently reworked. Point coordinates are **no longer computed automatically**. To evaluate a point, you must explicitly invoke the `calculate_point` operation.
>See the examples under **`benchmark/IMO/2010-2014`** to know how to use the operations.

### Workflow

1. **Setup**: Register points and add geometric constraints
2. **Construct**: Use high-level methods for centers, reflections, etc.
3. **Query**: Retrieve symbolic expressions for coordinates, distances, or constraint polynomials
4. **Verify**: Check if constraints are satisfied (polynomials evaluate to zero)

### Key Invariants

- **Squared Distance**: `|P - Q|² = (z_P - z_Q)(zb_P - zb_Q)`
- **Cross Ratio**: `(z_A - z_C)(z_B - z_D) / ((z_A - z_D)(z_B - z_C))`
- **Angle Relations**: Complex exponential forms for angle constraints


## 2. Quick Start

## Quick Start

### Basic Setup

```python
from geometry_engine import GeometryEngine

# Create engine instance
eng = GeometryEngine()

# Point O is automatically registered at origin (0, 0)
# Add your points
eng.add_point("A")
eng.add_point("B")
eng.add_point("C")

# Constrain points to unit circle (optional)
eng.add_unit_circle("A")
eng.add_unit_circle("B")
eng.add_unit_circle("C")
```

### Retrieving Coordinates

```python
# Get symbolic expressions (with all substitutions applied)
z_A = eng.z("A")   # z-coordinate of A
zb_A = eng.zb("A") # conjugate coordinate of A

# Format for display
print("A =", eng.format_expr(z_A, style="text"))
```

### Triangle Centers Example

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()

# Set up equilateral triangle on unit circle
for vertex in ["A", "B", "C"]:
    eng.add_point(vertex)
    eng.add_unit_circle(vertex)

# Construct key centers
eng.circumcenter("A", "B", "C", "U")      # Circumcenter
eng.orthocenter_via_altitudes("A", "B", "C", "H")  # Orthocenter
eng.centroid("A", "B", "C", "G")          # Centroid
eng.euler_center("A", "B", "C", "N")      # Euler center (midpoint U-H)

# Display results
centers = ["U", "H", "G", "N"]
for center in centers:
    coord = eng.format_expr(eng.z(center), style="text")
    print(f"{center} = {coord}")
```

### Verifying Geometric Properties

```python
# Check if points are collinear (Euler line property)
num, den = eng.constraint_conjugate_free("collinear", ["U", "G", "H"])[0]
print("Euler line polynomial:", eng.format_expr(num, style="text"))
# Should equal zero for unit circle triangles
```


## Core API Reference

This section provides a comprehensive reference of all public methods in `GeometryEngine`. Methods are grouped by functionality.

### Points and Basic Utilities

#### Point Registration
- `add_point(name: str)` - Register a new point with given label
- `add_unit_circle(name: str)` - Constrain point to lie on unit circle (`z * zb = 1`)
- `add_distinct_points(P: str, Q: str)` - Ensure points P and Q are distinct

#### Coordinate Access
- `z(name: str)` - Get simplified z-coordinate expression
- `zb(name: str)` - Get simplified zb-coordinate expression
- `squared_distance(P: str, Q: str)` - Compute `|P - Q|²` symbolically
- `cross_ratio(A: str, B: str, C: str, D: str)` - Compute cross ratio `(A,B;C,D)`

#### Display and Introspection
- `format_expr(expr, style="latex"|"text")` - Format symbolic expression
- `format_symbol(symbol, style="latex"|"text")` - Format individual symbol
- `point_summary(names=None, style="latex"|"text")` - Summary of point coordinates
- `constraint_strings(style="text")` - List all stored constraints
- `learned_rules()` - Dictionary of auto-learned conjugate substitutions
- `display_learned_rules(style="latex"|"text")` - Formatted learned rules


### Geometric Constraints

#### Basic Relations
- `add_collinear(A, B, C)` - Points A, B, C are collinear
- `add_perpendicular(A, B, C, D)` - Lines AB ⟂ CD
- `add_concyclic(A, B, C, D)` - Points A, B, C, D are concyclic
- `add_angle_value(A, B, C, angle_radians)` - ∠ABC = specified angle (radians)

#### Triangle Relations
- `add_triangle_similarity(A, B, C, D, E, F, directed=True)` - Triangles ABC ~ DEF
- `add_triangle_congruence(A, B, C, D, E, F, directed=True)` - Triangles ABC ≅ DEF

#### Point Constructions (as Constraints)
- `add_circumcenter(A, B, C, U)` - U is circumcenter of △ABC
- `add_midpoint(P, Q, M)` - M is midpoint of PQ
- `add_centroid_constraint(A, B, C, G)` - G is centroid of △ABC

#### Reflections and Projections
- `add_point_reflection(P, O, Q)` - Q is reflection of P over point O
- `add_line_reflection(P, A, B, Q)` - Q is reflection of P over line AB
- `add_projection_to_line(P, A, B, H)` - H is projection of P onto line AB

#### Advanced Triangle Properties
- `add_isogonal_reflection(A, B, C, D, E)` - AD, AE are isogonal wrt ∠BAC
- `add_isogonal_conjugate(A, B, C, P, Q)` - Q is isogonal conjugate of P wrt △ABC

#### Low-level Polynomial Helpers
For symbolic inspection without adding constraints:
- `collinear_poly(A, B, C)` - Raw collinearity polynomial
- `perpendicular_poly(A, B, C, D)` - Raw perpendicularity polynomial
- `concyclic_poly(A, B, C, D)` - Raw concyclicity polynomial
- `angle_value_poly(A, B, C, angle)` - Raw angle constraint polynomial
- `constraint_conjugate_free(constraint_type, args)` - Get conjugate-free form


### High-Level Constructions

These methods construct new points by solving geometric constraints and return their z-coordinates.

#### Triangle Centers
- `circumcenter(A, B, C, U)` - Construct circumcenter U of △ABC
- `orthocenter_via_altitudes(A, B, C, H)` - Construct orthocenter H of △ABC
- `centroid(A, B, C, G)` - Construct centroid G of △ABC
- `euler_center(A, B, C, N, circumcenter_name=None, orthocenter_name=None)` - Construct Euler center N
- `lemoine_point(A, B, C, L, circumcenter_name=None)` - Construct symmedian point L

#### Basic Operations
- `midpoint(P, Q, M)` - Construct midpoint M of segment PQ
- `reflect_point_over_point(P, O, Q)` - Reflect P over point O to get Q
- `reflect_point_over_line(P, A, B, Q)` - Reflect P over line AB to get Q
- `project_point_to_line(P, A, B, H)` - Project P onto line AB to get H

#### Special Points
- `isogonal_conjugate_point(A, B, C, P, Q)` - Isogonal conjugate of P wrt △ABC
- `add_fermat_points(A, B, C, F1, F2)` - Construct first and second Fermat points

#### Lines and Circles
- `line_through_points(P, Q)` - Create Line object through P and Q
- `symbolic_line(name)` - Create parameterized line with given name
- `line_intersection(line1, line2, name)` - Intersection point of two Line objects
- `circle_from_three_points(name, A, B, C)` - Circle through A, B, C

#### Intersections and Tangents
- `line_circle_intersection(line_pt1, line_pt2, center, radius_pt, name, avoid=None)` - Line-circle intersection
- `circle_intersection(center1, radius_pt1, center2, radius_pt2, name, avoid=None)` - Circle-circle intersection
- `line_circle_object_intersections(line_name, circle_name, point_names, avoid=None)` - Multiple line-circle intersections
- `circle_object_intersections(circle1_name, circle2_name, point_names, avoid=None)` - Multiple circle-circle intersections
- `tangent_lines_from_point_to_circle(point, circle_name, line_names, tangent_point_names=None)` - Tangent lines from point to circle

#### Advanced Constructions
- `radical_line_from_circles(circle1_name, circle2_name, line_name)` - Radical axis of two circles


### Unit Triangle Helpers

For triangles with all vertices on the unit circle, the engine provides specialized methods that use classical algebraic formulas.

#### Setup
```python
# First set up the main unit triangle
eng.set_main_unit_triangle("A", "B", "C")  # All vertices must be on unit circle
```

#### Incenter and Excenters
- `main_triangle_incenter(name)` - Construct incenter using formula -(xy + yz + zx)
- `main_triangle_excenter(which, name)` - Construct excenter opposite vertex `which` ∈ {"A", "B", "C"}

#### Arc Midpoints
- `main_triangle_arc_midpoint(which, name, containing_vertex=False)` - Midpoint of arc opposite vertex `which`
  - `containing_vertex=False`: Midpoint of arc not containing the vertex
  - `containing_vertex=True`: Midpoint of arc containing the vertex

These methods use square-root formulations and maintain consistent internal auxiliary labels for reproducible results.


### Introspection and Debugging

#### Learned Rules
- `learned_rules()` - Dictionary of auto-learned conjugate substitutions
- `display_learned_rules(style="latex"|"text")` - Formatted display of learned rules

#### Point Information
- `point_summary(names=None, style="latex"|"text")` - Coordinate summary for points

#### Constraints
- `constraint_strings(style="text")` - List of all stored constraint polynomials

#### Formatting Utilities
- `format_expr(expr, style="latex"|"text")` - Format symbolic expressions
- `format_symbol(symbol, style="latex"|"text")` - Format individual symbols


## Common Usage Recipes

### Checking Concyclicity

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()
for name in ["A", "B", "C", "D"]:
    eng.add_point(name)

# Add your constraints (e.g., unit circle, collinearities)
eng.add_unit_circle("A")
eng.add_unit_circle("B")
eng.add_unit_circle("C")
eng.add_unit_circle("D")

# Check if A,B,C,D are concyclic
num, den = eng.constraint_conjugate_free("concyclic", ["A", "B", "C", "D"])[0]
print("Concyclic polynomial numerator:", eng.format_expr(num, style="text"))
# Should equal zero if points are concyclic
```

### Verifying Euler Line

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()

# Unit circle triangle
for vertex in ["A", "B", "C"]:
    eng.add_point(vertex)
    eng.add_unit_circle(vertex)

# Construct centers
eng.circumcenter("A", "B", "C", "U")
eng.orthocenter_via_altitudes("A", "B", "C", "H")
eng.centroid("A", "B", "C", "G")

# Check Euler line: U, G, H should be collinear
num, den = eng.constraint_conjugate_free("collinear", ["U", "G", "H"])[0]
print("Euler line check:", eng.format_expr(num, style="text"))
# Should be zero for unit circle triangles
```

### Isogonal Conjugates

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()

# Triangle vertices on unit circle
for name in ["A", "B", "C"]:
    eng.add_point(name)
    eng.add_unit_circle(name)

# Arbitrary point P
eng.add_point("P")

# Construct isogonal conjugate Q of P
eng.isogonal_conjugate_point("A", "B", "C", "P", "Q")
print("Isogonal conjugate Q =", eng.format_expr(eng.z("Q"), style="text"))
```

### Fermat Points

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()

# Triangle vertices
for name in ["A", "B", "C"]:
    eng.add_point(name)

# Construct both Fermat points
eng.add_fermat_points("A", "B", "C", "F1", "F2")
print("First Fermat point F1 =", eng.format_expr(eng.z("F1"), style="text"))
print("Second Fermat point F2 =", eng.format_expr(eng.z("F2"), style="text"))
```

### Line and Circle Intersections

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()

# Define line through A,B and circle with center O through C
eng.add_point("A")
eng.add_point("B")
eng.add_point("O")
eng.add_point("C")

# Create geometric objects
line = eng.line_through_points("A", "B")
circle = eng.circle_from_three_points("circ", "O", "C", "D")  # Need third point

# Find intersections
eng.line_circle_object_intersections("line", "circ", ["P1", "P2"])
print("Intersection P1 =", eng.format_expr(eng.z("P1"), style="text"))
print("Intersection P2 =", eng.format_expr(eng.z("P2"), style="text"))
```


## Error Handling

All geometry operations raise `GeometryError` for invalid configurations:

```python
from geometry_engine import GeometryEngine, GeometryError

eng = GeometryEngine()
try:
    eng.add_collinear("A", "B", "C")  # Points not registered
except GeometryError as exc:
    print("Error:", exc)
```

Common error scenarios:
- **Missing points**: Referencing unregistered point labels
- **Invalid operations**: Degenerate inputs (identical points for line definition)
- **Inconsistent constraints**: Over-constrained or impossible geometric configurations
- **Duplicate labels**: Using same label where distinct points required

## Import and Configuration

### Basic Import

```python
from geometry_engine import GeometryEngine, GeometryError
```

### PYTHONPATH Setup

For direct repository usage:

```bash
export PYTHONPATH=src/lib:src
```

### Project Integration

1. **Create one engine per problem**:
   ```python
   eng = GeometryEngine()  # Fresh instance for each geometric configuration
   ```

2. **Use only documented public methods** - the API is designed to be stable

3. **SymPy expressions** for angles and parameters:
   ```python
   import sympy as sp
   eng.add_angle_value("A", "B", "C", sp.pi/2)  # 90-degree angle
   ```

## Method Reference

### Core Construction Methods
- **Triangle centers**: `circumcenter`, `orthocenter_via_altitudes`, `centroid`, `euler_center`, `lemoine_point`
- **Basic operations**: `midpoint`, `reflect_point_over_point`, `reflect_point_over_line`, `project_point_to_line`
- **Special points**: `isogonal_conjugate_point`, `add_fermat_points`
- **Geometric objects**: `line_through_points`, `circle_from_three_points`
- **Intersections**: `line_intersection`, `line_circle_intersection`, `circle_intersection`

### Constraint Methods
- **Basic relations**: `add_collinear`, `add_perpendicular`, `add_concyclic`, `add_angle_value`
- **Triangle properties**: `add_triangle_similarity`, `add_triangle_congruence`
- **Transformations**: `add_point_reflection`, `add_line_reflection`, `add_projection_to_line`
- **Advanced**: `add_isogonal_conjugate`, `add_distinct_points`

### Query and Display Methods
- **Coordinates**: `z`, `zb`, `squared_distance`, `cross_ratio`
- **Inspection**: `learned_rules`, `point_summary`, `constraint_strings`
- **Formatting**: `format_expr`, `display_learned_rules`

### Unit Triangle Methods
- **Setup**: `set_main_unit_triangle`
- **Centers**: `main_triangle_incenter`, `main_triangle_excenter`
- **Arc points**: `main_triangle_arc_midpoint`

Use this document as your primary reference when scripting with the engine.