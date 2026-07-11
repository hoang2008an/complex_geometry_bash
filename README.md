# Complex Geometry Bash

A powerful symbolic geometry engine for complex-number computations in the circumcircle model, built with Python and SymPy. Perfect for olympiad geometry problems, automated theorem proving, and educational explorations of advanced geometric constructions.

## Overview

This project provides a comprehensive toolkit for symbolic geometry:

- **Core Library**: [`GeometryEngine`](src/lib/geometry_engine.py:70) - A black-box API for complex plane geometry with automatic constraint solving and conjugate substitutions.
- **CLI Tool**: [`geometry_cli`](src/lib/geometry_cli.py:577) - JSON-driven command-line interface for reproducible geometric workflows.
- **Rich Features**: Supports points, lines, circles, reflections, projections, triangle centers, similarity checks, and more.
- **Educational Focus**: Designed for olympiad-style problems with algebraic verification of geometric relations.

The library handles all symbolic complexity internally, allowing users to describe geometric configurations and query results without deep SymPy knowledge.

## Key Features

- **Symbolic Coordinate System**: Each point has independent complex coordinates `(z, zb)` with automatic conjugate learning.
- **Geometric Primitives**: Points, lines, circles with intersection and tangent computations.
- **Triangle Centers**: Circumcenter, orthocenter, centroid, Euler center, Lemoine point, Fermat points.
- **Transformations**: Reflections over points/lines, projections, isogonal conjugates.
- **Constraint System**: Collinearity, perpendicularity, concyclicity, angle values, triangle similarity/congruence.
- **Unit Circle Support**: Specialized helpers for unit-circle triangles (incenters, excenters, arc midpoints).
- **Extensible Design**: Easy to add new constructions while maintaining API stability.

For detailed API documentation and usage examples, see:
- [`docs/LIBRARY_USAGE.md`](docs/LIBRARY_USAGE.md)


## Installation and Setup

### Quick Start

```bash
# Clone the repository
git clone <repository-url>
cd complex_geometry_bash

# Install dependencies
pip install -r requirements.txt
```

### Using the Library
Please check out benchmark/IMO for example usage
```bash
python src/lib/geometry_cli.py run path/to/scripts.json
```
<!-- The engine requires `src/lib` on your `PYTHONPATH`. Two options:

**Option 1: Direct import (recommended for development)**
```bash
export PYTHONPATH=src/lib:src
```

**Option 2: Use the wrapper script**
```bash
# The bin/geometry script sets PYTHONPATH automatically
./bin/geometry demo
```  -->

### Runtime Dependencies

- `sympy>=1.12` (see [`requirements.txt`](requirements.txt))
- Python 3.8+

### Development Setup

For contributors:

```bash
# Install in editable mode with dev dependencies
pip install -e .
pip install pytest  # For testing

# Run tests
pytest

# Run demo
./bin/geometry demo
```


## Using the Library

### Basic Usage

```python
from geometry_engine import GeometryEngine, GeometryError

# Create engine instance
eng = GeometryEngine()

# Set up a unit-circle triangle
for name in ["A", "B", "C"]:
    eng.add_point(name)
    eng.add_unit_circle(name)

# Compute triangle centers
eng.circumcenter("A", "B", "C", "U")  # Circumcenter
eng.orthocenter_via_altitudes("A", "B", "C", "H")  # Orthocenter
eng.centroid("A", "B", "C", "G")  # Centroid

# Display results
print("Circumcenter U =", eng.format_expr(eng.z("U"), style="text"))
print("Orthocenter H =", eng.format_expr(eng.z("H"), style="text"))
print("Centroid G =", eng.format_expr(eng.z("G"), style="text"))
```

### Core Concepts

- **Points**: Each point `X` has complex coordinates `(z_X, zb_X)` with automatic conjugate handling.
- **Constraints**: Geometric relations are stored as polynomial equations.
- **Constructions**: High-level methods for centers, reflections, projections, etc.
- **Auto-Learning**: The engine automatically discovers conjugate substitutions from constraints.
- **Black-Box API**: Focus on geometric description, let the engine handle symbolic algebra.

### Advanced Example: Verifying Euler Line

```python
from geometry_engine import GeometryEngine

eng = GeometryEngine()

# Unit circle triangle
for name in ["A", "B", "C"]:
    eng.add_point(name)
    eng.add_unit_circle(name)

# Triangle centers
eng.circumcenter("A", "B", "C", "U")
eng.orthocenter_via_altitudes("A", "B", "C", "H")
eng.centroid("A", "B", "C", "G")

# Check if U, G, H are collinear (Euler line property)
num, den = eng.constraint_conjugate_free("collinear", ["U", "G", "H"])[0]
print("Euler line polynomial:", eng.format_expr(num, style="text"))
# Should be zero for unit circle triangles
```

For comprehensive API documentation, see [`docs/LIBRARY_USAGE.md`](docs/LIBRARY_USAGE.md).


## CLI Usage

The command-line interface provides a JSON-based workflow system for reproducible geometric computations. Perfect for batch processing, sharing configurations, and automated testing.

### Getting Started

```bash
# Solve IMO 2011 Problem 6
python src/lib/geometry_cli.py run benchmark/IMO/2010_2014/imo_2011_p6.json
```
The solver may take a few minutes to finish. It will produce a large amount of output, but you only need to verify that the final result contains:

N = 0
D = 1

If these values are reported, the run completed successfully.
### JSON Script Format

Scripts are JSON files with a `steps` array. Each step is an object with an `op` field and parameters:

```json
{
  "steps": [
    {"op": "add_point", "name": "A", "unit_circle": true},
    {"op": "add_point", "name": "B", "unit_circle": true},
    {"op": "add_point", "name": "C", "unit_circle": true},
    {"op": "circumcenter", "args": ["A", "B", "C", "U"]},
    {"op": "orthocenter", "args": ["A", "B", "C", "H"]},
    {"op": "print_points", "names": ["U", "H"]},
    {"op": "constraint_check", "constraint": "collinear", "args": ["A", "U", "H"]}
  ]
}
```

### Available Operations

#### Point Management
- `add_point`, `add_unit_circle`

#### Constraints & Predicates
- `add_collinear`, `add_perpendicular`, `add_concyclic`, `add_angle_value`
- `add_triangle_similarity`, `add_triangle_congruence`

#### Constructions
- **Centers**: `circumcenter`, `orthocenter`, `centroid`, `euler_center`, `lemoine_point`
- **Reflections**: `reflect_point_over_point`, `reflect_point_over_line`
- **Projections**: `project_point_to_line`
- **Special Points**: `add_fermat_points`, `isogonal_conjugate`
- **Lines & Circles**: `line_from_points`, `circle_from_three_points`, intersections, tangents

#### Unit Triangle Features
- `set_main_unit_triangle`, `unit_triangle_incenter`, `unit_triangle_excenter`, `unit_triangle_arc_midpoint`

#### Inspection & Debugging
- `print_points`, `constraint_check`, `learned_rules`, `print_constraints`, `squared_distance`

See [`geometry_cli.py`](src/lib/geometry_cli.py:61) for the complete operation mapping.


## Testing and Examples

### Running Tests

The project includes comprehensive tests covering all major functionality:

```bash
# Activate virtual environment
source .venv/bin/activate

# Run all tests
pytest

# Run specific test file
pytest tests/test_special_points.py

# Run with verbose output
pytest -v
```

### Example Scripts

The [`scripts/`](scripts/) directory contains JSON scripts demonstrating various geometric constructions:

```bash
# Run an example script
./bin/geometry run scripts/incenter_arc_demo.json

# List available scripts
ls scripts/*.json
```

Key examples include:
- `incenter_arc_demo.json` - Unit triangle incenter and arc constructions
- `perp_bisector_OH.json` - Perpendicular bisector and orthocenter relations
- Various olympiad-style problems (P3, P4 series)


## Contributing

### For Library Users

If you're using the library, the key resources are:
- This README for quick start and overview
- [`docs/LIBRARY_USAGE.md`](docs/LIBRARY_USAGE.md) for complete API reference
- Example scripts in [`scripts/`](scripts/) for common patterns

### For Contributors

The project welcomes contributions! See [`AGENTS.md`](AGENTS.md) for detailed development guidelines.

#### Architecture Overview

- **`src/lib/geometry_engine.py`** - Core symbolic geometry engine with constraint solving
- **`src/lib/geometry_cli.py`** - CLI interface and JSON operation dispatch
- **`bin/geometry`** - Thin wrapper script setting PYTHONPATH
- **`tests/`** - Comprehensive test suite
- **`scripts/`** - Example JSON workflows

#### Development Workflow

1. **Setup development environment**:
   ```bash
   git clone <repository-url>
   cd complex_geometry_bash
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   pip install pytest  # For testing
   ```

2. **Follow coding standards**:
   - PEP 8 style with 4-space indentation (2 spaces allowed in SymPy-heavy blocks)
   - Type hints for public methods
   - SymPy-native expressions throughout
   - Shared `*_poly` and `add_*` helper patterns

3. **Add new features**:
   - Implement symbolic logic in `geometry_engine.py`
   - Add CLI operations in `geometry_cli.py`
   - Write comprehensive tests
   - Update documentation

4. **Testing**:
   ```bash
   pytest  # Run all tests
   pytest tests/test_your_feature.py  # Run specific tests
   ```

5. **Commit guidelines**:
   - Use Conventional Commits: `feat:`, `fix:`, `docs:`, `chore:`
   - Keep subject lines under 72 characters
   - Reference issues with `Refs #123`

#### Project Structure

```
complex_geometry_bash/
├── src/lib/
│   ├── geometry_engine.py    # Core engine
│   ├── geometry_cli.py       # CLI interface
│   └── __init__.py
├── bin/
│   └── geometry              # Launcher script
├── tests/                    # Test suite
├── scripts/                  # Example JSON scripts
├── docs/                     # Documentation
├── AGENTS.md                 # Development guidelines
└── README.md                 # This file
```

The modular design ensures the `GeometryEngine` remains a stable, black-box API for users while allowing easy extension by contributors.
