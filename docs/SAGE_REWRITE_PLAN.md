# Isolated Sage Math Mode V1

Sage mode is an experimental backend for exact rational-function workflows. It
is intentionally separate from the existing SymPy engine:

- Engine: `src/lib/geometry_sage_engine.py`
- CLI: `src/lib/geometry_sage_cli.py`
- Launcher: `bin/geometry-sage`

The default `bin/geometry` command and `src/lib/geometry_engine.py` remain the
normal SymPy path.

## Backend Model

`SageGeometryEngine` builds a Sage fraction field with coordinate variables
`z_X` and `zb_X` for registered points. The default coefficient field is
`QQbar`, but purely rational symbolic scripts can use `QQ`:

```bash
bin/geometry-sage run --field QQ scripts/perp_bisector_OH.json
```

The CLI pre-scans JSON scripts for point labels so the field is created before
execution.

Sage mode v1 does not use Sage's symbolic ring as a fallback.

Use `QQbar` when exact algebraic constants are required, for example Fermat
points or fixed algebraic angle constants. Use `QQ` for unit-circle/rational
workflows that do not need algebraic coefficients; it can be much faster because
the fraction-field coefficients are simpler.

## Supported Workflows

The v1 surface supports rational/linear exact operations:

- point registration
- fixed exact point values via `set_point_value`
- unit-circle declarations and learned conjugate rules
- collinear, perpendicular, fixed-angle, and angle-bisector constraints
- midpoint, centroid, circumcenter, orthocenter, projection
- line intersection
- line-circle intersection when one supplied line endpoint is already known to
  lie on the circle; the other intersection is selected with `avoid`
- Lemoine point
- Fermat points using fixed algebraic third-root constants
- point summaries, learned rules, constraints, and constraint checks

Run a script with:

```bash
bin/geometry-sage run scripts/example.json
```

Run the small Sage demo with:

```bash
bin/geometry-sage demo
```

## Square-Root Policy

General symbolic square roots are out of scope for Sage mode v1:

- no general radical towers
- no default symbolic-ring fallback
- no circle intersections
- no general line-circle intersections
- no tangent or radical-line workflows that require branch selection

Unsupported branch-producing operations fail with a clear message mentioning the
square-root branch limitation.

The current line-circle support is only the rational secant case where one
intersection is already known. If general circle intersections are added later,
they should stay isolated inside `geometry_sage_engine.py`, return explicit
candidate branches, and support `avoid` without changing the SymPy engine.

## Benchmark Policy

`P4TQH` must not be used as a required Sage-mode benchmark or regression test in
v1. Sage mode should first grow coverage around small rational workflows before
branch-heavy scripts are considered.
