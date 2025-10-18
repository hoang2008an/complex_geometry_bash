# Repository Guidelines

Complex Geometry Bash contains a symbolic complex-geometry engine powered by Python/SymPy. Use these practices when extending the system.

## Project Structure & Module Organization
- Core engine lives in `src/lib/geometry_engine.py`; keep all symbolic logic and auto-learning rules here.
- CLI orchestration sits in `src/lib/geometry_cli.py`; add new JSON ops or parser subcommands in this module.
- `bin/geometry` is the thin launcher; keep shell logic minimal and delegate to Python.
- Declare dependencies in `requirements.txt`; create scripts or fixtures under `scripts/` and `tests/` as the project grows.
- A local `.venv/` is expected but untracked; document any additional setup in README snippets.

## Build, Test, and Development Commands
- `python3 -m venv .venv` then `.venv/bin/pip install -r requirements.txt` — bootstrap tooling (SymPy).
- `bin/geometry demo` — smoke-test the canonical workflow and ensure learned substitutions print as expected.
- `bin/geometry run scripts/<scenario>.json` — execute custom problem setups; capture failing scripts under `scripts/`.
- `bin/geometry poly concyclic A B C D` / `bin/geometry poly circumcenter A B C U` / `bin/geometry poly angle A B C <θ>` — inspect raw constraint polynomials.
- `pytest` (placeholders today) — target new tests to `tests/` covering engine methods and CLI behaviors.

## Coding Style & Naming Conventions
- Follow PEP 8 with two-space indents inside SymPy-heavy blocks only when readability demands; otherwise default to four spaces.
- Type-hint public methods, prefer `dataclass` records for point state, and keep expressions SymPy-native (`sp.*`).
- When adding constraints, expose a `*_poly` helper and an `add_*` guard so CLI and scripts share one implementation path.
- Bin scripts stay POSIX-compatible; expand shell usage only when unavoidable.

## Testing Guidelines
- Organize tests as `tests/test_<area>.py`; favor `pytest` parametrization to exercise multiple geometric configurations.
- Assert that auto-learned substitutions propagate by inspecting `GeometryEngine.learned_rules()`.
- For CLI checks, run `subprocess.run(["bin/geometry", ...], env={"PYTHONPATH": "src"})` and compare stdout against curated baselines.
- Add regression JSON scripts for previously failing loci and run them inside tests to guard against symbolic regressions.

## Commit & Pull Request Guidelines
- Follow Conventional Commits (`feat:`, `fix:`, `docs:`, `chore:`). Keep the subject under 72 characters.
- Describe the "why" in the body, link related issues (`Refs #12`), and note test commands executed.
- PRs must include a brief summary, screenshots or sample output when behavior changes, and a checklist of updated docs/tests.
- Rebase onto the main branch before requesting review and ensure CI is green; avoid merge commits within feature branches.
