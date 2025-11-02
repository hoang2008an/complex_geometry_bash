from geometry_engine import GeometryEngine


def test_point_summary_uses_point_labels() -> None:
    engine = GeometryEngine()
    engine.add_point("A")
    engine.add_point("B")
    engine.add_point("M")
    engine.midpoint("A", "B", "M")

    summary = engine.point_summary(["M"], style="latex")
    value = summary["M"]

    assert "z_" not in value
    assert "zb_" not in value
    assert "A" in value and "B" in value


def test_constraint_strings_are_sanitized() -> None:
    engine = GeometryEngine()
    engine.add_point("A")
    engine.add_unit_circle("A")

    constraints = engine.constraint_strings(style="text")
    assert constraints, "Expected at least one constraint after adding A on the unit circle."

    for constraint in constraints:
        assert "z_" not in constraint
        assert "zb_" not in constraint
        assert "conjugate(A)" in constraint or "A" in constraint


def test_learned_rules_display_without_prefixes() -> None:
    engine = GeometryEngine()
    engine.add_point("A")
    engine.add_point("B")

    engine.learned_subs[engine.zb_symbol("A")] = engine.z_symbol("B")

    latex_rules = engine.display_learned_rules(style="latex")
    text_rules = engine.display_learned_rules(style="text")

    assert latex_rules == {"\\overline{A}": "B"}
    assert text_rules == {"conjugate(A)": "B"}
