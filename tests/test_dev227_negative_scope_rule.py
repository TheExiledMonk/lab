import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_frozen_condition_rule_preserves_negative_levels():
    data = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/frozen_condition_closure_rule.json").read_text())
    assert data["FROZEN_CONDITION_CLOSURE_RULE_ADDED"] is True
    assert data["formula"] == "R-=R-(M,P,G,B,O,Δx,Δt,T)"
    assert data["NO_RESULT_MOTIVATED_REOPENING"] is True
    assert data["negative_levels"] == ["NEGATIVE_OBSERVABLE", "NEGATIVE_MECHANISM_UNDER_FROZEN_CONDITIONS", "NEGATIVE_PHYSICAL_MECHANISM"]
