import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_lab_failure_is_inconclusive_not_absence():
    data = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/two_body_lab_representation_gate.json").read_text())
    assert data["TWO_BODY_LAB_REPRESENTATION_GATE"] == "BLOCKED_SOURCE_STATE"
    assert data["negative_pattern_result_if_run"] == "INCONCLUSIVE_LAB_REPRESENTATION"
