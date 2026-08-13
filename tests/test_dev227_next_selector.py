import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_next_selector_prioritizes_source_state():
    data = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/dev228_test_selection.json").read_text())
    assert data == {"DEV228_TEST_SELECTION": "TWO_BODY_SOURCE_STATE_VALIDITY_GATE", "DEV228_TEST_SELECTION_FROZEN": True, "priority_reason": "MAGNET_LIKE_SOURCE_STATE_VALIDITY=NOT_DERIVED takes priority over pattern mapping."}
