import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_source_and_passive_material_are_not_assumed():
    source = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/magnet_like_source_state_validity.json").read_text())
    iron = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/magnet_iron_representation_status.json").read_text())
    assert source["MAGNET_LIKE_SOURCE_STATE_VALIDITY"] == "NOT_DERIVED"
    assert iron["MAGNET_IRON_REPRESENTATION_STATUS"] == "NOT_DERIVED"
    assert iron["NO_FAKE_IRON_AS_SECOND_MAGNET"] is True
