import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "runs/dev228_two_body_source_state_validity"
def test_inventory_has_all_source_classes():
    data=json.loads((OUT / "native_source_state_inventory.json").read_text())
    assert data["NATIVE_SOURCE_STATE_INVENTORY_COMPLETE"]
    assert [x["id"] for x in data["classes"]] == ["S1", "S2", "S3", "S4", "S5"]
