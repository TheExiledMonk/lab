import json
from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_dev221_canonical_records_agree():
    contract = json.loads((ROOT / "runs/dev222_dev221_reconciliation/final_contract.json").read_text())
    registry = json.loads((ROOT / "docs/PBUF_MECHANISM_REGISTRY.json").read_text())
    handoff = (ROOT / "runs/dev221_extended_relational_geometry/discussion_handoff.md").read_text()
    target = next(item for item in registry["targets"] if item["target_id"] == "native_extended_directional_geometry")
    assert contract["NATIVE_EXTENDED_DIRECTIONAL_GEOMETRY"] == "DERIVED_STRONG"
    assert target["current_status"] == "CANONICAL"
    assert "nonzero longitudinal odd component" in handoff
    assert "end-for-end symmetric" not in (ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md").read_text()
