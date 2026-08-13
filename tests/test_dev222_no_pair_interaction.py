import json
from pathlib import Path


def test_dev222_did_not_run_pair_interaction():
    data = json.loads((Path(__file__).parents[1] / "runs/dev222_dev221_reconciliation/final_contract.json").read_text())
    assert data["NO_NEW_PAIR_INTERACTION"]
    assert data["PAIR_ORIENTATION_INTERACTION_GATE"] == "REMAINS_BLOCKED"
