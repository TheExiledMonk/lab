import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
OUT = ROOT / "runs/dev222_dev221_reconciliation"


def test_frozen_numeric_replay_is_exact_and_materially_nonzero():
    data = json.loads((OUT / "odd_geometry_replay.json").read_text())
    assert data["stored_profile_matches_replay_exactly"]
    assert data["LONGITUDINAL_ODD_GEOMETRY_CONTENT"] == "PRESENT"
    assert min(data["O_q"]) > 1e-2
