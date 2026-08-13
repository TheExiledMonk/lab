import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_required_history_has_scoped_entries():
    data = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/historical_negative_scope_matrix.json").read_text())
    entries = {entry["dev"]: entry for entry in data["entries"]}
    for dev in ["DEV159", "DEV195", "DEV206", "DEV207", "DEV211", "DEV214", "DEV215", "DEV217", "DEV218", "DEV223", "DEV226"]:
        assert dev in entries
        assert entries[dev]["valid_closure_scope"]
        assert entries[dev]["prohibited_broader_claim"]
    assert entries["DEV226"]["result"] == "ALIGNED_DOMINANT"
    assert "two-body interstitial stress patterns" in entries["DEV226"]["prohibited_broader_claim"]
