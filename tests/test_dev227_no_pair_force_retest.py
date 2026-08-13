import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_audit_did_not_reopen_pair_force_or_run_pattern():
    final = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/final_contract.json").read_text())
    assert final["PAIR_ORIENTATION_INTERACTION_GATE"] == "REMAINS_BLOCKED"
    assert final["NO_REOPENING_DEV218_FORCE_SIGN"] is True
    assert final["NO_INTERSTITIAL_PATTERN_RESULT_RUN"] is True
    assert final["NO_DYNAMIC_INTERSTITIAL_TEST"] is True
