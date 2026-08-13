import json
from pathlib import Path


def test_canonical_n6_signature_mismatch_authorizes_only_representation_audit():
    root = Path(__file__).parents[1] / "runs/dev222_dev221_reconciliation"
    candidate = json.loads((root / "native_pattern_mismatch_observable_candidate.json").read_text())
    gate = json.loads((root / "pattern_boundary_audit_gate.json").read_text())
    assert candidate["NATIVE_PATTERN_MISMATCH_OBSERVABLE_CANDIDATE"] == "DERIVABLE"
    assert candidate["NO_NORM_SELECTED"] and candidate["NO_BOUNDARY_THRESHOLD"]
    assert gate["PATTERN_BOUNDARY_AUDIT_GATE"] == "AUTHORIZED"
