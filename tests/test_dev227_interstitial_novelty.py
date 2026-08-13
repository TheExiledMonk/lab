import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_interstitial_route_is_distinct_from_force_and_single_body_audits():
    candidate = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/two_body_interstitial_candidate_definition.json").read_text())
    equivalence = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion/two_body_prior_test_equivalence.json").read_text())
    assert candidate["candidate"] == "TWO_BODY_INTERSTITIAL_RELATIONAL_STRESS_PATTERN"
    assert equivalence["TWO_BODY_INTERSTITIAL_PATTERN_PRIOR_TEST_EQUIVALENCE"] == "PARTIAL_OVERLAP_DEV217_218"
    assert equivalence["INTERSTITIAL_PATTERN_EQUIVALENT_TO_DEV221"] is False
    for name in ["dev218_force_observer_scope.json", "dev223_single_structure_scope.json"]:
        payload = json.loads((ROOT / "runs/dev227_magnetic_candidate_exhaustion" / name).read_text())
        assert False in payload.values()
