import json
from pathlib import Path
def test_finite_x_is_aggregate_preparation_not_pairwise_reduction():
    d=json.loads((Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity/finite_x_composition_audit.json").read_text())
    assert d["FINITE_X_NATIVE_COMPOSITION_GATE"] == "DERIVED"
    assert d["NO_PAIRWISE_REDUCIBILITY_ASSUMPTION"]
