import json
from pathlib import Path
def test_dev213_is_reused_without_linear_superposition():
    d=json.loads((Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity/two_body_source_composition.json").read_text())
    assert d["DEV213_MULTI_STRUCTURE_COMPOSITION_REUSED"] and d["NO_POSTHOC_LINEAR_SUPERPOSITION"]
    assert d["TWO_BODY_SOURCE_COMPOSITION"] == "BLOCKED_PREPARATION"
