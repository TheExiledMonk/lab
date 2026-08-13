import json
from pathlib import Path
def test_no_interstitial_or_dynamic_pattern_was_run():
    d=json.loads((Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity/final_contract.json").read_text())
    assert d["NO_INTERSTITIAL_PATTERN_RESULT_RUN"] and d["NO_DYNAMIC_PATTERN_RESULT_RUN"]
