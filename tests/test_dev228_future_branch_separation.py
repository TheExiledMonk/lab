import json
from pathlib import Path
ROOT=Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity"
def test_collective_and_composite_lanes_remain_separate():
    collective=json.loads((ROOT/"collective_x_body_candidate.json").read_text())
    composite=json.loads((ROOT/"three_constituent_composite_candidate.json").read_text())
    assert collective["target_id"] != composite["target_id"]
    assert composite["NO_QUARK_IDENTIFICATION"] and composite["NO_MAJORITY_STATE_RULE"]
