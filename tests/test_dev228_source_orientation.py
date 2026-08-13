import json
from pathlib import Path
def test_orientation_is_dynamic_not_magnetic_poles():
    d=json.loads((Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity/source_orientation_state.json").read_text())
    assert d["SOURCE_ORIENTATION_STATE"] == "DYNAMIC_ONLY"
    assert d["NO_NORTH_SOUTH_PRIMITIVE"] and d["NO_MAGNETIC_POLE_LABEL"]
