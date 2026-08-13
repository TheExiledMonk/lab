import json
from pathlib import Path
def test_no_persistent_source_is_assumed():
    d=json.loads((Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity/native_source_persistence.json").read_text())
    assert d["NATIVE_SOURCE_PERSISTENCE"] == "NOT_DERIVED"
