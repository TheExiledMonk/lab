import json
from pathlib import Path
def test_next_selector_returns_to_source_derivation():
    d=json.loads((Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity/dev229_test_selection.json").read_text())
    assert d["DEV229_TEST_SELECTION"] == "PERSISTENT_NATIVE_SOURCE_DERIVATION_GATE"
    assert d["DEV229_TEST_SELECTION_FROZEN"]
