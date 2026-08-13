import json
from pathlib import Path
def test_gaussian_packet_is_not_canonically_compact_source():
    d=json.loads((Path(__file__).parents[1]/"runs/dev228_two_body_source_state_validity/native_source_localization.json").read_text())
    assert d["NATIVE_SOURCE_LOCALIZATION"] == "NONUNIQUE"
