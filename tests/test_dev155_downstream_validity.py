from pathlib import Path
import json
def test_downstream_matrix_reopens_transport_and_dev153():
    p=Path(__file__).parents[1]/"runs/n6_native_excitation_restoration001/downstream_validity_matrix.json"
    if p.exists():
      d=json.loads(p.read_text()); assert d["DEV148_STATE_RESULT_SURVIVES"] and d["DEV148_TRANSPORT_RESULT_REOPENED"] and d["DEV153_CROSS_COUPLING_NULL_REOPENED"]
