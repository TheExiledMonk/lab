import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/dev187_physical_native_observer'
def test_dev187_transfer_detector_contract():
 f=json.loads((OUT/'final_contract.json').read_text());t=json.loads((OUT/'transfer_function_status.json').read_text())
 assert f['OBSERVATIONAL_COMPARISON_GATE']=='CLOSED'
 assert t['NATIVE_TRANSFER_OPERATOR_DERIVED'] and t['DIRECT_INTENSITY_NOT_YET_DERIVED']
 assert json.loads((OUT/'native_shape_tensor_gate.json').read_text())['NATIVE_SHAPE_TENSOR_GATE']=='BLOCKED_TRANSFER_FUNCTION_SEMANTICS'
