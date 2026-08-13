import json
from pathlib import Path
def test_order_values_were_not_inspected():
 d=json.loads((Path('runs/dev225_local_handedness_representation_gate')/'native_local_handedness_relation_gate.json').read_text())
 assert d['NO_STAGGERED_ORDER_RESULT_INSPECTION'] and d['NO_RESULT_SELECTED_TIMESTEP']
