import json
from pathlib import Path
def test_exact_dev203_tensor_is_reused():
 d=json.loads((Path('runs/dev225_local_handedness_representation_gate')/'dev203_antisymmetric_tensor_definition.json').read_text())
 assert d['DEV203_ANTISYMMETRIC_TENSOR_REUSED'] and d['definition'].startswith('A_ij=(M_ij-M_ji)/2')
