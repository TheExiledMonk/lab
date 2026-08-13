import json
from pathlib import Path
def test_tensor_does_not_authorize_scalar_handedness():
 d=json.loads((Path('runs/dev224_magnetic_minimality_reassessment')/'native_local_handedness_gate.json').read_text())
 assert d['NATIVE_LOCAL_HANDEDNESS_GATE'] == 'BLOCKED_NONUNIQUE'
