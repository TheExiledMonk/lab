import json
from pathlib import Path
def test_next_test_is_representation_gate():
 d=json.loads((Path('runs/dev224_magnetic_minimality_reassessment')/'dev225_test_selection.json').read_text())
 assert d == {'DEV225_TEST_SELECTION':'LOCAL_HANDEDNESS_REPRESENTATION_GATE','DEV225_TEST_SELECTION_FROZEN':True}
