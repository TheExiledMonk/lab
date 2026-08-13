import json
from pathlib import Path
def test_no_component_fishing():
 d=json.loads((Path('runs/dev225_local_handedness_representation_gate')/'least_reduced_representation_selection.json').read_text())
 assert all(d[k] for k in ('NO_SELECTED_TENSOR_COMPONENT','NO_COMPONENT_MAGNITUDE_RANKING','NO_COMPONENT_SWEEP'))
