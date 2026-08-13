import json
from pathlib import Path
R=Path(__file__).parents[1]/'runs/dev219_magnetic_mechanism_wide_net'
def test_gate():
 d=json.loads((R/'dev220_test_selection.json').read_text()); assert d['DEV220_TEST_SELECTION']=='SPATIAL_WINDING_AUDIT' and d['one_observable']
