import json
from pathlib import Path
def test_unloaded_control(): assert json.loads((Path(__file__).parents[1]/'runs/dev221_extended_relational_geometry/unloaded_extended_geometry.json').read_text())['q_zero']
