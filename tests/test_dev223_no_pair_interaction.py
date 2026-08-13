import json
from pathlib import Path
def test_pair_interaction_remains_blocked():
 d=json.loads((Path('runs/dev223_pattern_boundary_interface')/'final_contract.json').read_text())
 assert d['NO_PAIR_INTERACTION'] and d['PAIR_ORIENTATION_INTERACTION_GATE']=='REMAINS_BLOCKED'
