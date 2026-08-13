import json
from pathlib import Path
R=Path(__file__).parents[1]/'runs/dev219_magnetic_mechanism_wide_net'
def test_matrix():
 d=json.loads((R/'magnetic_candidate_family_matrix.json').read_text()); assert len(d['families'])==10
 for x in d['families']: assert set(x['scores']) >= {'REQUIRES_NEW_LAW','REQUIRES_NEW_DOF','CONSISTENT_WITH_DEV218'}
def test_external_inventory(): assert len(json.loads((R/'magnetic_phenomenology_matrix.json').read_text())['rows'])>=31
