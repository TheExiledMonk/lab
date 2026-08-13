import json
from pathlib import Path
def test_no_pair_route():
 d=json.loads((Path(__file__).parents[1]/'runs/dev221_extended_relational_geometry/final_contract.json').read_text()); assert d['NO_PAIR_INTERACTION_TEST'] and d['NO_NEW_FORCE']
