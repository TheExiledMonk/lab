import json
from pathlib import Path
def test_bipartite_reference_is_geometric_not_a_state_flip():
 d=json.loads((Path('runs/dev224_magnetic_minimality_reassessment')/'bipartite_n6_geometry_contract.json').read_text())
 assert d['N6_BIPARTITE'] and d['NO_MANUAL_CHECKERBOARD_FLIP']
