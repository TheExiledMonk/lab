import json
from pathlib import Path
def test_prior_closures_remain_preserved():
 d=json.loads((Path('runs/dev223_pattern_boundary_interface')/'final_contract.json').read_text())
 assert d['DEV220_SPATIAL_WINDING_CLOSURE_PRESERVED'] and d['DEV222_GEOMETRY_INTERACTION_TRANSFER_BLOCK_PRESERVED']
