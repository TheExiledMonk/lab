import json
from pathlib import Path
def test_no_scalar_norm_or_threshold_selects_boundary():
 r=Path('runs/dev223_pattern_boundary_interface'); a=json.loads((r/'pattern_mismatch_representation.json').read_text()); b=json.loads((r/'native_relational_pattern_boundary.json').read_text())
 assert a['SCALAR_MISMATCH_REPRESENTATION']=='NONUNIQUE' and b['NO_BOUNDARY_THRESHOLD']
