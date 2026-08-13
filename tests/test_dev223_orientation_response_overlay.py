import json
from pathlib import Path
def test_overlay_follows_frozen_boundary_result():
 d=json.loads((Path('runs/dev223_pattern_boundary_interface')/'pattern_boundary_orientation_response_relation.json').read_text())
 assert d['PRIMARY_BOUNDARY_RESULT_FROZEN_BEFORE_MECHANICAL_OVERLAY'] is True
