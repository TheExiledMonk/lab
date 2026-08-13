import json
from pathlib import Path
def test_reflection_contract_is_classified():
 d=json.loads((Path('runs/dev223_pattern_boundary_interface')/'pattern_mismatch_reflection_covariance.json').read_text())
 assert d['PATTERN_MISMATCH_REFLECTION_COVARIANCE'] in {'EXACT','ROUND_OFF'}
