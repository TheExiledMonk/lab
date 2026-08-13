import json
from pathlib import Path
def test_translation_covariance_is_classified():
 assert json.loads((Path('runs/dev223_pattern_boundary_interface')/'pattern_mismatch_translation_covariance.json').read_text())['PATTERN_MISMATCH_TRANSLATION_COVARIANCE'] in {'EXACT','ROUND_OFF'}
