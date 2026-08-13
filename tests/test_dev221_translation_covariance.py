import json
from pathlib import Path
def test_translation_contract():
 d=json.loads((Path(__file__).parents[1]/'runs/dev221_extended_relational_geometry/extended_geometry_translation_covariance.json').read_text()); assert d['EXTENDED_GEOMETRY_TRANSLATION_COVARIANCE'] in {'EXACT','ROUND_OFF'}
