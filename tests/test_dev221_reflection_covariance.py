import json
from pathlib import Path
def test_reflection_contract():
 d=json.loads((Path(__file__).parents[1]/'runs/dev221_extended_relational_geometry/extended_geometry_reflection_covariance.json').read_text()); assert d['EXTENDED_GEOMETRY_REFLECTION_COVARIANCE'] in {'EXACT','ROUND_OFF'}
