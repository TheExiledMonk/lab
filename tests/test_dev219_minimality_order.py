import json
from pathlib import Path
R=Path(__file__).parents[1]/'runs/dev219_magnetic_mechanism_wide_net'
def test_nonunique_pareto():
 p=json.loads((R/'magnetic_mechanism_pareto_set.json').read_text()); r=json.loads((R/'magnetic_mechanism_minimality_ranking.json').read_text()); assert len(p['MAGNETIC_MECHANISM_PARETO_SET'])>1 and r['MINIMAL_NATIVE_MAGNETIC_CANDIDATE']=='NONUNIQUE'
