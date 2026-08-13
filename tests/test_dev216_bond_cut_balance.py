import json
from pathlib import Path
RUN=Path(__file__).parents[1]/'runs/dev216_bond_cut_dynamic_polarity'
def test_closed_regions_balance_to_roundoff():
 c=json.loads((RUN/'region_momentum_balance.json').read_text()); assert c['REGION_MOMENTUM_BALANCE'] in ('EXACT','ROUND_OFF')
