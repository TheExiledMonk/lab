import json
from pathlib import Path

def test_final_contract_preserves_closed_routes():
    x=json.loads((Path(__file__).parents[1]/'runs/dev226_staggered_local_order/final_contract.json').read_text())
    assert x['NO_RL_STATE'] and x['NO_PAIR_INTERACTION'] and x['MAGNETIC_IDENTITY_NOT_DERIVED']
