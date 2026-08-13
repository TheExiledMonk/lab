import json
from pathlib import Path

def test_closed_routes_and_nonmagnetic_boundary_remain_explicit():
    data=json.loads((Path(__file__).parents[1]/'runs/dev220_native_spatial_winding/final_contract.json').read_text())
    assert data['DEV215_TEMPORAL_CYCLE_RESULT_PRESERVED'] and data['DEV218_MOMENTUM_POLARITY_CLOSURE_PRESERVED']
    assert data['NO_FORCE_TEST'] and data['NO_B_FIELD']
