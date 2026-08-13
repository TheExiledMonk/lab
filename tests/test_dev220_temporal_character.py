import json
from pathlib import Path

def test_full_frozen_trajectory_is_used_without_cycle_interpretation():
    p=Path(__file__).parents[1]/'runs/dev220_native_spatial_winding'
    assert json.loads((p/'dev203_trajectory_contract.json').read_text())['FULL_FROZEN_TRAJECTORY_USED']
    assert json.loads((p/'spatial_winding_temporal_character.json').read_text())['NO_CYCLE_RATE_DERIVED']
