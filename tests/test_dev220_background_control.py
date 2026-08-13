import json
from pathlib import Path

def test_existing_unloaded_control_is_not_subtracted_or_thresholded():
    data=json.loads((Path(__file__).parents[1]/'runs/dev220_native_spatial_winding/unloaded_spatial_winding.json').read_text())
    assert data['UNLOADED_CONTROL_REUSED'] and data['UNLOADED_SPATIAL_WINDING'] == 'ZERO'
