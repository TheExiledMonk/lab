import json
from pathlib import Path

def test_dev141_known_depth_contract_is_compatible():
    p=Path("runs/wl_spatial_wave_emergent_time_closure001/final_source_reconstruction_contract.json")
    d=json.loads(p.read_text())
    assert d["known_depth_reconstruction_established"] and d["source_size_recoverable"] and d["source_layout_recoverable"]
    assert not d["physical_distance_required"] and not d["time_required"]
