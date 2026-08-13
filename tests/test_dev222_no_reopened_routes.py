import json
from pathlib import Path


def test_closed_routes_remain_closed():
    data = json.loads((Path(__file__).parents[1] / "runs/dev222_dev221_reconciliation/final_contract.json").read_text())
    for key in ("DEV211_STATIC_ROUTE_PRESERVED", "DEV215_TEMPORAL_CYCLE_CLOSURE_PRESERVED", "DEV218_MOMENTUM_POLARITY_CLOSURE_PRESERVED", "DEV220_SPATIAL_WINDING_CLOSURE_PRESERVED"):
        assert data[key]
