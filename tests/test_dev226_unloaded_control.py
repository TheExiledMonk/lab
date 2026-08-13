import json
from pathlib import Path

def test_unloaded_control_is_zero():
    x=json.loads((Path(__file__).parents[1]/'runs/dev226_staggered_local_order/unloaded_staggered_local_order.json').read_text())
    assert x['UNLOADED_STAGGERED_LOCAL_ORDER'] == 'ZERO'
