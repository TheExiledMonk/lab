import json
from pathlib import Path
def test_global_winding_and_staggered_order_are_distinct():
 d=json.loads((Path('runs/dev224_magnetic_minimality_reassessment')/'global_winding_vs_staggered_order.json').read_text())
 assert d['GLOBAL_WINDING_VS_STAGGERED_ORDER_DISTINGUISHED'] and d['DEV220_ABSENCE_COMPATIBLE_WITH_STAGGERED_ORDER'] == 'YES'
