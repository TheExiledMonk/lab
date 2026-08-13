import json
from pathlib import Path
def test_next_test_is_frozen_staggered_order_audit():
 d=json.loads((Path('runs/dev225_local_handedness_representation_gate')/'dev226_test_selection.json').read_text())
 assert d['DEV226_TEST_SELECTION']=='STAGGERED_LOCAL_ORDER_AUDIT' and d['DEV226_TEST_SELECTION_FROZEN']
