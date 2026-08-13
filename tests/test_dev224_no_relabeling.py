import json
from pathlib import Path
def test_no_new_sign_state_or_failed_route_relabeling():
 d=json.loads((Path('runs/dev224_magnetic_minimality_reassessment')/'final_contract.json').read_text())
 assert d['NO_MANUAL_PLUS_MINUS_STATE'] and d['NO_MANUAL_CHECKERBOARD_FLIP'] and d['NO_RELABELED_FAILED_MECHANISM']
