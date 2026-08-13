import json
from pathlib import Path
def test_required_history_is_recorded():
 d=json.loads((Path('runs/dev224_magnetic_minimality_reassessment')/'final_contract.json').read_text())
 assert all(d[f'{x}_READ'] for x in ('DEV203','DEV204','DEV206','DEV212','DEV214','DEV215','DEV218','DEV219','DEV220','DEV221','DEV222','DEV223'))
