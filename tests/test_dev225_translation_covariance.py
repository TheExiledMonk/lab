import json
from pathlib import Path
def test_translation_contract_is_exact():
 d=json.loads((Path('runs/dev225_local_handedness_representation_gate')/'handedness_translation_covariance.json').read_text())
 assert d['LOCAL_HANDEDNESS_TRANSLATION_COVARIANCE']=='EXACT'
