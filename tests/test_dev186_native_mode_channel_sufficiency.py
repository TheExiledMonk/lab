import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/dev186_native_mode_channel_sufficiency'
def test_dev186_exact_content_sufficiency():
 f=json.loads((OUT/'final_contract.json').read_text()); c=json.loads((OUT/'native_channel_sufficiency.json').read_text())
 assert f['NATIVE_CHANNEL_SUFFICIENCY']=='REDUCED_EXACT_SUFFICIENT_SET_FOUND'
 assert c['removed_storage_redundancy']==['W01','W03','W04']
 assert json.loads((OUT/'physical_observer_derivation_gate.json').read_text())['PHYSICAL_OBSERVER_DERIVATION_GATE']=='AUTHORIZED_EXACT_REDUCED_STATE'
