import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/dev185_received_state_distribution_sufficiency'

def test_dev185_c100_distribution_contract():
    f=json.loads((OUT/'final_contract.json').read_text())
    assert f['DEV185_COMPLETE'] and f['C100_INPUT_HASHES_VERIFIED']
    assert f['ALL_EIGHT_REALIZATIONS_INCLUDED'] and f['ALL_121_LAUNCHES_PER_REALIZATION_INCLUDED']
    assert f['NATIVE_MODE_CHANNEL_GATE']=='AUTHORIZED'
    assert f['PHYSICAL_OBSERVER_GATE_BLOCKED_PENDING_MODE_CHANNEL']
    assert not f['NO_CHANNEL_SELECTION'] is False
