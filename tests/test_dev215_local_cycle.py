import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; RUN=ROOT/'runs/dev215_lattice_state_cycle'
def test_full_descriptor_and_threshold_free_nonclosure():
    c=json.loads((RUN/'final_contract.json').read_text()); z=np.load(RUN/'local_state_histories.npz')
    assert c['LOCAL_NATIVE_STATE_DESCRIPTOR_DEFINED'] and c['NO_PHASE_ASSUMED']
    assert z['relation'].shape[-2:] == (6,3) and z['force'].shape[-2:] == (6,3)
    assert c['LOCAL_NATIVE_STATE_CYCLE'] == 'OSCILLATORY_NONCLOSED'
