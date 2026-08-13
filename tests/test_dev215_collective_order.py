import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; RUN=ROOT/'runs/dev215_lattice_state_cycle'
def test_fixed_native_region_and_no_phase_mapping():
    r=json.loads((RUN/'lattice_region_contract.json').read_text()); p=json.loads((RUN/'native_phase_definition.json').read_text()); z=np.load(RUN/'neighbor_state_cycle_correlation.npz')
    assert r['node_count'] == 25 and r['NO_REGION_SIZE_SWEEP_BEFORE_PRIMARY_CLASSIFICATION']
    assert p['NATIVE_PHASE_DEFINITION'] == 'NOT_DERIVED' and len(z['pairs']) > 0
