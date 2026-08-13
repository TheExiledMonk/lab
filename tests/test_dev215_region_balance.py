import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; RUN=ROOT/'runs/dev215_lattice_state_cycle'
def test_exact_bond_cut_closes_region_momentum():
    c=json.loads((RUN/'closed_region_momentum_balance.json').read_text()); z=np.load(RUN/'native_bond_cut_accounting.npz')
    assert c['CLOSED_REGION_MOMENTUM_BALANCE'] in ('EXACT','ROUND_OFF')
    assert np.max(np.abs(z['residual'])) < 1e-12
