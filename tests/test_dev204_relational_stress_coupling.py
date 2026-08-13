import json
from pathlib import Path
import numpy as np

from pbuf.observer.relational_stress_coupling import force_response_tensors


def test_force_response_reuses_dev203_native_tensor_definition():
    d = np.zeros((1, 1, 1, 6, 3)); d[..., 0, 1] = 2; d[..., 1, 1] = -2
    s, a = force_response_tensors(d)
    assert s.shape == a.shape == (1, 1, 1, 3, 3)


def test_dev204_artifacts_preserve_native_boundaries():
    out = Path(__file__).resolve().parents[1] / 'runs/dev204_relational_stress_coupling'
    contract = json.loads((out / 'final_contract.json').read_text())
    assert contract['NO_ADDED_FIELD'] and contract['NO_ADDED_DOF']
    assert contract['FINITE_STEP_CROSS_TERM_PRESERVED']
    assert contract['ELECTRIC_FIELD_MAPPING'] == 'UNESTABLISHED'
    assert contract['MAGNETIC_FIELD_MAPPING'] == 'UNESTABLISHED'
    arrays = np.load(out / 'bond_force_change.npz')
    assert np.allclose(arrays['delta_force'], arrays['magnitude'] + arrays['orientation'] + arrays['cross'])
