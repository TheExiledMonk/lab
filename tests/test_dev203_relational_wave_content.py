import numpy as np
import json
from pathlib import Path
from pbuf.observer.relational_wave_content import decompose, native_symmetric_antisymmetric, opposite_components


def test_transverse_and_native_tensor_decomposition():
    v = np.array([[[[1., 2., 3.]]]])
    para, perp = decompose(v, np.array([1., 0., 0.]))
    assert np.allclose(para[..., 0], 1) and np.allclose(perp[..., 0], 0)
    d = np.zeros((1, 1, 1, 6, 3)); d[..., 0, 1] = 1; d[..., 3, 0] = 2
    s, a = native_symmetric_antisymmetric(d); pair_s, pair_a = opposite_components(d)
    assert s.shape == a.shape == (1, 1, 1, 3, 3)
    assert pair_s.shape == pair_a.shape == (1, 1, 1, 3, 3)


def test_dev203_frozen_artifact_contract():
    root = Path(__file__).resolve().parents[1]
    out = root / 'runs/dev203_relational_wave'
    contract = json.loads((out / 'final_contract.json').read_text())
    assert contract['NO_TANGENT_PROBE'] and contract['NO_ADDED_PERTURBATION']
    assert contract['UNLOADED_RELATIONAL_MOTION'] == 'ZERO_OR_MACHINE_FLOOR'
    assert contract['PROPAGATION_AS_RELATIONAL_MOTION'] == 'DERIVED'
    assert (out / 'node_n6_motion_pattern.npz').exists()
