import numpy as np
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState
from pbuf.observer.native_internal_state_inventory import reverse_momentum, state_summary


def test_momentum_reversal_preserves_shape_and_energy():
    u = np.zeros((3, 3, 3, 3)); p = np.ones_like(u) * 0.01
    a = VectorPairState(u, p); b = reverse_momentum(a)
    assert np.array_equal(a.displacement, b.displacement)
    assert np.array_equal(a.momentum, -b.momentum)
    assert state_summary(a)["total_energy"] == state_summary(b)["total_energy"]
