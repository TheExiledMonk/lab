import numpy as np
from pbuf.foundation.loaded_native_tangent import bond_tangent, tangent_net_force


def test_loaded_bond_contains_geometric_term_and_unloaded_does_not():
    u = np.zeros((3, 3, 3, 3)); u[1:, :, :, 0] = .01
    b = bond_tangent(u)
    assert np.max(np.abs(b['k_perp'])) > 0
    z = bond_tangent(np.zeros_like(u))
    assert np.allclose(z['k_perp'], 0)
    d = np.zeros_like(u); d[..., 1] = 1
    assert tangent_net_force(u, d).shape == u.shape
