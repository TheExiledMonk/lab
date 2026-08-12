import numpy as np
from pbuf.observer.relative_neighbor_motion import motion_state, directed_momentum_difference


def test_directed_relative_motion_and_momentum_are_exact():
    u = np.zeros((3, 3, 3, 3)); p = np.zeros_like(u); u1 = u.copy(); u1[1, 0, 0, 1] = .25
    w = motion_state(u, p, u1, p)
    assert w['delta_relation'].shape == (3, 3, 3, 6, 3)
    assert np.allclose(w['delta_strain'], np.linalg.norm(w['next_relation'], axis=-1) - 1)
    assert np.allclose(directed_momentum_difference(p), 0)
