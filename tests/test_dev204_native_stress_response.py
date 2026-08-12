import numpy as np

from pbuf.observer.native_stress_response import finite_step_force_response


def test_exact_finite_step_force_split_includes_cross_term():
    e0 = np.array([0.1]); e1 = np.array([0.2])
    u0 = np.array([[1.0, 0.0, 0.0]]); u1 = np.array([[0.0, 1.0, 0.0]])
    out = finite_step_force_response(e0, e1, u0, u1)
    s0 = e0 / (1 - e0**2); s1 = e1 / (1 - e1**2)
    assert np.allclose(out['delta_force'], s1[:, None] * u1 - s0[:, None] * u0)
    assert not np.allclose(out['cross'], 0.0)
