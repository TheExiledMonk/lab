from pbuf.observer.native_multistate_interaction import blocked_matrices


def test_closed_gate_has_no_synthetic_pair_matrix():
    m = blocked_matrices()
    assert m["radial_force"].shape == (0, 0)
    assert m["torque"].shape == (0, 0, 3)
