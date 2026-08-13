import numpy as np

from pbuf.observer.native_stress_relay import directed_forces, finite_components


def test_directed_force_uses_only_reciprocal_dev167_bonds():
    u = np.zeros((3, 3, 3, 3))
    f = directed_forces(u)
    assert f.shape == (3, 3, 3, 6, 3)
    assert np.array_equal(f, np.zeros_like(f))


def test_dev204_finite_split_is_retained_by_relay_observer():
    r0 = np.broadcast_to(np.array([1.0, 0.0, 0.0]), (2, 6, 3)).copy()
    r1 = r0.copy(); r1[..., 0] += 0.01
    parts = finite_components(r0, r1)
    assert np.allclose(parts["delta_force"], parts["magnitude"] + parts["orientation"] + parts["cross"])
