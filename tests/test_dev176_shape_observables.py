import numpy as np

from pbuf.direct_shape_observables import (
    local_deformation_tensor, project_to_screen, quadrupole_tensor,
    spin2_from_tensor, spin2_rotate, weighted_second_moment_tensor,
)


def test_tensor_symmetry_trace_normalization_and_zero_support():
    tensor, _, _ = weighted_second_moment_tensor([[0, 0], [2, 0], [0, 1]])
    assert np.allclose(tensor, tensor.T)
    e = spin2_from_tensor(tensor)
    assert np.isclose(np.hypot(*e), np.hypot(tensor[0, 0]-tensor[1, 1], 2*tensor[0, 1])/np.trace(tensor))
    empty, _, n = weighted_second_moment_tensor([], [])
    assert n == 0 and np.isnan(empty).all() and np.isnan(spin2_from_tensor(np.zeros((2, 2)))).all()


def test_translation_invariance_and_spin2_rotation_covariance():
    p = np.array([[0., 0.], [2., 0.], [0., 1.], [3., 2.]])
    a, _, _ = weighted_second_moment_tensor(p)
    b, _, _ = weighted_second_moment_tensor(p + [17., -9.])
    np.testing.assert_allclose(a, b)
    phi = .37
    rot = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    rr, _, _ = weighted_second_moment_tensor(p @ rot.T)
    # Rotating vectors is an active transform; DEV176's published convention
    # is passive screen-basis rotation, hence the opposite sign here.
    np.testing.assert_allclose(spin2_from_tensor(rr), spin2_rotate(spin2_from_tensor(a), -phi), atol=1e-12)


def test_projection_and_deformation_are_deterministic():
    xyz = np.array([[1., 2., 3.], [4., 5., 6.]])
    np.testing.assert_allclose(project_to_screen(xyz, [1, 0, 0], [0, 1, 0]), [[2, 3], [5, 6]])
    source = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]])
    received = source @ np.array([[1.2, .1], [.3, .9]]).T + [5, 7]
    j, stf = local_deformation_tensor(source, received)
    np.testing.assert_allclose(j, [[1.2, .1], [.3, .9]])
    assert np.isclose(np.trace(stf), 0)
    assert local_deformation_tensor(source[:2], received[:2]) is None


def test_direction_tensor_is_symmetric():
    tensor, n = quadrupole_tensor([[1, 0], [0, 1], [1, 1]])
    assert n == 3 and np.allclose(tensor, tensor.T)
