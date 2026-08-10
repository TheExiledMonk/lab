import numpy as np
import pytest

from pbuf.wl.spin2_basis_trace import (apply_spin2, rotation_matrix, spin2_matrix,
    synthetic_state, tensor_from_components, components_from_tensor)


@pytest.mark.parametrize("degrees", [0, 22.5, 45, 90])
def test_rotations(degrees):
    q = np.array([0.37, -0.81])
    got = spin2_matrix(rotation_matrix(degrees)) @ q
    phi = np.deg2rad(2 * degrees)
    expected = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]]) @ q
    assert np.allclose(got, expected, atol=1e-12)


def test_rotation_90_arbitrary():
    assert np.allclose(apply_spin2(.37, -.81, rotation_matrix(90)), (-.37, .81), atol=1e-12)


def test_rotation_180_spin2_equivalence():
    assert np.allclose(synthetic_state(17), synthetic_state(197), atol=1e-12)


@pytest.mark.parametrize("a,expected", [
    (np.diag([-1., 1.]), np.diag([1., -1.])),
    (np.diag([1., -1.]), np.diag([1., -1.])),
    (np.array([[0., 1.], [1., 0.]]), np.diag([-1., 1.])),
    (rotation_matrix(30) @ np.diag([-1., 1.]), None),
])
def test_reflections_swaps(a, expected):
    s = spin2_matrix(a)
    if expected is not None:
        assert np.allclose(s, expected, atol=1e-12)
    assert np.linalg.det(a) < 0
    assert np.allclose(np.linalg.inv(s) @ s, np.eye(2), atol=1e-12)


def test_tensor_vs_complex_representation():
    for theta in (0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5):
        q = synthetic_state(theta)
        for phi in (0, 22.5, 45, 90):
            a = rotation_matrix(phi)
            tensor = components_from_tensor(a @ tensor_from_components(*q) @ a.T)
            complex_result = (q[0] + 1j*q[1]) * np.exp(2j*np.deg2rad(phi))
            assert np.allclose(tensor, [complex_result.real, complex_result.imag], atol=1e-12)
