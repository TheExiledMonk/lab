import numpy as np
from pbuf.foundation.native_linear_modes import allowed_wavevectors, tangent_symbol, update_matrix


def test_unloaded_symbol_and_exact_k_grid():
    indices, ks = allowed_wavevectors((11, 11, 11))
    assert len(ks) == 11**3
    k = np.array([2*np.pi/11, 0., 0.])
    assert np.allclose(tangent_symbol(k), np.diag([-4*np.sin(np.pi/11)**2, 0, 0]))
    assert np.all(indices[:, 0] >= -5) and np.all(indices[:, 0] <= 5)


def test_kick_drift_eigenpairs_reconstruct_and_axial_transverse_is_zero():
    M = update_matrix(np.array([2*np.pi/11, 0., 0.]), .04)
    w, v = np.linalg.eig(M)
    assert np.max(np.abs(M @ v - v*w)) < 1e-12
    assert np.allclose(np.abs(w), 1.0, atol=1e-12)
    # Zero-force transverse displacements are fixed; their momenta generate
    # the expected free linear drift (a unit-eigenvalue Jordan block).
    assert np.allclose(M @ np.array([0.,1.,0.,0.,0.,0.]), np.array([0.,1.,0.,0.,0.,0.]))
    assert np.allclose(M[4:, 1:3], 0.0)
