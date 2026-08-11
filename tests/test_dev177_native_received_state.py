import numpy as np

from pbuf.labs.foundation.native_channel_information_geometry_dev177 import information_geometry
from pbuf.labs.foundation.native_received_j3_dev177 import fit_j3


def test_information_geometry_preserves_missing_support_without_imputation():
    metric = information_geometry(np.array([[0., 1.], [np.nan, 3.], [2., 5.]]))
    assert metric["n_rows"] == 2
    assert metric["missing_rows_dropped_not_imputed"] is True


def test_j3_affine_fixture_and_intrinsic_metric_identity():
    source = np.array([[0., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0., 1., 1.]])
    expected = np.array([[2., 3.], [1., -1.], [.5, 2.]])
    received = source[:, 1:] @ expected.T + [4., -2., 7.]
    fitted = fit_j3(source, received)
    assert fitted["J3_STATUS"] == "DEFINED"
    np.testing.assert_allclose(fitted["J3"], expected)
    np.testing.assert_allclose(fitted["G3"], fitted["J3"].T @ fitted["J3"])


def test_j3_refuses_rank_deficient_source_support():
    source = np.array([[0., 0., 0.], [0., 1., 0.], [0., 2., 0.]])
    assert fit_j3(source, source)["J3_STATUS"] == "UNDEFINED_INSUFFICIENT_SOURCE_SUPPORT"
