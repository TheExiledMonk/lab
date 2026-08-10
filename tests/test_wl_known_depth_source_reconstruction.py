import numpy as np
from pbuf.wl.native_source_controls import source_cloud
from pbuf.wl.native_source_reconstruction import source_size_metrics

def test_supplied_depth_affine_reverse_recovers_layout():
    truth=source_cloud("two_component_source",1,64); A=np.array([[1.1,.05],[.02,.9]]); shift=np.array([.3,-.2])
    received=truth@A.T+shift; recovered=(received-shift)@np.linalg.inv(A).T
    assert np.allclose(recovered,truth)
    assert np.isclose(source_size_metrics(recovered)["R_rms_native"],source_size_metrics(truth)["R_rms_native"])
