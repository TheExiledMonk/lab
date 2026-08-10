import numpy as np
from pbuf.wl.native_scalar_transport_audit import reverse_path

def test_reverse_path_closes():
    assert reverse_path(np.array([1.2,.9,1.4,.7]))["pass"]
