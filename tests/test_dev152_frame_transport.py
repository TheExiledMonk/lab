import numpy as np
from pbuf.foundation.native_neighbor_mixed_state import construct_case
from pbuf.foundation.native_neighbor_frame_transport import transport_map

def test_f04_is_norm_preserving():
    f=construct_case(5,0)["frames"]
    for a,b in zip(f[:-1],f[1:]):
        q=transport_map(a,b,"F04"); assert np.allclose(q.T@q,np.eye(2))

