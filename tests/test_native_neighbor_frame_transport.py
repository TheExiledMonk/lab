import numpy as np
from pbuf.foundation.native_neighbor_state import local_link_frame,frame_overlap
def test_overlap_transport_is_orthogonal():
    q=frame_overlap(local_link_frame([1,0,0]),local_link_frame([1,.2,.1])); assert np.allclose(q.T@q,np.eye(2))
