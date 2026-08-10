import numpy as np
from pbuf.foundation.native_transverse_transfer_capacity import orthogonal_transport,compose,ordering_audit
def test_identity_capacity_commutes_and_conserves():
    R=orthogonal_transport(.3); x=np.array([.4,.7]); y=compose(np.eye(2),R,x)
    assert np.isclose(y@y,x@x) and ordering_audit(np.eye(2),R)["classification"]=="COMMUTING_EQUIVALENT"
