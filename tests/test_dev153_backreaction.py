import numpy as np
from pbuf.foundation.native_link_response_invariants import audit
def test_no_free_wave_persistent_loading():
    x=np.ones((8,2)); l=np.zeros(8); result=audit(x,x,l,l); assert result["backreaction"]=="NO_BACKREACTION" and result["J05"]["established"]
