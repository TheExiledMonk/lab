import numpy as np
from pbuf.foundation.native_neighbor_invariants import basis_invariance
def test_transverse_norm_is_basis_invariant(): assert basis_invariance(np.arange(24).reshape(12,2))['status']=='ESTABLISHED'
