import numpy as np
from pbuf.foundation.native_neighbor_state import NativeNeighborState
from pbuf.foundation.native_neighbor_invariants import joint_invariant_audit
def test_all_joint_families_audited():
    r=joint_invariant_audit(NativeNeighborState(np.zeros(4),np.ones((4,2)))); assert {f'J{i:02d}' for i in range(1,11)}==set(r)
