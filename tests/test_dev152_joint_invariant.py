from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case
from pbuf.foundation.native_neighbor_mixed_invariants import audit

def test_joint_invariant_survives_mixed_progression():
    rows=audit(progress_case(construct_case(3,4))["history"])
    assert "J01" in rows["conserved_candidates"]

