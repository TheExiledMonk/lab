from pbuf.foundation.native_neighbor_mixed_state import construct_case
from pbuf.foundation.native_neighbor_frame_transport import audit_frame_candidates

def test_orthogonal_transports_covariant():
    rows={r["candidate"]:r for r in audit_frame_candidates(construct_case(2,0)["frames"])}
    assert not rows["F01"]["basis_covariant"]
    assert all(rows[f"F0{i}"]["basis_covariant"] for i in range(2,7))

