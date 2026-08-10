from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case
from pbuf.foundation.native_neighbor_mixed_observer import observe

def test_no_unsupported_composite_claim():
    assert observe(progress_case(construct_case(6,7))["history"])["loaded_composite"] is False

