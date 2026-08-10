from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case
from pbuf.foundation.native_neighbor_mixed_observer import observe

def test_no_unsupported_localization_claim():
    assert observe(progress_case(construct_case(2,5))["history"])["localization"]=="NO_LOCALIZATION"

