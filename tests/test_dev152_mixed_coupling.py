from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case

def test_mixed_progression_has_no_new_coefficient():
    assert progress_case(construct_case(8,0))["new_interaction_coefficients"]==0

