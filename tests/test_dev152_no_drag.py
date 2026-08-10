from pbuf.foundation.native_neighbor_mixed_state import construct_case, progress_case

def test_homogeneous_progression_has_no_norm_drag():
    assert progress_case(construct_case(0,6),steps=48)["relative_norm_drift"]<1e-12

