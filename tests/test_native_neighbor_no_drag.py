from pbuf.foundation.native_neighbor_loaded_excitation import run_case
def test_homogeneous_progression_has_no_norm_drag():
    r=run_case(0,0,steps=31); assert abs(r['norm_after']-r['norm_before']) < 1e-12
