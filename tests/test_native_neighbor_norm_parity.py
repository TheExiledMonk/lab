from pbuf.foundation.native_neighbor_loaded_excitation import run_case
def test_loaded_frame_transport_preserves_norm():
    r=run_case(5,2); assert abs(r['norm_after']-r['norm_before']) < 1e-10
