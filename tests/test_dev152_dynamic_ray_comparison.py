from pbuf.foundation.native_neighbor_mixed_observer import ray_comparison

def test_ray_comparison_exact_control():
    assert ray_comparison([0,1,2],[0,1,2])["status"]=="RAY_MORPHOLOGY_MATCH"

