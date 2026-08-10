from pbuf.wl.native_spatial_scale import scale_invariance_control

def test_dimensionless_synthetic_structure_is_scale_invariant():
    audit=scale_invariance_control()
    assert audit["passed"]
    assert audit["outcome"]=="GLOBAL_SCALE_DEGENERACY_SUPPORTED"
    assert audit["propagation_physics_modified"] is False
