from pbuf.wl.native_scalar_transport_audit import refinement_cv, uniform_identity

def test_uniform_identity_tolerance(): assert uniform_identity()["pass"]
def test_endpoint_ratio_refinement_stable():
    assert refinement_cv(lambda density: 1.25)["classification"] == "STRONG"
