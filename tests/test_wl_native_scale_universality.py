from pbuf.wl.native_scale_universality import resolution_stability,stability

def test_resolution_artifact_rejected():
    n=[32,48,64,96,128]
    a=resolution_stability(n,n)
    assert a["rejection_reason"]=="REJECT_GRID_RESOLUTION_ARTIFACT"

def test_universal_candidate_stable():
    assert stability([2,2,2,2,2],[.25,.5,1,2,4])["classification"]=="STABLE"

def test_mass_dependent_false_ruler_fails():
    m=[.25,.5,1,2,4]
    assert stability([x**.5 for x in m],m)["classification"]=="FAIL"
