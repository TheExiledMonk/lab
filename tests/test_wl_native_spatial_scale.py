from pbuf.wl.native_spatial_scale import PhysicalScaleCandidate, classify_candidate, recover_scale

def test_known_physical_spacing_recovered_exactly():
    candidate=recover_scale(2.5,"kpc",[1.0])
    assert candidate.status=="AUTHORITATIVE" and candidate.value==2.5

def test_normalized_grid_scale_recovery():
    candidate=recover_scale(10.0,"m",[4.0])
    assert candidate.status=="AUTHORITATIVE" and candidate.value==2.5

def test_numerical_and_target_candidates_rejected():
    numerical=classify_candidate(PhysicalScaleCandidate("e","launch","extent=8",8,"native","configured",(),False,False,True))
    target=classify_candidate(PhysicalScaleCandidate("t","image","field match",1,"arcsec","image matching",(),True,True,True))
    assert numerical.status=="NUMERICAL_ONLY"
    assert target.status=="TARGET_CONTAMINATED"

def test_forbidden_scale_sources_rejected():
    for source in ("Rmax", "strength=0.18", "undeclared Planck length", "LCDM distance"):
        c=classify_candidate(PhysicalScaleCandidate("x","legacy",source,1,"m",source,(),False,True,True))
        assert c.status=="REJECTED" and c.rejection_reason.startswith("FORBIDDEN_")
