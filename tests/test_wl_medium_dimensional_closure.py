from pbuf.wl.medium_dimensional_closure import DimensionalSystem, default_dimensional_system

def test_exact_rank_and_nullity():
    a=DimensionalSystem(("L0","U0"),((1,-1),(0,1))).audit()
    assert (a["matrix_rank"],a["nullity"])==(2,0)
    assert a["L0_identifiability"]=="UNIQUELY_IDENTIFIABLE_INTERNALLY"

def test_one_free_length_scale():
    a=DimensionalSystem(("L0",),()).audit()
    assert a["L0_identifiability"]=="IDENTIFIABLE_WITH_ONE_EXTERNAL_ANCHOR"

def test_response_length_codegeneracy():
    a=DimensionalSystem(("L0","U0"),((-1,1),)).audit()
    assert a["L0_identifiability"]=="CO_DEGENERATE_WITH_OTHER_UNIT_SCALE"
    assert a["nullspace_basis"]==[["1","1"]]

def test_current_system_is_not_closed():
    assert default_dimensional_system().audit()["L0_identifiability"]=="NON_IDENTIFIABLE_FROM_CURRENT_PHYSICS"
