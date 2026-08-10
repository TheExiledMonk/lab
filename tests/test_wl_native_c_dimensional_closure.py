from pbuf.wl.medium_dimensional_closure import default_dimensional_system

def test_c_relation_is_unavailable_without_native_time():
    a=default_dimensional_system().audit()
    assert a["nullity"] == 5 and a["L0_identifiability"] == "NON_IDENTIFIABLE_FROM_CURRENT_PHYSICS"
