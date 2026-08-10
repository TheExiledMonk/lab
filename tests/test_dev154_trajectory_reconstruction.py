from pbuf.audit.native_system_inventory import domains
def test_trajectory_is_centroid_diagnostic_and_direction_absent():
    d={x["id"]:x for x in domains()}; assert r"\Gamma=(\bar r_0" in d["D26"]["Mathematical definition"]
    assert d["D27"]["Status"]=="UNRESOLVED" and d["D28"]["Status"]=="UNRESOLVED"
