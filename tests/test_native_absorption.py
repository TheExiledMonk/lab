from pbuf.quantum.native_emission_absorption import absorption_audit
def test_absorption_waits_for_interaction_law():
    assert absorption_audit()['status']=='MISSING_INTERACTION_LAW'
