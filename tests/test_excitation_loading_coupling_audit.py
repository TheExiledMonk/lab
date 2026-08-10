from pbuf.foundation.excitation_loading_coupling_audit import coupling_audit, loading_contract, progression_comparison

def test_no_coefficient_free_coupling_is_claimed():
    assert len(coupling_audit()['families'])==12
    assert not progression_comparison()['native_beta_measurable']
    c=loading_contract(); assert not c['existing_native_coupling_found'] and not c['zero_mass_propagation_changed']

