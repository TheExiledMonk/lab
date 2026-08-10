from pbuf.foundation.excitation_propagation_provenance import provenance_contract

def test_dev146_requirement_is_confirmed_by_contract():
    c=provenance_contract()
    assert c['new_dynamic_dof_required'] is True
    assert c['raw_magnitude_exists'] and not c['raw_magnitude_physical']
    assert not c['normalization_discards_physical_state']
