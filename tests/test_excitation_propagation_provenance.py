from pbuf.foundation.excitation_propagation_provenance import candidate_manifest, provenance_contract, state_inventory

def test_inventory_schema_and_gate():
    rows=state_inventory(); assert rows and all(len(r)==22 for r in rows)
    assert len(candidate_manifest())==20
    c=provenance_contract(); assert c['trajectory_semantics']=='GEOMETRIC_TRACER_ONLY'
    assert c['new_dynamic_dof_required'] and not c['dynamic_excitation_state_found']

