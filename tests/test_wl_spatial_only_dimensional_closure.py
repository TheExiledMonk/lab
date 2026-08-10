from pbuf.wl.medium_dimensional_closure import spatial_only_native_basis,dev137_time_ontology_reconciliation

def test_spatial_basis_has_no_time():
    x=spatial_only_native_basis(); assert x["base_dimensions"]==["M","L"] and "T0" not in x["unknowns"]
def test_t0_is_ontology_artifact():
    x=dev137_time_ontology_reconciliation(); assert x["T0_degeneracy_classification"]=="ONTOLOGY_ARTIFACT_REMOVED"
