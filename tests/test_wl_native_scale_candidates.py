from pbuf.wl.native_scale_candidates import NativeScaleEstimate,candidate_registry,execute_internal_candidates,validate_candidate
from pbuf.wl.medium_dimensional_closure import default_dimensional_system

def test_all_internal_families_attempted():
    assert len(execute_internal_candidates(default_dimensional_system().audit()))==21
    assert len(candidate_registry())==30

def test_target_and_lcdm_contamination_rejected():
    c=NativeScaleEstimate("X","test",L0_m_per_native=1,status="SCALE_ESTIMATE",independent_of_target=False)
    assert validate_candidate(c).status=="FORBIDDEN_INPUT_DEPENDENCE"
    c=NativeScaleEstimate("X","test",L0_m_per_native=1,status="SCALE_ESTIMATE",independent_of_lcdm=False)
    assert validate_candidate(c).rejection_reason=="FORBIDDEN_LCDM_DEPENDENCE"

def test_forbidden_false_rulers_rejected():
    for text in ("Rmax candidate","historical 0.18","Planck length by framework name"):
        c=NativeScaleEstimate("X",text,L0_m_per_native=1,status="SCALE_ESTIMATE")
        assert validate_candidate(c).status=="FORBIDDEN_INPUT_DEPENDENCE"
