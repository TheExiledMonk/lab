from pbuf.quantum.native_excitation_interaction import coefficient_free_coupling_audit
from pbuf.quantum.native_transition_quantization import quantization_location_registry,quantization_family_registry
def test_missing_coupling_is_explicit():
    assert not coefficient_free_coupling_audit()['coefficient_free_binding_law_found']
    assert len(quantization_location_registry())==10 and len(quantization_family_registry())==20
