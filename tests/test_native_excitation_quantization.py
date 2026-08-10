from pbuf.quantum.native_excitation_modes import carrier_mode
from pbuf.quantum.native_excitation_quantization import candidate_registry,divisibility_audit,boundary_modes
def test_free_propagation_is_continuous_and_boundary_modes_are_geometric():
    assert divisibility_audit(carrier_mode(128,16))['quantization_established'] is False; assert len(candidate_registry())==20; assert boundary_modes(8)[0]['wavelength']==16

