from pbuf.quantum.native_excitation_modes import carrier_mode
from pbuf.quantum.native_excitation_momentum import candidate_registry,momentum_audit
def test_directional_flux_is_structural_but_k_relation_unresolved():
    assert len(candidate_registry())==7; r=momentum_audit(carrier_mode(128,16),-1); assert r['direction']==-1 and r['momentum_like_established'] is False

