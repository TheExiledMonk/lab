from pbuf.quantum.native_excitation_interaction import norm_exchange
def test_transition_accounting_identity():
    assert norm_exchange(3,1,2)['conserved']
    assert not norm_exchange(3,1,1)['conserved']
