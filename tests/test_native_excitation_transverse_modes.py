import numpy as np
from pbuf.excitation.native_excitation_invariants import transverse_rank_audit
from pbuf.excitation.native_excitation_state import localized_packet

def test_two_independent_transverse_modes_and_combination():
    assert transverse_rank_audit()['physical_rank']==2
    a=localized_packet(polarization=(1,0)); b=localized_packet(polarization=(0,1)); c=localized_packet(polarization=(1,1))
    peak=int(np.argmax(np.linalg.norm(a,axis=1)))
    assert np.linalg.matrix_rank(np.stack([a[peak],b[peak]]))==2 and np.any(c[:,0]*c[:,1])
