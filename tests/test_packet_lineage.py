import numpy as np
import pytest
from pbuf.excitation.native_finite_receipt import NativeReceivedState
from pbuf.receipt.packet_lineage import PacketAwareReceiptCollection

def state(n=2):
 return NativeReceivedState(*(np.zeros((n,3)) for _ in range(3)),np.ones(n),np.zeros(n,int),np.zeros(n,int),*(np.zeros((n,3)) for _ in range(3)),np.zeros((n,4)),"BOND_FLUX")
def test_packet_lineage_rejects_unknown_and_misaligned_indices():
 good=({'launch_id':'L'},); packet=({'packet_id':'P'},); real=({'realization_id':'R'},)
 with pytest.raises(ValueError): PacketAwareReceiptCollection(state(),np.array([0]),np.zeros(2,int),np.zeros(2,int),good,packet,real)
 with pytest.raises(ValueError): PacketAwareReceiptCollection(state(),np.array([0,1]),np.zeros(2,int),np.zeros(2,int),good,packet,real)
