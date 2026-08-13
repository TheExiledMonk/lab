import numpy as np
from pbuf.observer.native_pair_partition import derive_partition, translate_partition

def test_exchange_and_exact_periodic_translation_covariance():
    p=derive_partition((11,11,11),(1,5,5),(1,9,5)); q=derive_partition((11,11,11),(1,9,5),(1,5,5)); r=translate_partition(p,(0,1,0))
    assert np.array_equal(q.omega_a,p.omega_b) and np.array_equal(q.omega_b,p.omega_a)
    assert np.array_equal(r.omega_a,np.roll(p.omega_a,(0,1,0),(0,1,2)))
    assert np.array_equal(np.roll(np.flip(p.omega_a,axis=1),4,axis=1),p.omega_b)
