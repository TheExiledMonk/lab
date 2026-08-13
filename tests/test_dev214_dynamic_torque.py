import numpy as np
from pbuf.observer.native_dynamic_pair_torque import support_torque
def test_fixed_provenance_torque_is_zero_for_zero_force():
 m=np.ones((3,3,3),bool); assert np.array_equal(support_torque(np.zeros((3,3,3,3)),m,np.zeros(3)),np.zeros(3))
