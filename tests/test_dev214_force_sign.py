import numpy as np
from pbuf.observer.native_dynamic_pair_force import classify
def test_action_reaction_gate_rejects_unbalanced_support_force():
 kind,*_=classify(np.array([1.,0,0]),np.zeros(3),np.array([1.,0,0]),1e-12)
 assert kind=='UNRESOLVED'
