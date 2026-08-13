import numpy as np
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pbuf.analysis.native_handedness_representation import axial_dual,frobenius_neighbor_relation
def test_axial_reflection_and_relation_invariance():
 q=np.diag([-1.,1.,1.]); a=np.array([[0.,2.,-3.],[-2.,0.,5.],[3.,-5.,0.]])
 ap=q@a@q.T
 assert np.array_equal(axial_dual(ap),np.linalg.det(q)*q@axial_dual(a))
 assert frobenius_neighbor_relation(ap,ap)==frobenius_neighbor_relation(a,a)
