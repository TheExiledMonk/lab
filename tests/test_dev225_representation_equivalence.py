import numpy as np
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pbuf.analysis.native_handedness_representation import axial_neighbor_relation,frobenius_neighbor_relation
def test_frobenius_is_positive_multiple_of_axial_dot():
 a=np.array([[0.,2.,-3.],[-2.,0.,5.],[3.,-5.,0.]])
 b=np.array([[0.,7.,11.],[-7.,0.,13.],[-11.,-13.,0.]])
 assert frobenius_neighbor_relation(a,b)==2*axial_neighbor_relation(a,b)
