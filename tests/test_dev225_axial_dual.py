import numpy as np
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pbuf.analysis.native_handedness_representation import axial_dual
def test_fixed_right_handed_axial_dual():
 a=np.array([[0.,2.,3.],[-2.,0.,5.],[-3.,-5.,0.]])
 assert np.array_equal(axial_dual(a),[5.,-3.,2.])
