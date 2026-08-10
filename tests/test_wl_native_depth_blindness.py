from dataclasses import fields
import numpy as np
from pbuf.wl.native_source_controls import *
def test_truth_erased_from_blind_package():
 t=SourceTruth("x",2.,np.zeros((3,2)),1.,1.)
 b=blind_package(t,lambda p,z:(p,np.ones_like(p)))
 assert not {"depth_native","positions_native","lens_depth_native"}&{x.name for x in fields(b)}
