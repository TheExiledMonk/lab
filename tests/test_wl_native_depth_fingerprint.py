import numpy as np
from pbuf.wl.native_depth_fingerprint import *
from pbuf.wl.native_source_controls import deterministic_split
def test_monotonic_and_training_only():
 z=np.arange(10.);assert monotonicity(z*z,z)["classification"]=="STRICTLY_MONOTONIC"
 ids=[f"x{i}" for i in range(10)];b=build_fingerprint_bank(z,z[:,None],ids,deterministic_split)
 assert nearest_fingerprint([4.],b) in b["depths"][b["split"]=="TRAIN"]
