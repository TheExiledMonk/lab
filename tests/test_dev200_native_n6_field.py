from pathlib import Path
import numpy as np
from pbuf.observer.native_n6_field import n6_field
ROOT=Path(__file__).resolve().parents[1]
def test_n6_exact_topology_and_artifacts():
 x=np.zeros((3,3,3,3)); f=n6_field(x,x)
 assert f['force'].shape == (3,3,3,6,3)
 assert np.allclose(f['strain'],0)
 out=ROOT/'runs/dev200_native_n6_field'; assert (out/'final_contract.json').exists()
 z=np.load(out/'canonical_packet_n6_field.npz'); assert z['force'].shape[-2:]==(6,3)
