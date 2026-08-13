import json
from pathlib import Path
from pbuf.observer.weak_event_regime import tangent_force,positive_force
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_tangent_and_contract_are_native_and_unthresholded():
 x=np.zeros((3,3,3,3)); d=np.zeros_like(x); d[0,0,0,0]=1e-6
 assert tangent_force(x,d).shape==positive_force(x).shape
 c=json.loads((ROOT/'runs/dev200_native_n6_field/final_contract.json').read_text())
 assert c['NO_PHYSICAL_THRESHOLD'] and c['NO_PACKET_AMPLITUDE_SWEEP']
 assert c['POLARIZATION_MAPPING']=='UNESTABLISHED'
