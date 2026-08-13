import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pbuf.observer.field_strength_cross_event import force_orientation,ratio
OUT=ROOT/'runs/dev198_field_strength_cross_event';IN=ROOT/'runs/dev196_sequential_event_independence'
def test_exact_zero_denominators_remain_undefined():
 r,d=ratio(np.array([[1.,0,0]]),np.zeros((1,3)));assert not d[0] and np.isnan(r[0])
 c,d=force_orientation(np.array([[1.,0,0]]),np.zeros((1,3)),np.array([True]));assert not d[0] and np.isnan(c[0])
def test_hashes_and_complete_fixed_time_sequence():
 m=json.loads((OUT/'dev197_input_manifest.json').read_text())['sha256']
 for n,h in m.items(): assert hashlib.sha256((IN/n).read_bytes()).hexdigest()==h
 t=json.loads((OUT/'field_strength_time_series.json').read_text())['records'];assert [x['time'] for x in t]==list(range(181))
def test_contract_forbids_new_physics_and_thresholds():
 c=json.loads((OUT/'final_contract.json').read_text())
 for key in ('NO_NEW_FORCE_LAW','NO_PACKET_AMPLITUDE_SWEEP','NO_PHYSICAL_THRESHOLD','NO_SUPPORT_THRESHOLD','NO_EPSILON_REGULARIZATION','NO_QED_FIT'):assert c[key]
def test_frozen_force_law_and_packet_amplitude_are_not_edited():
 law=(ROOT/'pbuf/excitation/native_vector_pair_dynamics.py').read_text();packet=(ROOT/'tools/generate_dev169_raw_abell_native_observer.py').read_text()
 assert 'return e / (1.0 - e*e)' in law and '.006*env' in packet
