import json, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from pbuf.observer.independent_event_transport import expected_accumulation, repeated_response_accumulation
OUT=ROOT/'runs/dev194_independent_event_wave_transport'
def load(n): return json.loads((OUT/n).read_text())
def test_contract_preserves_event_boundary():
 c=load('final_contract.json'); assert c['INDEPENDENT_EVENT_TRANSPORT']=='BLOCKED_EVENT_INDEPENDENCE'; assert c['NATIVE_EVENT_IMAGE_TRANSPORT']=='BLOCKED_EVENT_INDEPENDENCE'; assert c['DOUBLE_SLIT_LANE']=='FOUNDATIONAL_DIAGNOSTIC_ONLY'
def test_detector_side_accumulation_is_exact_algebra():
 z=np.load(ROOT/'runs/dev188_native_source_distribution_pushforward/transfer_kernel_weight.npz',allow_pickle=False); k=z['R00_weight']; n=np.zeros(121);n[[2,70,119]]=[3,2,5]
 assert np.array_equal(expected_accumulation(k,n), repeated_response_accumulation(k,n))
def test_required_controls_are_classified():
 assert load('two_path_nonadditivity.json')['status']=='NOT_REPRESENTABLE'
 assert load('post_event_medium_memory.json')['POST_EVENT_MEDIUM_MEMORY']=='PERSISTENT'
