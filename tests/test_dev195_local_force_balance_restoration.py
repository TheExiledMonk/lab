import json, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from pbuf.observer.local_force_balance_restoration import shell_partition
from pbuf.excitation.native_vector_pair_dynamics import pair_reciprocity_error
OUT=ROOT/'runs/dev195_local_force_balance_restoration'
def load(n): return json.loads((OUT/n).read_text())
def test_shells_are_an_exact_partition():
 d,s=shell_partition((11,11,11),(1,5,5)); assert sum(x.sum() for x in s)==d.size and np.all(sum(s)==1)
def test_dev195_preserves_frozen_contract_and_time_matching():
 c=load('final_contract.json'); assert c['NO_DAMPING'] and c['NO_RECOVERY_TERM'] and c['NO_ABSORPTION_TERM'] and c['TIME_MATCHED_COMPARISON']
 z=np.load(OUT/'time_matched_difference.npz',allow_pickle=False); assert z['displacement'].shape==z['momentum'].shape==z['net_force'].shape
def test_dev194_endpoint_is_spatially_resolved():
 d=load('dev194_endpoint_spatial_decomposition.json'); assert d['global_displacement_l2']>0 and 'dominant_displacement_shell' in d
def test_force_diagnostic_does_not_change_pair_reciprocity():
 z=np.load(OUT/'excited_trajectory.npz',allow_pickle=False)
 assert pair_reciprocity_error(z['displacement'][0])==0.0
