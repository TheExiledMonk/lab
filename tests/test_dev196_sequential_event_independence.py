import json, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from pbuf.observer.sequential_event_independence import inject, support_relation
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState
OUT=ROOT/'runs/dev196_sequential_event_independence'
def load(n): return json.loads((OUT/n).read_text())
def test_contract_and_semantics():
 c=load('final_contract.json'); assert c['NO_NEW_INJECTION_LAW'] and c['NO_RESET_BETWEEN_A_AND_B'] and c['NO_A_REMOVAL']
 assert load('second_event_injection_semantics.json')['SECOND_EVENT_INJECTION_SEMANTICS']=='EXISTING_SEMANTICS_APPLICABLE_WITH_STRUCTURAL_REUSE'
def test_injection_is_additive_and_deterministic():
 u=np.zeros((3,3,3,3)); p=u.copy(); du=np.ones_like(u)*.01; dp=np.ones_like(p)*.02
 s=inject(VectorPairState(u,p),du,dp); assert np.array_equal(s.displacement,u+du) and np.array_equal(s.momentum,p+dp)
def test_matched_subtractions_and_overlap_are_serialized():
 for label in ('T1','T2','T3'):
  a=np.load(OUT/f'delta_B_fresh_{label}.npz');b=np.load(OUT/f'delta_B_after_A_{label}.npz');r=np.load(OUT/f'interaction_residual_{label}.npz')
  for key in a.files: assert np.array_equal(r[key],b[key]-a[key])
  o=np.load(OUT/f'support_overlap_{label}.npz'); assert o['A_support'].shape==o['B_support'].shape==o['overlap'].shape
def test_overlap_classification_deterministic():
 a=np.zeros((2,2,2),bool);b=a.copy(); assert support_relation(a,b,False)=='DISJOINT'
