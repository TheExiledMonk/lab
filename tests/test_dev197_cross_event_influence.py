import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pbuf.observer.cross_event_influence import ratios,parallel_transverse
OUT=ROOT/'runs/dev197_cross_event_influence'; IN=ROOT/'runs/dev196_sequential_event_independence'
def test_zero_denominator_is_not_regularized():
 r,d=ratios(np.array([[1.,0.,0.]]),np.zeros((1,3)));assert not d[0] and np.isnan(r[0])
def test_parallel_transverse_is_exact_decomposition():
 e=np.array([[3.,4.,0.]]);b=np.array([[2.,0.,0.]]);p,t,d=parallel_transverse(e,b);assert d[0] and p[0]==3 and np.array_equal(t,np.array([[0.,4.,0.]]))
def test_frozen_dev196_residuals_unchanged_and_reused():
 m=json.loads((OUT/'dev196_input_manifest.json').read_text());import hashlib
 for n,h in m['sha256'].items():assert hashlib.sha256((IN/n).read_bytes()).hexdigest()==h
def test_ratios_and_receipt_are_serialized_without_composite_score():
 for t in ('T1','T2','T3'):
  x=json.loads((OUT/f'relative_influence_ratios_{t}.json').read_text());assert set(x)=={'displacement','momentum','force','flux'}
  r=json.loads((OUT/f'receipt_influence_{t}.json').read_text());assert 'receipt_cell_difference' in r and 'receipt_weight_difference' in r
 c=json.loads((OUT/'final_contract.json').read_text());assert c['EXACT_SUPPORT_ROLE']=='DIAGNOSTIC_ONLY' and c['NO_COMPOSITE_SCORE']
