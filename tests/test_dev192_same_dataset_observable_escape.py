import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from pbuf.observer.same_dataset_escape import REQUIRED_ESCAPE_FIELDS,validate_escape
OUT=ROOT/'runs/dev192_same_dataset_observable_escape'
def load(n):return json.loads((OUT/n).read_text())
def test_matrix_and_result_are_complete():
 rows=load('expanded_escape_candidate_matrix.json')['candidates'];assert len(rows)==22
 for x in rows: assert set(REQUIRED_ESCAPE_FIELDS)<=set(x);validate_escape(x)
 assert load('same_dataset_observable_escape_result.json')['SAME_DATASET_OBSERVABLE_ESCAPE']=='NONE'
def test_contract_preserves_locks_and_closed_gates():
 c=load('final_contract.json')
 for k in ('DEVELOPMENT_DATASET_LOCKED','NO_NEW_LENS_DATA_OPENED','ALL_ESCAPE_LANES_TESTED_BEFORE_SELECTION','NO_INTRINSIC_SOURCE_RECONSTRUCTION','NO_SOURCE_PRIOR_FIT','NO_LENSING_SCORE','NO_CHI2','PIPELINE_DETERMINISTIC'):assert c[k] is True
 assert c['OTHER_LENS_VALIDATION_GATE']=='CLOSED_PENDING_OBSERVER_FREEZE';assert c['OBSERVATIONAL_COMPARISON_GATE']=='CLOSED'
