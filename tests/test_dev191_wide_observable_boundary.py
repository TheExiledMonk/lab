import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev191_wide_observable_boundary'
sys.path.insert(0, str(ROOT))
from pbuf.observer.observable_boundary import REQUIRED_CANDIDATE_FIELDS, validate_candidate
def load(name): return json.loads((OUT/name).read_text())
def test_dev191_artifacts_and_matrix_are_complete():
    matrix=load('observable_candidate_matrix.json')['candidates']
    assert len(matrix) >= 29
    assert {x['candidate_id'] for x in matrix} >= {f'O{i:02d}' for i in range(1,30)}
    for row in matrix:
        assert set(REQUIRED_CANDIDATE_FIELDS) <= set(row); validate_candidate(row)
    assert load('observable_boundary_result.json')['OBSERVABLE_BOUNDARY_AUDIT']=='STRUCTURAL_EQUIVALENCE_ONLY'
def test_dev191_hard_gates_remain_closed():
    c=load('final_contract.json')
    for key in ('MECHANISM_REGISTRY_QUERIED','HISTORICAL_OBSERVER_WORK_INSPECTED','HISTORICAL_MOMENT_WORK_INSPECTED','EXTERNAL_METHODOLOGY_RESEARCH_COMPLETED','OBSERVED_E1_E2_VALUES_BLINDED','NO_LENSING_SCORE','NO_PBUF_OBSERVATION_RESIDUAL','NO_INTRINSIC_SOURCE_RECONSTRUCTION','NO_GR_LENS_INVERSION','NO_LCDM_DISTANCE','NO_FITTED_ROTATION','NO_FITTED_SCALE','NO_OBSERVABLE_FUSION','PIPELINE_DETERMINISTIC'):
        assert c[key] is True
    assert c['SPIN2_OBSERVABLE_GATE']=='CLOSED'; assert c['OBSERVATIONAL_COMPARISON_GATE']=='CLOSED'
