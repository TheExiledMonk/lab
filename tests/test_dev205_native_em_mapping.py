import json
from pathlib import Path
import numpy as np


def test_dev205_artifacts_keep_mapping_boundary_and_inventory():
    root=Path(__file__).resolve().parents[1]; out=root/'runs/dev205_native_em_mapping'
    c=json.loads((out/'final_contract.json').read_text())
    assert c['NO_NEW_FIELD'] and c['NO_NEW_DOF'] and c['NO_EM_COEFFICIENT_FITTING']
    assert c['NATIVE_RESULTS_FROZEN_BEFORE_EM_COMPARISON']
    inv=json.loads((out/'predeclared_candidate_inventory.json').read_text())
    assert inv['PREDECLARED_AXIAL_CANDIDATES_COMPLETE']
    a=np.load(out/'native_axial_candidates.npz')
    assert {'Q_R','Q_F'} <= set(a.files)


def test_dev205_does_not_promote_partial_geometry_to_full_em_identity():
    out=Path(__file__).resolve().parents[1]/'runs/dev205_native_em_mapping'
    assert json.loads((out/'electric_field_mapping_status.json').read_text())['ELECTRIC_FIELD_MAPPING']=='PARTIAL'
    assert json.loads((out/'magnetic_field_mapping_status.json').read_text())['MAGNETIC_FIELD_MAPPING']=='PARTIAL'
