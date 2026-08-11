import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/dev171_independent_3d_abell001'

def test_dev171_frozen_catalog_ensemble_contract():
    required=['spectroscopic_source_inventory.json','cluster_member_phase_space.npy','phase_space_component_model.json','component_depth_constraints.json','source_3d_ensemble_manifest.json','source_3d_freeze_contract.json','constrained_3d_observer_spread.json','depth_uncertainty_reduction.json','final_contract.json']
    assert all((OUT/name).exists() for name in required)
    final=json.loads((OUT/'final_contract.json').read_text())
    assert final['SOURCE_3D_ENSEMBLE_FROZEN'] and final['FULL_NATIVE_ENSEMBLE_EXECUTED']
    assert not final['SPECTROSCOPIC_REDSHIFT_DIRECTLY_USED_AS_GEOMETRIC_DEPTH']
    assert not final['LENSING_DERIVED_SOURCE_INFORMATION_USED']
    assert not final['DEV167_PAIR_LAW_MODIFIED'] and not final['DEV168_RECEIPT_MODIFIED']
    assert not final['OBSERVER_PHYSICS_MODIFIED'] and final['3D_NATIVE_SMOKE_PASS']
    manifest=json.loads((OUT/'source_3d_ensemble_manifest.json').read_text())
    assert manifest['ensemble_count']==8 and len(manifest['realizations'])==8
    q=np.load(OUT/'cluster_member_phase_space.npy')
    assert q.ndim==2 and q.shape[1]==3 and len(q)>0
