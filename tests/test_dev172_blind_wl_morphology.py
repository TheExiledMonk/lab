import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev172_blind_wl_morphology001"
DEV171 = ROOT / "runs/dev171_independent_3d_abell001"


def test_dev172_blocks_without_substituting_or_aligning_a_wl_map():
    final = json.loads((OUT / "final_contract.json").read_text())
    grid = json.loads((OUT / "comparison_grid_contract.json").read_text())
    provenance = json.loads((OUT / "observational_wl_provenance.json").read_text())
    assert final["DEV172_COMPLETE"] is False
    assert final["OUTCOME"] == "OUTCOME_D"
    assert final["COMPARISON_REQUIRES_MISSING_PHYSICAL_BRIDGE"] is True
    assert provenance["WL_DATA_PROVENANCE_FROZEN"] is False
    assert grid["common_grid_established"] is False
    assert all(grid[k] is False for k in ("TRANSLATION_FITTED_TO_WL", "ROTATION_FITTED_TO_WL", "SCALE_FITTED_TO_WL", "MIRROR_TRANSFORM_FITTED_TO_WL"))


def test_dev172_reconciles_stale_ledger_without_changing_dev171_outputs():
    reconciliation = json.loads((OUT / "dev171_metric_reconciliation.json").read_text())
    repo = json.loads((OUT / "repository_provenance.json").read_text())
    assert reconciliation["classification"] == "SAME_METRIC_STALE_LEDGER"
    assert reconciliation["correct_values"]["V_constrained"] == 0.10624528343935943
    assert reconciliation["correct_values"]["R_3D"] == 1.832080315334409
    for name, expected in repo["dev171_output_sha256"].items():
        assert hashlib.sha256((DEV171 / name).read_bytes()).hexdigest() == expected
