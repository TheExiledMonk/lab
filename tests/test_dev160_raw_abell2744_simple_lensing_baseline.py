import json
from pathlib import Path

from pbuf.labs.foundation import raw_abell2744_simple_lensing_baseline001 as dev160


def test_dev160_fail_closed_raw_audit(tmp_path, monkeypatch):
    monkeypatch.setattr(dev160, "OUT", tmp_path)
    assert dev160.main() == 0
    required = {
        "report.txt", "raw_input_inventory.json", "active_pipeline_inventory.json",
        "native_source_support.json", "native_lens_support.json", "propagation_object_contract.json",
        "lens_sampling_contract.json", "null_vs_loaded_response.json", "full_received_state_inventory.json",
        "projection_information_loss.json", "reversal_stage_audit.json", "blocker_classification.json",
        "dev159_relevance_contract.json", "downstream_validity_matrix.json", "final_raw_a2744_lensing_contract.json"}
    assert required <= {p.name for p in tmp_path.iterdir()}
    contract = json.loads((tmp_path / "final_raw_a2744_lensing_contract.json").read_text())
    assert contract["DEV160_AUDIT_COMPLETE"] is True
    assert contract["RAW_ABELL2744_PIPELINE_LOCATED"] == "PARTIAL"
    assert contract["ACTIVE_RAW_LENSING_RUNNER_IDENTIFIED"] is False
    assert contract["FIVE_CLUSTER_BENCHMARK_USED_AS_DEVELOPMENT_BASELINE"] is False
    assert contract["LENS_PHYSICS_MODIFIED"] is False
    assert contract["PROPAGATION_PHYSICS_MODIFIED"] is False
    assert contract["DEV159_STATE_SUBSTITUTED_INTO_LENSING"] is False


def test_current_production_contract_is_recorded_without_execution(tmp_path, monkeypatch):
    monkeypatch.setattr(dev160, "OUT", tmp_path)
    dev160.main()
    propagation = json.loads((tmp_path / "propagation_object_contract.json").read_text())
    sampling = json.loads((tmp_path / "lens_sampling_contract.json").read_text())
    response = json.loads((tmp_path / "null_vs_loaded_response.json").read_text())
    assert propagation["current_production_propagation_object"] == "ZERO_WIDTH_RAY"
    assert sampling["current_production_method"] == "POINT_SAMPLE"
    assert response["status"] == "NOT_EXECUTED"
    assert "benchmark rho3" in response["why_fail_closed"]
