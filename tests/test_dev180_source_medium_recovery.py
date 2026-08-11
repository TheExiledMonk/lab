import json
from pathlib import Path

from tools import generate_dev180_source_medium_recovery as dev180


def test_dev180_recovery_contract(tmp_path, monkeypatch):
    monkeypatch.setattr(dev180, "OUT", tmp_path)
    dev180.main()
    required = {"starting_state.json", "repo_search_manifest.json", "historical_source_medium_map.json", "historical_source_medium_concept_recovery.json", "pr36_matter_loading_recovery.json", "pr37_metric_strain_recovery.json", "pr69_neighbor_source_recovery.json", "pr70_constitutive_recovery.json", "pr72_n6_network_recovery.json", "pr74_native_accumulation_recovery.json", "current_source_contact_audit.json", "current_packet_launch_audit.json", "current_receipt_multiplicity_tree.json", "historical_25pct_density_semantics.json", "density_terminology.json", "density_architecture_result.json", "source_medium_correspondence_matrix.json", "dev179_scope_correction.json", "source_work_dev170_reconciliation.json", "current_native_reuse_candidates.json", "viewer_extension_status.json", "final_contract.json", "discussion_handoff.md"}
    assert required <= {p.name for p in Path(tmp_path).iterdir()}
    final = json.loads((Path(tmp_path) / "final_contract.json").read_text())
    assert final["NO_NEW_SOURCE_MEDIUM_LAW"] and final["NO_NEW_PAIR_LAW"]
    density = json.loads((Path(tmp_path) / "historical_25pct_density_semantics.json").read_text())
    assert density["HISTORICAL_25PCT_RAY_DENSITY_NOT_SOURCE_DENSITY"]
