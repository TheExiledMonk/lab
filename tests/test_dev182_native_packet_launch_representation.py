import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev182_native_packet_launch_representation"


def test_dev182_representation_artifacts_and_closure():
    required = {"starting_state.json", "registry_lookup.json", "historical_launch_semantics.json", "historical_c25_count_reconciliation.json", "current_packet_semantics.json", "current_launch_multiplicity_tree.json", "packet_initial_condition_space.json", "independent_replay_test.json", "launch_order_independence.json", "translation_covariance.json", "reflection_covariance.json", "launch_representation_status.json", "density_authorization.json", "final_contract.json"}
    assert required <= {p.name for p in OUT.iterdir()}
    final = json.loads((OUT / "final_contract.json").read_text())
    assert final["OUTCOME"] == "OUTCOME_B"
    assert final["CONTINUOUS_LAUNCH_COORDINATES_NOT_ASSUMED"]
    assert json.loads((OUT / "independent_replay_test.json").read_text())["INITIAL_LOADED_MEDIUM_IDENTICAL_FOR_ALL_REPLAYS"]
    assert json.loads((OUT / "launch_order_independence.json").read_text())["passed"]
