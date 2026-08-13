import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import tools.generate_dev210_exact_local_em_relay as dev210


def test_dev210_preserves_the_negative_local_launch_boundary(tmp_path, monkeypatch):
    monkeypatch.setattr(dev210, "OUT", tmp_path)
    monkeypatch.setattr(dev210, "update_memory", lambda: None)
    dev210.main()
    required = {"existing_local_preparation_inventory.json", "source_release_semantics.json",
                "source_force_support.json", "source_induced_state_support.json", "final_contract.json"}
    assert required <= {p.name for p in Path(tmp_path).iterdir()}
    assert json.loads((tmp_path / "local_launch_eligibility.json").read_text())["EXACT_LOCAL_NATIVE_LAUNCH"] == "NOT_DERIVED"
    assert json.loads((tmp_path / "source_force_support.json").read_text())["support_count"] == 6
    assert json.loads((tmp_path / "source_release_semantics.json").read_text())["NATIVE_SOURCE_RELEASE_SEMANTICS"] == "NOT_DERIVED"
