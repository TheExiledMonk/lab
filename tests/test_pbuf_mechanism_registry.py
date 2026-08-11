import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pbuf.registry.search import search
from pbuf.registry.validate import validate


def registry():
    return json.loads((ROOT / "docs/PBUF_MECHANISM_REGISTRY.json").read_text())


def test_registry_is_valid_and_required_history_is_present():
    data = registry()
    assert not validate(data)
    ids = {a["attempt_id"] for a in data["attempts"]}
    for required in [
        "dev167_vector_relational_dynamics", "dev168_finite_loaded_receipt",
        "dev171_independent_3d_source", "dev177_full_received_state",
        "dev178_vulkan_viewer", "dev179_subcell_source_closure",
        "dev180_density_reconciliation", "pr36_effective_matter_loading",
        "pr37_metric_strain_normalization", "pr69_forcing_redistribution",
        "pr70_bounded_strain_constitutive", "pr72_local_n6_redistribution",
        "pr74_a8_cstate_accumulation", "dev159_local_source_forcing",
    ]:
        assert required in ids


def test_alias_search_and_equivalence_are_retrievable():
    data = registry()
    targets, attempts = search(data, "source medium")
    assert "source_medium_coupling" in {t["target_id"] for t in targets}
    assert "dev180_density_reconciliation" in {a["attempt_id"] for a in attempts}
    targets, attempts = search(data, "70,756")
    assert "dev180_density_reconciliation" in {a["attempt_id"] for a in attempts}
    assert any(e["source"] == "pr70_bounded_strain_constitutive" and e["target"] == "dev167_vector_relational_dynamics" for e in data["equivalences"])
