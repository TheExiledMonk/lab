import json
import numpy as np
from pbuf.labs.native_mechanisms.candidate_base import evaluate_candidates, transfer_metrics

def test_transfer_metrics_and_wide_inventory():
    rows,tx=evaluate_candidates(np.array([.1,-.2,.4,-.1,.05,-.3]))
    assert len(rows)==16 and {r["ID"] for r in rows}=={f"H{i:02d}" for i in range(16)}
    assert np.isclose(tx["H07"]["loaded"]["T_total"],1)
    assert any(abs(x)>0 for x in tx["H07"]["delta_T"])
    assert next(r for r in rows if r["ID"]=="H07")["STATUS"]=="PARTIAL"
    assert next(r for r in rows if r["ID"]=="H00")["LOADED_DIRECTIONAL_REDIRECTION"]=="FALSE"

def test_runner_emits_contract_and_artifacts(tmp_path,monkeypatch):
    from pbuf.labs.foundation import native_medium_interaction_wide_net001 as lab
    monkeypatch.setattr(lab,"OUT",tmp_path)
    c=lab.main()
    assert c["DEV165_AUDIT_COMPLETE"] and c["CANDIDATE_COUNT"]==16
    assert c["SURVIVING_MECHANISM_COUNT"]==0 and not c["OBSERVATIONAL_TARGET_USED"]
    required={"report.txt","frozen_input_contract.json","candidate_inventory.json","candidate_complexity_inventory.json",
      "unloaded_equilibrium_results.json","n6_symmetry_results.json","free_propagation_results.json",
      "dev157_dispersion_results.json","dev159_static_compatibility.json","finite_state_compatibility.json",
      "loaded_directional_transfer.json","reversibility_results.json","invariant_results.json",
      "magnetic_like_results.json","separation_results.json","memory_results.json","allocation_results.json",
      "candidate_equivalence_classes.json","candidate_dominance_matrix.json","surviving_candidate_contract.json",
      "downstream_validity_matrix.json","final_native_mechanism_contract.json"}
    assert required<={p.name for p in tmp_path.iterdir()}
    inv=json.loads((tmp_path/"candidate_inventory.json").read_text())
    assert all(r["NEW_FREE_COEFFICIENT_COUNT"]==0 for r in inv["candidates"])
