import json

import numpy as np

from pbuf.excitation.native_loaded_background_dynamics import linearization_audit, perturbation_step
from pbuf.excitation.native_relational_dynamics import f03_step
from pbuf.lens.native_stationary_lens_from_source import stationary_distributed_response
from pbuf.source.projected_source_3d_family import diagnostic_family


def test_loaded_equilibrium_linearization_is_exactly_the_free_operator():
    image = np.arange(1, 43, dtype=float).reshape(6, 7); image /= image.sum()
    realization = diagnostic_family(image)[4]
    q0 = stationary_distributed_response(realization.source)
    rng = np.random.default_rng(163)
    dq = rng.normal(size=q0.shape) * 1e-7; dr = rng.normal(size=q0.shape) * 1e-7
    loaded = perturbation_step(q0, realization.source, dq, dr)
    free = f03_step(dq, dr)
    assert np.allclose(loaded[0], free[0], rtol=0, atol=2e-16)
    assert np.allclose(loaded[1], free[1], rtol=0, atol=2e-16)
    assert linearization_audit(q0, realization.source, dq, dr)["loaded_operator_depends_on_background"] is False


def test_dev163_stops_at_failed_coupling_gate(tmp_path, monkeypatch):
    from pbuf.labs.foundation import raw_abell2744_finite_native_lensing_gate001 as lab
    monkeypatch.setattr(lab, "OUT", tmp_path)
    contract = lab.main()
    assert contract["DEV163_AUDIT_COMPLETE"]
    assert contract["LOADED_DYNAMIC_COUPLING_DERIVED"] == "FALSE"
    assert contract["FINITE_NATIVE_LOADED_RESPONSE"] == "NOT_TESTED_GATE_FAILED"
    assert not contract["ARBITRARY_LOADING_COUPLING_INTRODUCED"]
    required = {"loaded_dynamic_candidate_inventory.json", "loaded_equilibrium_linearization.json",
                "zero_load_recovery.json", "loaded_dynamic_invariant.json",
                "loaded_coupling_contract.json", "required_test_results.json", "report.txt"}
    assert required <= {p.name for p in tmp_path.iterdir()}
    candidates = json.loads((tmp_path / "loaded_dynamic_candidate_inventory.json").read_text())
    assert candidates["surviving_nontrivial_loaded_candidate"] is None
    assert not (tmp_path / "loaded_propagation_per_depth.json").exists()
