import json
import numpy as np

from pbuf.excitation.native_finite_receipt import (
    NativeReceivedState, crossing_bond_flux, flux_vectors,
    local_content_candidates, unit_directions,
)
from pbuf.excitation.native_vector_pair_dynamics import pair_power_flux


def test_local_content_candidates_are_native_local_and_nonnegative():
    rng=np.random.default_rng(168)
    u=rng.normal(0,1e-4,(5,6,7,3)); p=rng.normal(0,1e-4,u.shape)
    w=local_content_candidates(u,p)
    assert w.shape==(5,6,7,4)
    assert np.isfinite(w).all() and np.all(w>=0)
    np.testing.assert_allclose(w[...,2],w[...,0]+w[...,1])


def test_flux_vector_and_predeclared_plane_are_exact_dev167_derived_values():
    rng=np.random.default_rng(169)
    u=rng.normal(0,1e-4,(5,6,7,3)); p=rng.normal(0,1e-4,u.shape)
    j=pair_power_flux(u,p)
    assert flux_vectors(u,p).shape==u.shape
    np.testing.assert_array_equal(crossing_bond_flux(u,p,3),j[2,...,0])


def test_received_state_preserves_position_direction_lineage_and_spread():
    n=4; xyz=np.arange(12,dtype=float).reshape(n,3); d=unit_directions(np.ones((n,3)))
    state=NativeReceivedState(xyz,xyz+.5,d,np.ones(n),np.arange(n),np.arange(n),
        np.zeros((n,3)),np.ones((n,3)),np.ones((n,3)),np.ones((n,4)),"BOND_FLUX")
    assert set(state.arrays()) >= {"source_positions","received_positions","directions","weights"}
    assert np.ptp(state.received_positions,axis=0).sum()>0


def test_dev168_committed_artifact_contract():
    from tools import generate_dev168_finite_receipt as lab
    required={"report.txt","repository_contract.json","ledger_frontier_update.json","historical_crosscheck.json",
      "packet_definition.json","packet_free_propagation.json","packet_loaded_propagation.json",
      "persistent_source_lane.json","frozen_load_lane.json","receipt_surface_contract.json","node_receipt.json",
      "flux_receipt.json","receipt_representation_comparison.json","received_state_schema.json",
      "received_state_full.npz","received_position_audit.json","received_direction_audit.json",
      "momentum_flux_direction_comparison.json","receipt_weight_candidates.json","receipt_content_closure.json",
      "transverse_receipt_audit.json","symmetry_controls.json","step_convergence.json",
      "packet_amplitude_ladder.json","source_loading_ladder.json","offset_ladder.json","invariant_accounting.json",
      "boundary_contamination_audit.json","observer_primitive_mapping.json","observer_adapter_contract.json",
      "final_contract.json","discussion_handoff.md"}
    assert required <= {p.name for p in lab.OUT.iterdir()}
    contract=json.loads((lab.OUT/"final_contract.json").read_text())
    assert contract["DEV168_COMPLETE"] and contract["PRIMARY_TOPOLOGY"]=="N6"
    assert contract["DEV167_MECHANISM_MODIFIED"] is False
    assert contract["OBSERVER_MODIFIED"] is False and contract["OBSERVER_EXECUTED_PRIMARY"] is False
    assert contract["ADAPTER_INTRODUCES_NEW_PHYSICS"] is False
    assert contract["NEW_NATIVE_PROPAGATION_LAW_INTRODUCED"] is False
