import json

import numpy as np

from pbuf.geometry.native_path_geometry import path_diagnostics, straight_path
from pbuf.geometry.native_relational_geometry import (directional_scalar_information,
    reciprocity_error_for_scalar_edges, undeformed_bond_vectors)
from pbuf.geometry.static_deformation_embedding import embedding_derivability


def test_scalar_information_is_directional_but_not_spatial_geometry():
    q=np.arange(4*5*6,dtype=float).reshape(4,5,6)
    d=directional_scalar_information(q)
    assert d.shape==q.shape+(6,)
    assert reciprocity_error_for_scalar_edges(q)==0
    gate=embedding_derivability(q)
    assert gate["directed_scalar_differences_available"]
    assert not gate["deformed_bond_lengths_derivable"]
    assert not gate["deformed_bond_directions_derivable"]
    assert not gate["global_node_embedding_derivable"]


def test_zero_load_cartesian_n6_and_straight_path_control():
    bonds=undeformed_bond_vectors((2,3,4))
    assert bonds.shape==(2,3,4,6,3)
    assert np.all(np.linalg.norm(bonds,axis=-1)==1)
    diag=path_diagnostics(straight_path((1,2,3),2,8))
    assert diag["maximum_turning_radians"]==0
    assert set(diag["bond_lengths"])=={1.0}


def test_runner_emits_outcome_a_contract_and_required_artifacts(tmp_path,monkeypatch):
    from pbuf.labs.foundation import static_native_deformation_to_relational_geometry001 as lab
    monkeypatch.setattr(lab,"OUT",tmp_path)
    contract=lab.main()
    assert contract["DEV164_AUDIT_COMPLETE"]
    assert contract["DEV163_NULL_COUPLING_RESULT_PRESERVED"]
    assert contract["STATIC_STATE_CONTAINS_DIRECTIONAL_INFORMATION"]=="PARTIAL"
    assert contract["DEFORMED_BOND_DIRECTIONS_DERIVABLE"]=="FALSE"
    assert contract["GLOBAL_NODE_EMBEDDING_DERIVABLE"]=="FALSE"
    assert not contract["RELATIONAL_GEOMETRY_READY_FOR_FINITE_PROPAGATION"]
    assert not contract["FULL_FINITE_NATIVE_LENSING_EXECUTED"]
    required={"report.txt","static_state_information_inventory.json","geometry_candidate_inventory.json",
      "scalar_to_vector_information_audit.json","bond_geometry_contract.json","integrability_audit.json",
      "symmetry_covariance_audit.json","unloaded_geometry_control.json","loaded_path_geometry.json",
      "trajectory_change_audit.json","depth_realization_geometry_stability.json","amplitude_geometry_scaling.json",
      "propagation_geometry_handoff.json","final_relational_geometry_contract.json"}
    assert required<={p.name for p in tmp_path.iterdir()}
    handoff=json.loads((tmp_path/"propagation_geometry_handoff.json").read_text())
    assert handoff["node_positions_zyx"] is None and not handoff["ready_for_finite_propagation"]
