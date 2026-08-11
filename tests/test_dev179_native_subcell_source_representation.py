import json
import numpy as np
import pytest

from pbuf.excitation.native_vector_pair_dynamics import source_contact_force
from pbuf.labs.foundation.native_subcell_geometry_dev179 import (
    SubcellSourceRepresentationNotDerived, node_contact_loading,
)


def test_exact_node_recovery_and_zero_source_are_the_existing_contact_law():
    got = node_contact_loading((11, 11, 11), (5., 5., 5.), .02)
    assert np.array_equal(got, source_contact_force((11, 11, 11), (5, 5, 5), .02))
    assert not np.any(node_contact_loading((11, 11, 11), (5., 5., 5.), 0.0))


@pytest.mark.parametrize("position", [(5., 5.5, 5.), (5., 5.25, 5.75)])
def test_subcell_coordinate_has_no_frozen_source_to_medium_mapping(position):
    with pytest.raises(SubcellSourceRepresentationNotDerived):
        node_contact_loading((11, 11, 11), position, .02)


def test_dev179_artifact_contract(tmp_path, monkeypatch):
    from tools import generate_dev179_native_subcell_source_representation as lab
    monkeypatch.setattr(lab, "OUT", tmp_path)
    contract = lab.main()
    assert contract["OUTCOME"] == "OUTCOME_D"
    assert contract["SOURCE_POSITION_NATIVE_DISCRETE"]
    assert not contract["SUBCELL_SOURCE_REPRESENTATION_DERIVED"]
    required = {"starting_state.json", "source_lattice_semantics_audit.json", "source_position_semantics.json",
      "native_discretization_roles.json", "historical_launch_architecture.json", "current_launch_blocker.json",
      "candidate_inventory.json", "candidate_derivation_basis.json", "synthetic_fixture_manifest.json", "node_recovery.json",
      "conservation_tests.json", "translation_path.json", "node_crossing_test.json", "reflection_test.json",
      "lattice_symmetry_test.json", "source_lineage_test.json", "receipt_lineage_test.json", "subcell_information_test.json",
      "refinement_convergence.json", "j3_subcell_diagnostic.json", "g3_subcell_diagnostic.json", "cpu_vulkan_parity.json",
      "viewer_extension_status.json", "candidate_status_matrix.json", "final_contract.json", "discussion_handoff.md"}
    assert required <= {p.name for p in tmp_path.iterdir()}
    assert json.loads((tmp_path / "current_launch_blocker.json").read_text())["INTERPOLATION_IS_NOT_PHYSICS"]
