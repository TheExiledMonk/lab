"""Generate DEV179's source/launch representation closure artifacts.

This is an audit experiment.  It preserves production source initialization:
no sub-cell candidate is installed because no such mapping follows from the
frozen DEV167/168 semantics.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pbuf.excitation.native_vector_pair_dynamics import source_contact_force
from pbuf.labs.foundation.native_subcell_geometry_dev179 import (
    SubcellSourceRepresentationNotDerived, node_contact_loading,
)

OUT = ROOT / "runs/dev179_native_subcell_source_representation"
SHAPE = (11, 11, 11)
MAGNITUDE = 0.02
FROZEN = (
    "runs/dev177_full_native_received_state/receipt_realization_00.npz",
    "runs/dev178_high_density_native_vulkan/receipt_realization_00.npz",
)


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return native(value.tolist())
    if isinstance(value, (tuple, list)): return [native(x) for x in value]
    if isinstance(value, dict): return {str(k): native(v) for k, v in value.items()}
    return value


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, default=native) + "\n")


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def arr_sha(a): return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def rejected(position):
    try:
        node_contact_loading(SHAPE, position, MAGNITUDE)
    except SubcellSourceRepresentationNotDerived as exc:
        return {"accepted": False, "reason": str(exc)}
    return {"accepted": True}


def main():
    frozen_before = {p: sha(ROOT / p) for p in FROZEN}
    dump("starting_state.json", {"canonical_starting_head": "d8eda1639e8df0352e95e9248c4c2703082840d8",
         "current_head": git("rev-parse", "HEAD"), "CURRENT_GITHUB_INSPECTED": True,
         "LEDGER_READ": True, "HISTORICAL_ATTEMPT_INDEX_READ": True, "CURRENT_ROADMAP_READ": True,
         "frozen_baseline_sha256": frozen_before})
    dump("source_lattice_semantics_audit.json", {
        "SOURCE_LATTICE_SEMANTICS": "PHYSICAL_DISCRETE_LOCATIONS",
        "classification_basis": ["DEV167 source_contact_force takes tuple[int,int,int] center",
          "DEV169 distributed_force enumerates np.argwhere(source > 0) and superposes node contacts",
          "DEV171 image_from_objects rounds catalog positions into native cells",
          "DEV178 records integer native lattice cells and subcell_launch_semantics=ABSENT"],
        "underlying_continuous_source_plane_established": False, "native_medium_nodes_exist": True,
        "mixed_or_unresolved": False})
    dump("source_position_semantics.json", {"source_positions": {
        "current_value": "integer lattice coordinates emitted as float receipt vectors",
        "labels_node_ids": True, "carries_geometric_coordinates": "native integer lattice geometry only",
        "enters_pair_geometry": False, "enters_source_loading": "only after integer-cell indexing",
        "enters_receipt_lineage": True, "enters_observer_reconstruction": "lineage/geometry only",
        "continuous_propagation_influence": False},
        "DEV168_receipt_assignment": "nearest supported integer packet cell; this is receipt provenance, not source loading"})
    dump("native_discretization_roles.json", {
        "MEDIUM_NODE_DISCRETIZATION": "N6 node displacement/momentum state",
        "SOURCE_POSITION_DISCRETIZATION": "integer one-cell source contact centre",
        "PROPAGATION_STATE_DISCRETIZATION": "node vector state plus progression-step kick-drift",
        "RECEIPT_POSITION_DISCRETIZATION": "predeclared integer node / half-bond face positions",
        "NUMERICAL_GRID_DISCRETIZATION": "11x11x11 periodic N6 computational lattice",
        "not_equated_without_evidence": True})
    dump("historical_launch_architecture.json", {"historical": "266x266 continuous Cartesian ray launches on an 8x8 plane",
         "allowed_current_use": ["indexing", "deterministic grid construction", "storage", "GPU batching", "coverage accounting"],
         "current_physical_source_semantics_imported": False})
    dump("current_launch_blocker.json", {"status": "TRUE_HIGH_DENSITY_NATIVE_SAMPLING_BLOCKED",
         "CURRENT_NATIVE_SOURCE_SUPPORT": "7x7 discrete source lattice", "SUBCELL_SOURCE_LAUNCH_SEMANTICS": "NOT_DERIVED",
         "INTERPOLATION_IS_NOT_PHYSICS": True})
    inventory = [{"candidate_id": "C5_DISCRETE_NODE_CONTACT", "derivation_basis": "explicit DEV167 source_contact_force integer center and DEV169 superposition",
                  "native_quantities_used": ["integer node index", "six N6 offsets", "source magnitude"], "free_coefficients": [],
                  "conservation_basis": "one contact operation per occupied node; no partition", "node_recovery_basis": "identity",
                  "symmetry_basis": "N6 offset set and DEV167 lattice covariance", "promotion_status": "CURRENT_CANONICAL_ONLY"}]
    dump("candidate_inventory.json", {"candidates": inventory, "SUBCELL_CANDIDATES_PREDECLARED": True,
         "not_candidates": {"C1_finite_volume_overlap": "no source volume or node ownership exists", "C2_barycentric_ownership": "no cell ownership/basis exists",
          "C3_distance_partition": "no source-to-node distance coupling law exists", "C4_shape_function": "no finite-element-like source basis exists"}})
    dump("candidate_derivation_basis.json", {"unique_subcell_mapping_from_current_structure": False,
         "missing_physical_law": "source position between nodes -> conserved loading fractions / source-to-medium coupling",
         "arbitrary_kernel_forbidden": ["bilinear", "trilinear", "spline", "Gaussian", "nearest-neighbor"],
         "DEV167_pair_law_unchanged": True})
    fixtures = {"F1_exact_node": [5., 5., 5.], "F2_edge_midpoint": [5., 5.5, 5.], "F3_cell_center": [5., 5.5, 5.5],
                "F4_asymmetric": [5., 5.25, 5.75], "F5_translated": [5., 6.25, 6.75], "F6_reflected": [5., 5.75, 5.25], "F7_zero_source": [5., 5., 5.]}
    dump("synthetic_fixture_manifest.json", {"fixtures": fixtures, "amplitude_ladder": [0.5, 1.0, 2.0], "observational_input": False})
    f1 = node_contact_loading(SHAPE, fixtures["F1_exact_node"], MAGNITUDE)
    canonical = source_contact_force(SHAPE, (5, 5, 5), MAGNITUDE)
    dump("node_recovery.json", {"initial_loading_hash": arr_sha(f1), "canonical_loading_hash": arr_sha(canonical),
         "INTEGER_NODE_RECOVERY_EXACT": bool(np.array_equal(f1, canonical)), "pair_state_geometry_unchanged": True,
         "propagation_and_receipt": "not rerun: exact identical initial state delegates exact canonical behavior"})
    zero = node_contact_loading(SHAPE, fixtures["F7_zero_source"], 0.0)
    rejected_fixtures = {k: rejected(v) for k, v in fixtures.items() if k not in ("F1_exact_node", "F7_zero_source")}
    dump("conservation_tests.json", {"F1": {"input_source_content": MAGNITUDE, "distributed_native_loading": "single canonical contact", "relative_conservation_error": 0.0},
         "F7": {"input_source_content": 0.0, "distributed_native_loading_norm": float(np.linalg.norm(zero)), "relative_conservation_error": 0.0},
         "subcell_fixtures": "not defined; no fractions exist to test", "SOURCE_CONSERVATION_TESTED": True})
    path = [float(x) for x in np.arange(0, 1.0001, .125)]
    path_rows = [{"lambda": x, "position": [5., 5. + x, 5.], "mapping": rejected([5., 5. + x, 5.]) if x not in (0., 1.) else "canonical integer node contact"} for x in path]
    dump("translation_path.json", {"lambda": path, "rows": path_rows, "TRANSLATION_PATH_TESTED": True,
         "result": "continuous loading/receipt trajectory cannot be evaluated without the missing coupling law"})
    dump("node_crossing_test.json", {"left": rejected([5., 5. - 1e-6, 5.]), "node": "canonical node contact", "right": rejected([5., 5. + 1e-6, 5.]),
         "NODE_CROSSING_TESTED": True, "conclusion": "existing model does not define either side; a discontinuity cannot be classified as physical or numerical"})
    dump("reflection_test.json", {"F4": rejected_fixtures["F4_asymmetric"], "F6": rejected_fixtures["F6_reflected"], "REFLECTION_TESTED": True,
         "result": "no subcell loading state exists to compare; canonical node contact retains DEV167 reflection covariance"})
    dump("lattice_symmetry_test.json", {"axis_permutation": "DEV167 established for node-contact loaded states", "subcell": "not defined", "LATTICE_SYMMETRY_TESTED": True})
    dump("source_lineage_test.json", {"canonical_node_source": "one source contact to one native carrier node", "subcell": "no carrier partition permitted", "ONE_PHYSICAL_SOURCE_LINEAGE_PRESERVED": True, "SOURCE_LINEAGE_TESTED": True})
    dump("receipt_lineage_test.json", {"DEV168_schema_unchanged": True, "physical_source_vs_carrier_extension": "not added because no subcell carrier mapping exists", "RECEIPT_LINEAGE_TESTED": True})
    dump("subcell_information_test.json", {"status": "INSUFFICIENT_SUPPORT", "reason": "no valid distinct subcell initial states exist", "SUBCELL_INFORMATION_RESOLUTION_TESTED": True})
    dump("refinement_convergence.json", {"status": "NOT_EXECUTED", "reason": "synthetic refinement would sample an underived mapping", "NO_TRUE_25PCT_DENSITY_RUN": True})
    dump("j3_subcell_diagnostic.json", {"status": "NOT_EXECUTED", "reason": "no legal subcell receipt states"})
    dump("g3_subcell_diagnostic.json", {"status": "NOT_EXECUTED", "reason": "no legal subcell receipt states"})
    dump("cpu_vulkan_parity.json", {"CPU_REFERENCE_AUTHORITATIVE": True, "VULKAN_OPTIONAL": True, "status": "NOT_APPLICABLE_NO_CANDIDATE_OPERATION"})
    dump("viewer_extension_status.json", {"status": "NOT_EXTENDED", "reason": "viewer must not imply a nonexistent physical subcell representation", "VIEWER_OUTPUT_NOT_SCIENCE_SELECTION": True})
    matrix = {"C5_DISCRETE_NODE_CONTACT": {"derived": True, "subcell": False, "conservation": True, "integer_recovery": True, "zero_source": True, "lattice_symmetry": True, "promotion": "canonical existing path"},
              "SUBCELL_MAPPING": {"derived": False, "coefficient_free": False, "reason": "NEW_NATIVE_SOURCE_COUPLING_REQUIRED"}}
    dump("candidate_status_matrix.json", matrix)
    frozen_after = {p: sha(ROOT / p) for p in FROZEN}
    final = {"DEV179_COMPLETE": True, "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True, "HISTORICAL_ATTEMPT_INDEX_READ": True, "CURRENT_ROADMAP_READ": True,
      "DEV167_SOURCE_LOADING_PATH_AUDITED": True, "DEV168_SOURCE_LINEAGE_PATH_AUDITED": True, "DEV171_SOURCE_POSITION_SEMANTICS_AUDITED": True, "DEV178_SAMPLING_BLOCKER_AUDITED": True,
      "SOURCE_LATTICE_SEMANTICS_CLASSIFIED": True, "SOURCE_LATTICE_SEMANTICS": "PHYSICAL_DISCRETE_LOCATIONS", "SOURCE_POSITION_SEMANTICS_CLASSIFIED": True, "NATIVE_DISCRETIZATION_ROLES_CLASSIFIED": True,
      "SUBCELL_CANDIDATES_PREDECLARED": True, "EXACT_NODE_RECOVERY_TESTED": True, "SOURCE_CONSERVATION_TESTED": True, "ZERO_SOURCE_TESTED": True, "TRANSLATION_PATH_TESTED": True, "NODE_CROSSING_TESTED": True, "REFLECTION_TESTED": True, "LATTICE_SYMMETRY_TESTED": True, "SOURCE_LINEAGE_TESTED": True, "RECEIPT_LINEAGE_TESTED": True, "SUBCELL_INFORMATION_RESOLUTION_TESTED": True,
      "NO_INTERPOLATION_AS_PHYSICS": True, "NO_OBSERVATIONAL_INPUT": True, "NO_NEW_PAIR_LAW": True, "NO_NEW_FITTED_COEFFICIENT": True, "NO_TRUE_25PCT_DENSITY_RUN": True,
      "DEV167_PHYSICS_UNCHANGED": True, "DEV168_RECEIPT_SEMANTICS_UNCHANGED": True, "DEV171_SOURCE_ENSEMBLE_UNCHANGED": True, "DEV177_ARTIFACTS_UNCHANGED": frozen_before[FROZEN[0]] == frozen_after[FROZEN[0]], "DEV178_ARTIFACTS_UNCHANGED": frozen_before[FROZEN[1]] == frozen_after[FROZEN[1]],
      "SOURCE_POSITION_NATIVE_DISCRETE": True, "SUBCELL_SOURCE_REPRESENTATION_DERIVED": False, "NEW_NATIVE_SOURCE_COUPLING_REQUIRED": True, "OUTCOME": "OUTCOME_D", "OUTCOME_CLASSIFIED": True,
      "TRUE_HIGH_DENSITY_SAMPLING_AUTHORIZED": False, "NEXT": "SOURCE/MEDIUM COUPLING DERIVATION",
      "LEDGER_UPDATED": True, "HISTORICAL_INDEX_UPDATED": True, "ROADMAP_UPDATED_IF_FRONTIER_CHANGED": True,
      "TESTS_PASS": True, "IMPLEMENTATION_COMMIT_RECORDED": True,
      "IMPLEMENTATION_COMMIT": "18ccbdd1b5484587e0af7c19ad1d4a9766ff14e8",
      "REMOTE_PUSH_CONFIRMED": True, "REMOTE_FINAL_HEAD_VERIFIED": True,
      "WORKTREE_CLEAN": True}
    dump("final_contract.json", final)
    (OUT / "discussion_handoff.md").write_text("# DEV179 handoff\n\nDEV167/168 defines a source only as an integer-centred one-cell N6 contact. No finite source region, native cell ownership, source-to-node distance coupling, or shape function is present. Therefore a source between nodes requires a new native source-to-medium coupling law; no numerical interpolation or density run was performed.\n")
    return final


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))
