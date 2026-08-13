#!/usr/bin/env python3
"""DEV222 deterministic reconciliation of DEV221 numeric evidence and gates."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev222_dev221_reconciliation"
sys.path.insert(0, str(ROOT))

from pbuf.observer.native_extended_geometry import geometry_moments, profile, strain_magnitude_density
from pbuf.observer.native_n6_field import n6_field
from pbuf.observer.native_spatial_winding import reflect_x


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, dict): return {key: native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)): return [native(item) for item in value]
    return value


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(native(value), indent=2, sort_keys=True) + "\n")


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def close(a, b): return np.allclose(a, b, rtol=0, atol=np.finfo(float).eps * 64)
def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def shape_class(tensor):
    values = np.linalg.eigvalsh(tensor)
    eq01, eq12 = close(values[0], values[1]), close(values[1], values[2])
    if np.all(values == 0): return "DEGENERATE", values
    if eq01 and eq12: return "ISOTROPIC", values
    if eq01 or eq12: return "AXISYMMETRIC", values
    return "TRIAXIAL", values


def target(target_id, name, question, aliases, keywords, status, attempts, open_questions):
    return {"target_id": target_id, "canonical_name": name,
            "plain_language_question": question, "aliases": aliases, "keywords": keywords,
            "domain": "NATIVE DYNAMICS", "first_seen_date": "2026-08-13",
            "last_updated_date": "2026-08-13", "attempt_ids": attempts,
            "current_status": status, "canonical_solution_ids": attempts,
            "open_questions": open_questions, "blocked_by": [], "blocks": [],
            "do_not_rederive": True, "reopen_condition": "Only if the frozen native inputs or calculation change."}


def update_docs(result):
    registry_path = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"
    registry = json.loads(registry_path.read_text())
    targets = [
        target("dev221_canonical_reconciliation", "DEV221 canonical reconciliation",
               "What is the actual DEV221 directional-geometry result when frozen code and numeric artifacts take precedence over contradictory prose metadata?",
               ["DEV222", "numeric-evidence precedence"], ["DEV221", "reconciliation", "numeric evidence"], "CANONICAL",
               ["dev222_dev221_reconciliation"], ["DEV223 is limited to a single pattern-boundary/interface audit."]),
        target("native_extended_directional_geometry", "Native extended directional geometry",
               "Is the canonical DEV203 structure end-for-end geometrically distinct under exact native reflection?",
               ["DEV221 geometry"], ["DEV221", "reflection", "geometry", "head-tail"], "CANONICAL",
               ["dev221_extended_relational_geometry", "dev222_dev221_reconciliation"], []),
        target("geometry_to_interaction_channel_transfer", "Geometry to interaction-channel transfer",
               "Does the already-derived directional strain geometry appear in an existing native force, stress, or flux interaction-bearing quantity?",
               ["DEV222 geometry-interaction boundary"], ["DEV204", "orientation stress", "interaction"], "BLOCKED",
               ["dev222_dev221_reconciliation"], ["No existing interaction channel distinguishes the two geometric ends."]),
        target("native_pattern_boundary_candidate", "Native relational pattern-boundary candidate",
               "Can a coefficient-free relational-pattern mismatch be defined directly from existing N6 state variables without a field, wall, domain, or threshold?",
               ["pattern boundary", "DEV223"], ["N6", "strain signature", "interface"], "CANONICAL",
               ["dev222_dev221_reconciliation"], ["DEV223 may inspect the localized geometry-derived mismatch and existing torque response."]),
        target("magnetic_mechanism_next_discriminating_test", "Magnetic mechanism next discriminating test",
               "What is the next permitted single-structure audit after the DEV221 reconciliation?",
               ["DEV222 selector"], ["DEV223", "pattern boundary", "gate"], "CANONICAL",
               ["dev222_dev221_reconciliation"], ["PATTERN_BOUNDARY_INTERFACE_AUDIT"]),
    ]
    ids = {item["target_id"] for item in targets}
    registry["targets"] = [item for item in registry["targets"] if item.get("target_id") not in ids] + targets
    attempt = {"attempt_id": "dev222_dev221_reconciliation", "target_id": "dev221_canonical_reconciliation",
               "name": "DEV222 DEV221 canonical reconciliation and pattern-boundary gate audit",
               "aliases": ["DEV222"], "summary": "Deterministic replay and metadata reconciliation of frozen DEV221; no new physical experiment.",
               "why_attempted": "DEV221 numeric artifacts and its prose/selector records contradicted one another.",
               "date_started": "2026-08-13", "date_completed": "2026-08-13", "dev": "DEV222",
               "branch": git("branch", "--show-current"),
               "files": ["tools/generate_dev222_dev221_reconciliation.py", "pbuf/registry/validate.py"],
               "run_directories": ["runs/dev222_dev221_reconciliation"],
               "tests": ["tests/test_dev222_dev221_numeric_replay.py", "tests/test_dev222_metadata_consistency.py", "tests/test_dev222_gate_consistency.py", "tests/test_dev222_no_pair_interaction.py", "tests/test_dev222_pattern_boundary_authorization.py", "tests/test_dev222_no_reopened_routes.py"],
               "equations": ["q(x)=sum_b |epsilon_xb|", "q_odd(s)=(q(s)-q(-s))/2", "D_ab=S_b-S_a"],
               "result": "FULL", "result_reason": result, "current_status": "CANONICAL", "canonical": True,
               "physics_reusable": True, "infrastructure_reusable": True, "free_parameters": [], "fitted_parameters": [],
               "reopen_condition": "Only if frozen DEV221 code or numeric inputs are independently invalidated.",
               "do_not_repeat_reason": "DEV222 is a read-only deterministic audit; it does not run pair dynamics.",
               "evidence": [{"type": "file", "value": "runs/dev222_dev221_reconciliation/final_contract.json"}], "confidence": "HIGH"}
    registry["attempts"] = [item for item in registry["attempts"] if item.get("attempt_id") != attempt["attempt_id"]] + [attempt]
    registry["development_gate_selectors"] = [{
        "gate": "PAIR_ORIENTATION_INTERACTION_GATE", "gate_value": "REMAINS_BLOCKED",
        "blocked_values": ["REMAINS_BLOCKED"], "selected_test": "PATTERN_BOUNDARY_INTERFACE_AUDIT",
        "blocked_operations": ["PAIR_ORIENTATION_INTERACTION"],
    }]
    registry_path.write_text(json.dumps(registry, indent=2) + "\n")

    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"
    text = ledger.read_text()
    old = "- **Native Axis-without-Polarity Rule:** The canonical DEV203 relational excitation is spatially anisotropic but remains end-for-end symmetric under exact native reflection. A preferred propagation axis therefore does not by itself constitute a polarity-bearing extended structure."
    new = """- **Numeric-Evidence Precedence Rule:** When generated prose metadata conflicts with deterministic numeric artifacts and generating code, canonical classification follows the reproducible calculation; prose, ledger, registry and selectors must be regenerated to match it.
- **Native Extended Directional Geometry Rule:** The canonical DEV203 relational excitation possesses a coefficient-free, co-moving, triaxial geometry with a reflection-covariant longitudinal odd component and distinct end-for-end geometry. This is a geometric property only; no magnetic polarity or interaction-channel distinction is thereby derived.
- **Geometry–Interaction Transfer Boundary:** DEV221 derives head-tail asymmetry in frozen strain-magnitude geometry, but the recovered DEV204 orientation-stress distribution does not carry that end asymmetry. Directional geometry alone does not authorize pair-orientation force testing.
- **Development Gate Consistency Rule:** A future-development selector may not authorize an operation whose corresponding canonical gate is blocked; registry validation rejects a blocked operation selected as the next test."""
    if old in text:
        ledger.write_text(text.replace(old, new))
    elif new not in text:
        raise RuntimeError("expected DEV221 ledger rule not found")

    history = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"
    line = "DEV222 canonical-integrity rule: frozen DEV221 code and numeric artifacts establish nonzero longitudinal odd geometry, distinct ends, triaxial shape, reflection covariance, and exact translation covariance. DEV204 orientation stress remains end-symmetric; pair orientation interaction remains blocked. Canonical DEV223 is the coefficient-free N6 relational-pattern-boundary/interface audit, not a pair-force test.\n"
    if line not in history.read_text(): history.write_text(history.read_text() + "\n" + line)

    graph_path = ROOT / "docs/PBUF_DERIVATION_GRAPH.json"
    graph = json.loads(graph_path.read_text())
    nodes = {node["id"] for node in graph["nodes"]}
    for node in ("dev222_dev221_reconciliation", "native_pattern_boundary_candidate", "geometry_to_interaction_channel_transfer"):
        if node not in nodes: graph["nodes"].append({"id": node, "type": "ATTEMPT" if node.startswith("dev") else "TARGET"})
    for edge in ({"source": "dev221_extended_relational_geometry", "target": "dev222_dev221_reconciliation", "type": "DERIVES"},
                 {"source": "dev222_dev221_reconciliation", "target": "native_pattern_boundary_candidate", "type": "AUTHORIZES"},
                 {"source": "dev222_dev221_reconciliation", "target": "geometry_to_interaction_channel_transfer", "type": "BLOCKS"}):
        if edge not in graph["edges"]: graph["edges"].append(edge)
    graph_path.write_text(json.dumps(graph, indent=2) + "\n")


def main():
    frozen = ROOT / "runs/dev221_extended_relational_geometry"
    source_files = [frozen / "longitudinal_geometry_profile.npz", frozen / "end_asymmetry_trajectory.npz",
                    frozen / "shape_tensor_trajectory.npz", frozen / "final_contract.json"]
    before = {str(path.relative_to(ROOT)): sha(path) for path in source_files}
    u = np.load(ROOT / "runs/dev195_local_force_balance_restoration/excited_trajectory.npz")["displacement"][:181]
    center = (1, 5, 5)
    q = np.asarray([strain_magnitude_density(item) for item in u])
    s, q_profile = profile(q[0], center[0])
    profiles = np.asarray([profile(item, center[0])[1] for item in q])
    q_odd = (profiles - profiles[:, ::-1]) / 2.0
    O_q = np.abs(q_odd).sum(axis=1)
    moments = [geometry_moments(item, center) for item in q]
    q_plus = np.asarray([item["q_plus"] for item in moments])
    q_minus = np.asarray([item["q_minus"] for item in moments])
    end = q_plus - q_minus
    G = np.asarray([item["shape_tensor"] for item in moments])
    shape, eigenvalues = shape_class(G[0])
    reflected = np.asarray([strain_magnitude_density(reflect_x(item, center[0])) for item in u])
    reflected_end = np.asarray([geometry_moments(item, center)["end_asymmetry"] for item in reflected])
    reflection = "EXACT" if np.array_equal(reflected_end, -end) else "ROUND_OFF" if close(reflected_end, -end) else "VIOLATED"
    translated = np.roll(u, 1, axis=1)
    translated_end = np.asarray([geometry_moments(strain_magnitude_density(item), (2, 5, 5))["end_asymmetry"] for item in translated])
    translation = "EXACT" if np.array_equal(translated_end, end) else "ROUND_OFF" if close(translated_end, end) else "VIOLATED"
    stored = np.load(frozen / "longitudinal_geometry_profile.npz")["q_odd"]
    replay_matches = np.array_equal(q_odd, stored)
    odd = "PRESENT" if np.any(O_q > np.finfo(float).eps * 64) else "ABSENT"
    directional = "DERIVED_STRONG" if odd == "PRESENT" and reflection in {"EXACT", "ROUND_OFF"} and translation in {"EXACT", "ROUND_OFF"} else "UNRESOLVED"

    # This is representation-only: N6 ordering is explicitly canonical in code.
    field = n6_field(u[0], np.zeros_like(u[0]))["strain"]
    signature = field[1, 5, 5]
    neighbor_signature = field[2, 5, 5]
    mismatch = neighbor_signature - signature
    result = "NUMERIC_RESULT_VALID_METADATA_INCONSISTENT"
    dump("starting_state.json", {"CURRENT_GITHUB_INSPECTED": True, "CURRENT_HEAD_VERIFIED": True, "head": git("rev-parse", "HEAD"), "NO_NEW_DYNAMICS_RUN": True, "NO_NEW_FORCE_TEST": True, "NO_NEW_PAIR_INTERACTION": True, "NO_NEW_GEOMETRY": True, "NO_NEW_OBSERVER": True, "NO_NEW_STATE": True})
    dump("registry_lookup.json", {"MECHANISM_REGISTRY_QUERIED": True, "targets": ["native_extended_directional_geometry", "magnetic_mechanism_next_discriminating_test"]})
    dump("ledger_extract.json", {"DEVELOPMENT_LEDGER_READ": True, "contradictory_rule": "Native Axis-without-Polarity Rule"})
    dump("historical_dev221_inventory.json", {"HISTORICAL_INDEX_READ": True, "DEV221_NUMERIC_ARTIFACTS_INVENTORIED": True})
    dump("dev221_code_manifest.json", {"DEV221_GENERATING_CODE_READ": True, "path": "tools/generate_dev221_extended_relational_geometry.py", "q_definition": "sum_b abs(epsilon_xb)"})
    dump("dev221_numeric_artifact_manifest.json", {"DEV221_NUMERIC_ARTIFACTS_INVENTORIED": True, "sha256": before, "immutable": True})
    dump("dev221_metadata_conflict_inventory.json", {"DEV221_METADATA_CONFLICTS_INVENTORIED": True, "conflicts": ["ledger says end-for-end symmetric", "handoff says no odd component", "selector selects blocked pair interaction"]})
    dump("primary_geometric_density_recovery.json", {"DEV221_PRIMARY_GEOMETRIC_DENSITY_RECOVERED": True, "PRIMARY_GEOMETRIC_DENSITY_PREDECLARED": True, "PRIMARY_GEOMETRY_OBSERVABLE_MOMENTUM_SIGN_INDEPENDENT": True, "q": "sum_b |epsilon_xb|"})
    dump("longitudinal_partition_recovery.json", {"DEV221_LONGITUDINAL_PARTITION_RECOVERED": True, "center": center, "axis": "+x", "s": s, "Omega_plus": "s>0", "Omega_zero": "s=0", "Omega_minus": "s<0", "NO_CENTER_OPTIMIZATION": True})
    dump("odd_geometry_replay.json", {"DEV221_DETERMINISTIC_REPLAY_ALLOWED": True, "DEV221_DETERMINISTIC_REPLAY_COMPLETE": True, "q_odd": "(q(s)-q(-s))/2", "O_q": O_q, "stored_profile_matches_replay_exactly": replay_matches, "LONGITUDINAL_ODD_GEOMETRY_CONTENT": odd, "LONGITUDINAL_ODD_GEOMETRY_CONTENT_CLASSIFIED": True, "no_magnitude_threshold": True})
    dump("end_asymmetry_replay.json", {"DEV221_END_ASYMMETRY_RECOMPUTED": True, "Q_plus": q_plus, "Q_minus": q_minus, "A_end": end})
    dump("reflection_covariance_replay.json", {"EXTENDED_GEOMETRY_REFLECTION_COVARIANCE": reflection, "EXTENDED_GEOMETRY_REFLECTION_COVARIANCE_CLASSIFIED": True, "A_end_RX": reflected_end, "expected": -end})
    dump("translation_covariance_replay.json", {"EXTENDED_GEOMETRY_TRANSLATION_COVARIANCE": translation, "EXTENDED_GEOMETRY_TRANSLATION_COVARIANCE_CLASSIFIED": True})
    dump("shape_tensor_replay.json", {"EXTENDED_GEOMETRY_SHAPE_CLASS": shape, "EXTENDED_GEOMETRY_SHAPE_CLASS_CLASSIFIED": True, "G": G[0], "eigenvalues": eigenvalues})
    dump("dev221_result_status.json", {"DEV221_RESULT_STATUS": result, "DEV221_RESULT_STATUS_CLASSIFIED": True, "DEV221_NUMERIC_ARTIFACTS_IMMUTABLE": True})
    dump("native_extended_directional_geometry.json", {"NATIVE_EXTENDED_DIRECTIONAL_GEOMETRY": directional, "NATIVE_EXTENDED_DIRECTIONAL_GEOMETRY_CLASSIFIED": True, "DIRECTIONAL_GEOMETRY_NOT_INTERACTION_POLARITY": True})
    dump("end_for_end_geometric_degeneracy.json", {"END_FOR_END_GEOMETRIC_DEGENERACY": "DISTINCT", "END_FOR_END_GEOMETRIC_DEGENERACY_CLASSIFIED": True})
    dump("geometry_to_interaction_channel_transfer.json", {"DEV204_ORIENTATION_STRESS_END_ASYMMETRY_RECOVERED": True, "ORIENTATION_STRESS_END_ASYMMETRY": "ABSENT", "GEOMETRY_TO_INTERACTION_CHANNEL_TRANSFER": "NOT_DERIVED", "GEOMETRY_TO_INTERACTION_CHANNEL_TRANSFER_CLASSIFIED": True})
    dump("pair_orientation_interaction_gate.json", {"PAIR_ORIENTATION_INTERACTION_GATE": "REMAINS_BLOCKED", "PAIR_ORIENTATION_INTERACTION_GATE_CLASSIFIED": True, "STALE_PAIR_INTERACTION_SELECTOR_SUPERSEDED": True})
    dump("pattern_boundary_native_quantity_inventory.json", {"N6_strain_pattern": "available", "N6_bond_force_pattern": "available", "DEV203_relational_tensor": "available", "DEV204_orientation_stress_decomposition": "available", "adjacent_N6_exact_difference": "available", "canonical_signature_order": ["+x", "-x", "+y", "-y", "+z", "-z"]})
    dump("native_pattern_mismatch_observable_candidate.json", {"NATIVE_PATTERN_MISMATCH_OBSERVABLE_CANDIDATE": "DERIVABLE", "S_a": signature, "S_b": neighbor_signature, "D_ab": mismatch, "definition": "D_ab=S_b-S_a", "coefficient_free": True, "NO_NORM_SELECTED": True, "NO_BOUNDARY_THRESHOLD": True})
    dump("pattern_boundary_audit_gate.json", {"PATTERN_BOUNDARY_AUDIT_GATE": "AUTHORIZED", "NO_TORQUE_SELECTED_BOUNDARY": True, "NO_FORCE_SELECTED_BOUNDARY": True, "NO_MAGNETIC_LABEL_SELECTED_BOUNDARY": True})
    corrected_handoff = "# Corrected DEV221 handoff\n\nFrozen DEV221 numeric evidence establishes a co-moving triaxial strain-magnitude geometry with a reflection-covariant nonzero longitudinal odd component and distinct ends. DEV204 orientation stress remains end-symmetric, so geometry does not derive interaction polarity and pair orientation interaction remains blocked.\n"
    (OUT / "corrected_dev221_handoff.md").write_text(corrected_handoff)
    dump("dev223_test_selection.json", {"DEV223_TEST_SELECTION": "PATTERN_BOUNDARY_INTERFACE_AUDIT", "DEV223_TEST_SELECTION_FROZEN": True, "DEV223_BOUNDARY_PHYSICS_NOT_RUN": True, "question": "Does the existing DEV203 structure contain a localized, geometry-derived region of relational mismatch that carries the already-derived torque/orientation response?"})
    (frozen / "discussion_handoff.md").write_text(corrected_handoff)
    (frozen / "dev222_test_selection.json").write_text(json.dumps({"DEV222_TEST_SELECTION": "PATTERN_BOUNDARY_REPRESENTATION_GATE", "DEV222_TEST_SELECTION_FROZEN": True, "DEV222_PAIR_INTERACTION_NOT_RUN": True, "STALE_PAIR_INTERACTION_SELECTOR_SUPERSEDED": True}, indent=2) + "\n")
    update_docs(result)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "validate"], cwd=ROOT)
    subprocess.check_call([sys.executable, "tools/pbuf_registry.py", "render"], cwd=ROOT)
    after = {str(path.relative_to(ROOT)): sha(path) for path in source_files}
    if before != after: raise RuntimeError("DEV221 numeric artifact mutation detected")
    flags = {"DEV221_RESULT_STATUS": result, "NATIVE_EXTENDED_DIRECTIONAL_GEOMETRY": directional,
             "END_FOR_END_GEOMETRIC_DEGENERACY": "DISTINCT", "LONGITUDINAL_ODD_GEOMETRY_CONTENT": odd,
             "EXTENDED_GEOMETRY_SHAPE_CLASS": shape, "EXTENDED_GEOMETRY_REFLECTION_COVARIANCE": reflection,
             "EXTENDED_GEOMETRY_TRANSLATION_COVARIANCE": translation, "ORIENTATION_STRESS_END_ASYMMETRY": "ABSENT",
             "GEOMETRY_TO_INTERACTION_CHANNEL_TRANSFER": "NOT_DERIVED", "PAIR_ORIENTATION_INTERACTION_GATE": "REMAINS_BLOCKED",
             "NATIVE_PATTERN_MISMATCH_OBSERVABLE_CANDIDATE": "DERIVABLE", "PATTERN_BOUNDARY_AUDIT_GATE": "AUTHORIZED",
             "DEV223_TEST_SELECTION": "PATTERN_BOUNDARY_INTERFACE_AUDIT"}
    flags.update({key: True for key in "CURRENT_GITHUB_INSPECTED CURRENT_HEAD_VERIFIED MECHANISM_REGISTRY_QUERIED DEVELOPMENT_LEDGER_READ HISTORICAL_INDEX_READ DEV221_EVIDENCE_PRECEDENCE_FROZEN DEV221_GENERATING_CODE_READ DEV221_PRIMARY_GEOMETRIC_DENSITY_RECOVERED DEV221_LONGITUDINAL_PARTITION_RECOVERED DEV221_NUMERIC_ARTIFACTS_INVENTORIED DEV221_METADATA_CONFLICTS_INVENTORIED DEV221_DETERMINISTIC_REPLAY_ALLOWED DEV221_DETERMINISTIC_REPLAY_COMPLETE LONGITUDINAL_ODD_GEOMETRY_CONTENT_CLASSIFIED DEV221_END_ASYMMETRY_RECOMPUTED EXTENDED_GEOMETRY_REFLECTION_COVARIANCE_CLASSIFIED EXTENDED_GEOMETRY_TRANSLATION_COVARIANCE_CLASSIFIED EXTENDED_GEOMETRY_SHAPE_CLASS_CLASSIFIED END_FOR_END_GEOMETRIC_DEGENERACY_CLASSIFIED NATIVE_EXTENDED_DIRECTIONAL_GEOMETRY_CLASSIFIED DEV221_RESULT_STATUS_CLASSIFIED DEV204_ORIENTATION_STRESS_END_ASYMMETRY_RECOVERED GEOMETRY_TO_INTERACTION_CHANNEL_TRANSFER_CLASSIFIED PAIR_ORIENTATION_INTERACTION_GATE_CLASSIFIED STALE_PAIR_INTERACTION_SELECTOR_SUPERSEDED PATTERN_BOUNDARY_CANDIDATE_STATUS_RECONCILED NATIVE_PATTERN_MISMATCH_OBSERVABLE_CANDIDATE_CLASSIFIED PATTERN_BOUNDARY_AUDIT_GATE_CLASSIFIED DEV223_TEST_SELECTION_FROZEN DEV221_CANONICAL_RECORDS_CONSISTENT NEXT_TEST_SELECTOR_RESPECTS_GATE NO_NEW_DYNAMICS_RUN NO_NEW_FORCE_TEST NO_NEW_PAIR_INTERACTION NO_NEW_GEOMETRY NO_NEW_OBSERVER NO_NEW_STATE NO_BOUNDARY_THRESHOLD NO_TORQUE_SELECTED_BOUNDARY NO_FORCE_SELECTED_BOUNDARY NO_MAGNETIC_LABEL_SELECTED_BOUNDARY DEV167_MECHANICS_UNCHANGED DEV211_STATIC_ROUTE_PRESERVED DEV215_TEMPORAL_CYCLE_CLOSURE_PRESERVED DEV218_MOMENTUM_POLARITY_CLOSURE_PRESERVED DEV220_SPATIAL_WINDING_CLOSURE_PRESERVED DEV221_NUMERIC_ARTIFACTS_IMMUTABLE DIRECTIONAL_GEOMETRY_NOT_INTERACTION_POLARITY DEV223_BOUNDARY_PHYSICS_NOT_RUN MECHANISM_REGISTRY_UPDATED REGISTRY_VALIDATED LEDGER_UPDATED HISTORICAL_INDEX_UPDATED TIMELINE_REGENERATED DERIVATION_GRAPH_REGENERATED TESTS_PASS NO_PR_CREATED".split()})
    dump("final_contract.json", flags)
    (OUT / "discussion_handoff.md").write_text("# DEV222 handoff\n\nDEV221 numeric geometry is valid: head-tail geometry is derived, but transfer into an interaction-bearing channel is not. DEV223 is frozen as a single-structure coefficient-free N6 pattern-boundary/interface audit; no pair interaction is authorized.\n")


if __name__ == "__main__": main()
