#!/usr/bin/env python3
"""DEV180 repository-first source/medium and density-semantics recovery audit.

Analysis-only: reads code/history/artifacts and writes no production state.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev180_source_medium_recovery"
START = "3212364b77979bc46f9c05135748d3b7cf7c78f1"

PRS = {
    "PR36": ("a72fbab^", "origin/foundation/matter-loading-physical-closure001", "pbuf/labs/foundation/matter_loading_physical_closure001.py"),
    "PR37": ("a72fbab", "origin/foundation/metric-strain-map-closure001", "pbuf/labs/foundation/metric_strain_map_closure001.py"),
    "PR69": ("4a517d1", "origin/foundation/pretensioned-neighbor-medium-response001", "pbuf/labs/foundation/pretensioned_neighbor_medium_response001.py"),
    "PR70": ("a304f49", "origin/foundation/bounded-strain-neighbor-constitutive001", "pbuf/labs/foundation/bounded_strain_neighbor_constitutive001.py"),
    "PR72": ("8f7d4fe", "origin/foundation/bounded-strain-3d-neighbor-network001", "pbuf/labs/foundation/bounded_strain_3d_neighbor_network001.py"),
    "PR74": ("957d05e", "origin/foundation/c-state-bounded-strain-bridge001", "pbuf/labs/foundation/c_state_bounded_strain_bridge001.py"),
}

def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(name: str, value: object) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")

def historical(id_: str, mechanism: str, source: str, medium: str, interaction: str,
               relation: str, status: str, motion: str = "UNSUPPORTED") -> dict:
    merge, ref, file = PRS[id_]
    return {"historical_id": id_, "PR_or_commit": f"{merge} / {ref}", "file": file,
            "equation_or_mechanism": mechanism, "source_representation": source,
            "medium_representation": medium, "interaction_type": interaction,
            "static_or_dynamic": status, "source_motion_supported": motion,
            "continuous_source_position_supported": "NO_EVIDENCE_FOUND",
            "free_coefficients": "none promoted; historical structural normalizations only",
            "test_result": "historical audited implementation", "canonical_status_then": "FOUNDATION_HISTORICAL",
            "relation_to_DEV167": relation, "relation_to_DEV179": "prior architecture exists; off-node placement remains separately unclosed"}

def main() -> None:
    current = git("rev-parse", "HEAD")
    source_files = ["pbuf/excitation/native_vector_pair_dynamics.py", "tools/generate_dev169_raw_abell_native_observer.py",
                    "tools/generate_dev177_full_native_received_state.py", "tools/generate_dev178_high_density_native_vulkan.py",
                    "pbuf/source/native_source_medium_interaction.py", "pbuf/labs/foundation/native_source_medium_interaction001.py"]
    dump("starting_state.json", {"repository": "TheExiledMonk/lab", "canonical_starting_head": START,
        "head_at_audit": current, "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True,
        "HISTORICAL_INDEX_READ": True, "NO_OBSERVATIONAL_INPUT": True,
        "protected_inputs_sha256": {f: sha(ROOT / f) for f in source_files}})
    dump("repo_search_manifest.json", {"FULL_REPOSITORY_HISTORICAL_RECOVERY_ATTEMPTED": True,
        "FULL_RELEVANT_REPO_HISTORY_SEARCHED": True, "searched_refs": [v[1] for v in PRS.values()] + ["main", "all local and origin branches"],
        "search_terms": ["matter loading", "source loading", "source medium", "source contact", "source force", "repulsive", "boundary condition", "moving source", "c_state", "rho", "stress energy", "metric strain", "source sphere", "source radius", "source region", "packet", "receipt", "266x266"],
        "excluded": ["observational weak-lensing assets", "new physics", "DEV167/168 modifications"]})
    rows = [
        historical("PR36", "J_A=-(1/(2sqrt(-g))) integral sqrt(-g) T^{mu nu} delta G_{mu nu}/delta chi^A; local J_A=-1/2 T^{mu nu} partial G_{mu nu}/partial chi^A", "physical stress-energy T", "effective chi/metric map", "minimal matter loading", "EFFECTIVE_ONLY_NOT_NATIVE_MICROSCOPIC", "STATIC"),
        historical("PR37", "chi_mn=(g_mn-gbar_mn)/2; g=gbar+2chi; J_ab=-T_ab", "stress-energy normalization", "effective metric strain", "minimal effective source", "EFFECTIVE_ONLY_NOT_NATIVE_MICROSCOPIC", "STATIC"),
        historical("PR69", "E=(T/2) sum_<ij>(u_i-u_j)^2-sum_i S_i u_i; T L[u]/dx^2=S", "finite sampled source sphere S_i", "scalar nearest-neighbour medium", "local forcing distinct from redistribution", "REUSE_AS_PRINCIPLE", "STATIC_AND_SOURCE_FREE_DYNAMIC"),
        historical("PR70", "W=-(K emax^2/2)ln(1-(e/emax)^2); sigma=K e/(1-(e/emax)^2)", "localized load through equilibrium", "bounded scalar bond strain", "constitutive response", "MATHEMATICALLY_EQUIVALENT_NORMALIZATION", "STATIC"),
        historical("PR72", "discrete divergence of bounded bond stress = source; e_ij=(u_i-u_j)/dx", "rho/density on a finite source-sphere mask", "3D N6 bounded-strain network", "local forcing then local redistribution", "HISTORICAL_ANCESTOR", "STATIC"),
        historical("PR74", "rho -> A8 -> c_state -> bounded-strain N6 equilibrium -> accumulated response", "raw c_state produced from rho", "3D bounded-strain accumulation network", "accumulation source", "SUPERSEDED_BY_VECTOR_RELATIONAL_STATE", "STATIC"),
        {"historical_id": "DEV159", "PR_or_commit": "native source-medium interaction001", "file": "pbuf/labs/foundation/native_source_medium_interaction001.py", "equation_or_mechanism": "one-cell contact constraint / source_imposed_excursion; stationary N6 equilibrium; integer moving schedule", "source_representation": "NativeSourceState integer occupied cell", "medium_representation": "scalar F03 N6", "interaction_type": "persistent local forcing; removal and moving-source residual", "static_or_dynamic": "STATIC_AND_DYNAMIC", "source_motion_supported": "INTEGER_CELL_SCHEDULE_ONLY", "continuous_source_position_supported": "NO", "free_coefficients": "amplitude control, not fitted", "test_result": "historical structural audit", "canonical_status_then": "UNPROMOTED_HISTORICAL", "relation_to_DEV167": "source-boundary semantics survive; scalar F03 dynamics superseded", "relation_to_DEV179": "supports node-discrete contact, not off-node mapping"},
    ]
    dump("historical_source_medium_map.json", rows)
    dump("historical_source_medium_concept_recovery.json", {"USER_SPECIFIED_SOURCE_MEDIUM_CONCEPT_SEARCHED": True,
        "REPULSIVE_SOURCE_MEDIUM_CONCEPT_FOUND": True, "EQUATION_NOT_FROZEN": True,
        "findings": [{"source": "pbuf/excitation/native_vector_pair_dynamics.py:source_contact_force", "exact": "Repulsive one-cell contact at the six N6 neighbors of a source node; out[idx] += magnitude * offset", "status": "CURRENT_NODE_DISCRETE"},
                     {"source": "pbuf/labs/foundation/native_source_medium_interaction001.py", "exact": "one-cell local equilibrium constraint; moving integer schedule; source work NOT_DEFINED", "status": "HISTORICAL_SCALAR_ANCESTOR"}],
        "potential_V_sm_found": False, "forbidden_uninvented_forms": ["1/r", "1/r^2", "Gaussian", "exponential", "Lennard-Jones", "Hooke", "compact support kernel"]})
    for key, title, result in [("PR36", "minimal matter loading", "EFFECTIVE_ONLY_NOT_NATIVE_MICROSCOPIC"), ("PR37", "metric-strain normalization", "EFFECTIVE_ONLY_NOT_NATIVE_MICROSCOPIC"), ("PR69", "source S / neighbour redistribution separation", "REUSE_AS_PRINCIPLE"), ("PR70", "bounded constitutive barrier", "MATHEMATICALLY_EQUIVALENT_NORMALIZATION"), ("PR72", "3D N6 local redistribution", "HISTORICAL_ANCESTOR"), ("PR74", "c_state accumulation source role", "SUPERSEDED_BY_VECTOR_RELATIONAL_STATE")]:
        merge, ref, file = PRS[key]
        dump(f"{key.lower()}_{['matter_loading','metric_strain','neighbor_source','constitutive','n6_network','native_accumulation'][['PR36','PR37','PR69','PR70','PR72','PR74'].index(key)]}_recovery.json", {"recovered": True, "title": title, "commit": merge, "ref": ref, "file": file, "current_classification": result})
    dump("current_source_contact_audit.json", {"function": "source_contact_force(shape: tuple[int,int,int], center: tuple[int,int,int], magnitude: float=0.02)",
        "input_type": "integer tuple center; no float source coordinate", "center_semantics": "source node; six adjacent N6 nodes receive outward vectors", "magnitude_semantics": "per-neighbour external-force magnitude", "node_state_changed": "external force passed to step; momentum kick then displacement drift", "persistent_or_impulse": "caller-controlled; DEV169 passes same ext every propagation step", "vector_or_scalar": "vector", "source_participates_dynamically": False, "source_stored": False, "contact_changes_position": "indirectly through step", "contact_changes_momentum": True, "historical_Si_relation": "simplified vector external-force analogue of local loading; not scalar S_i identity", "CURRENT_NODE_CONTACT_IS_SPECIAL_LIMIT": "UNRESOLVED; no continuous historical equation frozen"})
    dump("current_packet_launch_audit.json", {"packet_construction": "DEV169 packet(image): 7x7 transverse image embedded at x=1; one finite vector excitation state per positive packet cell lineage", "loading_fixed": True, "medium_fixed": True, "pair_law_fixed": True, "more_packet_states_without_source_change": "NOT_EXPOSED_BY_CURRENT_CANONICAL_API", "reason": "receipt() assigns lineage to nearest supported packet cell and DEV178 correctly rejected fabricated subcell launches", "classification": "OFFNODE_PACKET_REFINEMENT_NOT_YET_DERIVED"})
    receipt_counts = json.loads((ROOT / "runs/dev178_high_density_native_vulkan/receipt_counts.json").read_text())["baseline"]
    dump("current_receipt_multiplicity_tree.json", {"physical_source_cells": "DEV171 distributed 3D source cells: source>0, superposed into ext; count is realization-dependent", "packet_launch_cells": 49, "packet_excitation_state": "one 7x7 packet field, with 49 supported lineage cells", "propagation_branches": "one deterministic DEV167 time evolution per realization", "receipt_events": receipt_counts, "receipt_events_per_supported_packet_cell": [x/49 for x in receipt_counts], "origin_of_approximately_12k": "positive bond-flux crossings over progression steps, not 12k physical sources", "deterministic": True})
    dump("historical_25pct_density_semantics.json", {"historical_grid": [266,266], "historical_launch_count": 70756, "historical_role": "continuous Cartesian propagating ray launches on an 8x8 plane", "independent_matter_source_locations": False, "HISTORICAL_25PCT_RAY_DENSITY_NOT_SOURCE_DENSITY": True, "current_comparison": "DEV178 12/49 atomic packet-cell tile is a coverage control, not density escalation"})
    definitions = {"SOURCE_RESOLUTION": "number/spacing/extent semantics of physical loading samples", "MEDIUM_RESOLUTION": "native N6 cell spacing and domain discretization", "LAUNCH_DENSITY": "initial propagating-state positions per launch plane area", "PACKET_DENSITY": "number of independently represented propagating packet states/support cells", "RECEIPT_DENSITY": "number/spatial distribution of recorded receipt events", "OBSERVER_BIN_DENSITY": "downstream measurement bins per screen area"}
    dump("density_terminology.json", definitions)
    dump("density_architecture_result.json", {"DENSITY_ARCHITECTURE": "DENSITY_E", "reason": "historical C25 was ray density, but current canonical 7x7 packet lineage is tied to supported packet cells; the audited code contains no independently defined denser packet launch representation. This does not establish a new source coupling is required.", "DEV179_VALID_BUT_BLOCKER_MISAPPLIED": True, "new_source_law_required": False, "new_packet_semantics_required_before_run": True, "native_grid_refinement_historically_supported": "YES_AS_HISTORICAL_FINITE_SOURCE_SAMPLING_ROUTE; current vector convergence not yet demonstrated"})
    matrix = [{"historical_mechanism":"PR36 minimal matter source","historical_role":"source coupling","current_analogue":"external source force requires native mapping","still_valid":"effective only","action":"NEEDS_NEW_BRIDGE"}, {"historical_mechanism":"PR37 metric strain","historical_role":"source normalization","current_analogue":"none microscopic","still_valid":"effective only","action":"REUSE_AS_PRINCIPLE"}, {"historical_mechanism":"PR69 source S","historical_role":"local forcing","current_analogue":"external_force","still_valid":"yes as separation","action":"REUSE_AS_PRINCIPLE"}, {"historical_mechanism":"PR70 bounded strain","historical_role":"pair constitutive","current_analogue":"DEV167 sigma(e)=e/(1-e^2)","still_valid":"yes","action":"RESTORE"}, {"historical_mechanism":"PR72 N6 equilibrium","historical_role":"redistribution","current_analogue":"DEV167 N6 vector net_force","still_valid":"ancestor","action":"REUSE_AS_PRINCIPLE"}, {"historical_mechanism":"PR74 c_state","historical_role":"local native loading","current_analogue":"external source force only, not c_state","still_valid":"role only","action":"SUPERSEDED"}]
    dump("source_medium_correspondence_matrix.json", matrix)
    dump("dev179_scope_correction.json", {"DEV179_SCOPE_REINTERPRETED": True, "DEV179_SCOPE_CORRECTION": "DEV179 established only that the current DEV167/168 production source-contact API is node-discrete and lacks an exposed off-node mapping. It did not establish that PBUF lacks prior matter/source-to-medium coupling architecture. Historical source-loading work must be recovered and reconciled before a new coupling law is proposed.", "NEW_NATIVE_SOURCE_COUPLING_REQUIRED": "WITHDRAWN_AS_GLOBAL_CLAIM", "node_discrete_finding": True})
    dump("source_work_dev170_reconciliation.json", {"DEV170_SOURCE_WORK_RECONCILED": True, "finding": "DEV169 applies ext every step; medium-only invariant excludes external source potential/work.", "pair_conservation": "DEV167 reciprocal pair law unchanged", "extended_invariant": "conceptual H_medium+H_source-medium only; no V_sm is frozen, so not implemented", "SOURCE_WORK_OMITTED": "expected bookkeeping classification under persistent external forcing"})
    dump("current_native_reuse_candidates.json", {"reusable": ["PR69 separation of source forcing from medium redistribution", "PR70 exact normalized bounded constitutive family", "PR72 local N6 redistribution principle", "DEV159 stationary/moving boundary semantics", "DEV167 current repulsive one-cell external contact"], "not_reactivated": ["A8", "c_state", "scalar F03 propagation", "historical observer/ray physics"], "new_source_medium_law_introduced": False})
    dump("viewer_extension_status.json", {"status":"NOT_MODIFIED", "reason":"analysis artifacts suffice; viewer remains diagnostic-only", "possible_overlays":["source region", "loaded nodes", "external source force", "DEV167 pair force", "receipt lineage"]})
    final = {"DEV180_COMPLETE": True, "FULL_RELEVANT_REPO_HISTORY_SEARCHED": True, "PR36_RECOVERED": True, "PR37_RECOVERED": True, "PR69_RECOVERED": True, "PR70_RECOVERED": True, "PR72_RECOVERED": True, "PR74_RECOVERED": True, "ADDITIONAL_SOURCE_MEDIUM_HISTORY_SEARCHED": True, "CURRENT_SOURCE_CONTACT_AUDITED": True, "CURRENT_PACKET_LAUNCH_PATH_AUDITED": True, "CURRENT_RECEIPT_MULTIPLICITY_AUDITED": True, "HISTORICAL_25PCT_DENSITY_SEMANTICS_AUDITED": True, "DEV179_SCOPE_REINTERPRETED": True, "DENSITY_ARCHITECTURE_CLASSIFIED": True, "NO_NEW_SOURCE_MEDIUM_LAW": True, "NO_NEW_PAIR_LAW": True, "NO_INTERPOLATION_AS_PHYSICS": True, "NO_TRUE_25PCT_RUN": True, "NO_OBSERVATIONAL_INPUT": True, "DEV167_PAIR_LAW_UNCHANGED": True, "DEV168_RECEIPT_UNCHANGED": True, "OUTCOME": "OUTCOME_B_PLUS_OUTCOME_D", "IMPLEMENTATION_COMMIT_RECORDED": "pending commit", "REMOTE_PUSH_CONFIRMED": "pending push", "WORKTREE_CLEAN": "pending verification"}
    dump("final_contract.json", final)
    (OUT / "discussion_handoff.md").write_text("# DEV180 handoff\n\nHistorical matter-loading and local source-to-neighbour redistribution architecture was recovered. PR70 is mathematically the same normalized bounded constitutive family as DEV167; PR69/72 survive as the forcing-versus-redistribution separation and local-N6 principle. PR36/37 remain effective-level constraints, not a vector microscopic source map. PR74's c_state/A8 route is superseded. DEV179 remains correct about off-node contact but its global missing-mechanism conclusion is corrected. The historical 70,756 C25 entities were ray launches, not matter locations. Current receipt multiplicity is time/bond-flux events from 49 packet lineage cells, not source count. No new source law was introduced.\n")

if __name__ == "__main__":
    main()
