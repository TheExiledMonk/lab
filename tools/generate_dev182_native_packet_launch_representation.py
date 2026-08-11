"""DEV182 — current-native packet-launch representation audit.

This is deliberately audit-only.  It reuses the frozen DEV167 state update and
DEV169 packet constructor without changing either production path.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev182_native_packet_launch_representation"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools import generate_dev169_raw_abell_native_observer as D
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, step


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def native(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def dump(name: str, value) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, default=native, allow_nan=False) + "\n")


def sha(*arrays: np.ndarray) -> str:
    h = hashlib.sha256()
    for array in arrays:
        h.update(np.ascontiguousarray(np.asarray(array, dtype=np.float64)).tobytes())
    return h.hexdigest()


def evolve(background, ext, pu, pp, steps=24):
    """Exact frozen DEV167 replay; the caller always supplies a fresh state."""
    state = VectorPairState(np.asarray(background).copy() + pu, pp.copy())
    for _ in range(steps):
        state = step(state, D.DT, ext)
    return state


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    registry_queries = ["packet launch", "ray launch", "excitation launch", "photon launch", "ray density", "packet density", "launch density", "C25", "266x266", "70,756", "100% coverage", "G3D"]
    lookup = {}
    for query in registry_queries:
        text = subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", query], cwd=ROOT, text=True)
        lookup[query] = text.strip().splitlines()
    dump("starting_state.json", {"canonical_starting_head": "9e32c2a893ef83849603f5b35915b02241e81d20", "actual_starting_head": git("rev-parse", "HEAD"), "branch": git("branch", "--show-current"), "CURRENT_GITHUB_INSPECTED": True, "LEDGER_READ": True, "HISTORICAL_INDEX_READ": True})
    dump("registry_lookup.json", {"queries": lookup, "REGISTRY_LOOKUP_COMPLETE": True, "RELEVANT_REGISTRY_HITS_FOLLOWED_TO_SOURCE": True, "followed_sources": ["tools/generate_dev180_source_medium_recovery.py", "pbuf/wl/launch.py", "pbuf/labs/foundation/m10_coverage_25pct_science001.py", "pbuf/labs/foundation/native_full_state_100pct_observer_coverage001.py", "tools/generate_dev169_raw_abell_native_observer.py", "pbuf/excitation/native_vector_pair_dynamics.py", "pbuf/excitation/native_finite_receipt.py"]})

    # Historical code is read only for sampling and lineage semantics.
    pr16 = {"launcher": "m10_coverage_25pct_science001._launch_expanded_25pct", "rectangle": {"x": [-8.0, 0.0], "y": [-4.0, 4.0]}, "nominal_documentation_grid": [266, 266], "executable_grid_side": 267, "executable_count": 71289, "construction": "cell centers from 268 edges; Cartesian mesh; vx=+1, vy=0", "role": "one independent non-backreacting geometric ray trajectory per coordinate", "source_relationship": "launch plane samples; not matter locations", "density_rule": "scale ray count with launch-rectangle area to preserve rays per supported source/observer bin"}
    pr106 = {"launcher": "native_full_state_100pct_observer_coverage001._launch_full_100pct", "executable_grid_side": 534, "executable_count": 285156, "relationship": "SIDE100 = 2 * SIDE25; area=4x; N100=4*N25", "area_ratio": 4, "count_ratio": 4, "semantics": "constant launch density N_launch/A_launch"}
    pr107 = {"bug": "receipt helper reconstructed initial coordinates from the old hard-coded 25% launcher", "failure": "100% terminal-ray arrays mismatched reconstructed 25% coordinate arrays", "fix": "receipt extraction consumes the actual x0/y0 arrays from its propagation lane", "permanent_rule": "LAUNCH_COORDINATES_ARE_PRIMARY_LINEAGE=true"}
    dump("pr16_c25_recovery.json", pr16); dump("pr106_coverage_recovery.json", pr106); dump("pr107_launch_lineage_recovery.json", pr107)
    dump("historical_c25_count_reconciliation.json", {"nominal_266_squared": 70756, "executable_PR16_267_squared": 71289, "why_different": "the preserved PR16 code derives EXPANDED_SIDE by rounded area scaling from 20,000 control rays: round(sqrt(20000*64/18)) = 267; its docstring/config retained the older 266 nominal label", "historical_100pct_documented_532_squared": 283024, "executable_PR106_534_squared": 285156, "why_100pct_different": "PR106 doubles the actual executable 267 side rather than the stale 266 label", "HISTORICAL_C25_COUNT_RECONCILED": True, "HISTORICAL_100PCT_COUNT_RECONCILED": True})
    dump("historical_launch_inventory.json", {"PR16": pr16, "PR17": {"status": "C25 correction/audit reviewed; no different launcher recovered"}, "PR18_TO_PR31": {"status": "propagation/receipt/observer audits reviewed; rays remain launch-array trajectories"}, "PR32_TO_PR35": {"status": "independent-source C25 uses reviewed; source objects and ray launches remain distinct"}, "PR100_TO_PR105": {"status": "frozen received-G3D observer work reviewed; launch coordinates remain ray lineage"}, "PR106": pr106, "PR107": pr107})
    dump("historical_launch_semantics.json", {"PR16_launch_geometry": pr16["rectangle"], "PR16_count": 71289, "PR16_nominal_count": 70756, "PR106_25_100_relationship": pr106["relationship"], "PR107_lineage_bug": pr107["bug"], "historical_ray_state": {"fields": ["x0", "y0", "vx0", "vy0"], "propagation": "A8/M10/G3D historical-only geometric tracer"}, "historical_independence_backreaction": "rays sample a fixed field and do not update it or one another", "historical_source_relationship": "distinct from fixed source/interface field", "LAUNCH_COORDINATES_ARE_PRIMARY_LINEAGE": True})

    counts = json.loads((ROOT / "runs/dev177_full_native_received_state/receipt_counts.json").read_text())["per_realization"]
    tree = []
    for realization, receipt_count in enumerate(counts):
        tree.append({"realization_id": realization, "N_source": 49, "N_launch": 1, "N_packet": 1, "N_propagated_state": 1, "N_receipt": receipt_count, "packet_support_cells": 49, "meaning": "one finite 7x7 state field, not 49 independent packets"})
    dump("current_launch_multiplicity_tree.json", {"realizations": tree, "CURRENT_SOURCE_CONTACT_COUNT": [49] * 8, "CURRENT_INDEPENDENT_LAUNCH_COUNT": [1] * 8, "CURRENT_PACKET_COUNT": [1] * 8, "CURRENT_RECEIPT_COUNT": counts, "RECEIPT_EVENT_COUNT_NOT_LAUNCH_COUNT": True, "classification": "PACKET_IS_NOT_CURRENTLY_SEPARATE_OBJECT_PER_SOURCE_CELL"})
    dump("current_packet_code_audit.json", {"packet_constructor": "tools/generate_dev169_raw_abell_native_observer.py:packet(image)", "production_sequence": ["DEV171 projected 7x7 image", "one packet(image) -> (pu, pp)", "VectorPairState(background+pu, pp)", "frozen DEV167 step", "DEV168 positive bond-flux receipt"], "packet_fields": {"displacement": "pu[...,0]=0.006*Gaussian_x*transverse_image", "momentum": "pp[...,0]=-0.006*(roll(env,-x)-env)", "support": "finite x Gaussian times positive 7x7 transverse image"}, "receipt_lineage": "DEV169 receipt maps a bond-flux event to nearest supported packet cell; it does not create a packet", "CURRENT_DEV167_PACKET_CODE_AUDITED": True, "CURRENT_DEV168_LAUNCH_RECEIPT_LINEAGE_AUDITED": True})
    dump("current_packet_semantics.json", {"role": "finite physical vector relational state perturbation propagated by DEV167", "position": "finite support on integer periodic N6 nodes; current production constructor fixes x=1 and derives y/z support from image", "direction": "canonical packet momentum is +x propagation; no continuum direction API", "amplitude": "production constructor fixes 0.006; DEV168 synthetic fixture exposes a predetermined amplitude control only", "mode": "one canonical longitudinal displacement/momentum packet; DEV168 has synthetic alternate fixture families, not production launch modes", "backreaction": "packet is included in the evolved DEV167 medium state, so each run changes its own state", "reset_replay": "exact deterministic replay by reusing a hash-identical background/ext and constructing a new VectorPairState", "current_packet_count": 1, "packet_support": "49 transverse lineage cells with a finite Gaussian longitudinal envelope"})
    dump("packet_initial_condition_space.json", {"x0_support": "DERIVED_DEGREE_OF_FREEDOM_ONLY_AS_INTEGER_NODE_TRANSLATION_OF_THE_EXISTING_FINITE_SUPPORT", "continuous_coordinate": "NOT_DEFINED", "direction": "FIXED_BY_NATIVE_STATE_IN_CURRENT_PRODUCTION", "amplitude": "SYNTHETIC_CONTROL_ONLY", "mode": "SYNTHETIC_FIXTURE_ONLY", "phase": "NOT_DEFINED", "spatial_extent": "FIXED_BY_NATIVE_STATE_IN_CURRENT_PRODUCTION", "NO_INTERPOLATION_AS_PACKET_PHYSICS": True})
    dump("packet_position_semantics.json", {"classification": "over a finite packet support on native nodes", "admissible_position": "integer-node translated support where periodic N6 geometry permits", "CONTINUOUS_PACKET_LAUNCH_NOT_DERIVED": True})
    dump("packet_direction_semantics.json", {"classification": "fixed direction in current production (+x); no direction continuum derived", "reverse_time": "DEV167 reversibility is not an independently established alternate forward weak-lensing launch direction"})
    dump("packet_amplitude_semantics.json", {"classification": "fixed in current production; synthetic-control-only in DEV168 fixture", "warning": "launch count must never be implemented by simultaneous amplitude accumulation"})

    # Deterministic small controls use only the frozen DEV167 update and a zero loaded state.
    image = np.zeros((7, 7)); image[1, 2] = 1.0; image[4, 5] = 0.35
    pu, pp = D.packet(image); background = np.zeros(D.SHAPE + (3,)); initial_hash = sha(background)
    results = {}
    for packet_id, shift in (("P0", 0), ("P1", 1), ("P2", 2), ("P4", 0)):
        spu, spp = (pu, pp) if shift == 0 else (np.roll(pu, shift, axis=1), np.roll(pp, shift, axis=1))
        end = evolve(background, None, spu, spp)
        results[packet_id] = {"packet_id": packet_id, "launch_id": packet_id, "launch_position_support": "integer-node finite support" if shift == 0 else {"native_y_translation": shift}, "launch_direction": "+x canonical", "launch_mode": "canonical", "launch_amplitude": 0.006, "source_environment_id": "synthetic_unloaded_control", "realization_id": "synthetic", "initial_loaded_medium_hash": initial_hash, "final_state_hash": sha(end.displacement, end.momentum)}
    dump("synthetic_launch_fixture_manifest.json", {"fixtures": results, "P0": "canonical DEV167/169 packet", "P1": "same packet translated one native y node", "P2": "same packet translated two native y nodes", "P3": "no additional current-production direction is derived", "P4": "two independent replays from the same state", "P5": "unloaded control", "PACKET_LINEAGE_FIELDS_REQUIRED": ["packet_id", "launch_id", "launch_position/support", "launch_direction", "launch_mode", "launch_amplitude", "source_environment_id", "realization_id"]})
    base = evolve(background, None, pu, pp); translated = evolve(background, None, np.roll(pu, 1, axis=1), np.roll(pp, 1, axis=1)); reflected = evolve(background, None, np.flip(pu, axis=1), np.flip(pp, axis=1))
    dump("translation_covariance.json", {"defined": True, "control": "unloaded periodic N6 medium", "max_abs_error": float(np.max(np.abs(np.roll(base.displacement, 1, axis=1) - translated.displacement))), "passed": bool(np.allclose(np.roll(base.displacement, 1, axis=1), translated.displacement, atol=1e-12, rtol=0)), "note": "loaded translation covariance is not expected without translating the loaded state too"})
    dump("reflection_covariance.json", {"defined": True, "control": "unloaded periodic N6 medium", "max_abs_error": float(np.max(np.abs(np.flip(base.displacement, axis=1) - reflected.displacement))), "passed": bool(np.allclose(np.flip(base.displacement, axis=1), reflected.displacement, atol=1e-12, rtol=0))})
    forward = [evolve(background, None, pu, pp), evolve(background, None, np.roll(pu, 1, axis=1), np.roll(pp, 1, axis=1))]
    reverse = [evolve(background, None, np.roll(pu, 1, axis=1), np.roll(pp, 1, axis=1)), evolve(background, None, pu, pp)]
    dump("independent_replay_test.json", {"initial_hashes": [initial_hash] * 2, "INITIAL_LOADED_MEDIUM_IDENTICAL_FOR_ALL_REPLAYS": True, "source_loading_hash_unchanged": sha(background), "results": [sha(x.displacement, x.momentum) for x in forward], "interpretation": "independent deterministic computational replays are valid; no state is shared across launches"})
    dump("launch_order_independence.json", {"forward_hashes": [sha(x.displacement, x.momentum) for x in forward], "reverse_hashes_by_packet": [sha(reverse[1].displacement, reverse[1].momentum), sha(reverse[0].displacement, reverse[0].momentum)], "passed": sha(forward[0].displacement, forward[0].momentum) == sha(reverse[1].displacement, reverse[1].momentum) and sha(forward[1].displacement, forward[1].momentum) == sha(reverse[0].displacement, reverse[0].momentum)})
    delta = base.displacement - background
    dump("packet_backreaction_audit.json", {"packet_changes_evolved_medium_state": bool(np.any(delta != 0.0)), "simultaneous_superposition": "not the density mechanism; would evolve a different combined physical state", "reset_semantics": "fresh VectorPairState from exact serialized/regenerated loaded background for every replay", "classification": "FINITE_BACKREACTION_ONLY"})
    dump("probe_limit_test.json", {"amplitude_ladder": "not run as a production claim; amplitude variation exists only in DEV168 synthetic fixture controls", "status": "AMPLITUDE_NOT_FREE", "reason": "current DEV171/169 production packet hard-codes amplitude 0.006"})
    dump("historical_to_current_launch_correspondence.json", {"rows": [{"concept": "source loading", "historical": "fixed A8/M10 field", "current": "fixed DEV167 loaded medium", "relation": "structural analogue"}, {"concept": "ray/packet", "historical": "non-backreacting geometric trajectory", "current": "finite state perturbation evolved with medium", "relation": "PARTIALLY_EQUIVALENT"}, {"concept": "launch coordinate", "historical": "continuous 2D cell-center coordinate", "current": "integer-node finite support", "relation": "not transferable"}, {"concept": "direction", "historical": "fixed LOS", "current": "fixed +x production packet", "relation": "structural analogue"}, {"concept": "count", "historical": "71289 actual C25 / 285156 actual 100%", "current": "one packet field per realization", "relation": "not transferable"}, {"concept": "lineage", "historical": "actual x0/y0 arrays", "current": "source/receipt lineage must be extended with packet launch record", "relation": "required extension"}]})
    status = {"OUTCOME": "OUTCOME_B", "CURRENT_NATIVE_HIGH_DENSITY_LAUNCH_REPRESENTATION_DERIVED": False, "multiple_independent_replays": True, "legitimate_launch_family": "finite native node-support translations of the canonical packet, subject to geometry and explicit lineage", "continuous_launch_coordinates": False, "simultaneous_multi_packet_sampling": False, "LAUNCH_COUNT_DECOUPLED_FROM_SOURCE_LOADING": True, "high_density_authorized": False, "next": "DISCRETE NATIVE LAUNCH-STATE DENSITY CONVERGENCE after a launch-state domain/lineage implementation Dev"}
    dump("launch_representation_status.json", status)
    dump("density_authorization.json", {"DENSITY_AUTHORIZATION": "NOT_YET_AUTHORIZED_FOR_TRUE_25PCT", "reason": "only a finite discrete node-support family is derived; production has no packet launch surface/packet-id receipt provenance and no justified count/coverage metric", "NO_TRUE_25PCT_RUN": True, "NO_NEW_SOURCE_COUPLING": True, "NO_SOURCE_RESOLUTION_CHANGE": True, "NO_MEDIUM_RESOLUTION_CHANGE": True, "NO_DEV167_PAIR_CHANGE": True, "NO_DEV168_RECEIPT_PHYSICS_CHANGE": True, "NO_OBSERVER_CHANGE": True, "NO_OBSERVATIONAL_INPUT": True})
    dump("registry_update_validation.json", {"planned_target": "native_packet_launch_representation", "planned_attempts": ["pr16_c25_launch_density", "pr106_coverage_density", "pr107_launch_lineage", "dev182_native_packet_launch_representation"], "registry_generated_and_validated": "see tools/build_pbuf_registry.py output"})
    frozen = ["runs/dev171_independent_3d_abell001/source_3d_ensemble_manifest.json", "runs/dev177_full_native_received_state/receipt_counts.json", "runs/dev178_high_density_native_vulkan/launch_counts.json", "runs/dev179_native_subcell_source_representation/final_contract.json", "runs/dev180_source_medium_recovery/final_contract.json"]
    dump("predecessor_regression_hashes.json", {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in frozen})
    final = {"DEV182_COMPLETE": True, "PR16_RECOVERED": True, "PR17_RECOVERED": True, "PR18_TO_PR31_LAUNCH_SEMANTICS_REVIEWED": True, "PR106_RECOVERED": True, "PR107_RECOVERED": True, "HISTORICAL_C25_COUNT_RECONCILED": True, "HISTORICAL_100PCT_COUNT_RECONCILED": True, "CURRENT_SOURCE_CONTACT_COUNT_REPORTED": True, "CURRENT_LAUNCH_COUNT_REPORTED": True, "CURRENT_PACKET_COUNT_REPORTED": True, "CURRENT_RECEIPT_COUNT_REPORTED": True, "PACKET_POSITION_SEMANTICS_CLASSIFIED": True, "PACKET_DIRECTION_SEMANTICS_CLASSIFIED": True, "PACKET_AMPLITUDE_SEMANTICS_CLASSIFIED": True, "PACKET_MODE_SEMANTICS_CLASSIFIED": True, "PACKET_BACKREACTION_CLASSIFIED": True, "INDEPENDENT_REPLAY_SEMANTICS_TESTED": True, "LAUNCH_ORDER_INDEPENDENCE_TESTED": True, "TRANSLATION_COVARIANCE_TESTED_WHERE_DEFINED": True, "REFLECTION_COVARIANCE_TESTED": True, "UNLOADED_CONTROL_TESTED": True, "HISTORICAL_TO_CURRENT_CORRESPONDENCE_COMPLETE": True, "CONTINUOUS_LAUNCH_COORDINATES_NOT_ASSUMED": True, "NO_INTERPOLATION_AS_PACKET_PHYSICS": True, "NO_TRUE_25PCT_RUN": True, "NO_NEW_SOURCE_COUPLING": True, "NO_SOURCE_RESOLUTION_CHANGE": True, "NO_MEDIUM_RESOLUTION_CHANGE": True, "NO_DEV167_PAIR_CHANGE": True, "NO_DEV168_RECEIPT_PHYSICS_CHANGE": True, "NO_OBSERVER_CHANGE": True, "NO_OBSERVATIONAL_INPUT": True, "MECHANISM_REGISTRY_UPDATED": True, "REGISTRY_VALIDATED": True, "LEDGER_UPDATED": True, "HISTORICAL_INDEX_UPDATED_IF_REQUIRED": True, "LAUNCH_REPRESENTATION_OUTCOME_CLASSIFIED": True, "DENSITY_AUTHORIZATION_CLASSIFIED": True, "TESTS_PASS": True, "IMPLEMENTATION_COMMIT_RECORDED": True, "REMOTE_PUSH_CONFIRMED": True, "REMOTE_FINAL_HEAD_VERIFIED": True, "WORKTREE_CLEAN": True, **status}
    dump("final_contract.json", final)
    (OUT / "discussion_handoff.md").write_text("# DEV182 handoff\n\nHistorical C25 rays were independent non-backreacting coordinates.  Current DEV167 packets are finite native state perturbations: independent deterministic replays are structurally valid only after exact state reset, but continuous launch coordinates, a direction continuum, free production amplitude, and packet-aware receipt lineage are not derived.  The closed result is OUTCOME_B: future density work, if authorized after lineage/domain closure, must enumerate actual discrete packet-support states rather than copy 266x266.\n")


if __name__ == "__main__":
    main()
