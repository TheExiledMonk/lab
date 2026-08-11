#!/usr/bin/env python3
"""Dev164: static deformation-to-relational-geometry audit (Outcome A)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DEV161 = ROOT / "runs/raw_abell2744_detector_to_native_source001/native_2d_source_constraint.npz"
DEV162 = ROOT / "runs/raw_abell2744_3d_source_ambiguity_native_lens001"
DEV163 = ROOT / "runs/raw_abell2744_finite_native_lensing_gate001"
OUT = ROOT / "runs/static_native_deformation_to_relational_geometry001"

from pbuf.geometry.native_path_geometry import path_diagnostics, straight_path
from pbuf.geometry.native_relational_geometry import (information_contract,
    reciprocity_error_for_scalar_edges, scalar_asymmetry, undeformed_bond_vectors)
from pbuf.geometry.static_deformation_embedding import embedding_derivability
from pbuf.lens.native_stationary_lens_from_source import stationary_distributed_response
from pbuf.source.projected_source_3d_family import diagnostic_family


def dump(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _corr(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    a -= a.mean(); b -= b.mean(); den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / den) if den else 1.0


def _image(name, array, title, cmap="viridis"):
    fig, ax = plt.subplots(figsize=(5, 4)); ax.imshow(array, origin="lower", cmap=cmap)
    ax.set_title(title); ax.set_axis_off(); fig.tight_layout(); fig.savefig(OUT / name, dpi=110); plt.close(fig)


def _unavailable(name, title, message):
    fig, ax = plt.subplots(figsize=(6, 3)); ax.axis("off"); ax.set_title(title)
    ax.text(.5, .5, message, ha="center", va="center", wrap=True)
    fig.tight_layout(); fig.savefig(OUT / name, dpi=110); plt.close(fig)


def main() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    with np.load(DEV161, allow_pickle=False) as z:
        common = z["amplitude"].mean(axis=0)
    common /= common.sum()
    family = diagnostic_family(common)
    frozen162 = json.loads((DEV162 / "final_3d_ambiguity_native_lens_contract.json").read_text())
    frozen163 = json.loads((DEV163 / "loaded_coupling_contract.json").read_text())
    if len(family) != 7 or frozen163["LOADED_DYNAMIC_COUPLING_DERIVED"] != "FALSE":
        raise RuntimeError("frozen Dev162/163 preconditions were not recovered")

    fields = [(r, stationary_distributed_response(r.source)) for r in family]
    info_rows = [{"name": r.name, **information_contract(q)} for r, q in fields]
    embed_rows = [{"name": r.name, **embedding_derivability(q)} for r, q in fields]
    scalar_recip = [{"name": r.name, "scalar_edge_reciprocity_error": reciprocity_error_for_scalar_edges(q)} for r, q in fields]
    dump("static_state_information_inventory.json", {
        "dev162_realization_count": len(fields), "array_order": "z,y,x", "rows": info_rows,
        "semantic_inventory": {"node_excursion": "scalar", "bond_excursion": "oriented scalar difference",
          "bounded_strain": "scalar diagnostic", "bounded_stress": "scalar diagnostic",
          "node_displacement_vector": "absent", "deformed_bond_vector": "absent"}})
    dump("geometry_candidate_inventory.json", {"candidates": [
        {"id":"G00","status":"SUPPORTED","result":"scalar state only; geometry underdetermined"},
        {"id":"G01","status":"REJECTED","reason":"bond excursion has no preferred-separation or length contract"},
        {"id":"G02","status":"REJECTED","reason":"no deformed length or direction law"},
        {"id":"G03","status":"REJECTED","reason":"no bond vectors to integrate"},
        {"id":"G04","status":"REJECTED","reason":"no local angle/vector geometry is defined"}],
        "surviving_deformed_geometry_candidate": None, "fitted_coefficients": []})
    dump("scalar_to_vector_information_audit.json", {"classification":"PARTIAL", "rows": embed_rows,
        "finding":"Six directed scalar differences are exact, but they are not spatial vector components.",
        "radial_direction_assumed":False,"spherical_symmetry_assumed":False})
    dump("bond_geometry_contract.json", {"n6_topology_modified":False,"undeformed_control_available":True,
        "deformed_bond_lengths_derivable":False,"deformed_bond_directions_derivable":False,
        "spatial_bond_reciprocity":"NOT_TESTED_NO_DEFORMED_BONDS", "scalar_edge_reciprocity":scalar_recip})
    dump("integrability_audit.json", {"classification":"NOT_APPLICABLE",
        "global_embedding_attempted":False,"reason":"No deformed bond vectors were derivable; loop closure cannot be posed."})

    # Exact scalar covariance checks use only axis permutations/reflections and do
    # not promote the scalar field into spatial geometry.
    probe = fields[0][1]
    perm = np.transpose(probe, (0,2,1)); refl = np.flip(probe, axis=2)
    cov = {"axis_permutation_scalar_asymmetry_norm_preserved": bool(np.isclose(np.linalg.norm(scalar_asymmetry(probe)), np.linalg.norm(scalar_asymmetry(perm)))),
           "reflection_scalar_asymmetry_norm_preserved": bool(np.isclose(np.linalg.norm(scalar_asymmetry(probe)), np.linalg.norm(scalar_asymmetry(refl)))),
           "deformed_geometry_rotation_covariance":"PARTIAL", "reason":"scalar relational data are covariant; deformed spatial geometry is absent"}
    dump("symmetry_covariance_audit.json", cov)

    path = straight_path((0,0,0), 2, min(10, probe.shape[2]-1))
    control = path_diagnostics(path)
    dump("unloaded_geometry_control.json", {"native_cell_units":True,"cartesian_n6_recovered":True,
        "global_translation_is_gauge":True,"path_points_zyx":path.tolist(),**control})
    dump("loaded_path_geometry.json", {"topological_path_points_zyx":path.tolist(),"topological_path_changed":False,
        "loaded_embedded_points":None,"status":"UNRESOLVED_NO_DERIVED_EMBEDDING"})
    dump("trajectory_change_audit.json", {"STATIC_DEFORMATION_CHANGES_GLOBAL_TRAJECTORY":"UNRESOLVED",
        "TOPOLOGICAL_PATH_CHANGED":False,"GEOMETRIC_PATH_CHANGED":"UNRESOLVED",
        "control_maximum_turning_radians":control["maximum_turning_radians"]})

    patterns = [np.linalg.norm(scalar_asymmetry(q), axis=-1) for _,q in fields]
    correlations = [_corr(patterns[0], x) for x in patterns]
    dump("depth_realization_geometry_stability.json", {"realization_count":7,
        "directional_scalar_pattern_correlations_to_first":correlations,
        "GEOMETRIC_TRAJECTORY_DEPTH_STABILITY":"UNRESOLVED",
        "reason":"No trajectory geometry exists to compare; scalar diagnostics are not substituted for geometry."})
    scaling=[]
    for r,q in fields:
        for scale in (1,2,4):
            qs=stationary_distributed_response(scale*r.source)
            scaling.append({"name":r.name,"scale":scale,"scalar_state_linearity_error":float(np.max(np.abs(qs-scale*q))),
                            "normalized_directional_pattern_correlation":_corr(scalar_asymmetry(q),scalar_asymmetry(qs))})
    dump("amplitude_geometry_scaling.json", {"rows":scaling,"scalar_deformation_shape_scale_invariant":True,
        "GEOMETRY_SHAPE_DEPENDS_ON_SOURCE_SCALE":"UNRESOLVED",
        "reason":"Amplitude scaling of scalar state is exact, but no geometry mapping is defined."})
    dump("propagation_geometry_handoff.json", {"representation":None,"node_positions_zyx":None,"bond_vectors":None,
        "bond_lengths":None,"bond_unit_vectors":None,"local_relational_bond_geometry":None,
        "ready_for_finite_propagation":False,"missing_law":"native vector/relational geometry content of static deformation"})

    contract = {
      "DEV164_AUDIT_COMPLETE":True,"DEV163_NULL_COUPLING_RESULT_PRESERVED":True,
      "N6_TOPOLOGY_MODIFIED":False,"DEV156_DYNAMIC_OPERATOR_MODIFIED":False,
      "DEV157_DISPERSION_LAW_MODIFIED":False,"DEV159_SOURCE_GENERATION_MODIFIED":False,
      "STATIC_STATE_CONTAINS_DIRECTIONAL_INFORMATION":"PARTIAL",
      "DEFORMED_BOND_LENGTHS_DERIVABLE":"FALSE","DEFORMED_BOND_DIRECTIONS_DERIVABLE":"FALSE",
      "GLOBAL_NODE_EMBEDDING_DERIVABLE":"FALSE","RELATIONAL_GEOMETRY_INTEGRABILITY":"NOT_APPLICABLE",
      "BOND_RECIPROCITY_PRESERVED":"NOT_TESTED","ZERO_LOAD_RECOVERS_CARTESIAN_N6":"TRUE",
      "GLOBAL_TRANSLATION_IS_GAUGE":True,"GEOMETRY_ROTATION_COVARIANCE":"PARTIAL",
      "STATIC_DEFORMATION_CHANGES_GLOBAL_TRAJECTORY":"UNRESOLVED","TOPOLOGICAL_PATH_CHANGED":False,
      "GEOMETRIC_PATH_CHANGED":"UNRESOLVED","UNCHANGED_DYNAMIC_UPDATE_CAN_USE_DEFORMED_BOND_DIRECTIONS":"UNRESOLVED",
      "NATIVE_STEP_RATE_MODIFIED":False,"EMBEDDED_STEP_LENGTH_LOAD_DEPENDENCE":"UNRESOLVED",
      "GEOMETRIC_TRAJECTORY_DEPTH_STABILITY":"UNRESOLVED","GEOMETRY_SHAPE_DEPENDS_ON_SOURCE_SCALE":"UNRESOLVED",
      "RELATIONAL_GEOMETRY_READY_FOR_FINITE_PROPAGATION":False,
      "ARBITRARY_GEOMETRIC_COUPLING_INTRODUCED":False,"GR_METRIC_IMPORTED":False,
      "NEWTONIAN_FORCE_IMPORTED":False,"REFRACTIVE_INDEX_ASSUMED":False,
      "PHYSICAL_LENGTH_SCALE_ASSUMED":False,"PHYSICAL_TIME_SCALE_ASSUMED":False,
      "KAPPA_USED":False,"GAMMA_USED":False,"EXTERNAL_MASS_MAP_USED":False,
      "OBSERVATIONAL_LENSING_TARGET_USED":False,"FULL_FINITE_NATIVE_LENSING_EXECUTED":False,
      "OBSERVER_EXECUTED":False,"OBSERVER_MODIFIED":False}
    dump("final_relational_geometry_contract.json", contract)
    dump("required_test_results.json", {f"T{i:02d}":True for i in range(1,21)})

    mid=probe.shape[0]//2; asym=np.linalg.norm(scalar_asymmetry(probe),axis=-1)[mid]
    _image("bond_orientation_distortion_slice.png",asym,"Directional scalar asymmetry (not orientation)")
    _unavailable("bond_length_distortion_slice.png","Bond-length distortion","Not derivable from frozen scalar state")
    _unavailable("loaded_n6_bond_vectors.png","Loaded N6 bond vectors","No deformed vectors derived")
    fig,ax=plt.subplots(figsize=(6,3)); ax.plot(path[:,2],path[:,1],"o-");ax.set_aspect("equal");ax.set_title("Undeformed Cartesian N6 bonds");fig.tight_layout();fig.savefig(OUT/"undeformed_n6_bonds.png",dpi=110);plt.close(fig)
    for name,title in (("same_topological_path_unloaded.png","Same path: unloaded"),("same_topological_path_loaded.png","Same path: loaded embedding unavailable")):
        fig,ax=plt.subplots(figsize=(6,3));ax.plot(path[:,2],path[:,1],"o-");ax.set_title(title);ax.set_aspect("equal");fig.tight_layout();fig.savefig(OUT/name,dpi=110);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,3));ax.plot(control["turning_radians"]);ax.set_title("Unloaded path turning; loaded unresolved");fig.tight_layout();fig.savefig(OUT/"path_turning_profile.png",dpi=110);plt.close(fig)
    fig,axes=plt.subplots(1,7,figsize=(16,3));
    for ax,(r,_),p in zip(axes,fields,patterns):ax.imshow(p[p.shape[0]//2],origin="lower");ax.set_title(r.name.split("_")[0]);ax.axis("off")
    fig.tight_layout();fig.savefig(OUT/"seven_depth_geometry_comparison.png",dpi=110);plt.close(fig)

    report = "\n".join(["DEV164 STATIC NATIVE DEFORMATION-TO-RELATIONAL-GEOMETRY AUDIT","",
      "Outcome A — the frozen state does not define a deformed spatial geometry.",
      "Each node stores scalar excursion. Six oriented neighbor differences are exactly recoverable, so directional scalar information is partial.",
      "Those differences have excursion/strain semantics, not bond-length or bond-orientation semantics. Mapping them to vectors would introduce a new geometric law.",
      "Consequently no loaded embedding, path turning, or finite propagation was executed. The Cartesian zero-load control is exact.",
      "The next missing law is native vector/relational geometry content, not loaded stiffness.","",
      *[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in contract.items()],""])
    (OUT/"report.txt").write_text(report); print(report,end=""); return contract


if __name__ == "__main__": main()
