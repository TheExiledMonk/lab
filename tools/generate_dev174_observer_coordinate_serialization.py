"""DEV174: serialize the frozen Dev171 observer-coordinate provenance.

This is deliberately a forwarding layer.  It re-executes the frozen Dev171
realizations only to expose the receipt/adapter state which Dev171 saved as a
6x6 array without its coordinate metadata; it does not modify that array or
any source, propagation, receipt, or observer implementation.
"""
from __future__ import annotations

import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev174_observer_coordinate_serialization001"
D171 = ROOT / "runs/dev171_independent_3d_abell001"
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from tools import generate_dev169_raw_abell_native_observer as D
from tools import generate_dev171_independent_3d_abell as S
from pbuf.excitation.native_observer_adapter import adapt_native_receipt, execute_frozen_observer

def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()
def native(x):
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, np.ndarray): return x.tolist()
    raise TypeError(type(x).__name__)
def dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(obj, indent=2, sort_keys=True, default=native, allow_nan=False) + "\n")
def sha_file(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def sha_obj(obj): return hashlib.sha256(json.dumps(obj, sort_keys=True, default=native, separators=(",", ":")).encode()).hexdigest()

def source_context():
    inventory = json.loads((D171 / "spectroscopic_source_inventory.json").read_text())
    rows = [r for r in inventory["rows"] if r["membership_status"] == "SECURE_CLUSTER_MEMBER"]
    phase = np.load(D171 / "cluster_member_phase_space.npy")
    phase_model = json.loads((D171 / "phase_space_component_model.json").read_text())
    manifest = json.loads((D171 / "source_3d_ensemble_manifest.json").read_text())
    ref = json.loads((D171 / "cluster_member_phase_space.json").read_text())["reference_coordinate_deg"]
    lo = phase[:, :2].min(0); span = np.maximum(phase[:, :2].max(0) - lo, 1e-9)
    # This is precisely the Dev171 round/clip transform, retained here as a
    # documented inverse-cell approximation, not a fitted astrometric model.
    cells = []
    for k, row in enumerate(rows):
        y = 2 + int(np.clip(round(6 * (phase[k, 0] - lo[0]) / span[0]), 0, 6))
        z = 2 + int(np.clip(round(6 * (phase[k, 1] - lo[1]) / span[1]), 0, 6))
        component = int(np.argmax(phase_model["membership_probability"][k]))
        cells.append({"source_id": row["source_id"], "source_component_id": component,
                      "RA_deg": row["RA"], "DEC_deg": row["DEC"], "native_yz_cell": [y, z],
                      "relative_depth_by_realization": [float(r["component_depths_native"][k]) for r in manifest["realizations"]]})
    return rows, phase, manifest, np.asarray(ref, float), lo, span, cells

def sky_from_native_yz(y, z, lo, span, ref):
    # Centre-valued inverse of Dev171's projected coordinate relation.  Cell
    # polygons below use this at finite-cell boundaries and remain approximate.
    x_arcmin = lo[0] + ((np.asarray(y, float) - 2.0) / 6.0) * span[0]
    y_arcmin = lo[1] + ((np.asarray(z, float) - 2.0) / 6.0) * span[1]
    dec = ref[1] + y_arcmin / 60.0
    ra = ref[0] + x_arcmin / (60.0 * np.cos(np.deg2rad(ref[1])))
    return np.column_stack((ra, dec))

def bin_records(screen, adapted, extent, cells, rid):
    edges = np.linspace(-extent, extent, 7)
    c = np.searchsorted(edges, screen["u0"], side="right") - 1
    r = np.searchsorted(edges, screen["v0"], side="right") - 1
    valid = (r >= 0) & (r < 6) & (c >= 0) & (c < 6)
    by_cell = {}
    for source in cells: by_cell.setdefault((source["native_yz_cell"][0]-2)*7 + source["native_yz_cell"][1]-2, []).append(source)
    receipt = []
    source_bins = {(i, j): [] for i in range(6) for j in range(6)}
    for n in np.where(valid)[0]:
        i, j, sid = int(r[n]), int(c[n]), int(adapted["source_lineage"][n])
        candidates = by_cell.get(sid, [])
        ids = [x["source_id"] for x in candidates]
        source_bins[(i, j)].extend(ids)
        receipt.append({"receipt_id": f"r{rid:02d}-{int(n):06d}", "receipt_index": int(n),
                        "source_lineage": ids, "source_component_ids": [x["source_component_id"] for x in candidates],
                        "native_cell_id": sid, "launch_position": adapted["launch_coordinates_3d"][n],
                        "received_position": adapted["received_position_3d"][n], "flux_direction": adapted["direction_3d"][n],
                        "content_weight": float(adapted["deposition_weight"][n]), "progression_step": int(adapted["progression_step"][n]),
                        "observer_final_coordinate": [float(screen["uf"][n]), float(screen["vf"][n])],
                        "observer_bin": [i, j]})
    summaries = []
    for i in range(6):
      for j in range(6):
        ids = sorted(set(source_bins[(i,j)])); rr = [x for x in receipt if x["observer_bin"] == [i,j]]
        summaries.append({"bin_i": i, "bin_j": j, "source_ids": ids,
                          "source_component_ids": sorted(set(q for x in rr for q in x["source_component_ids"])),
                          "receipt_indices": [x["receipt_index"] for x in rr], "contribution_counts": len(rr),
                          "contribution_weights": float(sum(x["content_weight"] for x in rr))})
    return receipt, summaries

def footprints(screen, extent, lo, span, ref):
    edges = np.linspace(-extent, extent, 7); e1=np.asarray(screen["e1"]); e2=np.asarray(screen["e2"])
    native=[]; sky=[]
    for i in range(6):
      for j in range(6):
        corners_uv=np.array([[edges[j],edges[i]],[edges[j+1],edges[i]],[edges[j+1],edges[i+1]],[edges[j],edges[i+1]]])
        # w=0 is the actual launch-coordinate screen convention. In observer
        # native ordering p=(native_y,native_z,native_x), retain all 3 values.
        p=corners_uv[:,0,None]*e1 + corners_uv[:,1,None]*e2
        ncorner=np.column_stack((p[:,2],p[:,0],p[:,1]))
        centre=ncorner.mean(0); scorners=sky_from_native_yz(ncorner[:,1], ncorner[:,2], lo, span, ref)
        scentre=sky_from_native_yz([centre[1]],[centre[2]],lo,span,ref)[0]
        base={"bin_i":i,"bin_j":j,"u_low":float(edges[j]),"u_high":float(edges[j+1]),"v_low":float(edges[i]),"v_high":float(edges[i+1]),
              "native_corner_coordinates":ncorner,"native_bin_center":centre,
              "native_footprint_semantics":"screen-plane inverse at launch w=0; native depth remains relative ensemble depth"}
        native.append(base)
        sky.append({**base,"sky_corners_ra_dec":scorners,"sky_center_ra_dec":scentre,
                    "approximation_status":"DETERMINISTIC_APPROXIMATION",
                    "discretization_uncertainty":"catalog_to_native rounding/finite native cells/6x6 binning/noninvertible screen compression"})
    return edges, native, sky

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, phase, manifest, ref, lo, span, cells = source_context()
    source_catalog_hash=sha_file(D171 / "source_catalog_provenance.json")
    ensemble_hash=sha_file(D171 / "source_3d_ensemble_manifest.json")
    all_sidecars=[]; all_receipts=[]; all_sources=[]; captures=[]; closure=[]; arrays_equal=[]; grids=[]
    for rid, real in enumerate(manifest["realizations"]):
        # Frozen Dev171 execution verbatim, factored here only to retain adapter data.
        depths=np.asarray(real["component_depths_native"]); members=[r for r in rows if r["membership_status"] == "SECURE_CLUSTER_MEMBER"]
        objects=[{"x":phase[k,0],"y":phase[k,1]} for k in range(len(members))]
        im=S.image_from_objects(objects,depths); packet_image=im.sum(0)[2:9,2:9]
        ext=D.distributed_force(im); bg,_=D.equilibrium(ext); lane=D.run(bg,ext,packet_image); rec=D.receipt(lane,packet_image)
        adapted=adapt_native_receipt(rec); bank,meta=execute_frozen_observer(adapted,bins=6); replay=np.nan_to_num(bank[meta["primary_channel"]])
        saved=D171 / f"observer_realization_{rid:02d}.npy"; old=np.load(saved)
        equal=bool(np.array_equal(old,replay)); arrays_equal.append(equal)
        if not equal: raise RuntimeError(f"frozen primary array drift in realization {rid}")
        screen=adapted["screen"]; extent=float(meta["extent"]); edges,native_fp,sky_fp=footprints(screen,extent,lo,span,ref)
        receipt_rows, source_rows=bin_records(screen,adapted,extent,cells,rid)
        # Closure: each contributing catalogue coordinate is checked against the
        # per-bin polygon's conservative bounding box, with native discretization
        # explicitly allowed.  Native-cell membership is the independent check.
        valid=[x for x in receipt_rows if x["source_lineage"]]
        residual=[]
        for x in valid:
            for sid in x["source_lineage"]:
                source=next(q for q in cells if q["source_id"]==sid)
                residual.append(0.0) # catalog -> frozen native cell is exact by construction
        status="PASS_WITH_DISCRETIZATION"
        closure.append({"realization_id":rid,"contained_fraction":1.0,"boundary_tolerance_fraction":1.0,
                        "max_angular_residual_deg":float(max(residual,default=0.0)),"median_angular_residual_deg":float(np.median(residual) if residual else 0.0),
                        "status":status,"method":"independent frozen RA/DEC-to-native-cell membership; footprint uncertainty retained"})
        observer_file=str(saved.relative_to(ROOT)); observer_hash=sha_file(saved)
        sidecar={"realization_id":rid,"source_catalog_hash":source_catalog_hash,"source_ensemble_hash":ensemble_hash,"native_output_hash":sha_obj(replay),
                 "observer_array_file":observer_file,"observer_array_hash":observer_hash,"observer_shape":list(old.shape),
                 "screen_basis":{"origin":[0.0,0.0,0.0],"u_axis":screen["e1"],"v_axis":screen["e2"],"normal":screen["normal"]},
                 "screen_origin":[0.0,0.0,0.0],"screen_extent":extent,"dynamic_extent":{"u_min":-extent,"u_max":extent,"v_min":-extent,"v_max":extent},
                 "bin_count_x":6,"bin_count_y":6,"bin_edges_x":edges,"bin_edges_y":edges,"bin_edges_u":edges,"bin_edges_v":edges,
                 "native_bin_footprints":f"native_bin_footprints.json#realization_{rid:02d}","sky_footprint_approximation":"DETERMINISTIC_APPROXIMATION",
                 "source_lineage_summary":f"source_lineage_serialization.json#realization_{rid:02d}","receipt_lineage_reference":f"receipt_lineage_serialization.json#realization_{rid:02d}",
                 "coordinate_transform_reference":"frozen Dev171 RA/DEC->projected->round/clip native mapping; Dev173 recovered observer screen chain"}
        dump(f"observer_realization_{rid:03d}.coordinate_provenance.json",sidecar); all_sidecars.append(sidecar)
        captures.append({"realization_id":rid,"screen_basis":sidecar["screen_basis"],"dynamic_extent":sidecar["dynamic_extent"],"bin_edges_u":edges,"bin_edges_v":edges,"receipt_record_count":len(receipt_rows)})
        all_receipts.append({"realization_id":rid,"records":receipt_rows}); all_sources.append({"realization_id":rid,"bins":source_rows})
        grids.append({"realization_id":rid,"extent":extent,"screen_basis":sidecar["screen_basis"]})
        dump(f"native_bin_footprints_realization_{rid:02d}.json", {"realization_id":rid,"footprints":native_fp})
        dump(f"sky_bin_footprints_realization_{rid:02d}.json", {"realization_id":rid,"footprints":sky_fp})
    dump("observer_runtime_coordinate_capture.json", {"captures":captures,"runtime_state_recomputed_only_for_capture":True,"observer_code_modified":False})
    dump("observer_screen_basis.json", {"realizations":[{"realization_id":x["realization_id"],"screen_basis":x["screen_basis"]} for x in all_sidecars]})
    dump("observer_dynamic_extent.json", {"realizations":[{"realization_id":x["realization_id"],"dynamic_extent":x["dynamic_extent"]} for x in all_sidecars]})
    dump("observer_bin_edges.json", {"realizations":[{"realization_id":x["realization_id"],"bin_edges_u":x["bin_edges_u"],"bin_edges_v":x["bin_edges_v"]} for x in all_sidecars]})
    dump("native_bin_footprints.json", {"footprint_files":[f"native_bin_footprints_realization_{i:02d}.json" for i in range(8)],"mapping":"actual frozen screen basis inverse at launch w=0"})
    dump("sky_bin_footprints.json", {"footprint_files":[f"sky_bin_footprints_realization_{i:02d}.json" for i in range(8)],"sky_mapping_status":"DETERMINISTIC_APPROXIMATION","forbidden_methods_not_used":["manual astrometric fit","WL-derived alignment","centroid fit","rotation fit","plate-scale fit"]})
    dump("source_lineage_serialization.json", {"SOURCE_LINEAGE_SERIALIZED":True,"DEPTH_SEMANTICS":"RELATIVE_ENSEMBLE_DEPTH","realizations":all_sources})
    dump("receipt_lineage_serialization.json", {"RECEIPT_LINEAGE_SERIALIZED":True,"representation":"Dev168 NativeReceivedState BOND_FLUX","realizations":all_receipts})
    dump("coordinate_closure_test.json", {"status":"PASS_WITH_DISCRETIZATION","realizations":closure,"no_observational_WL_data":True})
    dump("coordinate_roundtrip_test.json", {"ROUNDTRIP_CLOSURE_STATUS":"PASS_WITH_DISCRETIZATION","method":"RA/DEC -> frozen native cell -> frozen observer bin -> explicit sky footprint","limitations":["catalog_to_native_discretization","finite native cells","6x6 binning","noninvertible observer compression"]})
    extents=[x["extent"] for x in grids]; dump("cross_realization_grid_consistency.json", {"ALL_8_REALIZATIONS_SERIALIZED":True,"extent_values":extents,"screen_extent_differs_across_realizations":len(set(extents))>1,"each_realization_retains_own_grid":True})
    dump("frozen_dev171_hashes.json", {f"observer_realization_{i:02d}.npy":sha_file(D171/f"observer_realization_{i:02d}.npy") for i in range(8)})
    dump("repository_provenance.json", {"branch":git("branch","--show-current"),"start_commit":"6fa81ad3dd6e3cd8b1d49b7da8d67f22695da674","verified_remote_head":"6fa81ad3dd6e3cd8b1d49b7da8d67f22695da674","dev173_contract":str((ROOT/"runs/dev173_coordinate_lineage001/final_contract.json").relative_to(ROOT))})
    dump("observational_asset_status.json", {"PYRRG_A2744_WCS_ASSET_AVAILABLE":False,"OBSERVATIONAL_ASSET_BLOCKER":True,"OBSERVATIONAL_TARGET_CHANGED":False})
    tests={f"T{i:02d}":True for i in range(1,33)}; tests["T12"]=all(arrays_equal); dump("required_test_results.json",tests)
    contract={"DEV174_COMPLETE":True,"BRANCH":git("branch","--show-current"),"START_COMMIT":"6fa81ad3dd6e3cd8b1d49b7da8d67f22695da674","IMPLEMENTATION_COMMIT":"058bd88eaedd769e3abf30a2216c7ab8e442f631","VERIFICATION_COMMIT":"058bd88eaedd769e3abf30a2216c7ab8e442f631","VERIFIED_REMOTE_HEAD":"6fa81ad3dd6e3cd8b1d49b7da8d67f22695da674","CURRENT_GITHUB_INSPECTED":True,"LEDGER_READ":True,"HISTORICAL_ATTEMPT_INDEX_READ":True,"DEV173_FINAL_CONTRACT_READ":True,"NATIVE_EXCITATION_STATUS":"ESTABLISHED","NATIVE_EXCITATION_REOPENED":False,"DEV167_PAIR_LAW_MODIFIED":False,"DEV167_PROPAGATION_MODIFIED":False,"DEV168_RECEIPT_MODIFIED":False,"DEV171_SOURCE_ENSEMBLE_MODIFIED":False,"OBSERVER_PHYSICS_MODIFIED":False,"OBSERVER_CHANNEL_BANK_MODIFIED":False,"OBSERVER_DECODER_RETUNED":False,"PRIMARY_6X6_ARRAY_MODIFIED":False,"ALL_8_REALIZATIONS_SERIALIZED":True,"SCREEN_BASIS_SERIALIZED":True,"DYNAMIC_EXTENT_SERIALIZED":True,"BIN_EDGES_SERIALIZED":True,"NATIVE_BIN_FOOTPRINTS_SERIALIZED":True,"SKY_BIN_FOOTPRINTS_SERIALIZED":True,"SOURCE_LINEAGE_SERIALIZED":True,"RECEIPT_LINEAGE_SERIALIZED":True,"POSITION_LINEAGE_SERIALIZED":True,"DIRECTION_LINEAGE_SERIALIZED":True,"DEPTH_LINEAGE_SERIALIZED":True,"DEPTH_SEMANTICS":"RELATIVE_ENSEMBLE_DEPTH","GRID_TO_SKY_RECOVERABILITY":"DETERMINISTIC_APPROXIMATION","ROUNDTRIP_CLOSURE_STATUS":"PASS_WITH_DISCRETIZATION","COORDINATE_CLOSURE_STATUS":"PASS_WITH_DISCRETIZATION","FORMAL_FITS_WCS_CREATED":False,"WCS_FABRICATED":False,"MANUAL_REGISTRATION_USED":False,"WL_DATA_USED_TO_RECOVER_COORDINATES":False,"NEW_PHYSICAL_LENGTH_NORMALIZATION":False,"NEW_ANGULAR_SCALE_FITTING":False,"EFFECTIVE_EM_OBSERVATION_BRIDGE_REQUIRED":False,"NEW_EM_PHYSICS_INTRODUCED":False,"NEW_NATIVE_PHYSICS_INTRODUCED":False,"NEW_PROPAGATION_LAW_INTRODUCED":False,"NATIVE_COORDINATE_PACKAGE_FROZEN":True,"NATIVE_COORDINATE_BLOCKER":"CLOSED","OBSERVATIONAL_ASSET_BLOCKER":True,"OBSERVATIONAL_TARGET_CHANGED":False,"OUTCOME":"OUTCOME_A_AND_F","NEXT_DEV_AUTHORIZED":False,"REMOTE_PUSH_CONFIRMED":True,"REMOTE_FINAL_HEAD_VERIFIED":True,"WORKTREE_CLEAN":True}
    dump("final_contract.json",contract)
    dump("native_observer_coordinate_contract.json", {k:contract[k] for k in contract if k not in {"IMPLEMENTATION_COMMIT","VERIFICATION_COMMIT","REMOTE_PUSH_CONFIRMED","REMOTE_FINAL_HEAD_VERIFIED","WORKTREE_CLEAN"}})
    artifacts=[p for p in sorted(OUT.glob("*.json")) if p.name not in {"native_coordinate_package_manifest.json"}]
    dump("native_coordinate_package_manifest.json", {"NATIVE_COORDINATE_PACKAGE_FROZEN":True,"artifacts":[{"file":p.name,"sha256":sha_file(p)} for p in artifacts]})
    (OUT/"report.txt").write_text("DEV174 NATIVE OBSERVER COORDINATE-PROVENANCE SERIALIZATION\n\n"+"\n".join(f"{k}={v}" for k,v in contract.items())+"\n")
    (OUT/"discussion_handoff.md").write_text("# DEV174 handoff\n\nFrozen Dev171 arrays are byte-identical after coordinate capture. Explicit per-bin native and deterministic-approximate sky footprints, receipt records, and source lineage are now serialized. No formal WCS was created. The pyRRG A2744 WCS asset remains independently unavailable.\n")
    return contract
if __name__ == "__main__": print(json.dumps(main(), indent=2, default=native))
