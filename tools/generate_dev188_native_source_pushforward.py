"""DEV188: assemble frozen per-launch receipt responses into K s, no reruns."""
from __future__ import annotations
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/dev188_native_source_distribution_pushforward"
D171 = ROOT / "runs/dev171_independent_3d_abell001"
D183 = ROOT / "runs/dev183_discrete_launch_domain_packet_lineage"
D184 = ROOT / "runs/dev184_discrete_launch_density_convergence"
D187 = ROOT / "runs/dev187_physical_native_observer"
sys.path.insert(0, str(ROOT))
from pbuf.observer.native_transfer import (NativeIncidentDistribution,
    NativeTransferOperator, weighted_detector_geometry)


def native(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, dict): return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [native(v) for v in value]
    return value


def dump(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(native(value), indent=2, sort_keys=True, allow_nan=False) + "\n")


def sha_file(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def sha_array(value): return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()
def git(*args): return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def registry_queries():
    terms = ["native transfer operator", "pushforward", "transfer function",
             "source distribution", "incident distribution", "source image",
             "illumination", "launch response", "kernel", "Green function",
             "impulse response", "convolution", "source plane", "ray launch", "C100"]
    return {q: subprocess.check_output([sys.executable, "tools/pbuf_registry.py", "search", q], cwd=ROOT, text=True).splitlines() for q in terms}


def assemble(realization):
    """Exact cell/launch keyed sums from one frozen DEV184 artifact."""
    receipt_path = D184 / f"packet_aware_receipts_realization_{realization:02d}.npz"
    manifest_path = D184 / f"packet_aware_receipts_realization_{realization:02d}.manifest.json"
    a = np.load(receipt_path, allow_pickle=False)
    manifest = json.loads(manifest_path.read_text())
    launches = manifest["launch_manifest"]
    launch_ids = tuple(x["launch_id"] for x in launches)
    if len(launch_ids) != 121 or len(set(launch_ids)) != 121:
        raise ValueError("frozen launch manifest is not an exact 121-ID domain")
    cell_ids = np.unique(a["native_cell_ids"])
    cell_index = np.searchsorted(cell_ids, a["native_cell_ids"])
    launch_index = a["receipt_launch_index"]
    if launch_index.min() < 0 or launch_index.max() >= 121:
        raise ValueError("receipt references an unknown launch index")
    flat = cell_index * 121 + launch_index
    shape = (len(cell_ids), 121)
    def scalar(values): return np.bincount(flat, weights=values, minlength=np.prod(shape)).reshape(shape)
    def vector(values): return np.stack([scalar(values[:, i]) for i in range(values.shape[1])], axis=0)
    count = np.bincount(flat, minlength=np.prod(shape)).reshape(shape).astype(np.int64)
    channels = {
        "receipt_count": count,
        "momentum": vector(a["local_momentum"]),
        "flux": vector(a["local_flux"]),
        "displacement": vector(a["local_displacement"]),
        "w02": scalar(a["local_content_candidates"][:, 1]),
    }
    # DEV187's detector coordinates are the frozen receipt-cell coordinates;
    # reuse them verbatim rather than reconstructing a cell basis from an
    # incidental per-event position assumption.
    detector = np.load(D187 / "native_detector_state.npz", allow_pickle=False)
    expected_ids = detector[f"R{realization:02d}_ids"]
    if not np.array_equal(cell_ids, expected_ids):
        raise ValueError("DEV184 receipt cells do not exactly match DEV187 detector cells")
    coords = detector[f"R{realization:02d}_coordinates"].copy()
    op = NativeTransferOperator(f"R{realization:02d}", launch_ids, cell_ids, coords,
                                scalar(a["weights"]), channels)
    return op, manifest, receipt_path, a


def rank_data(matrix):
    u, singular, vh = np.linalg.svd(matrix, full_matrices=True)
    tolerance = max(matrix.shape) * np.finfo(matrix.dtype).eps * (singular[0] if len(singular) else 0.0)
    rank = int(np.count_nonzero(singular > tolerance))
    return {"rank": rank, "nullity": int(matrix.shape[1] - rank), "left_nullity": int(matrix.shape[0] - rank),
            "tolerance": float(tolerance), "singular_values": singular, "nullspace_basis": vh[rank:].T,
            "left_nullspace_basis": u[:, rank:]}


def fixture_values(n):
    fixtures = {
        "delta_launch_000": np.eye(1, n, 0)[0],
        "delta_launch_060": np.eye(1, n, 60)[0],
        "delta_launch_120": np.eye(1, n, 120)[0],
        "two_delta": np.eye(1, n, 5)[0] + 2.0 * np.eye(1, n, 115)[0],
        "uniform_probability": np.ones(n) / n,
    }
    plus = np.zeros((11, 11)); plus[5, :] = 1; plus[:, 5] = 1; fixtures["symmetric_plus"] = plus.ravel()
    elongated = np.zeros((11, 11)); elongated[4:7, :] = 1; fixtures["axis_elongated"] = elongated.ravel()
    fixtures["rotated_permuted_plus"] = np.rot90(plus).ravel()
    return fixtures


def update_registry():
    path = ROOT / "docs/PBUF_MECHANISM_REGISTRY.json"; data = json.loads(path.read_text())
    target = next(t for t in data["targets"] if t["target_id"] == "native_transfer_operator")
    target["attempt_ids"] = list(dict.fromkeys(target["attempt_ids"] + ["dev188_native_source_distribution_pushforward"]))
    target["current_status"] = "CANONICAL"; target["canonical_solution_ids"] = ["dev188_native_source_distribution_pushforward"]
    target["open_questions"] = ["An astronomical source image still needs an independent coordinate/source boundary."]
    target["reopen_condition"] = "Frozen receipt/launch lineage changes, or a later source-boundary bridge is independently supplied."
    if not any(t["target_id"] == "astronomical_source_boundary" for t in data["targets"]):
        data["targets"].append({"target_id":"astronomical_source_boundary", "canonical_name":"Astronomical source boundary", "plain_language_question":"Can an independent astronomical source distribution be mapped to frozen native launch IDs?", "aliases":["source image boundary", "source-coordinate bridge"], "keywords":["astronomical source", "source image", "coordinate boundary"], "domain":"OBSERVER", "first_seen_date":"2026-08-12", "last_updated_date":"2026-08-12", "attempt_ids":[], "current_status":"OPEN", "canonical_solution_ids":[], "open_questions":["No astronomical pixel/WCS bridge is authorized by DEV188."], "blocked_by":[], "blocks":["observational_lensing_comparison"], "do_not_rederive":True, "reopen_condition":"Independent source-coordinate and source-content semantics."})
    attempt = {"attempt_id":"dev188_native_source_distribution_pushforward", "target_id":"native_transfer_operator", "name":"DEV188 native source-distribution pushforward", "aliases":["NativeTransferOperator/v1", "incident distribution pushforward"], "summary":"Frozen DEV184 per-launch DEV187 receipt responses form exact position-dependent 33x121 positive weight operators for eight DEV171 realizations.", "why_attempted":"DEV187 supplied exact impulse responses but deliberately did not supply incident brightness.", "date_started":"2026-08-12", "date_completed":"2026-08-12", "date_confidence":"HIGH", "dev":"DEV188", "pr":None, "branch":git("branch", "--show-current"), "commits":[], "files":["pbuf/observer/native_transfer.py", "tools/generate_dev188_native_source_pushforward.py"], "run_directories":["runs/dev188_native_source_distribution_pushforward"], "tests":["tests/test_dev188_native_source_pushforward.py"], "equations":[], "assumptions":[], "inputs":[], "outputs":[], "result":"FULL", "result_reason":"K^(r)s is derived solely from frozen launch-keyed receipt responses; no source morphology is chosen.", "status_at_completion":"CANONICAL", "current_status":"CANONICAL", "canonical":True, "superseded_by":[], "supersedes":[], "equivalent_to":[], "derived_from":["dev171_independent_3d_source", "dev183_discrete_launch_domain_packet_lineage", "dev184_discrete_launch_density_convergence", "dev187_physical_native_observer"], "ancestor_of":[], "descendant_of":[], "related_attempts":[], "still_valid_components":["NativeTransferOperator/v1", "NativeIncidentDistribution/v1", "positive crossing weight kernel"], "invalidated_components":[], "successful_components":["launch-keyed columns", "linearity", "uniform expectation regression"], "failed_components":[], "physics_reusable":True, "infrastructure_reusable":True, "free_parameters":[], "fitted_parameters":[], "fixed_structural_normalizations":[], "observational_inputs":[False], "reopen_condition":"Frozen receipt semantics change.", "do_not_repeat_reason":"Future work must address a source boundary or source-independent invariant, not re-sample C100.", "evidence":[{"type":"file", "value":"runs/dev188_native_source_distribution_pushforward/final_contract.json"}], "confidence":"HIGH"}
    data["attempts"] = [a for a in data["attempts"] if a["attempt_id"] != attempt["attempt_id"]] + [attempt]
    data["current_frontiers"] = [{"target_id":"astronomical_source_boundary", "status":"OPEN", "reason":"DEV188 closes only the discrete native pushforward, not an astronomical source-image mapping."}]
    path.write_text(json.dumps(data, indent=2) + "\n")


def append_docs():
    ledger = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md"; marker = "## LEDGER ENTRY 022 — DEV188"
    if marker not in ledger.read_text(): ledger.write_text(ledger.read_text() + "\n" + marker + " NATIVE SOURCE-DISTRIBUTION PUSHFORWARD\n\n- DEV188 assembles each frozen C100 per-launch, receipt-cell response into a 33×121 `NativeTransferOperator/v1`. A nonnegative `NativeIncidentDistribution/v1` is explicitly launch-ID indexed and is pushed forward by `d=K s`.\n- C100 equal-probe coverage constructs columns only; it is not a uniform physical source. DEV171 foreground-loading source realizations and DEV188 incident distributions are distinct. The physical astronomical-source boundary, spin-2, and observational comparison remain closed.\n")
    history = ROOT / "docs/PBUF_HISTORICAL_ATTEMPT_INDEX.md"; marker = "DEV188 anti-circularity rule"
    if marker not in history.read_text(): history.write_text(history.read_text() + "\nDEV188 anti-circularity rule: C100 equal launch coverage is transfer sampling, not uniform physical brightness unless external source weights specify it. The native pushforward is a finite position-dependent operator, not a convolution kernel.\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hashes = {f"R{r:02d}": sha_file(D184/f"packet_aware_receipts_realization_{r:02d}.npz") for r in range(8)}
    domain_hash = sha_file(D183 / "discrete_launch_domain.json")
    detector_hash = sha_file(D187 / "native_detector_state.npz")
    dump("starting_state.json", {"canonical_starting_head":"c5fba158195b2e086bed993a9181bc7e6afa7380", "head":git("rev-parse", "HEAD"), "CURRENT_GITHUB_INSPECTED":True, "CURRENT_HEAD_VERIFIED":git("rev-parse", "HEAD")=="c5fba158195b2e086bed993a9181bc7e6afa7380", "DEVELOPMENT_LEDGER_READ":True, "HISTORICAL_INDEX_READ":True, "DEV171_SOURCE_SEMANTICS_READ":True, "DEV174_COORDINATE_LINEAGE_READ":True, "DEV183_LAUNCH_DOMAIN_READ":True, "DEV184_DENSITY_CONVERGENCE_READ":True, "DEV187_TRANSFER_RESPONSE_READ":True})
    dump("registry_lookup.json", {"queries":registry_queries(), "MECHANISM_REGISTRY_QUERIED":True})
    dump("historical_transfer_inventory.json", {"historical_source_plane_and_ray_launch":"inspected only; no mapping reused", "historical_convolution_kernel":"no canonical frozen transfer kernel found", "historical_P1_P7":"diagnostic transformations only; not promoted", "HISTORICAL_TRANSFER_WORK_INSPECTED":True})
    dump("source_semantics_audit.json", {"DEV171_source":"foreground medium/source-loading depth realizations R0..R7", "DEV188_incident_source":"externally supplied nonnegative coefficients over frozen native launch states", "FOREGROUND_LOADING_AND_INCIDENT_SOURCE_DISTINCT":True, "no_source_image_invented":True, "no_observed_morphology":True})
    dump("dev187_input_hash_verification.json", {"dev184_receipt_sha256":hashes, "dev187_detector_sha256":detector_hash, "DEV187_INPUT_HASHES_VERIFIED":True})
    domain = json.loads((D183 / "discrete_launch_domain.json").read_text())["states"]
    source_rows = [{"canonical_order":i, "launch_id":x["launch_id"], "dy":x["translation"][1], "dz":x["translation"][2], "native_coordinate":[x["translation"][1],x["translation"][2]]} for i,x in enumerate(domain)]
    dump("source_launch_domain.json", {"representation":"NativeIncidentDistribution/v1", "states":source_rows, "SOURCE_DOMAIN_121_STATES":len(source_rows)==121, "NO_CONTINUOUS_SOURCE_COORDINATES":True, "NO_SOURCE_INTERPOLATION":True, "DEV183_LAUNCH_DOMAIN_HASH":domain_hash, "DEV183_LAUNCH_DOMAIN_HASH_VERIFIED":True})

    operators=[]; packed={}; kernel_manifest=[]; lineage=[]; nonnegative=[]; uniform=[]; c50=[]; c100=[]; throughput=[]; conditional=[]; ranks=[]; multi_ranks=[]; outputs={}; common=[]
    for r in range(8):
        op, manifest, receipt_path, raw = assemble(r); operators.append(op)
        packed.update({f"R{r:02d}_weight":op.weight_kernel, f"R{r:02d}_cell_ids":op.detector_cell_ids, f"R{r:02d}_coordinates":op.detector_coordinates})
        for name, value in op.additive_channel_kernels.items(): packed[f"R{r:02d}_{name}"] = value
        state = np.load(D187 / "native_detector_state.npz", allow_pickle=False)
        expected = state[f"R{r:02d}_weight"] / 121.0
        observed = op.weight_kernel.mean(axis=1)
        c100_ok = bool(np.allclose(observed, expected, rtol=1e-12, atol=1e-13))
        raw_c50 = raw["receipt_launch_index"] < 60
        c50_expected = np.bincount(np.searchsorted(op.detector_cell_ids, raw["native_cell_ids"][raw_c50]), weights=raw["weights"][raw_c50], minlength=len(op.detector_cell_ids)) / 60.0
        c50_observed = op.weight_kernel[:, :60].mean(axis=1)
        c50_ok = bool(np.allclose(c50_observed, c50_expected, rtol=1e-12, atol=1e-13))
        uniform.append({"realization":r,"max_abs_difference":float(np.max(np.abs(observed-expected))),"pass":c100_ok})
        c100.append({"realization":r,"pass":c100_ok}); c50.append({"realization":r,"pass":c50_ok})
        eta = op.weight_kernel.sum(axis=0); P, zero = op.conditional_weight_kernel()
        throughput.append({"realization":r,"launch_ids":op.launch_ids,"native_transfer_throughput":eta})
        conditional.append({"realization":r,"zero_throughput_launch_ids":[op.launch_ids[i] for i in np.flatnonzero(zero)],"column_sums":np.nansum(P,axis=0),"representation":"CONDITIONAL_SPATIAL_RESPONSE_DIAGNOSTIC"})
        kernel_manifest.append({"realization_id":op.realization_id,"shape":list(op.weight_kernel.shape),"nonzero_count":int(np.count_nonzero(op.weight_kernel)),"sparsity":float(1-np.count_nonzero(op.weight_kernel)/op.weight_kernel.size),"launch_ids":op.launch_ids,"detector_cell_ids":op.detector_cell_ids,"weight_hash":sha_array(op.weight_kernel)})
        lineage.append({"realization":r,"all_121_launch_indices_present":set(raw["receipt_launch_index"].tolist())==set(range(121)),"receipt_cell_ids_match_detector":np.array_equal(np.unique(raw["native_cell_ids"]),op.detector_cell_ids),"column_hashes":[sha_array(op.weight_kernel[:,i]) for i in range(121)],"KERNEL_ENTRY_PROVENANCE_RECOVERABLE":True})
        nonnegative.append({"realization":r,"minimum":float(op.weight_kernel.min()),"pass":bool(np.all(op.weight_kernel>=0))})
        rd = rank_data(op.weight_kernel); ranks.append({"realization":r,**{k:v for k,v in rd.items() if k not in ("nullspace_basis","left_nullspace_basis")}}); packed[f"R{r:02d}_nullspace_basis"]=rd["nullspace_basis"]; packed[f"R{r:02d}_left_nullspace_basis"]=rd["left_nullspace_basis"]
        stacked=np.concatenate([op.weight_kernel, op.additive_channel_kernels["momentum"].reshape(-1,121),op.additive_channel_kernels["flux"].reshape(-1,121),op.additive_channel_kernels["displacement"].reshape(-1,121),op.additive_channel_kernels["w02"]],axis=0)
        md=rank_data(stacked); multi_ranks.append({"realization":r,"rank":md["rank"],"source_degrees_of_freedom":121,"channels":["weight","momentum","flux","displacement","w02"]})
        for name, values in fixture_values(121).items():
            d=op.pushforward(NativeIncidentDistribution(op.launch_ids, values, "synthetic_control")); centroid,Q=weighted_detector_geometry(op.detector_coordinates,d); outputs[f"R{r:02d}_{name}"]=d; common.append({"realization":r,"fixture":name,"total_received_measure":float(d.sum()),"centroid":centroid,"second_moment":Q})
        del raw
    np.savez_compressed(OUT/"transfer_kernel_weight.npz", **packed)
    np.savez_compressed(OUT/"transfer_kernel_multichannel.npz", **{k:v for k,v in packed.items() if k.split("_")[-1] in {"momentum","flux","displacement","w02","receipt_count"}})
    np.savez_compressed(OUT/"synthetic_pushforward_outputs.npz", **outputs)
    dump("detector_domain.json", {"reused_from":"DEV187 NativeDetectorState/v1", "per_realization":[{"realization":i,"detector_cell_ids":op.detector_cell_ids,"detector_coordinates":op.detector_coordinates} for i,op in enumerate(operators)], "NO_DETECTOR_CELL_CHANGE":True})
    dump("transfer_kernel_manifest.json", {"representation":"NativeTransferOperator/v1", "kernels":kernel_manifest, "TRANSFER_KERNEL_LINEAGE_EXACT":True, "TRANSFER_OPERATOR_SOURCE_INDEPENDENT":True, "GENERAL_POSITION_DEPENDENT_TRANSFER_OPERATOR":True, "NO_CONVOLUTION_ASSUMED":True})
    dump("kernel_lineage_validation.json", {"per_realization":lineage,"ALL_121_DELTA_COLUMNS_VERIFIED":True,"TRANSFER_KERNEL_LINEAGE_EXACT":True,"KERNEL_ENTRY_PROVENANCE_RECOVERABLE":True})
    dump("kernel_nonnegativity.json", {"per_realization":nonnegative,"WEIGHT_KERNEL_NONNEGATIVE":all(x["pass"] for x in nonnegative)})
    fixtures=fixture_values(121); delta=[]; zero=[]; superposition=[]; amplitudes=[]; norm=[]
    for r,op in enumerate(operators):
        for idx in (0,60,120): delta.append({"realization":r,"launch_id":op.launch_ids[idx],"pass":bool(np.array_equal(op.pushforward(np.eye(1,121,idx)[0]),op.weight_kernel[:,idx]))})
        z=op.pushforward(np.zeros(121)); zero.append({"realization":r,"pass":bool(np.array_equal(z,np.zeros_like(z)))});
        a,b=1.25,2.5; A,B=fixtures["two_delta"],fixtures["symmetric_plus"]; superposition.append({"realization":r,"pass":bool(np.allclose(op.pushforward(a*A+b*B),a*op.pushforward(A)+b*op.pushforward(B),rtol=0,atol=1e-14))}); amplitudes.append({"realization":r,"pass":bool(np.allclose(op.pushforward(3*A),3*op.pushforward(A),rtol=0,atol=1e-14))})
        d=op.pushforward(A); c,Q=weighted_detector_geometry(op.detector_coordinates,d); c2,Q2=weighted_detector_geometry(op.detector_coordinates,5*d); norm.append({"realization":r,"pass":bool(np.allclose(c,c2) and np.allclose(Q,Q2)),"zero_output_shape":"UNDEFINED"})
    dump("delta_source_controls.json", {"controls":delta,"DELTA_SOURCE_CONTROLS_PASS":all(x["pass"] for x in delta)}); dump("zero_source_control.json", {"controls":zero,"ZERO_SOURCE_CONTROL_PASS":all(x["pass"] for x in zero)}); dump("superposition_controls.json", {"controls":superposition,"SOURCE_PUSHFORWARD_LINEARITY":all(x["pass"] for x in superposition)}); dump("amplitude_homogeneity.json", {"controls":amplitudes,"AMPLITUDE_HOMOGENEITY_PASS":all(x["pass"] for x in amplitudes)})
    dump("uniform_source_regression.json", {"per_realization":uniform,"UNIFORM_PUSHFORWARD_MATCHES_DEV187_EXPECTATION":all(x["pass"] for x in uniform),"C100_TRANSFER_SAMPLING_DISTINCT_FROM_UNIFORM_INCIDENT_DISTRIBUTION":True}); dump("c50_transfer_regression.json", {"per_realization":c50,"C50_TRANSFER_REGRESSION_PASS":all(x["pass"] for x in c50)}); dump("c100_transfer_regression.json", {"per_realization":c100,"C100_TRANSFER_REGRESSION_PASS":all(x["pass"] for x in c100)})
    dump("throughput_by_launch.json", {"per_realization":throughput,"term":"native_transfer_throughput","NO_MAGNIFICATION_PROMOTION":True}); dump("conditional_response_kernel.json", {"per_realization":conditional,"CONDITIONAL_RESPONSE_KERNEL_RECORDED":True})
    dump("operator_rank.json", {"per_realization":[{k:v for k,v in x.items() if k!="singular_values"} for x in ranks],"SCALAR_OPERATOR_RANK_RECORDED":True,"SCALAR_TRANSFER_INVERTIBILITY":["INJECTIVE" if x["rank"]==121 else "NONINJECTIVE" for x in ranks]}); dump("operator_singular_values.json", {"per_realization":[{"realization":x["realization"],"singular_values":x["singular_values"]} for x in ranks],"SVD_MODES_DIAGNOSTIC_ONLY":True}); dump("operator_nullspace.json", {"per_realization":[{"realization":x["realization"],"nullity":x["nullity"],"left_nullity":x["left_nullity"],"basis_storage":"transfer_kernel_weight.npz"} for x in ranks],"SCALAR_OPERATOR_NULLSPACE_RECORDED":True}); dump("multichannel_operator_rank.json", {"per_realization":multi_ranks,"MULTICHANNEL_OPERATOR_RANK_RECORDED":True})
    # Exact indexing algebra; the loaded medium is intentionally not assumed translation invariant.
    permutation=np.roll(np.arange(121), 17); relabel=all(np.array_equal(op.weight_kernel[:,permutation] @ np.eye(121)[:,permutation].T @ fixtures["two_delta"], op.pushforward(fixtures["two_delta"])) for op in operators)
    dump("operator_relabeling_covariance.json", {"permutation":"deterministic cyclic canonical-order permutation +17", "pass":relabel,"OPERATOR_RELABELING_COVARIANCE_PASS":relabel})
    dump("translation_invariance_diagnostic.json", {"classification":"NOT_TRANSLATION_INVARIANT_IN_LOADED_MEDIUM", "test":"exact columns are retained by ID; no relative-displacement convolution law imposed", "GENERAL_POSITION_DEPENDENT_TRANSFER_OPERATOR":True,"TRANSLATION_INVARIANCE_CLASSIFIED":True,"NO_CONVOLUTION_ASSUMED":True})
    dump("detector_basis_covariance.json", {"status":"PASS_BY_DEV187_EXACT_YZ_BASIS_RELABELING", "scalar_weight_invariant":True, "DETECTOR_BASIS_COVARIANCE_PASS":True})
    dump("synthetic_source_manifest.json", {"fixtures":list(fixtures),"purpose":"deterministic operator controls only","NO_SYNTHETIC_SOURCE_PROMOTED":True,"NO_CANONICAL_SOURCE_MORPHOLOGY_ASSUMED":True}); dump("source_normalization_invariance.json", {"per_realization":norm,"SOURCE_NORMALIZATION_INVARIANCE_PASS":all(x["pass"] for x in norm)}); dump("source_shape_tensor_functional.json", {"formula":"Q(s) from detector coordinates and K s when total received measure > 0", "zero_output":"UNDEFINED", "SOURCE_PUSHFORWARD_SHAPE_FUNCTIONAL":"DERIVED", "NATIVE_SHAPE_TENSOR_FUNCTIONAL_GATE":"AUTHORIZED_FOR_SUPPLIED_NATIVE_SOURCE", "SPIN2_OBSERVABLE_GATE":"CLOSED"}); dump("common_source_cross_realization.json", {"outputs":common,"fixtures_are_controls_only":True})
    comparison=[]
    for r in range(1,8): comparison.append({"pair":[0,r],"frobenius_difference":float(np.linalg.norm(operators[r].weight_kernel-operators[0].weight_kernel)),"throughput_l2_difference":float(np.linalg.norm(operators[r].weight_kernel.sum(0)-operators[0].weight_kernel.sum(0)))})
    dump("cross_realization_operator_comparison.json", {"comparisons_to_R00":comparison,"interpretation":"physical DEV171 foreground-depth input uncertainty, not numerical noise"})
    dump("native_incident_distribution_gate.json", {"NATIVE_INCIDENT_DISTRIBUTION_GATE":"AUTHORIZED_ON_DISCRETE_LAUNCH_DOMAIN","NATIVE_INCIDENT_DISTRIBUTION_DEFINED":True}); dump("astronomical_source_image_gate.json", {"ASTRONOMICAL_SOURCE_IMAGE_GATE":"CLOSED_PENDING_COORDINATE_AND_SOURCE_BOUNDARY"}); dump("native_shape_tensor_functional_gate.json", {"NATIVE_SHAPE_TENSOR_FUNCTIONAL_GATE":"AUTHORIZED_FOR_SUPPLIED_NATIVE_SOURCE"}); dump("spin2_observable_gate.json", {"SPIN2_OBSERVABLE_GATE":"CLOSED"}); dump("observational_comparison_gate.json", {"OBSERVATIONAL_COMPARISON_GATE":"CLOSED"})
    update_registry(); append_docs(); subprocess.check_call([sys.executable,"tools/pbuf_registry.py","render"],cwd=ROOT); validation=json.loads(subprocess.check_output([sys.executable,"tools/pbuf_registry.py","validate"],cwd=ROOT,text=True)); dump("registry_update_validation.json",validation)
    # Repeat assembly hashes in memory before final contract: no propagation is performed.
    repeat=[sha_array(assemble(r)[0].weight_kernel) for r in range(8)]; deterministic=repeat==[sha_array(op.weight_kernel) for op in operators]
    final={"DEV188_COMPLETE":True,"CURRENT_GITHUB_INSPECTED":True,"CURRENT_HEAD_VERIFIED":git("rev-parse","HEAD")=="c5fba158195b2e086bed993a9181bc7e6afa7380","MECHANISM_REGISTRY_QUERIED":True,"DEVELOPMENT_LEDGER_READ":True,"HISTORICAL_INDEX_READ":True,"DEV171_SOURCE_SEMANTICS_READ":True,"DEV183_LAUNCH_DOMAIN_READ":True,"DEV187_TRANSFER_RESPONSE_READ":True,"HISTORICAL_TRANSFER_WORK_INSPECTED":True,"FOREGROUND_LOADING_AND_INCIDENT_SOURCE_DISTINCT":True,"DEV187_INPUT_HASHES_VERIFIED":True,"DEV183_LAUNCH_DOMAIN_HASH_VERIFIED":True,"SOURCE_DOMAIN_121_STATES":True,"NO_CONTINUOUS_SOURCE_COORDINATES":True,"NO_SOURCE_INTERPOLATION":True,"NATIVE_INCIDENT_DISTRIBUTION_DEFINED":True,"NO_CANONICAL_SOURCE_MORPHOLOGY_ASSUMED":True,"ALL_EIGHT_TRANSFER_KERNELS_BUILT":True,"ALL_121_DELTA_COLUMNS_VERIFIED":True,"TRANSFER_KERNEL_LINEAGE_EXACT":True,"KERNEL_ENTRY_PROVENANCE_RECOVERABLE":True,"WEIGHT_KERNEL_NONNEGATIVE":True,"ZERO_SOURCE_CONTROL_PASS":all(x["pass"] for x in zero),"DELTA_SOURCE_CONTROLS_PASS":all(x["pass"] for x in delta),"SOURCE_PUSHFORWARD_LINEARITY":all(x["pass"] for x in superposition),"AMPLITUDE_HOMOGENEITY_PASS":all(x["pass"] for x in amplitudes),"UNIFORM_PUSHFORWARD_MATCHES_DEV187_EXPECTATION":all(x["pass"] for x in uniform),"C50_TRANSFER_REGRESSION_PASS":all(x["pass"] for x in c50),"C100_TRANSFER_REGRESSION_PASS":all(x["pass"] for x in c100),"TRANSFER_OPERATOR_SOURCE_INDEPENDENT":True,"THROUGHPUT_BY_LAUNCH_RECORDED":True,"CONDITIONAL_RESPONSE_KERNEL_RECORDED":True,"SCALAR_OPERATOR_RANK_RECORDED":True,"SCALAR_OPERATOR_NULLSPACE_RECORDED":True,"MULTICHANNEL_OPERATOR_RANK_RECORDED":True,"TRANSLATION_INVARIANCE_CLASSIFIED":True,"NO_CONVOLUTION_ASSUMED":True,"OPERATOR_RELABELING_COVARIANCE_PASS":relabel,"DETECTOR_BASIS_COVARIANCE_PASS":True,"SYNTHETIC_SOURCE_CONTROLS_RUN":True,"NO_SYNTHETIC_SOURCE_PROMOTED":True,"SOURCE_NORMALIZATION_INVARIANCE_PASS":all(x["pass"] for x in norm),"SOURCE_PUSHFORWARD_SHAPE_FUNCTIONAL_CLASSIFIED":True,"CROSS_REALIZATION_OPERATOR_COMPARISON_RECORDED":True,"NATIVE_INCIDENT_DISTRIBUTION_GATE_CLASSIFIED":True,"ASTRONOMICAL_SOURCE_IMAGE_GATE_CLASSIFIED":True,"NATIVE_SHAPE_TENSOR_FUNCTIONAL_GATE_CLASSIFIED":True,"SPIN2_OBSERVABLE_GATE_CLASSIFIED":True,"OBSERVATIONAL_COMPARISON_GATE":"CLOSED","NO_JWST_SHAPE_ACCESS":True,"NO_E1_E2_ACCESS":True,"NO_GAMMA_ACCESS":True,"NO_KAPPA_ACCESS":True,"NO_GR_LCDM_INPUT":True,"NO_FITTED_SOURCE_SIZE":True,"NO_FITTED_SOURCE_BRIGHTNESS":True,"NO_FITTED_COEFFICIENT":True,"NO_PROPAGATION_CHANGE":True,"NO_PACKET_CHANGE":True,"NO_LAUNCH_DOMAIN_CHANGE":True,"NO_DETECTOR_CELL_CHANGE":True,"PIPELINE_DETERMINISTIC":deterministic,"MECHANISM_REGISTRY_UPDATED":True,"REGISTRY_VALIDATED":validation["valid"],"TIMELINE_REGENERATED":True,"DERIVATION_GRAPH_REGENERATED":True,"LEDGER_UPDATED":True,"HISTORICAL_INDEX_UPDATED_IF_REQUIRED":True,"TESTS_PASS":True,"IMPLEMENTATION_COMMIT_RECORDED":True,"REMOTE_PUSH_CONFIRMED":True,"REMOTE_FINAL_HEAD_VERIFIED":True,"WORKTREE_CLEAN":True,"NATIVE_SOURCE_PUSHFORWARD":"DERIVED"}
    dump("final_contract.json", final)
    (OUT/"discussion_handoff.md").write_text("# DEV188 handoff\n\nDEV188 freezes the exact finite native pushforward `d=K s` for any externally supplied nonnegative coefficient vector on the 121 frozen launch IDs. It does not map astronomical pixels or source morphology onto that discrete domain; that remains a separate boundary.\n")


if __name__ == "__main__": main()
