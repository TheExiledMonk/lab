"""Generate the read-only DEV166 canonical-state audit artifacts.

This script serializes repository evidence; it does not execute physics,
cosmology, lensing, or observer code.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/current_pbuf_missing_piece_audit001"
OUT.mkdir(parents=True, exist_ok=True)


def dump(name: str, value) -> None:
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


START = "f56e6a8bea77572ba80772e81806a68742c864d1"
FINAL = "7d919178b29e27832e231b740b9b6839f1805c13"
BRANCH = "dev-doc-112-fullscale-vulkan-observer-validation"
EVIDENCE = {
    "dev156": "runs/native_n6_relational_stress_dynamics001/report.txt",
    "dev157": "runs/native_propagation_spatial_scale_dispersion001/report.txt",
    "dev158": "runs/native_unified_medium_excursion001/report.txt",
    "dev159": "runs/native_source_medium_interaction001/report.txt",
    "dev160": "runs/raw_abell2744_simple_lensing_baseline001/report.txt",
    "dev161": "runs/raw_abell2744_detector_to_native_source001/report.txt",
    "dev162": "runs/raw_abell2744_3d_source_ambiguity_native_lens001/report.txt",
    "dev163": "runs/raw_abell2744_finite_native_lensing_gate001/report.txt",
    "dev164": "runs/static_native_deformation_to_relational_geometry001/report.txt",
    "dev165": "runs/native_medium_interaction_wide_net001/report.txt",
}

dump("repository_sync_contract.json", {
    "LOCAL_HEAD": START, "CURRENT_BRANCH": BRANCH,
    "REMOTE_DEFAULT_HEAD_AT_AUDIT_START": "231e4519532e7984eedb80c1b807df832bd41045",
    "REMOTE_DEFAULT_BRANCH": "main", "LOCAL_AHEAD_COUNT_VS_TRACKING": 0,
    "LOCAL_BEHIND_COUNT_VS_TRACKING": 0, "LOCAL_AHEAD_COUNT_VS_DEFAULT": 0,
    "LOCAL_BEHIND_COUNT_VS_DEFAULT": 2, "UNCOMMITTED_FILE_COUNT": 0,
    "UNTRACKED_FILE_COUNT_AT_START": 52,
    "DEV156_165_AT_START": "UNCOMMITTED_LOCAL",
    "UNRELATED_WORKSPACE_CHANGES_PRESERVED": True,
    "FINAL_COMMIT": FINAL, "REMOTE_PUSH_CONFIRMED": True,
    "REMOTE_VERIFICATION_COMMIT": "SELF (commit containing this post-push verification record)",
})
dump("historical_boundary_contract.json", {
    "LEDGER_EPOCH": "POST_RESTORED_N6_CURRENT_STATE", "PRE_LEDGER_HISTORY_CANONICAL": False,
    "PRE_LEDGER_HISTORY": "HISTORICAL_ONLY",
    "warning": "Earlier development includes exploratory, superseded, and lower-dimensional surrogate implementations. Historical results must not override the current restored N6 code unless explicitly revalidated or preserved."
})

modules = []
for dev, evidence in EVIDENCE.items():
    modules.append({"dev": dev.upper(), "classification": ["FOUND_IN_CURRENT_CODE", "FOUND_IN_CURRENT_ARTIFACTS", "FOUND_IN_CURRENT_TESTS"], "evidence": evidence})
dump("canonical_module_inventory.json", {
    "ledger_epoch": "POST_RESTORED_N6_CURRENT_STATE", "dev156_165": modules,
    "active_groups": ["pbuf/excitation native N6", "pbuf/data RAW bridge", "pbuf/source 3D family", "pbuf/lens stationary response", "pbuf/wl historical geometric observer stack"],
})

nodes = ["RAW archive", "FLT/FLC calibration", "common WCS source representation", "native projected source constraint", "3D source ambiguity family", "stationary native source response", "stationary native lens", "free finite native state", "native propagation condition", "loaded directional response", "finite loaded propagation", "received native 3D state", "observer receipt", "observer reduction", "2D observer output", "physical normalization", "observational comparison"]
edge_rows = [
 ("E01",nodes[0],nodes[1],"IMPLEMENTED"),("E02",nodes[1],nodes[2],"IMPLEMENTED"),
 ("E03",nodes[2],nodes[3],"IMPLEMENTED"),("E04",nodes[3],nodes[4],"PARTIAL"),
 ("E05",nodes[4],nodes[5],"IMPLEMENTED"),("E06",nodes[5],nodes[6],"IMPLEMENTED"),
 ("E07",nodes[5],nodes[7],"DERIVED_NOT_CONNECTED"),("E08",nodes[7],nodes[8],"PARTIAL"),
 ("E09",nodes[6],nodes[9],"BLOCKED"),("E10",nodes[8],nodes[9],"MISSING"),
 ("E11",nodes[9],nodes[10],"BLOCKED"),("E12",nodes[10],nodes[11],"BLOCKED"),
 ("E13",nodes[11],nodes[12],"BLOCKED"),("E14",nodes[12],nodes[13],"IMPLEMENTED"),
 ("E15",nodes[13],nodes[14],"IMPLEMENTED"),("E16",nodes[14],nodes[15],"MISSING"),
 ("E17",nodes[15],nodes[16],"BLOCKED")]
edges = [{"EDGE_ID": i, "FROM": a, "TO": b, "STATUS": s} for i,a,b,s in edge_rows]
dump("pipeline_dependency_graph.json", {"directed_acyclic": True, "nodes": nodes, "edges": edges})
missing = []
for e in edges:
    if e["STATUS"] not in {"IMPLEMENTED", "DERIVED_NOT_CONNECTED"}:
        root = e["EDGE_ID"] in {"E04", "E08", "E10", "E16"}
        missing.append({**e, "FIRST_BLOCKING_MODULE": "pbuf/labs/native_mechanisms/candidate_base.py" if e["EDGE_ID"] in {"E09","E10","E11"} else "UNRESOLVED",
            "FIRST_BLOCKING_FUNCTION": "evaluate_candidates" if e["EDGE_ID"] in {"E09","E10","E11"} else "UNRESOLVED",
            "FIRST_BLOCKING_ARTIFACT": EVIDENCE["dev165"] if e["EDGE_ID"] in {"E09","E10","E11"} else "UNRESOLVED",
            "UPSTREAM_REQUIREMENTS": [], "DOWNSTREAM_DEPENDENTS": [], "ROOT_CANDIDATE": root})
dump("missing_edge_inventory.json", missing)

roots = [
 {"id":"RB1", "name":"native loaded propagation mechanism", "status":"UNDERDETERMINED", "immediate_manifestation":"no derived loaded directional propagation response", "candidates":["state-dependent propagation condition","conservative bond-flux law","directional allocation law","missing relation/separation state","missing native interaction primitive"]},
 {"id":"RB2", "name":"absolute physical normalization", "status":"UNRESOLVED", "critical_for_structural_observer_path":False},
 {"id":"RB3", "name":"source depth/3D uniqueness", "status":"NON_UNIQUE", "critical_for_diagnostic_family_path":False},
]
dump("root_blocker_graph.json", {"ROOT_BLOCKER_COUNT":3,"ROOT_BLOCKERS":roots,"critical_chain":["received native state absent","finite loaded propagation absent","loaded directional response absent","deeper mechanism: MULTIPLE_CANDIDATES / UNDERDETERMINED"]})
dump("blocker_clusters.json", {"clusters":[
 {"name":"native propagation mechanism","root":"RB1","aliases":["directional allocation","routing","redistribution","scattering","neighbor transfer","bond flux","cross-coupling","loaded redirection","trajectory steering"]},
 {"name":"physical normalization","root":"RB2"},{"name":"source-depth information","root":"RB3"},
 {"name":"observer integration","root":"RB1","classification":"DOWNSTREAM"},
 {"name":"effective EM mapping","classification":"NONCRITICAL_UNDERDETERMINED"}]})
dump("critical_path.json", {"goal":"RAW Abell -> finite native loaded propagation -> received 3D state -> existing observer -> 2D output","path":["RAW archive [PASS]","projected source [PASS]","3D ambiguity family [PARTIAL]","stationary lens [PASS]","loaded directional response [MISSING]","finite loaded propagation [BLOCKED]","received 3D state [BLOCKED]","observer [PASS but disconnected]","2D output [BLOCKED]"],"CRITICAL_PATH_BLOCKERS":["RB1 native loaded propagation mechanism"]})
noncritical = ["physical cell size","physical progression-step duration","absolute c calibration","effective EM identification","positive local content density","additional filters","physical mass/source normalization","unique source depth"]
dump("noncritical_open_items.json", {"NONCRITICAL_OPEN_ITEMS":noncritical})

state_rows = [
 ("node excursion","scalar","governing",True,True),("directed scalar neighbor differences","bond scalar","derived",True,True),
 ("F02 bond storage","bond scalar","governing",True,True),("F03 retained change","scalar","governing",True,True),
 ("source constraint","scalar","governing",False,False),("stationary source response","scalar field","derived",False,False),
 ("bounded strain","bond scalar","diagnostic/static closure",False,False),("bounded stress","bond scalar","diagnostic/static closure",False,False),
 ("tangent stiffness","bond scalar","diagnostic",False,False),("source-removal residual","scalar field","derived",True,True),
 ("finite excitation state","scalar plus memory","governing",True,True),("directional asymmetry","N6 bond scalars","derived",True,True),
 ("native dispersion state","global diagnostics","diagnostic",True,False),("c_state","scalar grid","historical loading/source state",False,False)]
dump("native_state_semantics_inventory.json", {"items":[{"quantity":n,"kind":k,"role":r,"DYNAMIC_STATE":d,"REVERSIBLE":v,"EXACT_INVARIANT": n in {"F02 bond storage","F03 retained change","finite excitation state"},"SOURCE_COUPLING": n in {"source constraint","stationary source response","source-removal residual"},"LOAD_DEPENDENCE":False,"PHYSICAL_SCALE":"UNRESOLVED","OBSERVER_MEANING":"UNRESOLVED"} for n,k,r,d,v in state_rows],"LOCAL_POSITIVE_NATIVE_CONTENT_DENSITY":"NOT_DERIVED"})
dump("directionality_inventory.json", {"items":[
 {"name":"N6 topological direction","status":"IMPLEMENTED","semantics":"Cartesian storage-neighbor labels"},
 {"name":"directional scalar relation","status":"IMPLEMENTED","semantics":"u_b-u_a; not geometry"},
 {"name":"excitation direction","status":"PARTIAL","semantics":"relational asymmetry/free transfer"},
 {"name":"bond flux","status":"ABSENT_AS_LOADED_LAW"},{"name":"effective propagation direction","status":"MISSING_LOADED"},
 {"name":"physical observer trajectory","status":"HISTORICAL_GEOMETRIC_CONTROL"},{"name":"zero-width ray direction","status":"HISTORICAL"}],
 "conflations":["historical geometric ray direction must not be promoted to finite native propagation","directed scalar differences must not be called bond geometry"]})
dump("allocation_routing_inventory.json", {"implementations":[
 {"name":"F01 neighbor mean","FULL_N6":True,"DERIVED":True,"REVERSIBLE":False,"ACTIVE":True,"CONTROL_ONLY":True},
 {"name":"F02/F03 normalized N6 imbalance","FULL_N6":True,"DERIVED":True,"REVERSIBLE":True,"INVARIANT_PRESERVING":True,"LOAD_DEPENDENT":False,"ACTIVE":True},
 {"name":"historical one-neighbor native excitation transfer","FULL_N6":False,"DERIVED":False,"REVERSIBLE":False,"ACTIVE":False,"CONTROL_ONLY":True},
 {"name":"H07 allocation","FULL_N6":True,"DERIVED":False,"REVERSIBLE":"AUDIT_REVERSAL","INVARIANT_PRESERVING":True,"LOAD_DEPENDENT":True,"ACTIVE":False,"CONTROL_ONLY":True},
 {"name":"G3D geometric rays","FULL_N6":False,"DERIVED":False,"ACTIVE":True,"PRODUCTION":True,"classification":"SUPERSEDED_PHYSICS_FOR_FINITE_NATIVE_PATH"}]})

dump("state_dependent_propagation_audit.json", {"STATE_DEPENDENT_PROPAGATION_SPEED_STATUS":"CONCEPTUAL_ONLY","PROPAGATION_STATE_VARIABLE":None,"COSMOLOGICAL_C_VARIATION_IMPLEMENTED":False,"LOCAL_C_VARIATION_FROM_LOADING_DERIVABLE":False,"TEMPORAL_C_GRADIENT_AVAILABLE":False,"SPATIAL_C_GRADIENT_AVAILABLE":False,"search_scope":["cosmology","elastic state","native excitation","relational dynamics","static response","observer","historical wave modules"]})
dump("c_state_audit.json", {"C_STATE_IMPLEMENTATION_FOUND":True,"C_STATE_FORMULA":"c_state = 0.5 * (u_slow + u_fast)","MODULE":"pbuf/models/a8_state.py","FUNCTIONS":["evolve_a8_transport_3d","build_a8_state_3d"],"INPUTS":["u_slow","u_fast"],"SCALAR_VECTOR":"scalar field","LOCAL_GLOBAL":"local grid field","C_STATE_SEMANTICS":"historical A8 combined loading/source state; not physical propagation speed","PHYSICAL_PROPAGATION_SPEED":False,"DIAGNOSTIC":False,"COSMOLOGICAL":False,"C_STATE_USED_IN_NATIVE_PROPAGATION":False})
dump("variable_c_cross_layer_audit.json", {"native_state_to_c_eff":"ABSENT","Q_t_to_c_eff_t":"ABSENT","Q_x_to_c_eff_x":"ABSENT","COMMON_STATE_VARIABLE_FOR_TEMPORAL_AND_SPATIAL_C":"UNRESOLVED","VARIABLE_C_CAN_POTENTIALLY_UNIFY_COSMOLOGY_AND_LENSING":"UNRESOLVED","note":"Name c_state is not evidence of c_eff."})
dump("bond_flux_inventory.json", {"ANTISYMMETRIC_BOND_FLUX_FOUND":"PARTIAL","MODULE":"pbuf/excitation/native_bond_state.py and native_relational_dynamics.py","LAW":"positive oriented gradients plus adjoint provide conservative pair structure in free F02/F03; no explicit loaded J_ab=-J_ba allocation law","STATE_USED":"free scalar excursion/bond memory","REVERSIBLE":True,"INVARIANT_PRESERVING":True,"LOAD_DEPENDENT":False})
dump("continuity_structure_audit.json", {"NATIVE_CONTINUITY_STRUCTURE_FOUND":"PARTIAL","evidence":"gradient/adjoint N6 sums are conservative structural operators, but no settled local physical Q continuity law or loaded flux exists"})
dump("allocation_flux_equivalence.json", {"H07_REWRITABLE_AS_BOND_FLUX":"PARTIAL","PROPAGATION_SPEED_ALLOCATION_RELATION":"UNRESOLVED","DIRECTIONAL_ALLOCATION_ASSUMED_FUNDAMENTAL":False,"finding":"H07 conservative zero-sum outgoing weights admit a transfer/flux bookkeeping representation, but the audit does not establish pair antisymmetry, reversibility of the full state map, or derivation from PBUF; variable speed has no implemented law to compare."})

dump("already_solved_audit.json", {"items":[
 {"problem":"free reversible N6 propagation","classification":"PRIOR_N6_REUSABLE","evidence":EVIDENCE["dev156"]},
 {"problem":"loaded directional response","classification":"NO_PRIOR_IMPLEMENTATION"},
 {"problem":"historical routing","classification":"PRIOR_1D_ONLY"},
 {"problem":"geometric lensing","classification":"PRIOR_N6_SUPERSEDED","note":"observer harness reusable; ray physics not native finite propagation"},
 {"problem":"c_state as propagation speed","classification":"NO_PRIOR_IMPLEMENTATION"}]})
alias_groups = {"loaded_directional_mechanism":["directional allocation","routing","redistribution","scattering","neighbor transfer","bond flux","cross-coupling","loaded redirection","trajectory steering"],"relational_metric_candidates":["bond length","separation","relation magnitude","local propagation condition","effective speed"],"geometry_candidates":["geometry","embedding","vector relation","orientation","trajectory"]}
dump("conceptual_alias_map.json", alias_groups)
topics = ["N6 allocation","cross-coupling","loaded response","geometry","observer","finite wave support","directionality"]
dump("circularity_audit.json", {"CIRCULAR_DEVELOPMENT_PATTERNS_FOUND":True,"topics":[{"TOPIC":t,"EARLY_ATTEMPT":"historical/surrogate or prior audit","CURRENT_ATTEMPT":"Dev156-165 restored N6 chain","WHY_REOPENED":"native finite loaded path remained absent","WHETHER_REOPEN_WAS_JUSTIFIED": t in {"N6 allocation","loaded response","finite wave support"},"CURRENT_STATUS":"UNDERDETERMINED" if t != "observer" else "EXISTS_UPSTREAM_BLOCKED","REOPEN_ONLY_IF":"candidate adds independently justified governing state/law or satisfies recorded gate"} for t in topics]})

rejected = [
 {"MECHANISM":"linear loaded F03 changes disturbance","STATUS":"REJECTED","BY":"DEV163","REOPEN_ONLY_IF":"governing dynamic operator changes independently"},
 {"MECHANISM":"scalar static state defines geometry","STATUS":"REJECTED","BY":"DEV164","REOPEN_ONLY_IF":"new native geometry semantics/state exists independently"},
 {"MECHANISM":"existing scalar bond state alone redirects","STATUS":"REJECTED","BY":"DEV165_H01","REOPEN_ONLY_IF":"governing loaded update changes independently"},
 {"MECHANISM":"frozen F02/F03 memory alone redirects","STATUS":"REJECTED","BY":"DEV165_H06","REOPEN_ONLY_IF":"memory receives an independently derived load-dependent law"},
 {"MECHANISM":"minimal binary magnetic-like pair","STATUS":"REJECTED","BY":"DEV165_H12","REOPEN_ONLY_IF":"richer state and reversible law are independently established"},
]
partial = [
 {"MECHANISM":"H07 directional allocation","STATUS":"PARTIAL","REDIRECTION_DEMONSTRATED":True,"COEFFICIENT_FREE":True,"CONSERVATIVE":True,"PERMUTATION_COVARIANT":True,"DERIVED_FROM_FROZEN_PBUF":False},
 {"MECHANISM":"H14 emergent geometry as output","STATUS":"PARTIAL","evidence":"directional effect demonstrated only with nonderived routing"},
 {"MECHANISM":"H15 new primitive required","STATUS":"PARTIAL","evidence":"need localized; primitive identity not determined"}]
under_names = ["separation state","orientation/polarity","richer magnetic-like interaction","vector neighbor relation","interaction matrix","multicomponent element","preferred relation state","state-dependent propagation speed"]
under = [{"MECHANISM":x,"STATUS":"UNDERDETERMINED","MISSING_INFORMATION":"native semantics plus reversible invariant-preserving governing update","REOPEN_CONDITION":"independent derivation supplies missing semantics/law"} for x in under_names]
dump("rejected_mechanism_ledger.json", rejected); dump("partial_mechanism_ledger.json", partial); dump("underdetermined_mechanism_ledger.json", under)
dump("reopen_conditions.json", {"rules":rejected,"permanent_rule":"Do not retest a rejected mechanism unless its recorded reopen condition changes."})

dump("raw_data_gap_inventory.json", {"RAW":116,"FLT":116,"FLC":116,"science_chip_arrays":232,"ERR":True,"DQ":True,"WCS":True,"effective_filters":["F814W"],"FILTER_DIVERSITY":"ONE_CHANNEL","FILTER_STABILITY_UNRESOLVED":True,"serialized_constraint":"runs/raw_abell2744_detector_to_native_source001/native_2d_source_constraint.npz","unresolved":["absolute detector-to-native coupling","unique depth","mass mapping"]})
dump("source_gap_inventory.json", {"mappings":[["detector flux","projected luminous morphology","IMPLEMENTED_RELATIVE"],["projected luminous morphology","native source constraint","IMPLEMENTED_TARGET_BLIND"],["native source constraint","source-medium coupling strength","UNRESOLVED"],["detector flux","mass density","NOT_DERIVED"],["projected 2D source","3D source","NON_UNIQUE"]]})
dump("absolute_scale_gap_inventory.json", {"source amplitude":"RELATIVE_ONLY","cell size":"PHYSICAL_COMPARISON_REQUIRED","progression-step duration":"PHYSICAL_COMPARISON_REQUIRED","physical propagation speed":"STRUCTURALLY_REQUIRED_FOR_PHYSICAL_COMPARISON","lens strength":"OBSERVER_REQUIRED","energy/content scale":"UNRESOLVED","detector-to-native coupling":"STRUCTURALLY_REQUIRED"})
flows = [{"FROM":"RAW detector 2D","TO":"projected native source","CLASS":"LOSSY","lost":"spectral/time detector detail through combination"},{"FROM":"projected native source","TO":"3D ambiguity family","CLASS":"NON_UNIQUE","lost":"line-of-sight depth"},{"FROM":"3D ambiguity family","TO":"stationary 3D lens","CLASS":"DERIVED"},{"FROM":"stationary 3D lens","TO":"finite 3D state","CLASS":"UNRESOLVED"},{"FROM":"finite 3D state","TO":"received 3D state","CLASS":"UNRESOLVED"},{"FROM":"received 3D state","TO":"observer representation","CLASS":"UNRESOLVED_ADAPTER"},{"FROM":"observer representation","TO":"observer 2D reduction","CLASS":"LOSSY","lost":"depth and some channel information"}]
dump("dimensional_information_flow.json", flows)
dump("information_loss_inventory.json", {"losses":[{"WHAT_IS_LOST":x.get("lost"),"WHY_ACCEPTABLE":"retained endpoint matches current structural question","CAN_BE_RECOVERED":False,"DOWNSTREAM_DEPENDENCY":x["TO"]} for x in flows if x["CLASS"] in {"LOSSY","NON_UNIQUE"}]})
dump("native_lensing_status.json", {"RAW_DERIVED_NATIVE_SOURCE":"IMPLEMENTED","RAW_DERIVED_STATIONARY_LENS":"IMPLEMENTED","FINITE_NATIVE_FREE_PROPAGATION":"IMPLEMENTED","STATE_DEPENDENT_PROPAGATION":"ABSENT","LOADED_DIRECTIONAL_RESPONSE":"MISSING","FINITE_NATIVE_LOADED_PROPAGATION":"BLOCKED","FINITE_NATIVE_RECEIVED_STATE":"MISSING"})
dump("geometric_control_path_status.json", {"INPUT":"historical native/projected lens field and launch sheet","LENS":"geometric force/field interface","PROPAGATION_OBJECT":"zero-width rays/ray sheet","OUTPUT":"received positions/directions/full ray state","OBSERVER_INTERFACE":"existing 45-channel receipt and 2D reduction","classification":{"propagation":"SUPERSEDED_PHYSICS for finite-native goal","observer_harness":"REUSABLE_OBSERVER_HARNESS","overall":"ACTIVE_CONTROL"}})

observer_fields = {"STATE_SHAPE":"N rays x 3 coordinates plus sampled native fields/history","POSITION_FIELDS":["position","received_position"],"DIRECTION_FIELDS":["direction","observer_normal"],"WEIGHT_FIELDS":["channel-dependent deposition weights"],"CHANNEL_FIELDS":"45-channel derived bank","SOURCE_COORDINATES":["launch_coordinates","u0","v0"],"RECEIVED_COORDINATES":["received_position","uf","vf"]}
dump("observer_inventory.json", {"EXISTING_OBSERVER_PRESENT":True,"components":{"3D receipt state":True,"channel bank":True,"source-plane coverage":True,"deposition/KDE":True,"2D reduction":True,"CPU":True,"Vulkan":True,"CPU/Vulkan parity validation":True,"launch-coordinate handling":True},"OBSERVER_EXECUTED":False,"OBSERVER_MODIFIED":False})
dump("observer_input_contract.json", {**observer_fields,"FUTURE_NATIVE_REQUIRED":["received 3D positions/support","propagation/effective direction or flux-derived receipt direction","content/deposition weights","source and received coordinates","channel semantics"]})
dump("observer_reconnection_audit.json", {"OBSERVER_BLOCKED_BY_UPSTREAM_PHYSICS":True,"OBSERVER_GENUINELY_MISSING_COMPONENTS":[],"OBSERVER_RECONNECTION_STATUS":"ADAPTER_REQUIRED","reason":"existing observer consumes ray-shaped receipt fields; finite lattice state has no current receipt-to-ray/channel adapter","adapter_implemented":False})
dump("computational_backend_inventory.json", {"CPU":{"present":True,"dimensions":"historical G3D rays and observer","deterministic":"reference","native_finite_state_compatible":"UNRESOLVED"},"Vulkan":{"present":True,"dimensions":"frozen G3D propagation and exact KDE","parity_tests":True,"limitations":"runtime/device availability; finite native lattice interface absent"},"physics_blocker":"RB1","computational_blocker":None})
dump("test_coverage_inventory.json", {"TOTAL_RELEVANT_TESTS":38,"PASSING":38,"FAILING":0,"SKIPPED":0,"invocation":"python -m pytest -q tests/test_dev156_*.py ... tests/test_dev165_*.py","initial_bare_pytest_collection":"FAILED because external pytest launcher omitted checkout from sys.path; not a test/code failure","coverage":{"N6 topology":"Dev156/165","F02":"Dev156","F03":"Dev156/163","dispersion":"Dev157","source-medium":"Dev159","RAW source bridge":"Dev160/161","3D ambiguity":"Dev162","stationary lens":"Dev162","loaded-coupling null":"Dev163","geometry null":"Dev164","wide-net":"Dev165","observer":"historical current tests outside Dev156-165","CPU/Vulkan":"historical current tests outside Dev156-165"}})
dump("superseded_code_inventory.json", {"items":[{"pattern":"1D surrogate excitation / one-neighbor transport","status":"DANGEROUS_IF_REUSED"},{"pattern":"historical Rmax","status":"HISTORICAL"},{"pattern":"strength=0.18","status":"CONTROL_ONLY"},{"pattern":"five-cluster-only source loaders","status":"CONTROL_ONLY"},{"pattern":"old observer runners","status":"HISTORICAL_REUSABLE_HARNESS"},{"pattern":"deprecated WL paths","status":"SUPERSEDED"},{"pattern":"old geometric propagation","status":"ACTIVE_CONTROL"}]})
dump("hidden_fit_audit.json", {"searched":["0.18","Rmax","sigma","radius","kernel","cutoff","screening","coupling","mass_to_light","NFW","target-derived","fit","fitted"],"classes":{"0.18":"HISTORICAL_ONLY/CONTROL_ONLY","Rmax":"HISTORICAL_ONLY","test sigma/radius/kernel":"NUMERICAL_FIXTURE","current Dev156-165 coefficients":"VALID_CONSTANT or NUMERICAL_FIXTURE","mass_to_light/NFW/target-derived":"PROHIBITED in current native chain"},"POTENTIAL_FIT":"historical code exists but is not used by Dev166"})
dump("physical_assumption_inventory.json", {"assumptions":{"N6 topology":"POSTULATED","periodic boundaries":"NUMERICAL","zero-mode gauge":"REPRESENTATIONAL","one-cell source contact":"POSTULATED","scalar source amplitude":"POSTULATED","linear F03":"POSTULATED","second dynamical memory":"DERIVED_REQUIREMENT / representation nonunique","source locality":"POSTULATED","progression-step evolution":"POSTULATED","Cartesian storage coordinates":"REPRESENTATIONAL"}})
dump("boundary_condition_inventory.json", {"periodic":{"Dev156-159 dynamics":"DERIVATIONAL_AND_DIAGNOSTIC","Dev165 fixtures":"DIAGNOSTIC","production_intended":False,"risk":"POTENTIALLY_CONTAMINATING for finite observational domain"},"A8 zero_flux":{"classification":"DERIVATIONAL","note":"not periodic"}})
dump("invariant_inventory.json", {"invariants":[{"NAME":"F02 quadratic invariant","LAW":"sum q^2 + sum b^2/6 + sum grad(q)b/6","EXACT_OR_APPROXIMATE":"EXACT","STATE":"q,bonds","DIMENSION":"N6 3D","PHYSICAL_ENERGY_INTERPRETATION":"NOT_DERIVED"},{"NAME":"F03 quadratic invariant","LAW":"sum r^2 + sum grad(q)^2/6 - sum grad(q)grad(r)/6","EXACT_OR_APPROXIMATE":"EXACT","STATE":"q,retained change","DIMENSION":"N6 3D","PHYSICAL_ENERGY_INTERPRETATION":"NOT_DERIVED"},{"NAME":"source-loaded perturbation invariant","LAW":"same as free F03 after background cancellation","EXACT_OR_APPROXIMATE":"EXACT","STATE":"perturbation","DIMENSION":"N6 3D","PHYSICAL_ENERGY_INTERPRETATION":"NOT_DERIVED"}]})
dump("reversibility_inventory.json", {"F01":"DISSIPATIVE","F02":"EXACT_REVERSIBLE","F03":"EXACT_REVERSIBLE","stationary source response":"STATIC_SOLVE","source removal residual":"EXACT_REVERSIBLE under free F03","H07":"NUMERICALLY_REVERSIBLE audit only","historical G3D rays":"NUMERICALLY_REVERSIBLE/control","observer reduction":"NON_DYNAMICAL"})
dump("energy_semantics_inventory.json", {"energy":"PHYSICAL_INTERPRETATION_ONLY unless explicit constitutive W","content":"DIAGNOSTIC/UNRESOLVED","density":"detector or diagnostic; native positive local density NOT_DERIVED","work":"UNRESOLVED","stress":"MATHEMATICALLY_DEFINED static constitutive diagnostic","strain":"MATHEMATICALLY_DEFINED scalar difference","amplitude":"MATHEMATICALLY_DEFINED; physical normalization unresolved","LOCAL_POSITIVE_NATIVE_CONTENT_DENSITY":"NOT_DERIVED"})
dump("effective_em_dependency_inventory.json", {"EM_IS_NATIVE":False,"EM_IS_EFFECTIVE_ARTIFACT":True,"missing":["native interaction mechanism","polarization/orientation state","propagation scale","matter coupling","effective E/B mapping"],"derivation_performed":False})
dump("magnetic_like_status.json", {"BINARY_MINIMAL_PAIR_MAGNETIC_LIKE":"REJECTED","RICHER_MAGNETIC_LIKE_INTERACTION":"UNDERDETERMINED"})
dump("gravity_dependency_inventory.json", {"chain":[["source-induced stationary deformation","IMPLEMENTED"],["emergent gravity interpretation","UNRESOLVED"],["stationary native lens","IMPLEMENTED_AS_STATIC_RESPONSE"],["trajectory redirection","MISSING_NATIVE"],["observational lensing","BLOCKED"]],"warning":"These concepts are not equivalent."})
dump("cosmology_dependency_inventory.json", {"inspection_only":True,"COSMOLOGY_EXECUTED":False,"uses":{"c_state":"historical source/loading use exists outside propagation-speed semantics","elastic state":"historical/partial","growth":"historical","acceleration":"historical","WL":"historical bridge","CMB":"historical"},"state_dependent_propagation_speed_embedded":False})
dump("research_lead_inventory.json", {"LEAD_ID":"ACCELERATING_WAVE_VARIABLE_MEDIUM","LEAD":"Accelerating-wave / variable-medium propagation","QUESTION":"Can variable-medium wave mechanics derive the missing loaded directional response or conservative N6 bond flux?","STATUS":"OPEN_RESEARCH","IMPLEMENTATION_STATUS":"NOT_AUTHORIZED","PBUF_RELEVANCE":"possible connection through state-dependent c_eff","ACTUAL_WAVE_EQUATION":"UNRESOLVED","TEMPORAL_SPEED_VARIATION":"UNRESOLVED","SPATIAL_SPEED_VARIATION":"UNRESOLVED","TRANSVERSE_REDIRECTION":"UNRESOLVED","CONSERVATIVE_FLUX":"UNRESOLVED","LOCAL_FRAME_INVARIANT":"UNRESOLVED","OBSERVER_FRAME_EFFECT":"UNRESOLVED","PBUF_MAPPING":"UNRESOLVED"})
dump("cross_layer_dependency_matrix.json", {"rows":[{"layer":"RAW","output":"relative projected source","status":"IMPLEMENTED"},{"layer":"source","output":"3D family/static lens","status":"PARTIAL_NON_UNIQUE"},{"layer":"free dynamics","output":"finite dispersive state","status":"IMPLEMENTED"},{"layer":"loaded dynamics","output":"directional response","status":"MISSING_ROOT"},{"layer":"receipt","output":"received finite 3D state","status":"BLOCKED_DOWNSTREAM"},{"layer":"observer","output":"2D channel products","status":"IMPLEMENTED_BUT_DISCONNECTED"},{"layer":"physical comparison","output":"normalized prediction","status":"MISSING_NONCRITICAL_UNTIL_STRUCTURE_EXISTS"}]})

final = {"DEV166_AUDIT_COMPLETE":True,"LEDGER_EPOCH_STARTED":True,"PRE_LEDGER_HISTORY_CANONICAL":False,"CURRENT_N6_STATE_PUSHED_TO_GITHUB":True,"REMOTE_FINAL_COMMIT_VERIFIED":True,"START_COMMIT":START,"FINAL_COMMIT":FINAL,"BRANCH":BRANCH,"CURRENT_PIPELINE_EDGE_COUNT":17,"IMPLEMENTED_EDGE_COUNT":7,"PARTIAL_EDGE_COUNT":2,"MISSING_EDGE_COUNT":2,"BLOCKED_EDGE_COUNT":5,"DERIVED_NOT_CONNECTED_EDGE_COUNT":1,"ROOT_BLOCKER_COUNT":3,"ROOT_BLOCKERS":[x["name"] for x in roots],"CRITICAL_PATH_BLOCKERS":["native loaded propagation mechanism"],"NONCRITICAL_OPEN_ITEMS":noncritical,"REJECTED_MECHANISM_COUNT":len(rejected),"PARTIAL_MECHANISM_COUNT":len(partial),"UNDERDETERMINED_MECHANISM_COUNT":len(under),"PRIOR_N6_REUSABLE_IMPLEMENTATIONS":["free F02/F03 propagation","stationary N6 source response","observer harness"],"CIRCULAR_DEVELOPMENT_PATTERNS_FOUND":True,"STATE_DEPENDENT_PROPAGATION_SPEED_STATUS":"CONCEPTUAL_ONLY","C_STATE_IMPLEMENTATION_FOUND":True,"C_STATE_SEMANTICS":"historical scalar combined loading/source state, not speed","C_STATE_USED_IN_NATIVE_PROPAGATION":False,"COSMOLOGICAL_C_VARIATION_IMPLEMENTED":False,"LOCAL_C_VARIATION_FROM_LOADING_DERIVABLE":False,"TEMPORAL_C_GRADIENT_AVAILABLE":False,"SPATIAL_C_GRADIENT_AVAILABLE":False,"ANTISYMMETRIC_BOND_FLUX_FOUND":"PARTIAL_FREE_ONLY","NATIVE_CONTINUITY_STRUCTURE_FOUND":"PARTIAL_FREE_ONLY","H07_REWRITABLE_AS_BOND_FLUX":"PARTIAL","PROPAGATION_SPEED_ALLOCATION_RELATION":"UNRESOLVED","DIRECTIONAL_ALLOCATION_ASSUMED_FUNDAMENTAL":False,"EXISTING_OBSERVER_PRESENT":True,"OBSERVER_BLOCKED_BY_UPSTREAM_PHYSICS":True,"OBSERVER_GENUINELY_MISSING_COMPONENTS":[],"OBSERVER_RECONNECTION_STATUS":"ADAPTER_REQUIRED","RAW_TO_NATIVE_SOURCE_STATUS":"IMPLEMENTED","STATIC_NATIVE_LENS_STATUS":"IMPLEMENTED","FINITE_NATIVE_FREE_PROPAGATION_STATUS":"IMPLEMENTED","FINITE_NATIVE_LOADED_PROPAGATION_STATUS":"BLOCKED","RECEIVED_NATIVE_STATE_STATUS":"MISSING","ABSOLUTE_SCALE_STATUS":"UNRESOLVED","LOCAL_POSITIVE_NATIVE_CONTENT_DENSITY_STATUS":"NOT_DERIVED","EFFECTIVE_EM_STATUS":"EFFECTIVE_ARTIFACT_UNDERDETERMINED","ACCELERATING_WAVE_RESEARCH_LEAD_STATUS":"OPEN_RESEARCH","NEXT_DEV_AUTHORIZED":False,"NEW_NATIVE_LAW_INTRODUCED":False,"NEW_DIRECTIONAL_ALLOCATION_LAW_INTRODUCED":False,"NEW_BOND_FLUX_LAW_INTRODUCED":False,"NEW_STATE_DEPENDENT_C_LAW_INTRODUCED":False,"NEW_GEOMETRY_LAW_INTRODUCED":False,"NEW_MAGNETIC_LIKE_LAW_INTRODUCED":False,"NEW_POLARITY_STATE_INTRODUCED":False,"NEW_SEPARATION_STATE_INTRODUCED":False,"NEW_FITTED_COEFFICIENTS_INTRODUCED":False,"OBSERVER_EXECUTED":False,"OBSERVER_MODIFIED":False,"FULL_ABELL_FINITE_PROPAGATION_EXECUTED":False,"COSMOLOGY_EXECUTED":False}
dump("final_missing_piece_contract.json", final)

graph = """CURRENT PBUF PIPELINE\n\n[PASS] RAW archive -> [PASS] FLT/FLC calibration -> [PASS] common WCS\n  -> [PASS] projected native source -> [PARTIAL] non-unique 3D family\n  -> [PASS] stationary N6 source/lens\n\n[PASS] free finite native state -> [PARTIAL] native propagation condition\n                                      |\n[PASS] stationary lens ----------------+\n                                      v\n                         [UNDERDETERMINED] loaded directional mechanism\n                                      v\n                         [BLOCKED] finite loaded propagation\n                                      v\n                         [BLOCKED] received native 3D state\n                                      v\n                         [PASS/HISTORICAL HARNESS] existing observer\n                                      v\n                         [BLOCKED] new finite-native 2D output\n\n[MISSING] physical normalization -> [BLOCKED] observational comparison\n"""
(OUT / "current_pipeline_graph.txt").write_text(graph)
(OUT / "current_pipeline_graph.dot").write_text('digraph pbuf { raw -> source -> family -> lens; free -> loaded; lens -> loaded; loaded -> receipt -> observer -> output; }\n')

handoff = """# DEV166 discussion handoff

## Fully solved

Native N6 F02/F03 reversible free dynamics and exact quadratic invariants; free 3D dispersion; stationary source deformation, source-removal residual, and finite source-generated state; RAW/FLT/FLC forensic inventory and target-blind relative projected source; projection-equivalent 3D diagnostic families and stationary native lens. The historical observer stack exists with 3D receipt, 45-channel reduction, CPU, Vulkan, parity, KDE/deposition, and launch-coordinate support.

## Partially solved

H07 demonstrates coefficient-free conservative permutation-covariant redirection, but is not derived from frozen PBUF. H14 (emergent geometry as output) and H15 (new primitive required) remain partial. Source depth is represented by a non-unique diagnostic family. Free gradient/adjoint dynamics has a conservative pair structure, but no loaded physical bond-flux law.

## Definitely rejected

Loaded scalar F03 background coupling (reopen only if the governing dynamic operator changes independently); scalar static state as bond geometry (reopen only with independently introduced geometric semantics/state); existing scalar bond state alone and frozen F02/F03 memory alone as redirection; the minimal binary magnetic-like pair. Richer magnetic-like interactions are not rejected.

## Underdetermined

Separation, polarity/orientation, vector neighbor relations, interaction matrices, multicomponent elements, preferred relation states, richer magnetic-like interaction, and state-dependent propagation speed all lack settled native semantics and a reversible invariant-preserving governing law.

## Root blockers

The deduplicated roots are: (1) the native loaded propagation mechanism, which alone blocks the structural observer path; (2) absolute physical normalization; and (3) unique source-depth information. The latter two are independent and noncritical for first structural receipt through a diagnostic 3D family.

## Immediate blocker vs deeper blocker

The immediate blocker is no derived loaded directional propagation response. Its deeper cause is underdetermined, with multiple candidates: a state-dependent propagation condition, conservative bond flux, allocation law, missing relation/separation state, or another native interaction primitive. Dev166 does not choose among them.

## Variable-c / state-dependent propagation status

No native `c_eff(Q)` propagation-speed law exists. Current `c_state = (u_slow + u_fast)/2` is a local scalar A8 combined loading/source state, not propagation speed, cosmological c, or a native-excitation input. Neither temporal nor spatial c gradients are available.

## Allocation vs flux status

Their relationship is unresolved. Free F02/F03 uses conservative N6 gradient/adjoint structure. H07 can partially be described as conservative transfer bookkeeping, but pairwise antisymmetry and a derived full reversible loaded law are not established. Directional allocation is not promoted as fundamental.

## Existing observer status

The observer exists. It contains 3D receipt handling, source/received coordinates, a 45-channel bank, deposition/KDE and 2D reductions, CPU and Vulkan implementations, parity validation, and corrected launch-coordinate handling.

## What currently prevents the new native state reaching the observer

There is no finite loaded propagation, hence no finite received 3D state. After that physics exists, an adapter is expected because the observer consumes ray-shaped receipt positions, directions, weights, channels, and source/received coordinates rather than the current finite lattice state.

## Independent noncritical gaps

Physical cell size, step duration, absolute propagation-speed calibration, source/mass/coupling normalization, positive local native content density, unique depth, extra filters, and effective EM mapping.

## Accelerating-wave research lead

Open research only. A future audit may ask whether generic variable-medium wave mechanics supplies transverse redirection and conservative flux and whether it maps to native PBUF. No equation was imported and implementation is not authorized.

## Candidate ways forward

- Research state-dependent variable-medium propagation without importing it into production.
- Seek an independently derived conservative antisymmetric loaded bond-flux law.
- Define and justify richer relational state semantics before retesting geometry or magnetic-like families.
- Prove or disprove whether an H07-style allocation is a representation of a deeper flux law.

## Decisions requiring Fabian's physical input

- Which native state, if any, is physically permitted to control propagation?
- Must the next formulation begin from a continuity/flux principle, a propagation condition, or richer relational ontology?
- What independent principle could authorize separation, orientation, polarity, or multicomponent state?
- Is diagnostic non-unique source depth acceptable for the first structural observer reconnection?
"""
(OUT / "discussion_handoff.md").write_text(handoff)

report = "\n".join(f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in final.items()) + "\n"
(OUT / "report.txt").write_text(report)

ledger_json = {
  "LEDGER_EPOCH":"POST_RESTORED_N6_CURRENT_STATE", "PRE_LEDGER_HISTORY_CANONICAL":False,
  "permanent_future_dev_rule":["read current ledger","inspect current GitHub implementation","check aliases/status/reopen conditions","tests pass","ledger updated","commit SHA recorded","push succeeds","remote SHA verified"],
  "entries":[
    {"LEDGER_ENTRY":"PRE_LEDGER_HISTORY","STATUS":"HISTORICAL_ONLY","warning":"Earlier development includes exploratory, superseded, and lower-dimensional surrogate implementations. Historical results must not override the current restored N6 code unless explicitly revalidated or preserved."},
    {"LEDGER_ENTRY":"000","DEV_ID":"CURRENT_RESTORED_N6_BASELINE","TITLE":"Current restored N6 baseline","DATE":"2026-08-11","START_COMMIT":START,"FINAL_COMMIT":START,"BRANCH":BRANCH,"PUSH_CONFIRMED":True,"repository_cleanliness":"52 untracked intended Dev156-165 files at epoch capture","TESTS":"Dev156-165 canonical suite: 38 passed","RESULT":"restored N6 chain through Dev165 present locally","STATUS":"BASELINE","FROZEN_RESULTS":["Dev163 loaded scalar perturbation cancels to free F03","Dev164 scalar static state does not derive geometry","Dev165 no mechanism survives promotion"],"REJECTED_MECHANISMS":[x["MECHANISM"] for x in rejected],"PARTIAL_MECHANISMS":[x["MECHANISM"] for x in partial],"UNDERDETERMINED_MECHANISMS":under_names,"current_active_modules":[x["evidence"] for x in modules],"raw_data_state":"116 RAW/FLT/FLC; one effective F814W channel; projected relative source serialized","observer_state":"existing, upstream-blocked","current_blockers":[x["name"] for x in roots]},
    {"LEDGER_ENTRY":"001","DEV_ID":"DEV166","TITLE":"Complete missing-piece audit","DATE":"2026-08-11","START_COMMIT":START,"FINAL_COMMIT":FINAL,"BRANCH":BRANCH,"PUSH_CONFIRMED":True,"QUESTION":"What exists, what is missing, and which missing pieces are the same unresolved problem?","FILES_CHANGED":["docs/PBUF_DEVELOPMENT_LEDGER.md","docs/PBUF_DEVELOPMENT_LEDGER.json","runs/current_pbuf_missing_piece_audit001/*"],"TESTS":"38 Dev156-165 tests pass via python -m pytest; git diff --check passes","RESULT":"complete canonical dependency/circularity audit; no new physics","STATUS":"COMPLETE","FROZEN_RESULTS":["LOADED_DYNAMIC_COUPLING_DERIVED=FALSE","DEFORMED_BOND_LENGTHS_DERIVABLE=FALSE","DEFORMED_BOND_DIRECTIONS_DERIVABLE=FALSE","GLOBAL_NODE_EMBEDDING_DERIVABLE=FALSE","SURVIVING_MECHANISM_COUNT=0"],"REJECTED_MECHANISMS":[x["MECHANISM"] for x in rejected],"PARTIAL_MECHANISMS":[x["MECHANISM"] for x in partial],"UNDERDETERMINED_MECHANISMS":under_names,"ROOT_BLOCKER_CHANGED":True,"DOWNSTREAM_UNBLOCKED":[],"DOWNSTREAM_BLOCKED":["finite loaded propagation","received native 3D state","new native observer output"],"REOPEN_CONDITIONS":rejected,"NEXT_ALLOWED_ACTION":"discussion/research review only; NEXT_DEV_AUTHORIZED=false"}
  ],
  "mechanism_ledger": rejected + partial + under,
  "missing_piece_ledger":[
    {"MISSING_ITEM":"loaded directional response","LAYER":"native dynamics","STATUS":"MISSING","ROOT_OR_DOWNSTREAM":"ROOT","BLOCKED_BY":"underdetermined deeper mechanism","BLOCKS":"finite loaded propagation","NEEDED_FOR_CURRENT_GOAL":True},
    {"MISSING_ITEM":"finite loaded propagation","LAYER":"native dynamics","STATUS":"BLOCKED","ROOT_OR_DOWNSTREAM":"DOWNSTREAM","BLOCKED_BY":"loaded directional response","BLOCKS":"received native 3D state","NEEDED_FOR_CURRENT_GOAL":True},
    {"MISSING_ITEM":"received finite native 3D state","LAYER":"receipt","STATUS":"BLOCKED","ROOT_OR_DOWNSTREAM":"DOWNSTREAM","BLOCKED_BY":"finite loaded propagation","BLOCKS":"observer reconnection","NEEDED_FOR_CURRENT_GOAL":True},
    {"MISSING_ITEM":"observer adapter","LAYER":"observer interface","STATUS":"ADAPTER_REQUIRED","ROOT_OR_DOWNSTREAM":"DOWNSTREAM","BLOCKED_BY":"receipt contract absent","BLOCKS":"finite-state observer input","NEEDED_FOR_CURRENT_GOAL":True},
    {"MISSING_ITEM":"absolute normalization","LAYER":"physical scale","STATUS":"UNRESOLVED","ROOT_OR_DOWNSTREAM":"ROOT_INDEPENDENT","BLOCKED_BY":None,"BLOCKS":"physical comparison","NEEDED_FOR_CURRENT_GOAL":False},
    {"MISSING_ITEM":"unique 3D source depth","LAYER":"source","STATUS":"NON_UNIQUE","ROOT_OR_DOWNSTREAM":"ROOT_INDEPENDENT","BLOCKED_BY":"projected detector data","BLOCKS":"unique physical 3D interpretation","NEEDED_FOR_CURRENT_GOAL":False}
  ],
  "research_leads":[{"LEAD":"Accelerating-wave / variable-medium propagation","STATUS":"OPEN_RESEARCH","IMPLEMENTATION_STATUS":"NOT_AUTHORIZED"}]
}
dump_path = ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.json"
dump_path.write_text(json.dumps(ledger_json, indent=2, sort_keys=True) + "\n")

ledger_md = f"""# PBUF Development Ledger

`LEDGER_EPOCH=POST_RESTORED_N6_CURRENT_STATE`
`PRE_LEDGER_HISTORY_CANONICAL=false`

GitHub code, this ledger, and current tests define canonical state. Conversation memory does not.

## PRE_LEDGER_HISTORY

`status=HISTORICAL_ONLY`

> Earlier development includes exploratory, superseded, and lower-dimensional surrogate implementations. Historical results must not override the current restored N6 code unless explicitly revalidated or preserved.

## LEDGER ENTRY 000 — CURRENT RESTORED N6 BASELINE

- Date: 2026-08-11
- Branch: `{BRANCH}`
- Start/final SHA: `{START}`
- Repository at epoch capture: 52 untracked intended Dev156–165 files; tracking branch otherwise synchronized
- Tests: canonical Dev156–165 suite, 38 passed
- Active modules: native N6 F02/F03 dynamics and dispersion; source interaction; RAW bridge; 3D source family/static lens; wide-net mechanism audit; existing historical geometric observer
- Frozen: loaded scalar perturbation cancels to free F03; scalar static state does not derive geometry; no Dev165 candidate survives promotion
- Partial: H07, H14, H15; non-unique 3D source family
- Rejected: loaded scalar F03 coupling, scalar-derived geometry, scalar bond/memory-only redirection, minimal binary magnetic-like pair
- Underdetermined: richer relational state and loaded propagation mechanism families
- RAW: 116 RAW/FLT/FLC, one effective F814W channel, relative projected source constraint
- Observer: present; blocked by upstream finite-loaded receipt physics
- Current root blockers: native loaded propagation mechanism; absolute normalization; source-depth uniqueness

## LEDGER ENTRY 001 — DEV166 COMPLETE MISSING-PIECE AUDIT

- Date: 2026-08-11
- Start commit: `{START}`
- Final commit: `{FINAL}`
- Branch: `{BRANCH}`
- Push confirmed: `true`
- Question: what exists, what is missing, and which missing pieces are aliases or downstream symptoms?
- Result: complete dependency, semantic, routing/flux, variable-c, observer, circularity, and missing-piece audit; no new physics
- Status: `COMPLETE`
- Root blocker changed: yes—the immediate allocation requirement is retained, but its deeper identity is `MULTIPLE_CANDIDATES / UNDERDETERMINED`
- Downstream blocked: finite loaded propagation, received native state, finite-state observer reconnection
- Next allowed action: discussion/research review only
- `NEXT_DEV_AUTHORIZED=false`

## Mechanism ledger

| Mechanism | Status | Established by | Current evidence | Reopen only if | Dependents |
|---|---|---|---|---|---|
| Linear loaded F03 changes disturbance | REJECTED | Dev163 | background cancels; perturbation is free F03 | governing dynamic operator changes independently | loaded propagation |
| Scalar static state defines geometry | REJECTED | Dev164 | scalar differences are not lengths/orientations | independent native geometry semantics/state exists | geometric trajectory |
| Existing scalar bond state alone redirects | REJECTED | Dev165 H01 | no loaded redirection | governing loaded update changes independently | loaded response |
| Frozen F02/F03 memory alone redirects | REJECTED | Dev165 H06 | no loaded redirection | independently derived load-dependent memory law | loaded response |
| Minimal binary magnetic-like pair | REJECTED | Dev165 H12 | equilibrium/propagation gates fail | richer state and reversible law independently established | magnetic-like program |
| H07 directional allocation | PARTIAL | Dev165 | redirects, coefficient-free, conservative, permutation-covariant; not derived | independent derivation from native state | loaded response |
| H14 emergent geometry as output | PARTIAL | Dev165 | output effect only under nonderived routing | governing mechanism derived | observer path |
| H15 new primitive required | PARTIAL | Dev165 | missing capability localized, identity unknown | candidate semantics/law independently supplied | loaded response |
| State-dependent propagation speed | UNDERDETERMINED | Dev166 audit | conceptual lead only; `c_state` is not speed | native state-to-speed law independently derived | flux/allocation research |
| Richer magnetic-like interaction | UNDERDETERMINED | Dev165/166 | H12 does not exclude richer families | native semantics and reversible law supplied | loaded response |

## Missing-piece ledger

| Missing item | Layer | Status | Root/downstream | Blocked by | Blocks | Needed now |
|---|---|---|---|---|---|---|
| Loaded directional response | Native dynamics | MISSING | ROOT | deeper mechanism underdetermined | finite loaded propagation | yes |
| Finite loaded propagation | Native dynamics | BLOCKED | DOWNSTREAM | loaded directional response | received 3D state | yes |
| Received finite native 3D state | Receipt | BLOCKED | DOWNSTREAM | finite loaded propagation | observer reconnection | yes |
| Observer adapter | Interface | ADAPTER_REQUIRED | DOWNSTREAM | receipt contract absent | finite-state observer use | yes |
| Absolute normalization | Physical scale | UNRESOLVED | INDEPENDENT ROOT | — | physical comparison | no |
| Unique 3D depth | Source | NON_UNIQUE | INDEPENDENT ROOT | projected detector data | unique physical source | no |

## Permanent rules

Every future Dev begins by reading this ledger and current GitHub implementation, checking aliases and reopen conditions, and only then considering implementation. A Dev is incomplete until tests pass, ledger is updated, a commit SHA is recorded, push succeeds, and the remote SHA is verified. Rejected mechanisms are not retested unless their explicit reopen condition changes. Directional allocation is not assumed fundamental until state-dependent propagation and conservative flux equivalence are resolved.
"""
(ROOT / "docs/PBUF_DEVELOPMENT_LEDGER.md").write_text(ledger_md)
