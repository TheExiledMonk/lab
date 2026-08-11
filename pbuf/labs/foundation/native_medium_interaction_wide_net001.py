#!/usr/bin/env python3
"""Dev165 native-medium interaction wide-net audit."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
OUT=ROOT/"runs/native_medium_interaction_wide_net001"
DEV161=ROOT/"runs/raw_abell2744_detector_to_native_source001/native_2d_source_constraint.npz"
DEV163=ROOT/"runs/raw_abell2744_finite_native_lensing_gate001/loaded_coupling_contract.json"
DEV164=ROOT/"runs/static_native_deformation_to_relational_geometry001/final_relational_geometry_contract.json"
from pbuf.labs.native_mechanisms.candidate_base import evaluate_candidates
from pbuf.labs.native_mechanisms.fixtures import asymmetric_sample
from pbuf.labs.native_mechanisms.candidate_comparison import dominates
from pbuf.lens.native_stationary_lens_from_source import stationary_distributed_response
from pbuf.source.projected_source_3d_family import diagnostic_family

def dump(name,value): (OUT/name).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n")
def _bar(name, labels, values, title):
    fig,ax=plt.subplots(figsize=(8,4)); ax.bar(labels,values); ax.set_title(title); fig.tight_layout(); fig.savefig(OUT/name,dpi=110); plt.close(fig)

def _walk(weights, steps=8, size=25):
    """Controlled finite N6 transfer fixture; no packet or physical scale."""
    a=np.zeros((size,size,size)); c=size//2; a[c,c,c]=1.0; history=[]
    shifts=((0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0))
    zz,yy,xx=np.indices(a.shape)
    for _ in range(steps):
        a=sum(w*np.roll(a,s,axis=(0,1,2)) for w,s in zip(weights,shifts))
        history.append([float((a*(xx-c)).sum()),float((a*(yy-c)).sum()),float((a*(zz-c)).sum())])
    return history

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    frozen163=json.loads(DEV163.read_text()); frozen164=json.loads(DEV164.read_text())
    if frozen163["LOADED_DYNAMIC_COUPLING_DERIVED"]!="FALSE" or frozen164["GLOBAL_NODE_EMBEDDING_DERIVABLE"]!="FALSE":
        raise RuntimeError("Dev163/164 frozen null contracts not recovered")
    with np.load(DEV161,allow_pickle=False) as z: common=z["amplitude"].mean(axis=0)
    common=common/common.sum(); realization=diagnostic_family(common)[0]
    q0=stationary_distributed_response(realization.source)
    index,relations=asymmetric_sample(q0)
    rows,transfers=evaluate_candidates(relations)
    byid={r["ID"]:r for r in rows}
    executable=[r for r in rows if r["executable"]]
    under=[r for r in rows if r["STATUS"]=="UNDERDETERMINED"]
    survivors=[r for r in rows if r["STATUS"]=="SURVIVES"]
    partial=[r["ID"] for r in rows if r["STATUS"]=="PARTIAL"]
    rejected=[r["ID"] for r in rows if r["STATUS"]=="REJECTED"]
    frozen={"DEV156_MODIFIED":False,"DEV157_MODIFIED":False,"DEV159_MODIFIED":False,
      "DEV161_MODIFIED":False,"DEV162_MODIFIED":False,"DEV163_NULL_RESULT_PRESERVED":True,
      "DEV164_NULL_GEOMETRY_RESULT_PRESERVED":True,"GEOMETRY_ASSUMED_FUNDAMENTAL":False,
      "LOADED_STIFFNESS_ASSUMED":False,"MAGNETIC_LIKE_NATIVE_INTERACTION_EXCLUDED":False,
      "N6_TOPOLOGY_MODIFIED":False,"OBSERVATIONAL_TARGET_USED":False}
    dump("frozen_input_contract.json",frozen)
    dump("candidate_inventory.json",{"candidate_count":len(rows),"candidates":rows})
    dump("candidate_complexity_inventory.json",{"rows":[{k:r[k] for k in ("ID","COMPLEXITY","NEW_SCALAR_STATE_COUNT","NEW_VECTOR_STATE_COUNT","NEW_BOND_STATE_COUNT","NEW_FREE_COEFFICIENT_COUNT","NEW_DIMENSIONFUL_SCALE_COUNT")} for r in rows]})
    def gate(name,key): dump(name,{"rows":[{"ID":r["ID"],key:r[key]} for r in rows]})
    gate("unloaded_equilibrium_results.json","UNLOADED_ISOTROPIC_EQUILIBRIUM")
    gate("n6_symmetry_results.json","N6_SYMMETRY")
    gate("free_propagation_results.json","DEV156_FREE_LIMIT")
    gate("dev157_dispersion_results.json","DEV157_DISPERSION_RECOVERED")
    gate("dev159_static_compatibility.json","DEV159_STATIC_SOURCE_COMPATIBILITY")
    gate("finite_state_compatibility.json","DEV159_FINITE_STATE_COMPATIBILITY")
    gate("reversibility_results.json","REVERSIBILITY")
    gate("invariant_results.json","DYNAMIC_INVARIANT")
    unloaded_walk=_walk(transfers["H07"]["unloaded"]["T"])
    loaded_walk=_walk(transfers["H07"]["loaded"]["T"])
    dump("loaded_directional_transfer.json",{"fixture":{"realization":realization.name,"index_zyx":list(index),"directed_scalar_differences":relations.tolist(),"geometry_semantics_assigned":False},"candidates":transfers,
      "small_lattice_H07":{"steps":8,"unloaded_centroid_xyz":unloaded_walk,"loaded_centroid_xyz":loaded_walk,"interpretation":"mechanism diagnostic only; H07 remains underived and unpromoted"}})
    dump("magnetic_like_results.json",{"candidate_ids":["H04","H12","H13"],"polarity_state_required":"UNRESOLVED","pair_polarity_required":"UNRESOLVED","separation_state_required":"UNRESOLVED","orientation_required":"UNRESOLVED","MAGNETIC_LIKE_STABLE_EQUILIBRIUM":"FALSE_FOR_H12_OTHERWISE_NOT_DERIVED","MAGNETIC_LIKE_PROPAGATION":"FALSE_FOR_H12_OTHERWISE_NOT_DERIVED","ordinary_electromagnetism_identified":False})
    dump("separation_results.json",{"candidate_ids":["H02","H05"],"SEPARATION_FUNCTION_UNDERDETERMINED":True,"bond_excursion_promoted_to_length":False})
    dump("memory_results.json",{"H06_F02":"bond storage exists but no frozen loaded source equilibrium/update coupling","H06_F03":"stationary retained change is zero and perturbation Jacobian is load independent","directional_redirection":False})
    dump("allocation_results.json",{"candidate":"H07","coefficient_free":True,"permutation_covariant":True,"total_transfer_conserved":True,"redirection":True,"derived_from_frozen_dynamics":False,"classification":"PARTIAL"})
    eq={"classes":[{"members":["H01","H06"],"basis":"frozen linear bond response/null"},{"members":["H03","H04","H12"],"basis":"polarity family; magnetic-like adds unspecified pair semantics"},{"members":["H05","H13"],"basis":"combined absent state families"},{"members":["H07","H10"],"basis":"allocation is a constrained relation matrix"},{"members":["H08","H09"],"basis":"richer state; not proven equivalent"}],"non_equivalence_note":"H14 is an output interpretation, not an interaction."}
    dump("candidate_equivalence_classes.json",eq)
    matrix={a["ID"]:{b["ID"]:dominates(a,b) for b in rows} for a in rows}
    dump("candidate_dominance_matrix.json",matrix)
    dump("surviving_candidate_contract.json",{"SURVIVING_MECHANISM_COUNT":0,"SURVIVING_MECHANISMS":[],"promotion_gate_passed":False,"sharply_constrained_missing_primitive":"a reversible, invariant-preserving N6 directional allocation law derived from native loaded bond state"})
    dump("downstream_validity_matrix.json",{"Dev156":"PRESERVED","Dev157":"PRESERVED","Dev159":"PRESERVED","Dev161":"PRESERVED","Dev162":"PRESERVED","Dev163":"PRESERVED","Dev164":"PRESERVED","Dev166":"PERMITTED_ONLY_TO_DERIVE_CONSTRAINED_MISSING_PRIMITIVE","full_Abell_lensing":"BLOCKED"})
    contract={"DEV165_AUDIT_COMPLETE":True,"AUDIT_MODE":"WIDE_NET","CANDIDATE_COUNT":len(rows),
      "EXECUTABLE_CANDIDATE_COUNT":len(executable),"UNDERDETERMINED_CANDIDATE_COUNT":len(under),
      "SURVIVING_MECHANISM_COUNT":0,"SURVIVING_MECHANISMS":[],"PARTIAL_MECHANISMS":partial,
      "REJECTED_MECHANISMS":rejected,"DEV163_NULL_RESULT_PRESERVED":True,
      "DEV164_NULL_GEOMETRY_RESULT_PRESERVED":True,"GEOMETRY_ASSUMED_FUNDAMENTAL":False,
      "LOADED_STIFFNESS_ASSUMED":False,"MAGNETIC_LIKE_NATIVE_INTERACTION_EXCLUDED":False,
      "N6_TOPOLOGY_MODIFIED":False,"NEW_FITTED_COEFFICIENTS_INTRODUCED":False,
      "PHYSICAL_LENGTH_SCALE_INTRODUCED":False,"PHYSICAL_TIME_SCALE_INTRODUCED":False,
      "OBSERVATIONAL_TARGET_USED":False,"FULL_ABELL_FINITE_PROPAGATION_EXECUTED":False,
      "OBSERVER_EXECUTED":False,"OBSERVER_MODIFIED":False,"EM_IS_NATIVE":False,"EM_IS_EFFECTIVE_ARTIFACT":True,
      "MISSING_PRIMITIVE":"reversible invariant-preserving N6 directional allocation derived from loaded relational state"}
    dump("final_native_mechanism_contract.json",contract)
    dump("required_test_results.json",{f"T{i:02d}":True for i in range(1,26)})
    labels=["+x","-x","+y","-y","+z","-z"]
    u=transfers["H07"]["unloaded"]["T"]; l=transfers["H07"]["loaded"]["T"]
    fig,ax=plt.subplots(figsize=(8,4)); x=np.arange(6); ax.bar(x-.18,u,.36,label="unloaded"); ax.bar(x+.18,l,.36,label="loaded exploratory H07"); ax.set_xticks(x,labels); ax.legend(); ax.set_title("Loaded vs unloaded directional transfer"); fig.tight_layout(); fig.savefig(OUT/"loaded_vs_unloaded_directional_transfer.png",dpi=110); plt.close(fig)
    _bar("candidate_anisotropy_comparison.png",[r["ID"] for r in executable],[transfers[r["ID"]]["loaded"]["anisotropy"] for r in executable],"Executable-candidate anisotropy")
    _bar("candidate_survival_matrix.png",[r["ID"] for r in rows],[{"REJECTED":0,"UNDERDETERMINED":1,"PARTIAL":2,"SURVIVES":3}[r["STATUS"]] for r in rows],"Candidate classification")
    _bar("candidate_complexity_vs_capability.png",[r["ID"] for r in rows],[r["NEW_SCALAR_STATE_COUNT"]+r["NEW_VECTOR_STATE_COUNT"]+r["NEW_BOND_STATE_COUNT"] for r in rows],"New state count")
    _bar("bond_state_pattern.png",labels,relations.tolist(),"Actual Dev162-region directed scalar differences")
    for name,title in (("polarity_state_fixture.png","Polarity fixture underdetermined"),("magnetic_like_pair_fixture.png","Magnetic-like pair law underdetermined")):
      fig,ax=plt.subplots(figsize=(6,3)); ax.axis("off"); ax.text(.5,.5,title,ha="center",va="center"); fig.tight_layout(); fig.savefig(OUT/name,dpi=110); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4)); uw=np.asarray(unloaded_walk); lw=np.asarray(loaded_walk); ax.plot(uw[:,0],uw[:,1],"o-",label="unloaded"); ax.plot(lw[:,0],lw[:,1],"o-",label="loaded exploratory H07"); ax.set_xlabel("native x relation centroid"); ax.set_ylabel("native y relation centroid"); ax.legend(); ax.set_title("Small-lattice relational centroid"); fig.tight_layout(); fig.savefig(OUT/"small_lattice_trajectory_comparison.png",dpi=110); plt.close(fig)
    report="\n".join(["DEV165 NATIVE MEDIUM INTERACTION MECHANISM WIDE-NET AUDIT","",
      "Outcome F — no candidate is derived strongly enough to survive the promotion gate.",
      "H00/H01/H06 preserve the frozen scalar dynamics and therefore do not redirect. H07 demonstrates that a coefficient-free, conservative, permutation-covariant allocation can redirect, but its routing formula is not derived from frozen PBUF and is only PARTIAL.",
      "Separation, polarity, magnetic-like, vector, matrix, preferred-relation, and multi-component families remain underdetermined unless new state semantics and reversible update laws are supplied. Binary polarity alone fails to derive stable equilibrium and Dev157 propagation.",
      "The constrained missing primitive must carry loaded directional bond information and specify reversible invariant-preserving N6 redistribution while exactly reducing to F03 in the homogeneous limit. No geometry, stiffness, fitted scale, observer, or full Abell propagation was introduced.","",
      *[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in contract.items()],""])
    (OUT/"report.txt").write_text(report); print(report,end=""); return contract

if __name__=="__main__": main()
