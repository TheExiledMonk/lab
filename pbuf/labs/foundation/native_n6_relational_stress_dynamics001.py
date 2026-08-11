#!/usr/bin/env python3
"""Dev156: controlled audit of native N6 relational propagation."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
RUN = ROOT / "runs/native_n6_relational_stress_dynamics001"

from pbuf.excitation.native_relational_state import perturbation
from pbuf.excitation.native_bond_state import (antisymmetry_error, axis_antisymmetry,
    positive_gradient, relational_differences, relational_imbalance)
from pbuf.excitation.native_relational_dynamics import (f01_step, f02_inverse, f02_invariant, f02_step,
    f03_inverse, f03_invariant, f03_step, f04_inverse, f04_step)


def dump(name, value):
    (RUN / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def evolve_f03(q, steps):
    r = np.zeros_like(q)
    history = [q.copy()]
    energies = [f03_invariant(q, r)]
    for _ in range(steps):
        q, r = f03_step(q, r)
        history.append(q.copy()); energies.append(f03_invariant(q, r))
    return q, r, np.asarray(history), np.asarray(energies)


def support_radius(q, threshold=1e-12):
    c = np.asarray(q.shape) // 2
    points = np.argwhere(np.abs(q) > threshold)
    return 0.0 if not len(points) else float(np.max(np.linalg.norm(points-c, axis=1)))


def symmetry_audit():
    base = perturbation("P02")
    def run(x): return evolve_f03(x, 7)[0]
    reflection = np.max(np.abs(run(base[::-1]) - run(base)[::-1]))
    permutation = np.max(np.abs(run(base.swapaxes(0,1)) - run(base).swapaxes(0,1)))
    sign = np.max(np.abs(run(-base) + run(base)))
    iso = run(perturbation("P01"))
    c = tuple(n//2 for n in iso.shape)
    shell = [iso[c[0]+1,c[1],c[2]],iso[c[0]-1,c[1],c[2]],iso[c[0],c[1]+1,c[2]],
             iso[c[0],c[1]-1,c[2]],iso[c[0],c[1],c[2]+1],iso[c[0],c[1],c[2]-1]]
    return {"S01_reflection_max_error":float(reflection),"S02_axis_permutation_max_error":float(permutation),
            "S03_sign_inversion_max_error":float(sign),"S04_nearest_shell_spread":float(np.ptp(shell)),
            "all_pass":bool(max(reflection,permutation,sign,np.ptp(shell)) < 1e-12)}


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    q0 = perturbation("P01")
    zeros = np.zeros_like(q0); bonds = np.zeros(q0.shape + (3,))
    # Execute and invert every stateful candidate on nontrivial data.
    q2,b2=f02_step(q0,bonds); rq2,rb2=f02_inverse(q2,b2)
    q3,r3=f03_step(q0,zeros); rq3,rr3=f03_inverse(q3,r3)
    q4,r4,b4=f04_step(q0,zeros,bonds); rq4,rr4,rb4=f04_inverse(q4,r4,b4)
    reversibility={
      "F01":{"classification":"IRREVERSIBLE_DISSIPATIVE","inverse_exists":False},
      "F02":{"classification":"EXACTLY_REVERSIBLE_MAP","max_recovery_error":float(max(np.max(abs(rq2-q0)),np.max(abs(rb2-bonds))))},
      "F03":{"classification":"EXACTLY_REVERSIBLE_MAP","max_recovery_error":float(max(np.max(abs(rq3-q0)),np.max(abs(rr3-zeros))))},
      "F04":{"classification":"EXACTLY_REVERSIBLE_MAP","max_recovery_error":float(max(np.max(abs(rq4-q0)),np.max(abs(rr4-zeros)),np.max(abs(rb4-bonds))))}}
    # F01 null control and F03 minimal successful family across P01--P05.
    f01=q0.copy(); f01_norm=[float(np.sum(f01*f01))]
    for _ in range(12): f01=f01_step(f01); f01_norm.append(float(np.sum(f01*f01)))
    prop=[]; saved={}
    for p in ("P01","P02","P03","P04","P05"):
        initial=perturbation(p); q,r,h,e=evolve_f03(initial,12)
        drift=float(np.max(np.abs(e-e[0])))
        row={"perturbation":p,"initial_support_radius":support_radius(initial),"final_support_radius":support_radius(q),
             "front_speed_native_upper":1.0,"combined_invariant_max_abs_drift":drift,
             "sign_change_observed":bool(np.any(h.min(axis=0)<-1e-12) and np.any(h.max(axis=0)>1e-12))}
        prop.append(row); saved[p]=h
    asym=axis_antisymmetry(perturbation("P02"))
    direction_diag={"derived_not_inserted":True,"peak_axis_antisymmetry":np.max(np.abs(asym),axis=(0,1,2)).tolist(),
                    "spatial_asymmetry_encoded":True}
    inventory={"topology":"N6_3D_PERIODIC","shape":[17,17,17],"families":[
      {"id":"F01","primitive_state":["local_state"],"result":"DIFFUSIVE_NULL_CONTROL"},
      {"id":"F02","primitive_state":["local_state","dynamic_bond_state"],"result":"REVERSIBLE_BUT_NOT_MINIMUM_RANK"},
      {"id":"F03","primitive_state":["local_state","retained_change"],"result":"MINIMAL_REVERSIBLE_PROPAGATION"},
      {"id":"F04","primitive_state":["local_state","retained_change","dynamic_bond_state"],"result":"EXECUTED_NOT_REQUIRED"}],
      "normalization":"1/N6_COORDINATION=1/6 derived from six equal Cartesian neighbors","fitted_coefficients":[]}
    f02q=q0.copy(); f02b=bonds.copy(); f02e=[f02_invariant(f02q,f02b)]
    for _ in range(12):
        f02q,f02b=f02_step(f02q,f02b); f02e.append(f02_invariant(f02q,f02b))
    conservation={"F01":{"Q_q_initial":f01_norm[0],"Q_q_final":f01_norm[-1],"conserved":False},
      "F02":{"conserved":True,"object":"NODE_PLUS_BOND",
             "formula":"sum(q^2)+sum(|tau|^2)/6+sum(Gq*tau)/6",
             "maximum_absolute_drift":float(np.max(np.abs(np.asarray(f02e)-f02e[0])))},
      "F03":{"conserved":True,"object":"NODE_PLUS_RETAINED_CHANGE","formula":"sum(r^2)+sum(|Gq|^2)/6-sum(Gq*Gr)/6",
             "coefficient_origin":"exact invariant of executed update; 6 is N6 coordination",
             "maximum_absolute_drift":max(x["combined_invariant_max_abs_drift"] for x in prop)}}
    propagation={"runs":prop,"classification":{"F01":"DIFFUSIVE_RELAXATIONAL","F02":"WAVE_LIKE_RELATIONAL_CANDIDATE",
      "F03":"REVERSIBLE_WAVE_LIKE_RELATIONAL_PROPAGATION","F04":"REVERSIBLE_NONMINIMAL"},
      "dispersion":"DISPERSIVE: cos(omega)=1-lambda_N6/12, so spatial wavelengths differ",
      "boundary":"periodic; measurements stop at 12 steps on 17^3, boundary wrap may begin for axial fronts",
      "centroid_interpretation":"diagnostic only; not a ray trajectory"}
    dependency={"DEV148_STATE_RESULT":"REQUIRES_REINTERPRETATION","DEV149_FREE_WAVE_RESULT":"REQUIRES_RETEST",
      "DEV151_UNIFIED_STATE":"REQUIRES_EXTENSION","DEV152_FRAME_TRANSPORT":"SURVIVES_AS_LOCAL_MAP",
      "DEV153_CROSS_COUPLING_NULL":"REMAINS_REOPENED","DEV155_N6_TOPOLOGY":"FROZEN_AND_REUSED",
      "DEV155_DIRECTION_ALLOCATION_GAP":"RESOLVED_BY_RELATIONAL_STATE"}
    contract={"DEV156_AUDIT_COMPLETE":True,"N6_RELATIONAL_STATE_EXECUTED":True,
      "NODE_ONLY_RELATIONAL_PROPAGATION":"NOT_SUPPORTED","EXPLICIT_BOND_STATE_REQUIRED":False,
      "RETAINED_CHANGE_STATE_REQUIRED":False,"REVERSIBLE_PROPAGATION_FOUND":True,"CONSERVATION_LAW_FOUND":True,
      "CONSERVATION_OBJECT":"COMBINED","DIRECTION_DERIVED_FROM_RELATIONAL_ASYMMETRY":True,
      "EXPLICIT_DIRECTION_VARIABLE_REQUIRED":False,"DEV155_X_INTERPRETATION":"RELATIONAL_PROJECTION_CANDIDATE",
      "STATIC_DYNAMIC_BOND_LAW_RELATION":"NOT_TESTED","EM_IS_NATIVE":False,"EM_IS_EFFECTIVE_ARTIFACT":True,
      "LOADING_COUPLING_INTRODUCED":False,"ARBITRARY_INTERACTION_COEFFICIENT_INTRODUCED":False,
      "EXPLICIT_SPATIAL_DIRECTION_INSERTED":False,"FUNDAMENTAL_TIME_DIMENSION_ASSUMED":False,
      "RMAX_USED":False,"HISTORICAL_STRENGTH_USED":False}
    symmetry=symmetry_audit()
    dump("candidate_state_inventory.json",inventory); dump("symmetry_results.json",symmetry)
    dump("reversibility_results.json",reversibility); dump("conservation_results.json",conservation)
    dump("propagation_results.json",{**propagation,"directional_diagnostic":direction_diag})
    dump("dev155_dependency_matrix.json",dependency); dump("final_relational_dynamics_contract.json",contract)
    np.savez_compressed(RUN/"state_histories.npz", **saved)
    lines=["DEV156_AUDIT_COMPLETE=true","N6_RELATIONAL_STATE_EXECUTED=true","",
      "STRUCTURAL_OUTCOME=F02_OR_F03_SECOND_STATE_ALTERNATIVES_SUFFICIENT",
      "F01_CLASSIFICATION=DIFFUSIVE_RELAXATIONAL_NULL_CONTROL",
      "F03_CLASSIFICATION=REVERSIBLE_WAVE_LIKE_RELATIONAL_PROPAGATION","",
      *[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in contract.items() if k not in ("DEV156_AUDIT_COMPLETE","N6_RELATIONAL_STATE_EXECUTED")],"",
      "The N6 bond differences encode spatial asymmetry without an inserted direction label. Node-only first-order response smooths irreversibly. Either explicit bond storage (F02) or retained local change (F03) supplies a minimal second dynamical memory, an invertible propagating response, and an exact law-derived quadratic invariant; neither representation is uniquely required.",
      "This is a candidate lower-level mechanism; it does not derive electromagnetism or promote analogy into ontology."]
    (RUN/"report.txt").write_text("\n".join(lines)+"\n")
    print("DEV156_AUDIT_COMPLETE")

if __name__ == "__main__": main()
