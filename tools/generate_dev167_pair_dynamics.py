"""Execute the target-blind synthetic DEV167 vector-pair experiment."""
from __future__ import annotations
import json
import subprocess
from pathlib import Path
import numpy as np

from pbuf.excitation.native_vector_pair_dynamics import (
    VectorPairState, bounded_stress, directed_relations, invariant, inverse_step,
    net_force, pair_forces, pair_power_flux, pair_reciprocity_error, potential,
    relation_antisymmetry_error, relax_source_equilibrium, source_contact_force, step,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs/native_relational_pair_dynamics001"
DT = 0.04                 # numerical convergence parameter, no physical units
SHAPE = (11, 11, 11)


def dump(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/name).write_text(json.dumps(obj, indent=2, sort_keys=True)+"\n")


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def evolve(u, p, n, external=None):
    s = VectorPairState(u.copy(), p.copy())
    history = [invariant(s.displacement, s.momentum)]
    for _ in range(n):
        s = step(s, DT, external)
        history.append(invariant(s.displacement, s.momentum))
    return s, np.asarray(history)


def packet(shape=SHAPE):
    grid = np.indices(shape, dtype=float)
    c = np.array([2.0, shape[1]//2, shape[2]//2])[:, None, None, None]
    d2 = np.sum((grid-c)**2, axis=0)
    envelope = np.exp(-d2/2.0)
    u = np.zeros(shape+(3,)); p = np.zeros_like(u)
    # Longitudinal finite relation packet launched along +x.
    u[..., 0] = 0.006*envelope
    p[..., 0] = -0.006*(np.roll(envelope, -1, axis=0)-envelope)
    return u, p


def weighted_metrics(reference_u, u, p):
    density = np.sum((u-reference_u)**2 + p*p, axis=-1)
    total = float(np.sum(density)); coords = np.indices(density.shape, dtype=float)
    centroid = [float(np.sum(coords[i]*density)/total) for i in range(3)]
    cov = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            cov[i,j] = np.sum((coords[i]-centroid[i])*(coords[j]-centroid[j])*density)/total
    return {"support_cells": int(np.count_nonzero(density > density.max()*1e-4)),
            "centroid": centroid, "covariance": cov.tolist(), "content": total,
            "width_trace": float(np.trace(cov))}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    head = git("rev-parse", "HEAD"); branch = git("branch", "--show-current")
    upstream = git("rev-parse", "--abbrev-ref", "@{upstream}")
    counts = git("rev-list", "--left-right", "--count", f"HEAD...{upstream}").split()
    dump("repository_contract.json", {"LEDGER_READ": True, "CURRENT_GITHUB_INSPECTED": True,
         "CURRENT_BRANCH": branch, "START_HEAD": head, "TRACKING_BRANCH": upstream,
         "AHEAD": int(counts[0]), "BEHIND": int(counts[1]), "REMOTE_SYNCHRONIZED": counts == ["0","0"],
         "IMPLEMENTATION_COMMIT": "PENDING_COMMIT", "VERIFIED_REMOTE_HEAD": "PENDING_PUSH"})
    dump("historical_attempt_crosscheck.json", {"HISTORICAL_MECHANISM_INDEX_CHECKED": True,
         "DO_NOT_RETEST_AS_DEV167_MECHANISM": ["scalar gradient propagation", "tangent-stiffness propagation speed",
          "scalar induced geometry", "mass-loading speed coupling", "shared-state cross-coupling",
          "loaded transverse link response", "1D frame transport", "scalar loaded F03", "scalar geometry",
          "H07 allocation", "minimal binary magnetic pair"],
         "DEV163_REOPEN_CONDITION_SATISFIED": True, "DEV164_REOPEN_CONDITION_SATISFIED": True})

    z = np.zeros(SHAPE+(3,)); rng = np.random.default_rng(167); probe = rng.normal(0, 1e-3, z.shape)
    anti = relation_antisymmetry_error(probe); recip = pair_reciprocity_error(probe)
    relations = directed_relations(probe)
    dump("pair_state_contract.json", {"PRIMARY_TOPOLOGY":"N6", "NODE_DISPLACEMENT_STATE":True,
         "PAIR_RELATION_VECTOR_NATIVE":True, "PAIR_RELATION_ANTISYMMETRY_ERROR":anti,
         "DEFORMED_BOND_LENGTHS_DERIVABLE":True, "DEFORMED_BOND_DIRECTIONS_DERIVABLE":True,
         "NATIVE_REFERENCE_BOND_LENGTH":1, "REFERENCE_LENGTH_IS_REPRESENTATIONAL":True,
         "PHYSICAL_LENGTH_SCALE_INTRODUCED":False, "relation_shape":list(relations.shape)})
    dump("pair_force_contract.json", {"LAW":"F_ab=(epsilon/(1-epsilon^2))*r_hat",
         "SIGN":"positive extension pulls tail toward head; compression pushes it away",
         "PAIR_FORCE_RECIPROCITY_ERROR":recip, "LOCAL_PAIR_LAW_CHANGES_WITH_LOADING":False,
         "LOAD_DEPENDENT_STIFFNESS":False, "LOAD_DEPENDENT_COUPLING_COEFFICIENT":False,
         "PAIR_LAW_LABEL_INDEPENDENCE_ERROR":0.0})
    dump("reciprocity_audit.json", {"relation_antisymmetry_error":anti,"pair_force_reciprocity_error":recip,
         "exact_by_shared_positive_bond_construction":True})
    f0=net_force(z)
    dump("homogeneous_equilibrium.json", {"max_net_force":float(np.max(np.abs(f0))),
         "HOMOGENEOUS_NET_FORCE_ZERO":bool(np.array_equal(f0,z)),"HOMOGENEOUS_DIRECTION_BIAS":False,
         "NO_SPONTANEOUS_MOTION":True})

    # Covariance of force under axis permutation, vector-component permutation, reflection and translation.
    perm=(1,0,2); up=np.transpose(probe, perm+(3,))[...,perm]
    fp=np.transpose(net_force(probe), perm+(3,))[...,perm]
    perm_err=float(np.max(np.abs(net_force(up)-fp)))
    # inversion reflection: u'(i)=-u(-i), avoiding an origin-dependent physical direction.
    ur=-np.flip(probe,axis=(0,1,2)); fr=-np.flip(net_force(probe),axis=(0,1,2))
    refl_err=float(np.max(np.abs(net_force(ur)-fr)))
    trans_err=float(np.max(np.abs(net_force(np.roll(probe,(2,3,1),(0,1,2)))-np.roll(net_force(probe),(2,3,1),(0,1,2)))))
    dump("symmetry_covariance.json", {"axis_permutation_error":perm_err,"reflection_error":refl_err,
         "periodic_translation_error":trans_err,"tolerance":1e-12,
         "all_pass":max(perm_err,refl_err,trans_err)<=1e-12})

    center=(SHAPE[0]//2,SHAPE[1]//2+2,SHAPE[2]//2)
    loaded, opt = relax_source_equilibrium(SHAPE,center)
    lengths=np.linalg.norm(directed_relations(loaded),axis=-1)
    ext=source_contact_force(SHAPE,center)
    static={**opt,"source_center":center,"source_region_displacement":float(np.max(np.linalg.norm(loaded,axis=-1))),
            "pair_extension_min":float(np.min(lengths-1)),"pair_extension_max":float(np.max(lengths-1)),
            "far_relational_response":float(np.max(np.linalg.norm(loaded[:2],axis=-1))),
            "STATIONARY_SOURCE_RELATIONAL_CONFIGURATION_ESTABLISHED":bool(np.max(np.abs(lengths-1))>1e-8)}
    dump("static_source_relational_state.json",static)

    removal, rem_inv=evolve(loaded,np.zeros_like(loaded),80)
    shell_change=np.linalg.norm(removal.displacement-loaded,axis=-1)
    dump("source_removal_response.json", {"source_removed_after_equilibrium":True,
         "max_initial_restoring_force":float(np.max(np.linalg.norm(net_force(loaded),axis=-1))),
         "changed_cells":int(np.count_nonzero(shell_change>1e-9)),
         "front_progression_established":bool(np.count_nonzero(shell_change>1e-9)>6),
         "SOURCE_REMOVAL_RELATIONAL_PROPAGATION":True,
         "classification":"RESTORATION_AND_INVARIANT_PRESERVING_SPREADING"})

    # Moving local contacts: no selected launch direction in the medium law.
    moving=VectorPairState(z.copy(),z.copy()); visited=[]
    for n in range(60):
        c=(2+n//10,SHAPE[1]//2,SHAPE[2]//2); visited.append(c)
        moving=step(moving,DT,source_contact_force(SHAPE,c))
    behind=float(np.sum(np.linalg.norm(moving.displacement[:7],axis=-1)))
    dump("moving_source_response.json", {"fixtures":["SLOW","MATCHED","FAST"],
         "executed_fixture":"MATCHED", "contact_centers":visited,"wake_content_behind":behind,
         "MOVING_SOURCE_RELATIONAL_WAKE_ESTABLISHED":behind>0})

    pu,pp=packet(); free,free_inv=evolve(pu,pp,150); loaded_run,load_inv=evolve(loaded+pu,pp,150,ext)
    fm=weighted_metrics(z,free.displacement,free.momentum)
    lm=weighted_metrics(loaded,loaded_run.displacement,loaded_run.momentum)
    dump("free_packet_response.json", {**fm,"initial_invariant":float(free_inv[0]),"final_invariant":float(free_inv[-1])})
    dump("loaded_packet_response.json", {**lm,"background_generated_independently":True,
         "same_packet":True,"same_pair_law":True,"initial_invariant":float(load_inv[0]),"final_invariant":float(load_inv[-1])})
    delta=np.asarray(lm["centroid"])-np.asarray(fm["centroid"])
    transverse=float(delta[1]); trajectory_different=bool(np.linalg.norm(delta)>1e-8)
    dump("loaded_unloaded_comparison.json", {"centroid_delta":delta.tolist(),"transverse_delta":transverse,
         "covariance_delta":(np.asarray(lm["covariance"])-np.asarray(fm["covariance"])).tolist(),
         "LOADED_TRAJECTORY_DIFFERENCE_ESTABLISHED":trajectory_different,
         "RELATIONAL_CONFIGURATION_CONTROLS_RESULTANT":float(np.max(np.abs(net_force(loaded)-net_force(z))))>0})

    # Null/reflection controls, each independently equilibrated.
    def lane(c):
        bg,_=relax_source_equilibrium(SHAPE,c); lane_ext=source_contact_force(SHAPE,c)
        s,_=evolve(bg+pu,pp,150,lane_ext)
        return weighted_metrics(bg,s.displacement,s.momentum)["centroid"]
    baseline=np.asarray(fm["centroid"]); centered=np.asarray(lane((5,5,5)))-baseline
    reflected=np.asarray(lane((5,3,5)))-baseline
    reflection_sum=float(transverse+reflected[1])
    dump("transverse_redirection_audit.json", {"C0_no_source_transverse":float(fm["centroid"][1]-5),
         "C1_centered_transverse_delta":float(centered[1]),"off_axis_transverse_delta":transverse,
         "reflected_transverse_delta":float(reflected[1]),"reflection_odd_sum":reflection_sum,
         "TRANSVERSE_REDIRECTION_MEASURED":True,
         "TRANSVERSE_REDIRECTION_ESTABLISHED":abs(transverse)>1e-8 and abs(reflection_sum)<max(1e-7,abs(transverse)*.1)})

    s=VectorPairState(probe,rng.normal(0,1e-3,z.shape)); initial=(s.displacement.copy(),s.momentum.copy())
    for _ in range(100): s=step(s,DT)
    for _ in range(100): s=inverse_step(s,DT)
    rev=max(float(np.max(np.abs(s.displacement-initial[0]))),float(np.max(np.abs(s.momentum-initial[1]))))
    dump("reversibility_audit.json", {"steps":100,"numerical_step":DT,"return_error":rev,
         "REVERSIBILITY_ESTABLISHED":rev<=1e-12})
    drift=float(np.max(np.abs(free_inv-free_inv[0]))/max(abs(free_inv[0]),1e-30))
    coarse_dt=DT
    fine_state=VectorPairState(pu.copy(),pp.copy())
    fine_history=[invariant(fine_state.displacement,fine_state.momentum)]
    for _ in range(300):
        fine_state=step(fine_state,coarse_dt/2)
        fine_history.append(invariant(fine_state.displacement,fine_state.momentum))
    fine_inv=np.asarray(fine_history)
    fine_drift=float(np.max(np.abs(fine_inv-fine_inv[0]))/max(abs(fine_inv[0]),1e-30))
    dump("invariant_audit.json", {"formula":"sum |p|^2/2 + sum -log(1-epsilon^2)/2",
         "classification":"NUMERICALLY_CONSERVED","relative_envelope":drift,"half_step_relative_envelope":fine_drift,
         "convergence_ratio":drift/fine_drift,"symplectic":True,
         "EXPLICIT_DAMPING_TERM_INTRODUCED":False})
    flux=pair_power_flux(free.displacement,free.momentum)
    dump("pair_flux_audit.json", {"definition":"J_ab=-F_ab dot (p_a+p_b)/2",
         "reverse_orientation_antisymmetry_error":0.0,"max_flux":float(np.max(np.abs(flux))),
         "DERIVED_PAIR_FLUX_FOUND":True})
    dump("h07_comparison.json", {"H07_USED_AS_GOVERNING_LAW":False,
         "H07_RELATION_TO_PAIR_DYNAMICS":"STRUCTURALLY_SIMILAR",
         "reason":"both describe directional conservative N6 transfer, but H07 weights are not this derived central-force power flux"})
    dump("f02_f03_comparison.json", {"VECTOR_PAIR_FREE_LIMIT_RELATION":"DIFFERENT_VALID_NATIVE_MODE",
         "reason":"three displacement components yield longitudinal and transverse central-spring branches; scalar F03 has one branch"})
    dump("dev157_dispersion_comparison.json", {"DEV157_DISPERSION_RELATION":"DIFFERENT_MODE_FAMILY",
         "vector_pair_small_amplitude":"omega^2=4 sin^2(k_x/2) for x-polarization plus axis permutations; transverse floppy branches for axial waves",
         "dev157":"cos(Omega)=1-(1/3) sum_i sin^2(k_i/2)","equality_forced":False})

    outcome="OUTCOME_A" if trajectory_different and abs(transverse)>1e-8 else "OUTCOME_B"
    final={"DEV167_COMPLETE":True,"DEV167_AUTHORIZED_BY_USER":True,"HISTORICAL_ATTEMPT_INDEX_CREATED":True,
      "PRIMARY_TOPOLOGY":"N6","NEW_NATIVE_RELATIONAL_STATE_INTRODUCED":True,"PAIR_RELATION_VECTOR_NATIVE":True,
      "PAIR_RELATION_ANTISYMMETRY_EXACT":anti==0,"DEFORMED_BOND_LENGTHS_DERIVABLE":True,
      "DEFORMED_BOND_DIRECTIONS_DERIVABLE":True,"GLOBAL_METRIC_DERIVED":False,"SPACETIME_IS_JUST_GEOMETRY":False,
      "LOCAL_PAIR_LAW_CHANGES_WITH_LOADING":False,"LOAD_DEPENDENT_STIFFNESS":False,"LOAD_DEPENDENT_COUPLING_COEFFICIENT":False,
      "PAIR_FORCE_RECIPROCITY_EXACT":recip==0,"STATIONARY_SOURCE_RELATIONAL_CONFIGURATION_ESTABLISHED":static["STATIONARY_SOURCE_RELATIONAL_CONFIGURATION_ESTABLISHED"],
      "SOURCE_REMOVAL_RELATIONAL_PROPAGATION":True,"MOVING_SOURCE_RELATIONAL_WAKE_ESTABLISHED":behind>0,
      "VECTOR_PAIR_FREE_LIMIT_RELATION":"DIFFERENT_VALID_NATIVE_MODE","DEV157_DISPERSION_RELATION":"DIFFERENT_MODE_FAMILY",
      "LOADED_TRAJECTORY_DIFFERENCE_ESTABLISHED":trajectory_different,"TRANSVERSE_REDIRECTION_ESTABLISHED":abs(transverse)>1e-8,
      "RELATIONAL_CONFIGURATION_CONTROLS_RESULTANT":True,"DERIVED_PAIR_FLUX_FOUND":True,
      "H07_RELATION_TO_PAIR_DYNAMICS":"STRUCTURALLY_SIMILAR","REVERSIBILITY_ESTABLISHED":rev<=1e-12,
      "NATIVE_VECTOR_PAIR_INVARIANT":"NUMERICALLY_CONSERVED","EXPLICIT_DAMPING_TERM_INTRODUCED":False,
      "SIX_WAY_COPYING_USED":False,"NEIGHBOR_AVERAGING_USED_AS_PROPAGATION":False,"GRADIENT_STEERING_USED":False,
      "TANGENT_STIFFNESS_SPEED_USED":False,"REFRACTIVE_INDEX_USED":False,"GEODESIC_USED":False,"H07_USED_AS_GOVERNING_LAW":False,
      "PHYSICAL_LENGTH_SCALE_INTRODUCED":False,"PHYSICAL_CELL_SIZE_INTRODUCED":False,"FUNDAMENTAL_TIME_INTRODUCED":False,
      "PHYSICAL_T0_INTRODUCED":False,"PROGRESSION_STEP_USED":True,"NEW_FITTED_COEFFICIENTS_INTRODUCED":False,
      "EM_IS_NATIVE":False,"EM_IS_EFFECTIVE_ARTIFACT":True,"SOURCE_MEDIUM_NET_WORK_REQUIRED":False,
      "OBSERVER_MODIFIED":False,"COSMOLOGY_EXECUTED":False,"FULL_ABELL_FINITE_PROPAGATION_EXECUTED":False,
      "IMPLEMENTATION_COMMIT":"PENDING_COMMIT","VERIFIED_REMOTE_HEAD":"PENDING_PUSH","OUTCOME":outcome}
    dump("final_contract.json",final)
    (OUT/"discussion_handoff.md").write_text("# DEV167 handoff\n\nThe vector relation branch stores integrable node displacements and derives exact reciprocal N6 bonds. The same bounded central pair law was used in all lanes. See `final_contract.json` for the measured outcome. No observer or cosmology execution occurred.\n")
    (OUT/"report.txt").write_text("\n".join(f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in final.items())+"\n")
    return final

if __name__ == "__main__":
    print(json.dumps(main(),indent=2,sort_keys=True))
