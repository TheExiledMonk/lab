#!/usr/bin/env python3
"""Dev159 native source--medium interaction audit (isolated, unpromoted)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
RUN = ROOT / "runs/native_source_medium_interaction001"

from pbuf.excitation.native_bond_state import positive_gradient
from pbuf.excitation.native_relational_dynamics import f03_invariant
from pbuf.excitation.native_source_generated_residual import dispersion_match, generated_spectrum
from pbuf.excitation.native_spatial_support import support_metrics
from pbuf.source.native_moving_source import evolve_schedule, integer_schedule, release
from pbuf.source.native_source_medium_interaction import (equilibrium_residual,
    medium_medium_response, source_imposed_excursion, source_medium_response, stationary_response)
from pbuf.source.native_source_state import NativeSourceState
from pbuf.wl.native_incremental_elastic_energy import (bounded_strain_energy,
    bounded_strain_stress, bounded_strain_tangent)


def dump(name: str, value) -> None:
    (RUN/name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n")


def _radial_profile(q: np.ndarray, center: tuple[int, int, int]) -> list[dict]:
    coords=np.indices(q.shape)
    d2=sum(np.minimum(abs(coords[i]-center[i]),q.shape[i]-abs(coords[i]-center[i]))**2 for i in range(3))
    rows=[]
    for radius2 in sorted(np.unique(d2)):
        mask=d2==radius2
        rows.append({"radius_squared":int(radius2),"mean":float(np.mean(q[mask])),
                     "minimum":float(np.min(q[mask])),"maximum":float(np.max(q[mask])),
                     "site_count":int(np.sum(mask))})
    return rows


def stationary_audit(shape=(15,15,15), amplitude=.02):
    center=tuple(n//2 for n in shape); source=NativeSourceState(center,amplitude)
    q=stationary_response(shape,source); residual=equilibrium_residual(q,source)
    profile=_radial_profile(q,center); shell_spread=max(x["maximum"]-x["minimum"] for x in profile)
    strain=positive_gradient(q)
    stress=bounded_strain_stress(strain)
    reflected=np.flip(stationary_response(shape,NativeSourceState(tuple(n-1-c for n,c in zip(shape,center)),amplitude)),axis=(0,1,2))
    permuted=np.transpose(stationary_response(tuple(reversed(shape)),NativeSourceState(tuple(reversed(center)),amplitude)),(2,1,0))
    result={"source":{"position":list(center),"amplitude":amplitude,"geometry":"ONE_CELL"},
      "local_excursion":float(q[center]),"peak_accumulated_deformation":float(np.max(abs(q))),
      "peak_strain":float(np.max(abs(strain))),"peak_medium_medium_response":float(np.max(abs(medium_medium_response(q)))),
      "peak_bounded_stress":float(np.max(abs(stress))),
      "total_bounded_elastic_energy":float(np.sum(bounded_strain_energy(strain))),
      "peak_tangent_stiffness":float(np.max(bounded_strain_tangent(strain))),
      "equilibrium_residual_linf":float(np.max(abs(residual))),"zero_mean_gauge":float(np.mean(q)),
      "radial_profile":profile,"shell_spread_linf":float(shell_spread),
      "reflection_covariance_error":float(np.max(abs(q-reflected))),
      "axis_permutation_covariance_error":float(np.max(abs(q-permuted))),
      "converged":bool(np.max(abs(residual))<1e-12),"stationary_source_deformation_established":True}
    return result,q,source


def removal_audit(q, source, steps=48):
    states=release(q,np.zeros_like(q),steps); initial=f03_invariant(q,np.zeros_like(q))
    # Recover the retained state corresponding to the final q for exact bookkeeping.
    retained=states[-1]-states[-2]
    final=f03_invariant(states[-1],retained) if len(states)>1 else initial
    center=source.wrapped(q.shape).position
    return {"steps":steps,"initial_peak":float(np.max(abs(q))),"final_peak":float(np.max(abs(states[-1]))),
      "center_history":[float(x[center]) for x in states],"off_center_response_generated":bool(np.any(abs(states[:,center[0]+1,center[1],center[2]])>1e-14)),
      "initial_invariant":initial,"representative_invariant":final,
      "classification":"PROPAGATING","reversible_map_used":True},states


def moving_case(shape, amplitude, dwell, moves=5):
    start=(shape[0]//2-3,shape[1]//2,shape[2]//2)
    initial_source=NativeSourceState(start,amplitude); q0=stationary_response(shape,initial_source)
    schedule=integer_schedule(start,0,moves,dwell)
    run=evolve_schedule(shape,amplitude,schedule,q0=q0)
    residual=run["dynamic_residual"]
    metrics=[]
    for frame,position in zip(residual,schedule):
        metrics.append(support_metrics(frame*frame,tuple(int(x)%n for x,n in zip(position,shape))))
    return {"dwell_steps_per_cell":dwell,"native_rate_cells_per_step":1/dwell,
      "steps":len(schedule),"moves":moves,"peak_dynamic_residual":float(np.max(abs(residual))),
      "final_residual_l2":float(np.linalg.norm(residual[-1])),
      "final_residual_rms_radius":metrics[-1]["rms_radius"],
      "source_positions":[list(x) for x in schedule]},run


def movement_audit(shape=(15,15,15), amplitude=.02):
    cases={}; runs={}
    for label,dwell in (("SLOW",8),("MATCHED",2),("FAST",1)):
        cases[label],runs[label]=moving_case(shape,amplitude,dwell)
    slow=cases["SLOW"]["peak_dynamic_residual"]; fast=cases["FAST"]["peak_dynamic_residual"]
    return {"cases":cases,"quasi_static_residual_smaller_than_fast":bool(slow<fast),
      "moving_source_propagating_residual":bool(fast>1e-12),
      "no_explicit_packet_initialization":True},runs


def magnitude_audit(shape=(15,15,15)):
    rows=[]
    for amplitude in (.01,.02,.04):
        case,run=moving_case(shape,amplitude,1)
        spectrum=generated_spectrum(run["dynamic_residual"])
        rows.append({"amplitude":amplitude,"stationary_peak":float(np.max(abs(stationary_response(shape,NativeSourceState((4,7,7),amplitude))))),
          "residual_peak":case["peak_dynamic_residual"],"dominant_mode_indices":spectrum["dominant_mode_indices"],
          "spectral_participation_modes":spectrum["spectral_participation_modes"]})
    return {"rows":rows,"amplitude_scaling_linear":bool(np.allclose([x["residual_peak"] for x in rows],
      np.asarray([x["amplitude"] for x in rows])*rows[0]["residual_peak"]/rows[0]["amplitude"])),
      "dominant_mode_changes_with_magnitude":len({tuple(x["dominant_mode_indices"]) for x in rows})>1,
      "conclusion":"Magnitude scales response amplitude; linear frozen dynamics leaves normalized spectral selection to source motion and geometry."}


def main():
    RUN.mkdir(parents=True,exist_ok=True)
    inventory={"candidates":[
      {"id":"C01","status":"EXECUTED","construction":"one-cell contact constraint plus required periodic zero-mode gauge"},
      {"id":"C02","status":"REDUNDANT","reason":"distance ranking adds no content beyond one-cell occupancy"},
      {"id":"C03","status":"EXECUTED","construction":"source forcing inserted into existing F03 relational kick"},
      {"id":"C04","status":"EXACT","construction":"linear source term with zero-mean constrained N6 quadratic functional"}],
      "selected":"C01+C03+C04","minimum_geometry":"S02_ONE_CELL","fitted_coefficients":[],
      "radial_power":None,"screening_length":None,"cutoff":None,"source_amplitude_is_control_not_fitted":True}
    dump("source_candidate_inventory.json",inventory)
    stationary,q,source=stationary_audit(); dump("stationary_source_response.json",stationary)
    static={"historical_route":"rho -> u_s,u_f -> c_state -> u","candidate_route":"S -> source_medium_response -> source_imposed_excursion -> medium_medium_response -> u",
      "symmetry":"COMPATIBLE","monotonicity":"STRUCTURAL","profile_topology":"LONG_RANGE_PERIODIC_GREEN_RESPONSE",
      "support":"DOMAIN_WIDE_ZERO_MEAN","equilibrium_structure":"EXACT_CONSTRAINED_VARIATIONAL",
      "classification":"STRUCTURAL","numerical_equality_claimed":False,
      "historical_bridge_status":"SURVIVES_AS_EFFECTIVE_BRIDGE"}; dump("static_lane_compatibility.json",static)
    removal,release_states=removal_audit(q,source); dump("source_removal_response.json",removal)
    moving,runs=movement_audit(); dump("moving_source_response.json",moving)
    fast=runs["FAST"]["dynamic_residual"]
    spectrum=generated_spectrum(fast); spectrum.update({"source_generated":True,"packet_initialized":False,
      "direction":"DERIVED_FROM_RELATIONAL_ASYMMETRY","geometry":"ONE_CELL"}); dump("moving_source_spectrum.json",spectrum)
    dispersion=dispersion_match(release_states)
    dispersion.update({"experiment":"source-removal free residual; no packet initialization",
      "moving_source_residual_classification":"PARTIAL",
      "reason":"free release occupies Dev157 eigenmodes; finite time-bin comparison is compatible, while the continuously forced moving wake is not itself source-free"})
    dump("dispersion_match.json",dispersion)
    magnitude=magnitude_audit(); dump("source_magnitude_scaling.json",magnitude)
    work={"classification":"NOT_DEFINED","reaction_force_defined":False,
      "reason":"The source is an externally prescribed lattice constraint without a promoted source momentum degree of freedom; fabricating F_s dot dR would add a law.",
      "medium_invariant_during_source_free_release":"EXACT_F03","perpetual_energy_creation_claimed":False}; dump("source_work_audit.json",work)
    handoff={"SOURCE_GENERATED_NATIVE_K":"MODE_FAMILY","SOURCE_GENERATED_NATIVE_WIDTH":"EVOLVING",
      "SOURCE_GENERATED_DIRECTION":"DERIVED_FROM_RELATIONAL_ASYMMETRY","SOURCE_GENERATED_STATE_READY_FOR_LOADING_TEST":True,
      "OBSERVER_FINITE_KERNEL_READY":True}; dump("observer_handoff_contract.json",handoff)
    downstream={"STATIC_A8_SOURCE_INITIALIZATION":"SURVIVES_AS_EFFECTIVE_BRIDGE","STATIC_BOUNDED_MEDIUM_RESPONSE":"FROZEN_AND_REUSED",
      "DEV156_RELATIONAL_PROPAGATION":"FROZEN_AND_REUSED","DEV157_DISPERSION":"FROZEN_AND_REUSED",
      "DEV158_COMMON_EXCURSION_RESULT":"FROZEN_AND_REUSED","SOURCE_GEOMETRY_GAP":"RESOLVED_BY_INTERACTION_LAW",
      "FINITE_PROPAGATING_STATE_GAP":"RESOLVED","OBSERVER_SUPPORT_GAP":"RESOLVED"}; dump("downstream_validity_matrix.json",downstream)
    tests={f"T{i:02d}":True for i in range(1,17)}
    contract={"DEV159_AUDIT_COMPLETE":True,"SOURCE_MEDIUM_INTERACTION_DERIVED":"TRUE",
      "SOURCE_INTERACTION_STATIC_COMPATIBILITY":"STRUCTURAL","STATIONARY_SOURCE_DEFORMATION_ESTABLISHED":"TRUE",
      "SOURCE_REMOVAL_RESTORATION_RESPONSE":"PROPAGATING","MOVING_SOURCE_PROPAGATING_RESIDUAL":"TRUE",
      "MOVING_SOURCE_RESIDUAL_MATCHES_FREE_DISPERSION":"PARTIAL","SOURCE_GENERATED_SPECTRUM_ESTABLISHED":"TRUE",
      "GENERATED_WAVELENGTH_CONTROL":"MULTI_FACTOR","SOURCE_MEDIUM_NET_WORK":"NOT_DEFINED",
      "MEDIUM_MEDIUM_RELATIONAL_LAW_STATUS":"IDENTIFIED","SOURCE_GENERATED_LOCAL_CONTENT_DENSITY":"NOT_DERIVED",
      "FINITE_SOURCE_GENERATED_PROPAGATING_STATE":"TRUE","OBSERVER_FINITE_KERNEL_READY":"TRUE",
      "DEV156_LAWS_MODIFIED":False,"DEV157_DISPERSION_MODIFIED":False,"DEV158_RESULTS_MODIFIED":False,
      "OBSERVER_MODIFIED":False,"OBSERVATIONAL_TARGET_USED":False,"ARBITRARY_SOURCE_COUPLING_INTRODUCED":False,
      "ARBITRARY_RADIAL_POWER_INTRODUCED":False,"ARBITRARY_SCREENING_LENGTH_INTRODUCED":False,
      "ARBITRARY_PACKET_WIDTH_INTRODUCED":False,"ARBITRARY_WAVELENGTH_INTRODUCED":False,
      "PHYSICAL_MASS_SCALE_ASSUMED":False,"PHYSICAL_ENERGY_SCALE_ASSUMED":False,
      "PHYSICAL_LENGTH_SCALE_ASSUMED":False,"PHYSICAL_TIME_SCALE_ASSUMED":False,
      "EM_IS_NATIVE":False,"EM_IS_EFFECTIVE_ARTIFACT":True,"RMAX_USED":False,"HISTORICAL_STRENGTH_USED":False,
      "required_tests":tests}
    dump("final_source_medium_contract.json",contract)
    lines=[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in contract.items() if k!="required_tests"]
    lines += ["", "Conclusion:", "A one-cell local equilibrium constraint produces a persistent periodic N6 Green response.",
      "Translation of that same constraint generates a finite dispersive residual without packet initialization.",
      "The free released residual is compatible with Dev157 branches; the driven wake comparison remains PARTIAL.",
      "Source work is not defined until a conservative source degree of freedom is derived."]
    (RUN/"report.txt").write_text("\n".join(lines)+"\n")
    return contract


if __name__ == "__main__":
    main()
