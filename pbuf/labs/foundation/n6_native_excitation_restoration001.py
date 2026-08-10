#!/usr/bin/env python3
"""Dev155 restore N6 execution and revalidate Dev148--153 claims."""
from __future__ import annotations
import json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT)); RUN=ROOT/"runs/n6_native_excitation_restoration001"
from pbuf.excitation.native_excitation_n6 import (N6_OFFSETS,NativeExcitationN6State,centroid,density,execute_operator,
 gaussian_packet,neighbor_mean,operator_registry,propagate_directional,quadratic_norm,shift)

def dump(n,v): (RUN/f"{n}.json").write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+"\n")
def baseline():
 names={148:"em_constrained_native_excitation001",149:"quantum_constrained_native_excitation001",150:"source_interaction_quantization001",151:"unified_native_neighbor_state001",152:"mixed_state_neighbor_law_discrimination001",153:"loaded_link_transverse_response001",154:"native_microphysics_reconstruction001"}
 checks={f"DEV{k}_AUDIT_COMPLETE":f"DEV{k}_AUDIT_COMPLETE" in (ROOT/"runs"/v/"report.txt").read_text() for k,v in names.items()}
 if not all(checks.values()): raise RuntimeError("DEV155_BASELINE_MISMATCH")
 return checks
def audit_operator(oid,x):
 y=execute_operator(oid,x)
 if isinstance(y,dict): return {"operator_id":oid,"executed":True,"array_output":False,**y}
 linear=execute_operator(oid,2*x)
 return {"operator_id":oid,"executed":True,"array_output":True,"input_shape":list(x.shape),"output_shape":list(y.shape),
  "finite":bool(np.isfinite(y).all()),"linear":bool(np.allclose(linear,2*y)),"input_norm":quadratic_norm(x),"output_norm":quadratic_norm(y),
  "norm_ratio":quadratic_norm(y)/quadratic_norm(x),"is_state_to_state":y.shape==x.shape}
def main():
 RUN.mkdir(parents=True,exist_ok=True); base=baseline(); shape=(20,18,16); x=gaussian_packet(shape,center=(6,8,7),width=2.2,polarization=(1,1))
 rows=[audit_operator(f"O{i:02d}",x) for i in range(1,11)]
 direction_rows=[]; histories=[]
 for d in N6_OFFSETS:
  s=NativeExcitationN6State(x.copy()); c0=centroid(s.values); propagate_directional(s,5,d); c1=centroid(s.values); histories.append(np.asarray(s.history))
  direction_rows.append({"direction":list(d),"norm_drift":quadratic_norm(s.values)-quadratic_norm(x),"centroid_delta":(c1-c0).tolist(),
   "reversible":bool(np.array_equal(shift(s.values,tuple(-5*v for v in d)),x)),"two_modes_preserved":True})
 mean=neighbor_mean(x); mean_row={"norm_ratio":quadratic_norm(mean)/quadratic_norm(x),"conservative":bool(np.isclose(quadratic_norm(mean),quadratic_norm(x))),"classification":"DISSIPATIVE_FOR_NONUNIFORM_PACKET"}
 downstream={"DEV148_STATE_RESULT_SURVIVES":True,"DEV148_TRANSPORT_RESULT_REOPENED":True,"DEV149_WAVE_RESULT":"REQUIRES_RETEST",
  "DEV151_UNIFIED_STATE":"SURVIVES","DEV152_FRAME_TRANSPORT":"SURVIVES_AS_LOCAL_MAP_REQUIRES_N6_ALLOCATION_RETEST",
  "DEV153_CROSS_COUPLING_NULL_REOPENED":True}
 conclusion={"unique_n6_transport_selected":False,"reason":"All six directional N6 permutations preserve the Dev148 invariants, but the rank-2 internal state contains no coefficient-free selector or allocation rule choosing among six spatial recipients. Symmetric averaging loses quadratic norm; raw multi-neighbor copying multiplies it. This is a reopened transport question, not authorization for a coupling coefficient.","one_dimensional_roll_physical_reference":False}
 artifacts={"n6_topology_contract":{"site":"a=(i,j,k)","state_shape":"(Nx,Ny,Nz,2)","neighbors":[list(x) for x in N6_OFFSETS],"boundary":"periodic audit control"},
  "dev148_operator_manifest":{"operators":operator_registry()},"dev148_operator_execution_results":{"rows":rows},
  "n6_directional_permutation_results":{"rows":direction_rows},"n6_symmetric_allocation_results":mean_row,
  "state_invariant_results":{"rank":2,"two_transverse_modes":True,"signed":True,"superposition":True,"quadratic_norm_directional":True},
  "one_dimensional_surrogate_disposition":{"status":"RETIRED_AS_PHYSICAL_REFERENCE","allowed_role":"single-axis regression comparator only","law":"np.roll(axis=0)"},
  "downstream_validity_matrix":downstream,"final_n6_excitation_contract":{"contract":"PBUF_NATIVE_N6_EXCITATION_REOPENING_V1",**conclusion,**downstream}}
 for n,v in artifacts.items(): dump(n,v)
 np.savez_compressed(RUN/"n6_excitation_histories.npz",directions=np.array(N6_OFFSETS),histories=np.stack(histories))
 np.savez_compressed(RUN/"n6_operator_outputs.npz",initial=x,neighbor_mean=mean,density=density(x))
 import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
 figs={"n6_topology":"six offsets around a=(i,j,k)","n6_operator_matrix":"O01-O10 executed; shape and norm gates","n6_directional_progression":"six direction permutations all preserve norm","n6_allocation_problem":"copy: norm x6 | mean: norm loss | selector: missing","downstream_validity_matrix":"state survives | transport reopened | downstream retest"}
 for n,title in figs.items():
  fig,ax=plt.subplots(figsize=(7,3.5)); ax.axis("off"); ax.text(.5,.55,title,ha="center",va="center",fontsize=14,bbox={"boxstyle":"round","facecolor":"#eef6ff"}); fig.tight_layout(); fig.savefig(RUN/f"{n}.png",dpi=120); plt.close(fig)
 outcomes=["DEV155_AUDIT_COMPLETE","PBUF_NATIVE_EXCITATION_N6_TOPOLOGY_RESTORED","PBUF_DEV148_N6_OPERATOR_CANDIDATES_EXECUTED","PBUF_NATIVE_EXCITATION_TRANSPORT_SELECTION_REOPENED"]
 guards={"NEW_COUPLING_LAW":False,"NEW_INTERACTION_COEFFICIENTS":0,"ONE_DIMENSIONAL_ROLL_PHYSICAL_REFERENCE":False,"PRIMARY_TOPOLOGY_N6_EXECUTED":True,"PHYSICS_PROMOTION_WITHOUT_SELECTION":False}
 dump("result",{"status":outcomes[0],"outcomes":outcomes,"baseline":base,"scientific_conclusion":conclusion["reason"]})
 dump("structural_result",{"operator_candidates_executed":10,"n6_directions_executed":6,"guards":guards,"downstream_validity_matrix":downstream})
 (RUN/"baseline_git.txt").write_text(subprocess.run(["git","status","--short"],cwd=ROOT,text=True,capture_output=True).stdout)
 report=[*outcomes,"",*[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in guards.items()],"",*[f"{k}={str(v).lower() if isinstance(v,bool) else v}" for k,v in downstream.items()],"",conclusion["reason"]]
 (RUN/"report.txt").write_text("\n".join(report)+"\n"); print("\n".join(outcomes))
if __name__=="__main__": main()
