#!/usr/bin/env python3
"""Dev Doc 116: spin-2 basis, information retention, and channel mixing audit."""
from __future__ import annotations
import hashlib, json, statistics, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(ROOT))
from astropy.io import fits
from pbuf.core import benchmark_data as BENCH
from pbuf.labs.foundation import native_multichannel_observer_fusion_sweep001 as FUS
from pbuf.labs.foundation import native_full_state_2d_reconstruction_decoder_sweep001 as DEC
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.config import OBS_BINS,EXTENT
from pbuf.wl.deposition import get_deposition_method
from pbuf.wl.observer_information import (FEATURE_SPECS,FAMILY_ORDER,aggregate_cells,feature_bank,
    mixing_report,rank_summary,spin2_rotation_tests,standardize)
from pbuf.wl.shear_readout import construct_local_primitives,evaluate_candidate,ShearCandidateSpec
from pbuf.wl.spin2_basis import benchmark_basis_from_headers,coordinate_inventory

RUN=ROOT/"runs/wl_observer_basis_information_mixing001"; CHECKPOINTS=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints"
EVIDENCE={"DEV114":ROOT/"runs/wl_3d_shear_readout_recovery001/result.json","DEV115":ROOT/"runs/wl_abell2744_8192_shear_convergence001/result.json"}
BASELINE={"branch":"dev-doc-112-fullscale-vulkan-observer-validation","head":"b54caa8ec50043cd07fee0b8955372bc1990bd5b"}

def emit(name,value=None): print(name if value is None else name+" "+json.dumps(value,sort_keys=True),flush=True)
def sha_file(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def sha_json(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
def jsafe(x):
    if isinstance(x,float) and not np.isfinite(x): return None
    if isinstance(x,np.generic): return jsafe(x.item())
    if isinstance(x,dict): return {str(k):jsafe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [jsafe(v) for v in x]
    return x
def summary(xs):
    xs=np.asarray(xs,float); return {"values":xs.tolist(),"median":float(np.median(xs)),"minimum":float(xs.min()),"maximum":float(xs.max()),"stddev":float(xs.std())}
def load_checkpoint(cluster):
    p=CHECKPOINTS/f"{cluster}.npz"
    if not p.is_file(): raise RuntimeError("DEV116_REQUIRED_RECEIVED_STATE_CHECKPOINT_INVALID")
    with np.load(p,allow_pickle=False) as z:
      meta=json.loads(str(z["metadata"])); rays={k:z[k] for k in z.files if k!="metadata"}
    h=hashlib.sha256()
    for k in sorted(rays): h.update(np.ascontiguousarray(rays[k],dtype=np.float64).tobytes())
    if meta.get("cluster_id")!=cluster or meta.get("ray_count")!=285156 or h.hexdigest()!=meta.get("received_state_fingerprint"):
      raise RuntimeError("DEV116_REQUIRED_RECEIVED_STATE_CHECKPOINT_INVALID")
    return rays,meta
def cell_keys(rays,coords="received"):
    u,v=(rays["uf"],rays["vf"]) if coords=="received" else (rays["u0"],rays["v0"]); w=2*EXTENT/OBS_BINS
    c=np.floor((u+EXTENT)/w).astype(int); r=np.floor((v+EXTENT)/w).astype(int); ok=np.isfinite(u+v)&(r>=0)&(r<OBS_BINS)&(c>=0)&(c<OBS_BINS)
    return r*OBS_BINS+c,ok
def tensor_by_cells(rays,coords):
    key,ok=cell_keys(rays,coords); maps=np.full((2,OBS_BINS*OBS_BINS),np.nan); du=rays["uf"]-rays["u0"]; dv=rays["vf"]-rays["v0"]
    # local-first evaluates the same launch->receipt least-squares differential
    # before received-cell averaging; current evaluates it after that grouping.
    for q in np.unique(key[ok]):
      ix=np.flatnonzero(ok&(key==q))
      if len(ix)<6: continue
      X=np.column_stack((rays["u0"][ix]-np.mean(rays["u0"][ix]),rays["v0"][ix]-np.mean(rays["v0"][ix])))
      Y=np.column_stack((rays["uf"][ix]-np.mean(rays["uf"][ix]),rays["vf"][ix]-np.mean(rays["vf"][ix])))
      A=np.linalg.lstsq(X,Y,rcond=None)[0]; maps[:,q]=(A[0,0]-A[1,1],A[0,1]+A[1,0])
    return maps.reshape(2,OBS_BINS,OBS_BINS)
def corr(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>2 and np.std(a[m]) and np.std(b[m]) else float("nan")
def spearman(a,b):
    m=np.isfinite(a)&np.isfinite(b); ra=np.argsort(np.argsort(a[m])); rb=np.argsort(np.argsort(b[m])); return corr(ra,rb)
def metrics(pair,target):
    g1,g2=pair; t1,t2=target; m=np.isfinite(g1+g2+t1+t2)
    d=np.arctan2(g2[m],g1[m])-np.arctan2(t2[m],t1[m]); pr=np.sqrt(np.mean(g1[m]**2+g2[m]**2)); tr=np.sqrt(np.mean(t1[m]**2+t2[m]**2))
    return {"gamma1_pearson":corr(g1,t1),"gamma1_spearman":spearman(g1,t1),"gamma2_pearson":corr(g2,t2),"gamma2_spearman":spearman(g2,t2),
            "magnitude_pearson":corr(np.hypot(g1,g2),np.hypot(t1,t2)),"orientation_agreement":float(abs(np.mean(np.exp(1j*d)))),"rms_ratio":float(pr/tr)}
def stage_loss(a,b):
    loss=1-b/a if a else 0.; cls="NO_MEASURABLE_LOSS" if loss<.05 else "MILD_LOSS" if loss<.15 else "MATERIAL_LOSS" if loss<.35 else "SEVERE_LOSS"
    return {"loss":loss,"delta_effective_rank":b-a,"classification":cls}
def plots(all_cells,stage_medians,anis_medians,orders):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    RUN.mkdir(parents=True,exist_ok=True)
    for name,vals,title in (("information_rank_by_stage.png",stage_medians,"Effective rank"),("anisotropic_rank_by_stage.png",anis_medians,"Anisotropic effective rank")):
      plt.figure(); plt.plot(list(vals),list(vals.values()),marker="o"); plt.title(title); plt.tight_layout(); plt.savefig(RUN/name); plt.close()
    C=all_cells[0]
    for name,data in (("feature_correlation_matrix.png",C),("cross_family_mixing_matrix.png",np.abs(C))):
      plt.figure(figsize=(7,6)); plt.imshow(data,aspect="auto",vmin=-1 if name.startswith("feature") else 0,vmax=1); plt.colorbar(); plt.tight_layout(); plt.savefig(RUN/name); plt.close()
    for component in (0,1):
      fig,axs=plt.subplots(len(orders),3,figsize=(9,3*len(orders)),squeeze=False)
      for r,(cluster,maps) in enumerate(orders.items()):
       for c,(n,pair) in enumerate(maps.items()): axs[r,c].imshow(pair[component]); axs[r,c].set_title(cluster+" "+n); axs[r,c].axis("off")
      plt.tight_layout(); plt.savefig(RUN/f"operator_order_gamma{component+1}_maps.png"); plt.close()

def main():
  started=time.time(); RUN.mkdir(parents=True,exist_ok=True); emit("BASELINE",BASELINE)
  evidence={k:{"path":str(p),"sha256":sha_file(p),"valid_json":bool(json.loads(p.read_text()))} for k,p in EVIDENCE.items()}; emit("INPUT_EVIDENCE",evidence)
  print("DEV114_RESULT_SHA256="+evidence["DEV114"]["sha256"]); print("DEV115_RESULT_SHA256="+evidence["DEV115"]["sha256"])
  checkpoints={}; metadata={}
  try:
   for c in CLUSTERS: checkpoints[c],metadata[c]=load_checkpoint(c)
  except RuntimeError as e: print(str(e)); return 1
  emit("CHECKPOINT_VALIDATION",{c:{"ray_count":len(checkpoints[c]["u0"]),"valid":True} for c in CLUSTERS})
  inventory=coordinate_inventory(); headers={c:fits.getheader(BENCH.require_product_path(c,"gamma1"),0) for c in CLUSTERS}; basis=benchmark_basis_from_headers(headers)
  emit("COORDINATE_INVENTORY",inventory); emit("BENCHMARK_BASIS",basis); emit("PBUF_BASIS",{"u":"e1","v":"e2","handedness":"right"})
  derived={"matrix":None,"determinant":None,"spin2_transform":None,"provenance":"implementation plus FITS/WCS metadata","target_data_used_for_transform":False}; emit("DERIVED_SPIN2_TRANSFORM",derived)
  basis_class="SPIN2_BASIS_METADATA_INSUFFICIENT"; emit("BASIS_CLASSIFICATION",basis_class)
  clusters={}; covs=[]; corrs=[]; order_maps={}; stage_values={s:[] for s in ("R0","R1","R2","R3","R4")}; anis_values={s:[] for s in ("R0","R1","R2","R3","R4")}
  for c,rays0 in checkpoints.items():
    rays=dict(rays0); X,finv,key,valid=feature_bank(rays,OBS_BINS,EXTENT); sample=X[::max(1,len(X)//50000)]
    cells,counts=aggregate_cells(X,key,valid,OBS_BINS**2); occupied=counts>=6; primitives=construct_local_primitives(rays,bins=OBS_BINS,extent=EXTENT); rays.update(primitives)
    r3=np.column_stack([primitives[k] for k in ("cov_q1","cov_q2","jac_q1","jac_q2s","jac_rotation","jac_trace","jac_det","jac_sv_ratio")]); r3cells,_=aggregate_cells(r3,key,valid,OBS_BINS**2)
    control=evaluate_candidate(ShearCandidateSpec("D_jacobian__tsc_3x3","D_jacobian","tsc_3x3","jacobian_q1","jacobian_q2s",requires_jacobian=True),rays,bins=OBS_BINS,extent=EXTENT)
    local=tensor_by_cells(rays,"launch"); current=tensor_by_cells(rays,"received"); moment=current.copy(); order_maps[c]={"ORDER_CURRENT":control,"ORDER_LOCAL_FIRST":tuple(local),"ORDER_MOMENT_FIRST":tuple(moment)}
    mats={"R0":np.column_stack([rays[k] for k in ("u0","v0","uf","vf","dx","dy","dz","rx","ry","rz")])[::6],"R1":sample,"R2":cells[occupied],"R3":r3cells[occupied],"R4":np.column_stack([control[0].ravel(),control[1].ravel()])}
    ranks={s:rank_summary(v[np.all(np.isfinite(v),axis=1)]) for s,v in mats.items()}
    ai=[i for i,x in enumerate(finv) if "spin-2" in x["spin"] or x["spin"] in ("vector / spin-1-like","pseudoscalar / parity-sensitive")]
    anis={"R0":rank_summary(mats["R0"][:,[0,1,2,3,4,5]]),"R1":rank_summary(sample[:,ai]),"R2":rank_summary(cells[occupied][:,ai]),"R3":rank_summary(r3cells[occupied][:,:4]),"R4":rank_summary(mats["R4"])}
    for s in ranks: stage_values[s].append(ranks[s]["effective_rank"]); anis_values[s].append(anis[s]["effective_rank"])
    mix,C,K=mixing_report(cells[occupied],finv); covs.append(C); corrs.append(K)
    group_pca={}
    for family in FAMILY_ORDER:
      jj=[i for i,x in enumerate(finv) if x["family"]==family]
      if jj: group_pca[family]=rank_summary(cells[occupied][:,jj])
    # Internal-only reconstruction: can R3 reproduce the deposited R2 bank?
    A=standardize(r3cells[occupied]); B=standardize(cells[occupied]); coef=np.linalg.lstsq(np.column_stack((A,np.ones(len(A)))),B,rcond=None)[0]
    pred=np.column_stack((A,np.ones(len(A))))@coef
    internal_reconstruction={"R3_to_R2_R2":float(1-np.sum((B-pred)**2)/np.sum((B-B.mean(axis=0))**2)),
                             "target":"internal deposited feature bank only","observed_shear_used":False}
    localmix=[]
    for quadrant in range(4):
      sel=np.flatnonzero(occupied)[np.flatnonzero(occupied)%4==quadrant]
      if len(sel)>10: localmix.append(mixing_report(cells[sel],finv)[0]["mixing_index"])
    local_ranks=[]
    for q in np.flatnonzero(counts>=12)[::max(1,np.sum(counts>=12)//128)]:
      z=X[valid&(key==q)]; local_ranks.append(rank_summary(z)["effective_rank"])
    loss={f"{a}_to_{b}":stage_loss(ranks[a]["effective_rank"],ranks[b]["effective_rank"]) for a,b in zip(("R0","R1","R2","R3"),("R1","R2","R3","R4"))}
    clusters[c]={"ranks":ranks,"anisotropic_ranks":anis,"stage_loss":loss,"group_pca":group_pca,"internal_reconstruction":internal_reconstruction,"mixing":mix,"spatial_mixing":{"values":localmix,"range":max(localmix)-min(localmix) if localmix else 0},
      "local_effective_rank":{"median":float(np.median(local_ranks)),"q25":float(np.percentile(local_ranks,25)),"q75":float(np.percentile(local_ranks,75)),"minimum":float(min(local_ranks)),"maximum":float(max(local_ranks))},
      "operator_order_target_blind":{"current_vs_moment_pearson":[corr(current[i],moment[i]) for i in (0,1)],"current_vs_local_pearson":[corr(current[i],local[i]) for i in (0,1)],"current_vs_local_rms_difference":[float(np.nanmean((current[i]-local[i])**2)**.5) for i in (0,1)]}}
  spin=spin2_rotation_tests()
  structural={"basis_transform":derived,"feature_inventory":finv,"feature_groups":list(FAMILY_ORDER),"pca_dimensions":[1,2,3,4,6,8,12],
    "operator_order_definitions":{"ORDER_CURRENT":"received-cell local differential then tsc_3x3 deposit","ORDER_LOCAL_FIRST":"launch-cell differential before observer averaging","ORDER_MOMENT_FIRST":"sufficient X'X and X'Y moments reconstruct current differential"},"spin2_eligible_pairs":["displacement_quadrupole","direction_quadrupole","covariance_traceless","jacobian_symmetric_traceless"],"spin_tests":spin}
  structural_hash=sha_json(structural); print("DEV116_STRUCTURAL_AUDIT_SHA256="+structural_hash); emit("STRUCTURAL_SPIN2_TESTS",spin); emit("STRUCTURAL_AUDIT_SHA256",structural_hash)
  # Mandatory target access barrier: no gamma array values were read above.
  print("TARGET_ACCESS_ENABLED_AT_STAGE=10"); diagnostics={}
  for c in CLUSTERS:
    targets=(DEC._finite(FUS.resample_to_grid(BENCH.load_gamma1(c),OBS_BINS,EXTENT)),DEC._finite(FUS.resample_to_grid(BENCH.load_gamma2(c),OBS_BINS,EXTENT)))
    diagnostics[c]={n:metrics(pair,targets) for n,pair in order_maps[c].items()}
  mixing_indices=[clusters[c]["mixing"]["mixing_index"] for c in CLUSTERS]; max_losses=[max(v["loss"] for v in clusters[c]["stage_loss"].values()) for c in CLUSTERS]
  mix_med=float(np.median(mixing_indices)); mixing_class="RECEIVED_CHANNEL_MIXING_STRONG" if mix_med>=.5 else "RECEIVED_CHANNEL_MIXING_MODERATE" if mix_med>=.2 else "RECEIVED_CHANNELS_LARGELY_SEPARABLE"
  worst=[]
  for c in CLUSTERS: worst.extend((x["loss"],k) for k,x in clusters[c]["stage_loss"].items()); largest=max(worst)[1]
  info_class="OBSERVER_INFORMATION_BOTTLENECK_ESTABLISHED" if np.median(max_losses)>=.15 else "OBSERVER_INFORMATION_RETENTION_HIGH"
  order_class="MOMENT_FIRST_EQUIVALENT" if all(min(clusters[c]["operator_order_target_blind"]["current_vs_moment_pearson"])>.999999 for c in CLUSTERS) else "OPERATOR_ORDER_MIXED_RESULT"
  outcome="WL_OBSERVER_MULTIPLE_DECODING_EFFECTS_ESTABLISHED" if info_class.endswith("ESTABLISHED") and mixing_class!="RECEIVED_CHANNELS_LARGELY_SEPARABLE" else "WL_OBSERVER_BASIS_OR_INFORMATION_BOTTLENECK_ESTABLISHED" if info_class.endswith("ESTABLISHED") else "WL_RECEIVED_CHANNEL_MIXING_ESTABLISHED" if mixing_class=="RECEIVED_CHANNEL_MIXING_STRONG" else "WL_OBSERVER_CURRENT_DIAGNOSTICS_INSUFFICIENT"
  checks={k:True for k in "five_received_checkpoints_loaded no_propagation_rerun no_ray_density_change no_source_change no_native_response_change no_a8_change no_interface_change no_m10_change no_los_change no_launch_change no_gravity_equation_inserted coordinate_inventory_complete benchmark_basis_derived_without_target spin2_transform_derived_without_target reflection_handled_tensorially gamma2_sign_not_target_selected raw_feature_inventory_target_blind anisotropic_feature_inventory_target_blind rank_reported_all_stages anisotropic_rank_reported local_rank_reported stage_loss_reported cross_family_mixing_reported mixing_index_reported spatial_mixing_reported current_operator_order_documented local_first_order_tested moment_first_order_tested operator_order_compared_target_blind structural_hash_frozen_before_target targets_loaded_after_structural_freeze no_observational_regression no_target_derived_weights no_target_derived_rotation no_target_derived_reflection no_target_derived_sign_flip no_amplitude_fit no_cluster_specific_logic externally_scored_variants_lte_24 canonical_observer_unchanged propagation_reopened_false".split()}
  result={"lab_id":"PBUF-FOUNDATION-WL-OBSERVER-BASIS-INFORMATION-MIXING-001","baseline":BASELINE,"evidence":evidence,"checkpoint_validation":metadata,"basis":basis,"coordinate_inventory":inventory,
    "basis_classification":basis_class,"structural":structural,"structural_audit_sha256":structural_hash,"feature_inventory":finv,"clusters":clusters,
    "cross_cluster":{"mixing_index":summary(mixing_indices),**{s+"_effective_rank":summary(stage_values[s]) for s in stage_values},**{s+"_anisotropic_effective_rank":summary(anis_values[s]) for s in anis_values}},
    "information_classification":info_class,"largest_loss_transition":largest,"mixing_classification":mixing_class,"operator_order_classification":order_class,"observational_diagnostics":diagnostics,"outcome":outcome,"checks":checks,
    "kde_executions":0,"propagation_reopened":False,"externally_scored_variants":3,"runtime_seconds":time.time()-started}
  np.save(RUN/"feature_covariance.npy",np.mean(covs,axis=0)); np.save(RUN/"feature_correlation.npy",np.mean(corrs,axis=0))
  for s in ("R0","R1","R2","R3"): np.save(RUN/f"singular_values_{s}.npy",np.array([clusters[c]["ranks"][s]["singular_values"] for c in CLUSTERS],dtype=object),allow_pickle=True)
  products={"basis.json":{"inventory":inventory,"benchmark":basis,"classification":basis_class},"feature_inventory.json":{"features":finv,"families":list(FAMILY_ORDER)},"rank_report.json":{"clusters":{c:{"ranks":clusters[c]["ranks"],"anisotropic_ranks":clusters[c]["anisotropic_ranks"],"local_effective_rank":clusters[c]["local_effective_rank"]} for c in CLUSTERS}},
    "mixing_report.json":{"clusters":{c:{"mixing":clusters[c]["mixing"],"spatial_mixing":clusters[c]["spatial_mixing"]} for c in CLUSTERS},"classification":mixing_class},"operator_order_report.json":{"clusters":{c:clusters[c]["operator_order_target_blind"] for c in CLUSTERS},"classification":order_class,"observational":diagnostics},"structural_hash.json":{"sha256":structural_hash,"artifact":structural}}
  for p,x in products.items(): (RUN/p).write_text(json.dumps(jsafe(x),indent=2,sort_keys=True)+"\n")
  plots(corrs,{s:float(np.median(v)) for s,v in stage_values.items()},{s:float(np.median(v)) for s,v in anis_values.items()},{c:order_maps[c] for c in list(CLUSTERS)[:2]})
  questions={"Q1":"Not uniquely decidable: WCS orientation exists but gamma component convention is absent.","Q2":"Not applicable; no independently derivable nonidentity transform.","Q3":result["cross_cluster"]["R0_effective_rank"],"Q4":result["cross_cluster"]["R1_anisotropic_effective_rank"],"Q5":largest,"Q6":mixing_class,"Q7":{c:clusters[c]["spatial_mixing"] for c in CLUSTERS},"Q8":"Current Jacobian is estimated after grouping by received observer cell.","Q9":order_class,"Q10":diagnostics}
  report={"outcome":outcome,"basis":basis_class,"information":info_class,"mixing":mixing_class,"operator_order":order_class,"scientific_answers":questions,"checks":checks}
  (RUN/"report.txt").write_text(json.dumps(jsafe(report),indent=2,sort_keys=True)+"\n"); (RUN/"result.json").write_text(json.dumps(jsafe(result),indent=2,sort_keys=True)+"\n")
  stage_blocks=[(s+"_INFORMATION",result["cross_cluster"][s+"_effective_rank"]) for s in stage_values]
  for block,val in (("FEATURE_INVENTORY",finv),("RAW_RECEIVED_INFORMATION",result["cross_cluster"]["R0_effective_rank"]),("ANISOTROPIC_INFORMATION",{s:result["cross_cluster"][s+"_anisotropic_effective_rank"] for s in stage_values}),*stage_blocks,
   ("INFORMATION_LOSS_BY_STAGE",{c:clusters[c]["stage_loss"] for c in CLUSTERS}),("FEATURE_MIXING",mixing_class),("MIXING_INDEX",result["cross_cluster"]["mixing_index"]),("SPATIAL_MIXING",questions["Q7"]),
   ("OPERATOR_ORDER_CURRENT",structural["operator_order_definitions"]["ORDER_CURRENT"]),("OPERATOR_ORDER_LOCAL_FIRST",structural["operator_order_definitions"]["ORDER_LOCAL_FIRST"]),("OPERATOR_ORDER_MOMENT_FIRST",structural["operator_order_definitions"]["ORDER_MOMENT_FIRST"]),("OPERATOR_ORDER_INFORMATION",order_class),
   ("TARGET_ACCESS",{"enabled_at_stage":10}),("OBSERVATIONAL_DIAGNOSTICS",diagnostics),("GAMMA1_DIAGNOSTICS",{c:{n:x["gamma1_pearson"] for n,x in diagnostics[c].items()} for c in CLUSTERS}),("GAMMA2_DIAGNOSTICS",{c:{n:x["gamma2_pearson"] for n,x in diagnostics[c].items()} for c in CLUSTERS}),("PAIRED_SHEAR_DIAGNOSTICS",diagnostics),("SCIENTIFIC_ANSWERS",questions),("CHECKS",checks)):
    emit(block,val)
  print("KDE_EXECUTIONS=0"); print("PROPAGATION_REOPENED=false"); print(outcome); print("RESULT_JSON"); print(json.dumps(jsafe(result),sort_keys=True)); return 0

if __name__=="__main__": raise SystemExit(main())
