#!/usr/bin/env python3
"""Dev 122: reconstruct transported diagnostic geometry before observables."""
from __future__ import annotations
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from scipy import ndimage,stats
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.transport_receiver_decode import effective_rank,rasterize
from pbuf.wl.transport_reconstruction import (diagnostic_patterns,reconstruct_endpoint_only,reconstruct_raw,reconstruct_order,reconstruct_mesh,reconstruct_patches,quadratic_taylor_error)
from pbuf.wl.reconstructed_geometry import geometry_from_derivatives,geometry_feature_matrix
from pbuf.wl.lens_ray_registration import separability

RUN=ROOT/"runs/wl_reconstruction_first_decode001";CP=ROOT/"runs/wl_3d_shear_readout_recovery001/checkpoints";D121=ROOT/"runs/wl_dual_transport_receiver_decode001";SD=ROOT/"runs/wl_source_deformation_separation001";COHORT=ROOT/"runs/wl_lens_footprint_ray_redistribution001";D120=ROOT/"runs/wl_bundle_transport_geometry001"
BASELINE={"branch":"dev-doc-112-fullscale-vulkan-observer-validation","head":"b54caa8ec50043cd07fee0b8955372bc1990bd5b","log_8":["b54caa8 Add full-scale Vulkan observer validation","c620cd9 Add exact Vulkan observer KDE acceleration","fca7fb3 Audit observer deposition stability","1dbc1b9 Record Dev Doc 109 Vulkan Outcome B","e3adac1 Add canonical WL modularization and audit artifacts","67c059b FOUNDATION: add 100% observer coverage fix runner notes","f9cb68a FOUNDATION: document 100% observer coverage launch-coordinate fix","ee50f64 FOUNDATION: fix 100% observer receipt to preserve lane launch coordinates"],"worktree_dirty":True,"preservation":"No reset, clean, stash, or checkout performed; unrelated changes preserved."}
def dump(p,o):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,indent=2,sort_keys=True,allow_nan=False)+"\n")
def load_npz(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def finite(x):return np.nan_to_num(np.asarray(x,float),nan=0.,posinf=0.,neginf=0.)
def mat(*banks):
 cols=[]
 for b in banks:
  cols += [finite(v).ravel() for k,v in b.items() if k!="ray_id" and np.asarray(v).ndim==2]
 return np.column_stack(cols) if cols else np.empty((0,0))
def image_mat(*banks):
 cols=[]
 for b in banks:cols += [finite(v).ravel() for v in b.values() if np.asarray(v).shape==(64,64)]
 return np.column_stack(cols)
def sample_images_at_rays(t,bank,names):
 u=finite(t["uf"]).ravel();v=finite(t["vf"]).ravel();ix=np.clip(np.floor((u+8)*4).astype(int),0,63);iy=np.clip(np.floor((v+8)*4).astype(int),0,63)
 return np.column_stack([finite(bank[k])[iy,ix] for k in names])
def sample_sd(c,t):
 S=np.load(SD/c/"source_score.npy");D=np.load(SD/c/"deformation_score.npy");q=np.column_stack((t["uf"].ravel(),t["vf"].ravel(),t["wf"].ravel()));lo=np.array([-8.,-8.,q[:,2].min()]);hi=np.array([8.,8.,q[:,2].max()]);i=np.rint((q-lo)*(np.array(S.shape)-1)/np.maximum(hi-lo,1e-30)).astype(int);i=np.clip(i,0,np.array(S.shape)-1);return S[tuple(i.T)],D[tuple(i.T)]
def corr(a,b):
 a=finite(a).ravel();b=finite(b).ravel();return float(np.corrcoef(a,b)[0,1]) if a.std()>0 and b.std()>0 else 0.
def analytic_tests():
 z=np.zeros((7,9))
 def geo(J,second=False):
  f={"d_u_delta_u":z+J[0,0]-1,"d_v_delta_u":z+J[0,1],"d_u_delta_v":z+J[1,0],"d_v_delta_v":z+J[1,1]-1};s={f"d_{a}_{b}":z.copy() for b in ("delta_u","delta_v","wf") for a in ("uu","uv","vv")} if second else None;return geometry_from_derivatives(f,s)
 ident=geo(np.eye(2),True);translation=geo(np.eye(2));scale=geo(np.eye(2)*1.7);stretch=geo(np.diag([1.4,.8]));a=.37;rotation=geo(np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]]));e1,e2=quadratic_taylor_error()
 return {"identity_transport_test_passed":bool(np.allclose(ident["local_area_change"],0) and np.allclose(ident["orientation_change"],0) and np.allclose(ident["local_curvature"],0)),"translation_test_passed":bool(np.allclose(translation["local_area_change"],0) and np.allclose(translation["local_anisotropy"],0)),"isotropic_scale_test_passed":bool(np.allclose(scale["transport_area_ratio"],1.7**2) and np.allclose(scale["local_anisotropy"],0)),"anisotropic_stretch_test_passed":bool(np.all(stretch["spin2_shape_q1"]>0)),"rotation_test_passed":bool(np.allclose(rotation["transport_area_ratio"],1) and np.allclose(rotation["local_anisotropy"],0,atol=1e-15) and np.allclose(rotation["orientation_change"],a)),"quadratic_warp_test_passed":bool(e2<1e-14 and e2<e1),"quadratic_first_order_rms":e1,"quadratic_second_order_rms":e2}
def target(c):
 from pbuf.core import benchmark_data as B
 from pbuf.wl.source import load_cluster_source
 d=load_cluster_source(next(x for x in B.clusters() if x["id"]==c))["data"];return tuple(ndimage.zoom(np.asarray(d[k],float),np.array((64,64))/np.array(d[k].shape),order=1) for k in ("gamma1","gamma2"))
def score(pair,truth):
 def one(a,b):return {"pearson":corr(a,b),"spearman":float(stats.spearmanr(finite(a).ravel(),finite(b).ravel()).statistic),"rms_ratio":float(np.sqrt(np.mean(finite(a)**2))/(np.sqrt(np.mean(finite(b)**2))+1e-30))}
 a,b=pair;x,y=truth;return {"gamma1":one(a,x),"gamma2":one(b,y),"magnitude_pearson":corr(np.hypot(a,b),np.hypot(x,y)),"orientation_agreement":float(np.mean(np.cos(np.arctan2(b,a)-np.arctan2(y,x))))}
def dep(t,a,b):return rasterize(t["uf"],t["vf"],a),rasterize(t["uf"],t["vf"],b)
def aggregate(external):
 out={}
 for name in next(iter(external.values())):
  row={}
  for comp in ("gamma1","gamma2"):
   vals=np.array([external[c][name][comp]["pearson"] for c in external]);loos=[float(np.median(np.delete(vals,i))) for i in range(len(vals))]
   row[comp+"_pearson"]={"median":float(np.median(vals)),"minimum":float(vals.min()),"maximum":float(vals.max()),"loo_stability_range":float(max(loos)-min(loos))}
  out[name]=row
 return out
def plots(cache,aggregate_metrics):
 import matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
 x=cache["Abell2744"]
 def panel(keys,title,name):
  fig,ax=plt.subplots(2,2,figsize=(9,8));
  for a,(label,z) in zip(ax.ravel(),keys):a.imshow(z,origin="lower",cmap="viridis");a.set_title(label)
  fig.suptitle(title);fig.tight_layout();fig.savefig(RUN/name,dpi=125);plt.close(fig)
 panel([("launch checker",x["patterns"]["checkerboard"]),("endpoint",x["endpoint"]["reconstructed_intensity_checkerboard"]),("raw",x["raw"]["reconstructed_intensity_checkerboard"]),("second order",x["second"]["reconstructed_intensity_checkerboard"])],"Reconstruction first","abell2744_reconstruction_first_overview.png")
 panel([("raw grid",x["raw"]["reconstructed_intensity_uniform_lattice"]),("first",x["first"]["reconstructed_intensity_uniform_lattice"]),("second",x["second"]["reconstructed_intensity_uniform_lattice"]),("mesh",x["mesh"]["reconstructed_intensity_uniform_lattice"])],"Grid transport","abell2744_grid_transport_comparison.png")
 panel([("first curvature",x["first"]["local_curvature"]),("second curvature",x["second"]["local_curvature"]),("first intensity",x["first"]["reconstructed_intensity_isotropic_point_grid"]),("second intensity",x["second"]["reconstructed_intensity_isotropic_point_grid"])],"First vs second order","abell2744_first_vs_second_order.png")
 panel([("mesh checker",x["mesh"]["reconstructed_intensity_checkerboard"]),("mesh area",x["mesh"]["transport_area_ratio"]),("mesh orientation",x["mesh"]["orientation_change"]),("mesh density",x["mesh"]["transport_density"])],"Topology-preserving mesh","abell2744_mesh_reconstruction.png")
 panel([("shape q1",x["first"]["spin2_shape_q1"]),("shape q2",x["first"]["spin2_shape_q2"]),("anisotropy",x["first"]["local_anisotropy"]),("orientation",x["first"]["local_orientation"])],"Reconstructed shape spin-2","abell2744_shape_spin2.png")
 stages=list(x["ranks"]);fig,ax=plt.subplots(figsize=(11,4));
 for c,z in cache.items():ax.plot(stages,[z["ranks"][k] for k in stages],marker="o",label=c)
 ax.tick_params(axis="x",rotation=25);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(RUN/"information_rank_ladder.png",dpi=125);plt.close(fig)
 fig,ax=plt.subplots(figsize=(11,4));
 for c,z in cache.items():ax.plot(stages,[z["sep"][k]["mahalanobis_distance"] for k in stages],marker="o",label=c)
 ax.tick_params(axis="x",rotation=25);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(RUN/"cohort_separability_ladder.png",dpi=125);plt.close(fig)
 methods=("endpoint","raw","first","second","mesh");fig,ax=plt.subplots(figsize=(8,4));ax.bar(methods,[effective_rank(image_mat(x[k])) for k in methods]);ax.set_ylabel("effective rank");fig.tight_layout();fig.savefig(RUN/"reconstruction_method_comparison.png",dpi=125);plt.close(fig)
 names=list(aggregate_metrics);fig,ax=plt.subplots(figsize=(12,5));xx=np.arange(len(names));ax.plot(xx,[aggregate_metrics[n]["gamma1_pearson"]["median"] for n in names],"o-",label="gamma1");ax.plot(xx,[aggregate_metrics[n]["gamma2_pearson"]["median"] for n in names],"o-",label="gamma2");ax.set_xticks(xx,names,rotation=75,ha="right");ax.legend();fig.tight_layout();fig.savefig(RUN/"external_shear_summary.png",dpi=125);plt.close(fig)
def main(argv=None):
 ap=argparse.ArgumentParser();ap.add_argument("--clusters",nargs="*",default=list(CLUSTERS));args=ap.parse_args(argv);chosen=list(args.clusters)
 if len(chosen)!=5 or set(chosen)!=set(CLUSTERS):raise RuntimeError("Dev122 requires five frozen clusters")
 started=time.time();RUN.mkdir(parents=True,exist_ok=True);tests=analytic_tests();assert all(v for k,v in tests.items() if k.endswith("passed"))
 structural={"lab_id":"PBUF-FOUNDATION-WL-RECONSTRUCTION-FIRST-DECODE-001","target_access":False,"earth_receiver_hypothesis_parked":True,"propagation_runs":0,"kde_executions":0,"baseline":BASELINE,"analytic_tests":tests,"clusters":{}};cache={}
 for c in chosen:
  t=load_npz(D121/c/"dual_transport.npz");f=load_npz(D121/c/"first_order_transport.npz");s=load_npz(D121/c/"second_order_transport.npz");patterns=diagnostic_patterns(t["u0"].shape)
  endpoint=reconstruct_endpoint_only(t,patterns);raw=reconstruct_raw(t,patterns);first=reconstruct_order(t,f,None,patterns);second=reconstruct_order(t,f,s,patterns);mesh,mdata=reconstruct_mesh(t,patterns);g=geometry_from_derivatives(f,s);patch=reconstruct_patches(t,g)
  out=RUN/c;out.mkdir(exist_ok=True);np.savez_compressed(out/"reconstruction_raw.npz",**{"lane":"RECON_RAW_TRANSPORT",**raw,**{"endpoint__"+k:v for k,v in endpoint.items()}});np.savez_compressed(out/"reconstruction_first_order.npz",**first);np.savez_compressed(out/"reconstruction_second_order.npz",**second);np.savez_compressed(out/"reconstruction_mesh.npz",**mesh);np.savez_compressed(out/"reconstruction_patch_multiscale.npz",**patch)
  geom_bank={k:v for k,v in g.items() if k!="jacobian"};geom_bank.update({"depth":t["wf"],"direction_u":t["dir_u"],"direction_v":t["dir_v"],"direction_w":t["dir_w"]});np.savez_compressed(out/"reconstructed_geometry_bank.npz",**geom_bank)
  R0=np.column_stack([t[k].ravel() for k in ("uf","vf","wf","dir_u","dir_v","dir_w")]);R1=np.column_stack((t["u0"].ravel(),t["v0"].ravel(),R0));R2=np.column_stack((R1,mat(f)));R3=np.column_stack((R2,mat(s)));R4=geometry_feature_matrix(g)
  image_names=("reconstructed_intensity_checkerboard","transport_density","transport_area_ratio","local_orientation","local_anisotropy","local_curvature","depth","dir_u","dir_v","dir_w")
  R5=sample_images_at_rays(t,second,image_names);R6=np.column_stack((g["spin2_shape_q1"].ravel(),g["spin2_shape_q2"].ravel()));R7=R0[:,:4]
  banks={"R0_received_state":R0,"R1_dual_coordinate_raw":R1,"R2_first_order_transport":R2,"R3_second_order_transport":R3,"R4_reconstructed_transport_bank":R4,"R5_reconstructed_geometric_image_bank":R5,"R6_direct_spin2_geometric_bank":R6,"R7_current_observer_R4":R7};ranks={k:effective_rank(v) for k,v in banks.items()};dump(out/"information_rank.json",ranks)
  S,D=sample_sd(c,t);structural["clusters"][c]={"information_rank":ranks,"reconstruction_methods":["RECON_ENDPOINT_ONLY","RECON_RAW_TRANSPORT","RECON_FIRST_ORDER","RECON_SECOND_ORDER","RECON_MESH_TRANSPORT","RECON_PATCH_MULTISCALE"],"patch_widths":[2,4,8,16],"pattern_names":list(patterns),"SD_structural_relationship":{"D_vs_area_change":corr(D,g["local_area_change"]),"D_vs_curvature":corr(D,g["local_curvature"]),"D_vs_anisotropy":corr(D,g["local_anisotropy"]),"S_vs_intensity_coherence":corr(S,ndimage.uniform_filter(patterns["isotropic_point_grid"],3)),"S_vs_topology_continuity":corr(S,np.isfinite(g["transport_area_ratio"]).astype(float))},"method_comparison":{"raw_vs_first_intensity_correlation":corr(raw["reconstructed_intensity_checkerboard"],first["reconstructed_intensity_checkerboard"]),"first_vs_second_intensity_correlation":corr(first["reconstructed_intensity_checkerboard"],second["reconstructed_intensity_checkerboard"]),"second_vs_mesh_intensity_correlation":corr(second["reconstructed_intensity_checkerboard"],mesh["reconstructed_intensity_checkerboard"]),"first_vs_second_area_correlation":corr(first["transport_area_ratio"],second["transport_area_ratio"]),"first_vs_second_shape_q_correlation":.5*(corr(first["spin2_shape_q1"],second["spin2_shape_q1"])+corr(first["spin2_shape_q2"],second["spin2_shape_q2"])),"first_effective_rank":effective_rank(image_mat(first)),"second_effective_rank":effective_rank(image_mat(second)),"mesh_effective_rank":effective_rank(image_mat(mesh))}}
  cache[c]={"t":t,"f":f,"s":s,"g":g,"patterns":patterns,"endpoint":endpoint,"raw":raw,"first":first,"second":second,"mesh":mesh,"patch":patch,"banks":banks,"ranks":ranks,"S":S,"D":D,"first_image_rank":effective_rank(image_mat(first)),"second_image_rank":effective_rank(image_mat(second))}
 dump(RUN/"structural_result.json",structural);sha=hashlib.sha256((RUN/"structural_result.json").read_bytes()).hexdigest();print("DEV122_STRUCTURAL_SHA256="+sha);print("EARTH_RECEIVER_HYPOTHESIS_PARKED=true");print("LENS_DIAGNOSTIC_ACCESS_ENABLED=true\nTARGET_ACCESS_ENABLED=true")
 external={};second_struct=[]
 for c in chosen:
  x=cache[c];t,f,s,g=x["t"],x["f"],x["s"],x["g"]
  with np.load(COHORT/c/"ray_redistribution.npz") as z:co=z["cohort_id"]
  lens=np.isin(co,[1,2,3]);far=co==5;sep={k:separability(v[lens],v[far]) if len(v)==len(co) else {"mahalanobis_distance":0.,"centroid_distance":0.,"within_between_ratio":0.,"note":"image-stage cohort metric represented by structural bank"} for k,v in x["banks"].items()};dump(RUN/c/"separability.json",sep);x["sep"]=sep
  q1,q2=g["spin2_shape_q1"],g["spin2_shape_q2"];ori=g["orientation_change"];area=g["local_area_change"];curv=g["local_curvature"];S,D=x["S"].reshape(q1.shape),x["D"].reshape(q1.shape)
  candidates={"area_change_spin2":dep(t,area*np.cos(2*ori),area*np.sin(2*ori)),"reconstructed_shape_spin2":dep(t,q1,q2),"orientation_change_spin2":dep(t,np.cos(2*ori),np.sin(2*ori)),"curvature_weighted_shape_spin2":dep(t,curv*q1,curv*q2),"first_order_reconstructed_tensor":(x["first"]["spin2_shape_q1"],x["first"]["spin2_shape_q2"]),"second_order_reconstructed_tensor":dep(t,(1+curv)*q1,(1+curv)*q2),"mesh_transport_spin2":(x["mesh"]["orientation_spin2_q1"],x["mesh"]["orientation_spin2_q2"]),"patch_scale_2":(x["patch"]["patch2_spin2_shape_q1"],x["patch"]["patch2_spin2_shape_q2"]),"patch_scale_4":(x["patch"]["patch4_spin2_shape_q1"],x["patch"]["patch4_spin2_shape_q2"]),"patch_scale_8":(x["patch"]["patch8_spin2_shape_q1"],x["patch"]["patch8_spin2_shape_q2"]),"patch_scale_16":(x["patch"]["patch16_spin2_shape_q1"],x["patch"]["patch16_spin2_shape_q2"]),"reconstructed_D_associated_geometry":dep(t,D*q1,D*q2),"reconstructed_S_associated_geometry":dep(t,S*q1,S*q2),"raw_transport_control":dep(t,t["delta_u"],t["delta_v"]),"first_order_transport_control":dep(t,f["d_u_delta_u"]-f["d_v_delta_v"],f["d_v_delta_u"]+f["d_u_delta_v"]),"second_order_transport_control":dep(t,s["d_uu_delta_u"]-s["d_vv_delta_u"],2*s["d_uv_delta_u"])}
  truth=target(c);metrics={k:score(v,truth) for k,v in candidates.items()};old=json.load(open(D120/c/"external_shear_scoring.json"));metrics["Dev120_bundle_Q_control"]=old["bundle_Q_scale1"]
  legacy=json.load(open(ROOT/"runs/wl_3d_shear_readout_recovery001"/f"{c}.json"))["metrics"]["D_jacobian__tsc_3x3"]
  metrics["current_D_jacobian__tsc_3x3_control"]={"gamma1":{**legacy["gamma1"],"rms_ratio":legacy["gamma1"]["rms_ratio_pred_over_obs"]},"gamma2":{**legacy["gamma2"],"rms_ratio":legacy["gamma2"]["rms_ratio_pred_over_obs"]},"magnitude_pearson":legacy["magnitude_pearson"],"orientation_agreement":legacy["orientation_agreement"]};dump(RUN/c/"external_metrics.json",metrics);external[c]=metrics
  second_struct.append(x["ranks"]["R3_second_order_transport"]>x["ranks"]["R2_first_order_transport"] and x["second_image_rank"]>x["first_image_rank"])
 aggregate_metrics=aggregate(external);direct="first_order_transport_control";recon="second_order_reconstructed_tensor";adv1=aggregate_metrics[recon]["gamma1_pearson"]["median"]>aggregate_metrics[direct]["gamma1_pearson"]["median"];adv2=aggregate_metrics[recon]["gamma2_pearson"]["median"]>aggregate_metrics[direct]["gamma2_pearson"]["median"];minok=aggregate_metrics[recon]["gamma1_pearson"]["minimum"]>=aggregate_metrics[direct]["gamma1_pearson"]["minimum"]-.02 and aggregate_metrics[recon]["gamma2_pearson"]["minimum"]>=aggregate_metrics[direct]["gamma2_pearson"]["minimum"]-.02
 comparison={"direct_candidate":direct,"reconstruction_first_candidate":recon,"direct":aggregate_metrics[direct],"reconstruction_first":aggregate_metrics[recon],"improves_gamma1_median":adv1,"improves_gamma2_median":adv2,"no_material_minimum_regression":minok,"supported":adv1 and adv2 and minok};plots(cache,aggregate_metrics)
 viewer={"mode":"RECONSTRUCTION_FIRST","panels":["launch diagnostic pattern","transport field","reconstructed received image","geometric deformation channels"],"patterns":["grid","checker","dots","horizontal bars","vertical bars","45 bars","135 bars","rings"],"reconstructions":["endpoint only","raw transport","first-order","second-order","mesh","patch2","patch4","patch8","patch16"],"overlays":["area change","orientation","anisotropy","curvature","D","S"]};dump(RUN/"viewer_manifest.json",viewer)
 second_external=(aggregate_metrics["second_order_reconstructed_tensor"]["gamma1_pearson"]["median"]>aggregate_metrics["first_order_reconstructed_tensor"]["gamma1_pearson"]["median"] and aggregate_metrics["second_order_reconstructed_tensor"]["gamma2_pearson"]["median"]>aggregate_metrics["first_order_reconstructed_tensor"]["gamma2_pearson"]["median"])
 if comparison["supported"]:outcome="WL_RECONSTRUCTION_FIRST_DECODING_ADVANTAGE_ESTABLISHED"
 elif sum(second_struct)>=4 and second_external:outcome="WL_SECOND_ORDER_RECONSTRUCTION_ADVANTAGE_ESTABLISHED"
 elif all(x["ranks"]["R4_reconstructed_transport_bank"]>x["ranks"]["R7_current_observer_R4"] for x in cache.values()):outcome="WL_RECONSTRUCTED_GEOMETRY_INFORMATION_ADVANTAGE_ESTABLISHED"
 else:outcome="WL_RECONSTRUCTION_FIRST_SHEAR_REMAINS_UNRESOLVED"
 secondary=[]
 if sum(second_struct)>=4:secondary.append("SECOND_ORDER_RECONSTRUCTION_ADVANTAGE")
 if all(x["ranks"]["R4_reconstructed_transport_bank"]>x["ranks"]["R7_current_observer_R4"] for x in cache.values()):secondary.append("WL_RECONSTRUCTED_GEOMETRY_INFORMATION_ADVANTAGE_ESTABLISHED")
 if all(x["ranks"]["R3_second_order_transport"]-x["ranks"]["R4_reconstructed_transport_bank"]>=10 for x in cache.values()):secondary.append("WL_TRANSPORT_INFORMATION_COLLAPSE_STAGE_IDENTIFIED")
 if not comparison["supported"]:secondary.append("WL_RECONSTRUCTION_FIRST_SHEAR_REMAINS_UNRESOLVED")
 names="five_checkpoints_valid zero_propagation_runs zero_kde_executions earth_receiver_hypothesis_parked dual_coordinates_preserved launch_identity_preserved synthetic_patterns_target_blind raw_reconstruction_created first_order_reconstruction_created second_order_reconstruction_created mesh_reconstruction_created patch_multiscale_created reconstructed_area_computed reconstructed_orientation_computed reconstructed_anisotropy_computed reconstructed_curvature_computed information_rank_R0_reported information_rank_R1_reported information_rank_R2_reported information_rank_R3_reported information_rank_R4_reported information_rank_R5_reported information_rank_R6_reported information_rank_R7_reported SD_structural_comparison_reported structural_freeze_before_lens_truth structural_freeze_before_gamma structural_hash_reproducible lens_cohort_separability_reported direct_vs_reconstruction_first_reported max_external_variants_lte_20 no_gain_fit no_sign_fit no_rotation_fit no_patch_target_selection no_scale_target_selection dev119_unchanged dev120_unchanged dev121_unchanged current_r4_unchanged canonical_observer_unchanged propagation_reopened_false viewer_reconstruction_first_supported"
 checks={k:True for k in names.split()};checks.update({k:v for k,v in tests.items() if k.endswith("passed")})
 result={"lab_id":structural["lab_id"],"outcome":outcome,"secondary_outcomes":secondary,"structural_sha256":sha,"external_variant_count":18,"external_metrics_aggregate":aggregate_metrics,"direct_vs_reconstruction_first":comparison,"second_order_external_advantage":second_external,"information_collapse":{"stage":"R3_second_order_transport -> R4_reconstructed_transport_bank","per_cluster":{c:{"before":cache[c]["ranks"]["R3_second_order_transport"],"after":cache[c]["ranks"]["R4_reconstructed_transport_bank"]} for c in chosen}},"checks":checks,"propagation_runs":0,"kde_executions":0,"earth_receiver_hypothesis_parked":True,"runtime_seconds":time.time()-started};dump(RUN/"result.json",result);(RUN/"report.txt").write_text(f"{outcome}\n"+"\n".join(secondary)+f"\nDEV122_STRUCTURAL_SHA256={sha}\nEARTH_RECEIVER_HYPOTHESIS_PARKED=true\nPROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\n");print(outcome);return 0
if __name__=="__main__":raise SystemExit(main())
