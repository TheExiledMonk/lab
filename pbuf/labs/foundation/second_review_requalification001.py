#!/usr/bin/env python3
"""PBUF SECOND-REVIEW-REQUALIFICATION-001 — CORRECTION 001.

Supersedes the prior R2/R3 implementation. R2 now loads the real MACS0416
Frontier-Fields kappa FITS and runs the reviewed PL1_PM1_PS2 response through
the actual frozen photon propagation + ray-bundle Jacobian. R3 transforms the
real evolved scalar state, reruns the 3D candidate construction under RC0..RC6,
and inverse-transforms the resulting physical vector fields. No fitting or
synthetic substitution is allowed.
"""
from __future__ import annotations
import csv, hashlib, json, math, subprocess, sys, time
from pathlib import Path
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from a8_three_dimensional_projection_lab001 import CLUSTERS, PRODUCTION, construct_common_proxy, construct_rho_3d
from weak_lensing_observation001 import propagate as wl_propagate, resample_to_grid
import observable_lab001 as obs_lab
import source_plane_lab001 as src_lab
from pbuf.core import conventions as M01, coordinate_transforms as M02, vector_transforms as M03
from pbuf.core import pair_enumeration as M05, helmholtz_3d as M13, los_projection as M14
from pbuf.core import ray_interface as M15, observable_extraction as M16, pair_transfer as M08
from pbuf.core import midpoint_rasterization as M10
from pbuf.models import a8_state as M06_state, a8_pair_amplitude as M06, transverse_projector as M07

OUT = ROOT / 'runs' / 'verified_numerical_core_second_review_requalification001'
BENCHMARK = ROOT / 'PBUF_benchmark'
LAB_ID = 'PBUF-SECOND-REVIEW-REQUALIFICATION-001-CORRECTION-001'
CLUSTER_ID, CANDIDATE_ID, NZ, PROFILE, COV_TOL = 'MACS0416', 'PL1_PM1_PS2', 9, 'gaussian', 0.05
CFG = dict(PRODUCTION)


def wjson(name, obj):
    (OUT/name).write_text(json.dumps(obj, indent=2, default=lambda o: float(o) if isinstance(o,np.floating) else int(o) if isinstance(o,np.integer) else bool(o) if isinstance(o,np.bool_) else list(o) if isinstance(o,tuple) else str(o)))

def wcsv(name, rows):
    p=OUT/name
    if not rows: p.write_text(''); return
    fs=[]
    for r in rows:
        for k in r:
            if k not in fs: fs.append(k)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fs); w.writeheader(); w.writerows(rows)

def git(*args): return subprocess.check_output(['git',*args],cwd=str(ROOT),text=True).strip()
def sha_file(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def sha_arr(a): return hashlib.sha256(np.ascontiguousarray(np.asarray(a,dtype=np.float64)).tobytes()).hexdigest()
def vec_sha(v):
    h=hashlib.sha256()
    for a in v: h.update(np.ascontiguousarray(np.asarray(a,dtype=np.float64)).tobytes())
    return h.hexdigest()
def energy(v): return float(np.sum(v[0]**2+v[1]**2+v[2]**2))
def relerr(ref,test):
    d=float(np.sum((test[0]-ref[0])**2+(test[1]-ref[1])**2+(test[2]-ref[2])**2))
    return math.sqrt(d)/max(math.sqrt(energy(ref)),1e-300)
def rms(a): return float(np.sqrt(np.mean(np.asarray(a,dtype=np.float64)**2)))


def repo_state():
    tracked=git('diff','--name-only'); staged=git('diff','--name-only','--cached')
    return {'repository':'TheExiledMonk/lab','branch':git('rev-parse','--abbrev-ref','HEAD'),'head_sha':git('rev-parse','HEAD'),'tracked_changes':tracked,'staged_changes':staged,'working_tree_source_clean':tracked=='' and staged==''}


def load_real_input():
    c=[x for x in CLUSTERS if x['id']==CLUSTER_ID][0]
    p=BENCHMARK/c['directory']/f"hlsp_frontier_model_{c['slug']}_merten_v1_kappa.fits"
    if not p.exists(): raise FileNotFoundError(p)
    with fits.open(p) as h:
        k=np.asarray(h[0].data,dtype=np.float64); hdr=h[0].header
    rho2=construct_common_proxy(k,bins=CFG['bins'],extent=CFG['extent'])
    rho3=construct_rho_3d(rho2,NZ,profile=PROFILE)
    return {'cluster':c,'kappa':k,'rho2':rho2,'rho3':rho3,'prov':{'input_kind':'observed_frontier_fields_fits','cluster_id':CLUSTER_ID,'fits_path':str(p.relative_to(ROOT)),'fits_sha256':sha_file(p),'fits_shape':list(k.shape),'proxy_shape':list(rho2.shape),'rho3d_shape':list(rho3.shape),'proxy_sha256':sha_arr(rho2),'rho3d_sha256':sha_arr(rho3),'Z_L':float(hdr['Z_L']) if 'Z_L' in hdr else None,'Z_S':float(hdr['Z_S']) if 'Z_S' in hdr else None,'nz':NZ,'depth_profile':PROFILE}}


def candidate(state):
    shape=state['c_state'].shape
    pairs=M05.enumerate_internal_pairs(shape)
    ex,ey,ez,valid,g=M07.build_longitudinal_direction(state['c_state'])
    P=M07.build_transverse_projector(ex,ey,ez)
    amp=M06.compute_a8_pair_amplitudes(state['u_slow'],state['u_fast'],state['c_state'],pairs)
    pr=M08.build_pair_responses(pairs,amp,P,'PM1','PS2')
    end=M08.assemble_endpoint_field(pr,shape)
    iface=M10.rasterize_interface_field(pr,shape)
    return {'shape':shape,'pairs':pairs,'amp':amp,'pr':pr,'end':end,'iface':iface,'valid_count':int(np.count_nonzero(valid)),'gradient_rms':rms(g)}


def real_ray(Rx,Ry,real):
    meta={'candidate_id':CANDIDATE_ID,'cluster_id':CLUSTER_ID,'transform_id':'RC0','role':'los','source_artifact_ids':['real_endpoint_field']}
    art=M15.prepare_ray_input(Rx,Ry,meta,require_nontrivial=True)
    if Rx.shape != Ry.shape or Rx.ndim != 2 or Rx.shape[0] != Rx.shape[1]: raise RuntimeError('invalid ray image-plane shape')
    n=Rx.shape[0]; grid=np.linspace(-CFG['extent'],CFG['extent'],n)
    field={'xgrid':grid,'ygrid':grid,'rx':Rx,'ry':Ry}
    x0,y0,vx0,vy0=src_lab.launch_B_cartesian(CFG['nphotons'])
    ph=wl_propagate(field,CFG['step'],CFG['steps'],x0,y0,vx0,vy0); ph['x0']=x0; ph['y0']=y0
    jac=obs_lab.method_jacobian(x0,y0,ph['x'],ph['y'],CFG['extent'],CFG['bins'])
    obs_k=resample_to_grid(real['kappa'],CFG['bins'],CFG['extent'])
    obs=M16.package_lensing_observables(jac['convergence'],jac['shear_g1'],jac['shear_g2'],reference_kappa=obs_k)
    k=np.asarray(obs['kappa']); finite=np.isfinite(k); disp=np.hypot(ph['x']-x0,ph['y']-y0)
    th=hashlib.sha256()
    for nm in ('xs','ys','x','y','conservation'): th.update(np.ascontiguousarray(np.asarray(ph[nm],dtype=np.float64)).tobytes())
    return {'artifact':art,'obs':obs,'metrics':{'n_photons':int(len(x0)),'trajectory_sha256':th.hexdigest(),'mean_endpoint_displacement':float(np.mean(disp)),'max_endpoint_displacement':float(np.max(disp)),'conservation_max':float(np.max(ph['conservation'])),'kappa_finite_count':int(finite.sum()),'kappa_total_count':int(k.size),'kappa_finite_fraction':float(finite.mean()),'kappa_variance_finite':float(np.var(k[finite])) if finite.sum()>=2 else float('nan'),'kappa_rms_finite':rms(k[finite]) if finite.any() else float('nan'),'pearson_vs_observed':obs.get('pearson_vs_reference',float('nan')),'spearman_vs_observed':obs.get('spearman_vs_reference',float('nan'))}}


def run_R2():
    real=load_real_input(); state=M06_state.build_a8_state_3d(real['rho3'],strength=CFG['strength'],seed=12345); cand=candidate(state)
    end=cand['end']; v=(end['Rx_3d'],end['Ry_3d'],end['Rz_3d'])
    los=M14.project_vector_to_image_plane(*v,los_axis='z'); Rx,Ry=los['comp_1'],los['comp_2']
    ray=real_ray(Rx,Ry,real)
    hn=M13.helmholtz_decompose_3d(*v,padding='none'); hp=M13.helmholtz_decompose_3d(*v,padding='reflect_half')
    ee=float(end['statistics']['endpoint_energy']); ie=float(cand['iface']['statistics']['interface_energy']); cl=float(end['statistics']['global_vector_sum_norm']); km=ray['metrics']
    ok=bool(real['prov']['input_kind']=='observed_frontier_fields_fits' and ee>0 and ie>0 and ray['artifact'].statistics['ray_classification'] in ('structured_small','structured_normal','constant_nonzero') and km['kappa_finite_count']>=2 and np.isfinite(km['kappa_variance_finite']) and km['kappa_variance_finite']>0 and km['mean_endpoint_displacement']>0 and np.isfinite(km['conservation_max']))
    metrics={'cluster_id':CLUSTER_ID,'candidate_id':CANDIDATE_ID,'transform_id':'RC0','shape':list(cand['shape']),'n_pairs':int(len(cand['pairs'])),'endpoint_energy':ee,'interface_energy':ie,'endpoint_closure':cl,'los_rx_rms':rms(Rx),'los_ry_rms':rms(Ry),'ray_classification':ray['artifact'].statistics['ray_classification'],**km,'helmholtz_none':{k:hn[k] for k in ('field_reconstruction_error','energy_closure_error','orthogonality_error','f_irr_partition','f_sol_partition','f_irr_native','f_sol_native')},'helmholtz_padded':{k:hp[k] for k in ('field_reconstruction_error','energy_closure_error','orthogonality_error','f_irr_partition','f_sol_partition','f_irr_native','f_sol_native')}}
    lineage={'fits':real['prov']['fits_sha256'],'rho2':real['prov']['proxy_sha256'],'rho3':real['prov']['rho3d_sha256'],'u_slow':sha_arr(state['u_slow']),'u_fast':sha_arr(state['u_fast']),'c_state':sha_arr(state['c_state']),'endpoint':vec_sha(v),'los':vec_sha((Rx,Ry,np.zeros_like(Rx))),'ray_input':ray['artifact'].sha256,'trajectory':km['trajectory_sha256'],'kappa':sha_arr(np.nan_to_num(ray['obs']['kappa'],nan=0.0))}
    return {'passes':ok,'real':real,'state':state,'cand':cand,'metrics':metrics,'lineage':lineage}


def tstate(state,rc): return {k:M02.transform_scalar_field(state[k],rc) for k in ('rho_3d','u_slow','u_fast','c_state')}


def run_R3(r2):
    e0=r2['cand']['end']; i0=r2['cand']['iface']; refE=(e0['Rx_3d'],e0['Ry_3d'],e0['Rz_3d']); refI=(i0['Rx_3d_interface'],i0['Ry_3d_interface'],i0['Rz_3d_interface'])
    rows=[]
    for rc in M01.RC_TRANSFORMS:
        c=candidate(tstate(r2['state'],rc)); e=c['end']; i=c['iface']
        backE=M03.inverse_transform_vector_field(e['Rx_3d'],e['Ry_3d'],e['Rz_3d'],rc); backI=M03.inverse_transform_vector_field(i['Rx_3d_interface'],i['Ry_3d_interface'],i['Rz_3d_interface'],rc)
        wrong=M03.scalar_only_inverse_wrong_control(e['Rx_3d'],e['Ry_3d'],e['Rz_3d'],rc)
        ecE,ecI,ew=relerr(refE,backE),relerr(refI,backI),relerr(refE,wrong)
        ok=(ecE<1e-12 and ecI<1e-12 and ew<1e-12) if rc=='RC0' else (ecE<=COV_TOL and ecI<=COV_TOL and ew>0.3)
        rows.append({'transform':rc,'input_shape':list(c['shape']),'pair_count':len(c['pairs']),'E_cov_endpoint':ecE,'E_cov_interface':ecI,'E_cov_wrong_scalar_only':ew,'passes':bool(ok)})
    return {'rows':rows,'passes':all(r['passes'] for r in rows),'audit_kind':'real_MACS0416_state_full_3d_candidate_construction','covariance_tolerance':COV_TOL}


def invalidation():
    return {'prior_lab_id':'PBUF-SECOND-REVIEW-REQUALIFICATION-001','prior_merge_pr':8,'prior_outcome_A_valid':False,'module_level_results_retained':True,'R1_synthetic_results_retained':True,'R2_prior_result_invalid':True,'R2_reason':'synthetic density fixture labelled MACS0416; no cluster FITS loaded','R2_observable_prior_result_invalid':True,'R2_observable_reason':'kappa/gamma built algebraically from LOS response instead of photon propagation plus ray-bundle Jacobian','R3_prior_result_invalid':True,'R3_reason':'unrelated analytic vector field round-trip only; candidate path not rerun','full_candidate_rerun_allowed_before_corrected_run':False}


def main():
    t=time.perf_counter(); OUT.mkdir(parents=True,exist_ok=True); wjson('prior_outcome_invalidation.json',invalidation()); repo=repo_state(); wjson('repository_state.json',repo)
    if repo['branch']!='main' or not repo['working_tree_source_clean']:
        v={'lab_id':LAB_ID,'outcome':'Outcome B — REPOSITORY GATE FAILURE','second_review_status':'blocked','full_candidate_rerun_allowed':False}; wjson('validation.json',v); print(json.dumps(v,indent=2)); return 2
    print('[R2] real MACS0416 FITS -> reviewed pair response -> actual photons -> Jacobian')
    r2=run_R2(); wjson('real_cluster_input_provenance.json',r2['real']['prov']); wjson('field_lineage.json',r2['lineage']); wjson('restricted_recovery_statistics.json',r2['metrics'])
    if r2['passes']:
        print('[R3] transformed real evolved state -> rerun candidate RC0..RC6')
        r3=run_R3(r2); wcsv('covariance_revalidation.csv',r3['rows'])
    else:
        r3={'rows':[],'passes':False,'audit_kind':'not_run_due_to_R2_failure','covariance_tolerance':COV_TOL}; wcsv('covariance_revalidation.csv',[])
    outcome='Outcome A — CORE REQUALIFIED' if r2['passes'] and r3['passes'] else ('Outcome E — REAL-CANDIDATE COVARIANCE FAILURE' if r2['passes'] else 'Outcome D — REAL-CLUSTER RECOVERY FAILURE')
    allowed=outcome.startswith('Outcome A')
    v={'lab_id':LAB_ID,'outcome':outcome,'head_sha':repo['head_sha'],'module_core_status_from_prior_requalification':'retained; no core source module changed by this correction','R2_real_cluster_real_ray_pass':bool(r2['passes']),'R3_real_candidate_covariance_pass':bool(r3['passes']),'second_review_status':'accepted' if allowed else 'blocked','full_candidate_rerun_allowed':allowed,'next_permitted_experiment':'PBUF 3D PAIRWISE PRIMARY-CANDIDATE SCIENCE RE-RUN 001' if allowed else None,'duration_seconds':time.perf_counter()-t}; wjson('validation.json',v); wjson('run.json',{'lab_id':LAB_ID,'head_sha':repo['head_sha'],'cluster_id':CLUSTER_ID,'candidate_id':CANDIDATE_ID,'nz':NZ,'profile':PROFILE,'config':CFG,'duration_seconds':v['duration_seconds']}); print(json.dumps(v,indent=2)); return 0 if allowed else 1

if __name__=='__main__': raise SystemExit(main())
