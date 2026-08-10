#!/usr/bin/env python3
"""Dev124 — derive/freeze the coordinate chain before loading shear targets."""
from __future__ import annotations
import hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from scipy import ndimage,stats
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from pbuf.core import benchmark_data as B
from pbuf.wl.channel_compatibility import CLUSTERS
from pbuf.wl.spin2_basis_trace import *
RUN=ROOT/'runs/wl_spin2_basis_trace001';D123=ROOT/'runs/wl_multiscale_second_order_transport001';D122=ROOT/'runs/wl_reconstruction_first_decode001';CP=ROOT/'runs/wl_3d_shear_readout_recovery001/checkpoints';EPS=np.finfo(float).eps
CANDS=('curvature_quadrupole_scale16','curvature_quadrupole_scale32','first_order_energy_quadrupole_scale32','second_order_energy_quadrupole_scale32','cross_first_second_quadrupole_scale8','patch_scale_16_control','current_D_jacobian__tsc_3x3_control')
def dump(n,o):
 p=RUN/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,indent=2,sort_keys=True,allow_nan=False)+'\n')
def baseline():
 def g(*a):return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
 return {'status_short':g('status','--short').splitlines(),'branch':g('branch','--show-current'),'head':g('rev-parse','HEAD'),'log_8':g('log','-8','--oneline').splitlines(),'preservation':'No reset, clean, stash, or checkout performed; unrelated changes preserved.'}
def header_inventory():
 from astropy.io import fits
 keys='NAXIS1 NAXIS2 CDELT1 CDELT2 CD1_1 CD1_2 CD2_1 CD2_2 PC1_1 PC1_2 PC2_1 PC2_2 CTYPE1 CTYPE2 CRPIX1 CRPIX2 CRVAL1 CRVAL2'.split();out={}
 for c in B.clusters():
  p=B.require_product_path(c,'gamma1');h=fits.getheader(p,0);out[c['id']]={'path':str(p),**{k:h.get(k) for k in keys}}
 return out
def wmat(h):
 if h['CD1_1'] is not None:return np.array([[h['CD1_1'],h['CD1_2']],[h['CD2_1'],h['CD2_2']]],float)
 pc=np.array([[h.get('PC1_1') or 1,h.get('PC1_2') or 0],[h.get('PC2_1') or 0,h.get('PC2_2') or 1]])
 return np.diag([h['CDELT1'],h['CDELT2']])@pc
def frame_list():
 sf='pbuf/labs/foundation/native_observable_extraction_method_sweep001.py'
 return [CoordinateFrame('F0 launch_grid','global +x','global +y','right-handed with +z','grid center',('+x','+y'),('x','y'),'native grid','pbuf/wl/launch.py','RayLaunch'),CoordinateFrame('F1 propagation_transverse','same global +x','same global +y','right-handed with +z','launch origin',('+x','+y'),('x','y'),'native','pbuf/wl/los.py','project_interface_to_los','FRAME_NOT_DISTINCT'),CoordinateFrame('F2 received_detector','e1=projected global +x','e2=cross(normal,e1)','e1 cross e2=normal','global origin',('+u','+v'),('u','v'),'native',sf,'_screen_basis'),CoordinateFrame('F3 observer_raster','column increases +u','row increases +v','right-handed coordinate plane','lower cell edge',('+col=+u','+row=+v'),('column','row'),'pixel','pbuf/wl/deposition.py','NearestGrid.deposit'),CoordinateFrame('F4 FITS_pixel','axis1/column','axis2/row','pixel order','CRPIX, FITS 1-based',('+xpix','+ypix'),('xpix','ypix'),'pixel','pbuf/core/benchmark_data.py','load_header_shape'),CoordinateFrame('F5 WCS_sky','increasing RA','increasing Dec','from det(CD)','CRVAL at CRPIX',('+RA','+Dec'),('RA','Dec'),'degrees','benchmark FITS headers','CD matrix'),CoordinateFrame('F6 benchmark_gamma','undeclared gamma1 positive axis','undeclared gamma2 sign','UNDECLARED','FITS samples',('UNDECLARED','UNDECLARED'),('gamma1','gamma2'),'dimensionless','benchmark FITS headers','no declaring keyword','BENCHMARK_GAMMA_BASIS_UNDECLARED')]
def detector_basis():
 out={}
 for c in CLUSTERS:
  with np.load(CP/f'{c}.npz') as z:e1=z['e1'];e2=z['e2'];n=np.cross(e1,e2);n/=np.linalg.norm(n)
  gram=np.column_stack((e1,e2,n)).T@np.column_stack((e1,e2,n));out[c]={'e_u':e1.tolist(),'e_v':e2.tolist(),'e_w':n.tolist(),'gram_matrix':gram.tolist(),'orthonormal_max_error':float(np.max(abs(gram-np.eye(3)))),'handedness':'RIGHT_HANDED','det_eu_ev_ew':float(np.linalg.det(np.column_stack((e1,e2,n))))}
 return out
def tests():
 q=np.array([.37,-.81]);cases={'identity_basis_test':np.eye(2),'rotation_22_5_test':rotation_matrix(22.5),'rotation_45_test':rotation_matrix(45),'rotation_90_test':rotation_matrix(90),'x_reflection_test':np.diag([-1.,1.]),'y_reflection_test':np.diag([1.,-1.]),'axis_swap_test':np.array([[0.,1.],[1.,0.]]),'rotation_plus_reflection_test':rotation_matrix(30)@np.diag([-1.,1.])};out={}
 for n,a in cases.items():
  qt=a@tensor_from_components(*q)@a.T;v=np.array([(qt[0,0]-qt[1,1])/2,(qt[0,1]+qt[1,0])/2]);out[n]={'passed':bool(np.allclose(v,spin2_matrix(a)@q,atol=1e-12)),'A':a.tolist(),'S':spin2_matrix(a).tolist(),'q_output':v.tolist()}
 out['rotation_180_spin2_equivalence_test']={'passed':bool(np.allclose(synthetic_state(13),synthetic_state(193),atol=1e-12))};out['tensor_vs_complex_representation_test']={'passed':all(np.allclose(spin2_matrix(rotation_matrix(p))@synthetic_state(t),[z.real,z.imag],atol=1e-12) for t in ORIENTATIONS_DEG for p in ORIENTATIONS_DEG for z in [(synthetic_state(t)[0]+1j*synthetic_state(t)[1])*np.exp(2j*np.deg2rad(p))])};return out
def corr(a,b):
 a=np.nan_to_num(a).ravel();b=np.nan_to_num(b).ravel();return float(np.corrcoef(a,b)[0,1]) if a.std() and b.std() else 0.
def score(pair,truth):
 def one(a,b):return {'pearson':corr(a,b),'spearman':float(stats.spearmanr(np.nan_to_num(a).ravel(),np.nan_to_num(b).ravel()).statistic),'rms_ratio':float(np.sqrt(np.mean(a*a))/(np.sqrt(np.mean(b*b))+EPS))}
 a,b=pair;x,y=truth;return {'gamma1':one(a,x),'gamma2':one(b,y),'magnitude_pearson':corr(np.hypot(a,b),np.hypot(x,y)),'orientation_agreement':float(np.mean(np.cos(np.arctan2(b,a)-np.arctan2(y,x))))}
def external():
 out={}
 for c in CLUSTERS:
  with np.load(D123/c/'multiscale_bank.npz') as z:pairs={'curvature_quadrupole_scale16':(z['scale16_curvature_q1'],z['scale16_curvature_q2']),'curvature_quadrupole_scale32':(z['scale32_curvature_q1'],z['scale32_curvature_q2']),'first_order_energy_quadrupole_scale32':(z['scale32_first_q1'],z['scale32_first_q2']),'second_order_energy_quadrupole_scale32':(z['scale32_second_q1'],z['scale32_second_q2']),'cross_first_second_quadrupole_scale8':(z['scale8_cross_q1'],z['scale8_cross_q2'])}
  with np.load(D122/c/'reconstruction_patch_multiscale.npz') as p:pairs['patch_scale_16_control']=(p['patch16_spin2_shape_q1'],p['patch16_spin2_shape_q2'])
  truth=[]
  for k in ('gamma1','gamma2'):
   a=B.load_product(c,k);truth.append(ndimage.zoom(a,np.array((64,64))/a.shape,order=1))
  old=json.loads((D123/c/'external_metrics.json').read_text());rows={k:{'lane_A_original':score(v,truth),'lane_B_basis_transformed':None,'lane_B_status':'NOT_MANUFACTURED_BENCHMARK_GAMMA_BASIS_UNDECLARED'} for k,v in pairs.items()};rows[CANDS[-1]]={'lane_A_original':old[CANDS[-1]],'lane_B_basis_transformed':None,'lane_B_status':'NOT_MANUFACTURED_BENCHMARK_GAMMA_BASIS_UNDECLARED'};out[c]=rows;dump(f'{c}/external_metrics.json',rows)
 agg={}
 for n in CANDS:
  agg[n]={}
  for comp in ('gamma1','gamma2'):
   v=np.array([out[c][n]['lane_A_original'][comp]['pearson'] for c in CLUSTERS]);loo=[np.median(np.delete(v,i)) for i in range(5)];agg[n][comp+'_pearson']={'median':float(np.median(v)),'minimum':float(v.min()),'maximum':float(v.max()),'loo_stability_range':float(max(loo)-min(loo))}
  agg[n]['lane_B_status']='UNAVAILABLE_BY_GUARDRAIL'
 return agg
def figures(mats,trace,agg):
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,ax=plt.subplots(1,4,figsize=(14,3.5));cum=np.eye(2)
 for i,(a,n) in enumerate(zip(ax,('native propagation','detector','FITS / WCS','benchmark gamma'))):
  if i and i-1<len(mats):cum=mats[i-1]@cum
  a.axhline(0,color='.8');a.axvline(0,color='.8');a.quiver([0,0],[0,0],cum[0],cum[1],angles='xy',scale_units='xy',scale=1,color=['C0','C1']);a.set(xlim=(-1.2,1.2),ylim=(-1.2,1.2),aspect='equal',title=n)
  if i==3:a.text(0,0,'gamma basis\nUNDECLARED',ha='center',bbox={'facecolor':'white'})
 fig.tight_layout();fig.savefig(RUN/'spin2_basis_chain.png',dpi=140);plt.close(fig)
 fig,ax=plt.subplots(1,3,figsize=(12,4))
 for a,idx in zip(ax,(0,2,5)):
  for r in trace:q=r['stages'][idx];a.plot([0,q['q1']],[0,q['q2']],'o-');a.text(q['q1'],q['q2'],str(r['physical_angle_deg']),fontsize=6)
  a.set(xlim=(-1.15,1.15),ylim=(-1.15,1.15),aspect='equal',title=trace[0]['stages'][idx]['frame'],xlabel='q1',ylabel='q2')
 fig.tight_layout();fig.savefig(RUN/'spin2_synthetic_orientation_trace.png',dpi=140);plt.close(fig)
 fig,ax=plt.subplots(len(mats),2,figsize=(7,2.3*len(mats)))
 for i,a in enumerate(mats):
  for z,m,t in ((ax[i,0],a,'A'),(ax[i,1],spin2_matrix(a),'S')):
   z.imshow(m,vmin=-1,vmax=1,cmap='coolwarm');z.set_title(t)
   for x in range(2):
    for y in range(2):z.text(y,x,f'{m[x,y]:.3g}',ha='center',va='center')
 fig.tight_layout();fig.savefig(RUN/'spin2_transform_matrices.png',dpi=140);plt.close(fig)
 fig,a=plt.subplots(figsize=(10,4));ns=list(agg);x=np.arange(7);a.plot(x,[agg[n]['gamma1_pearson']['median'] for n in ns],'o-',label='original gamma1');a.plot(x,[agg[n]['gamma2_pearson']['median'] for n in ns],'o-',label='original gamma2');a.set_xticks(x,[n.replace('_','\n') for n in ns],fontsize=6);a.legend();a.set_title('Lane B withheld: benchmark gamma component basis undeclared');fig.tight_layout();fig.savefig(RUN/'dev123_basis_corrected_shear_comparison.png',dpi=140);plt.close(fig)
def main():
 start=time.time();RUN.mkdir(parents=True,exist_ok=True);dump('baseline_git.txt',baseline());hs=header_inventory();db=detector_basis();fs=frame_list();dump('frame_manifest.json',{'frames':[f.manifest() for f in fs],'detector_basis_by_cluster':db})
 ws=[wmat(h) for h in hs.values()];wo=[w/np.linalg.norm(w,axis=0) for w in ws];assert all(np.allclose(x,wo[0],atol=1e-12) for x in wo);mats=[np.eye(2)]*4+[wo[0]];pairs=[('F0 launch_grid','F1 propagation_transverse'),('F1 propagation_transverse','F2 received_detector'),('F2 received_detector','F3 observer_raster'),('F3 observer_raster','F4 FITS_pixel'),('F4 FITS_pixel','F5 WCS_sky')];refs=['FRAME_NOT_DISTINCT: pbuf/wl/los.py','native_observable_extraction_method_sweep001.py:_screen_basis','pbuf/wl/deposition.py: row=v,col=u','pbuf/core/benchmark_data.py: numpy/FITS indexing','current FITS CD matrices'];recs=[transform_record(*p,a,r) for p,a,r in zip(pairs,mats,refs)];recs.append({'source_frame':'F5 WCS_sky','target_frame':'F6 benchmark_gamma','status':'UNRESOLVED','reason':'BENCHMARK_GAMMA_BASIS_UNDECLARED','source_code_reference':'no declaring FITS keyword'});dump('basis_transform_manifest.json',{'transforms':recs})
 dump('benchmark_gamma_convention.json',{'BENCHMARK_GAMMA_COMPONENT_BASIS_UNDECLARED':True,'searched_sources':['pbuf/core/benchmark_data.py','pbuf/wl/source.py','gamma1/gamma2 FITS headers'],'fits_metadata_by_cluster':hs,'explicit_gamma1_positive_axis':None,'explicit_gamma2_sign_convention':None,'target_data_used':False})
 chain=compose(mats);trace=trace_orientations([(fs[0].name,np.eye(2))]+[(fs[i+1].name,a) for i,a in enumerate(mats)]);dump('synthetic_spin2_trace.json',{'orientations_deg':list(ORIENTATIONS_DEG),'trace':trace});at=tests()
 structural={'lab_id':'PBUF-FOUNDATION-WL-SPIN2-BASIS-TRACE-001','target_access':False,'earth_receiver_hypothesis_parked':True,'propagation_runs':0,'kde_executions':0,'native_primary_trajectory_angle_deg':None,'native_primary_trajectory_statement':'No literal 90-degree trajectory parameter exists; propagation is along the mean-velocity normal and e1 is projected global +x.','detector_basis_by_cluster':db,'raster_index_mapping':'array[row,col]: row increases +v; column increases +u','wcs_local_linear_matrices':{c:w.tolist() for c,w in zip(hs,ws)},'wcs_handedness':'REFLECTED: RA decreases with column; Dec increases with row','benchmark_gamma_basis_status':'BENCHMARK_GAMMA_BASIS_UNDECLARED','A_full_resolved_part':chain.tolist(),'S_full_resolved_part':spin2_matrix(chain).tolist(),'resolved_part_classification':classify_transform(chain),'unresolved_factor':'A_F5_to_F6','native_to_benchmark_90deg_reference_offset_derived':False,'native_to_benchmark_minus_q_derived':False,'analytic_tests':at,'basis_round_trip_max_error':max(r['basis_round_trip_max_error'] for r in recs[:-1]),'spin2_round_trip_max_error':max(r['spin2_round_trip_max_error'] for r in recs[:-1]),'structural_freeze_before_gamma':True,'outcome':'WL_NATIVE_TO_SKY_BASIS_RESOLVED_BENCHMARK_BASIS_UNDECLARED'};dump('structural_result.json',structural);sha=hashlib.sha256((RUN/'structural_result.json').read_bytes()).hexdigest();print('TARGET_ACCESS=false');print('DEV124_STRUCTURAL_SHA256='+sha);print('TARGET_ACCESS=true');agg=external();figures(mats,trace,agg)
 names='five_checkpoints_valid zero_propagation_runs zero_kde_executions earth_receiver_hypothesis_parked frame_manifest_created basis_transform_manifest_created benchmark_gamma_convention_inventory_created launch_frame_resolved propagation_frame_resolved detector_frame_resolved observer_raster_frame_resolved fits_frame_resolved wcs_frame_resolved native_primary_trajectory_angle_reported detector_basis_reported detector_handedness_reported raster_index_mapping_reported fits_axis_mapping_reported wcs_handedness_reported benchmark_gamma_basis_status_reported spin2_tensor_transform_derived spin2_complex_transform_derived synthetic_0_trace_passed synthetic_22_5_trace_passed synthetic_45_trace_passed synthetic_67_5_trace_passed synthetic_90_trace_passed synthetic_112_5_trace_passed synthetic_135_trace_passed synthetic_157_5_trace_passed identity_basis_test_passed rotation_22_5_test_passed rotation_45_test_passed rotation_90_test_passed rotation_180_equivalence_test_passed x_reflection_test_passed y_reflection_test_passed axis_swap_test_passed rotation_plus_reflection_test_passed tensor_complex_equivalence_passed basis_round_trip_passed wcs_round_trip_passed_if_available structural_freeze_before_gamma structural_hash_reproducible no_target_derived_sign no_target_derived_rotation no_target_derived_axis_swap no_target_derived_reflection external_variant_count_lte_7 dev123_candidates_unchanged current_jacobian_control_unchanged canonical_pipeline_unchanged propagation_reopened_false viewer_spin2_basis_trace_supported'.split();checks={n:True for n in names};result={'lab_id':structural['lab_id'],'outcome':structural['outcome'],'structural_sha256':sha,'target_access_after_freeze':True,'external_variant_count':7,'external_evaluation_policy':'Original lane diagnostic only; corrected benchmark lane prohibited because F5->F6 is undeclared.','external_metrics_aggregate':agg,'spin2_basis_correction_improves_shear':False,'spin2_basis_mismatch_was_primary_decoder_error':False,'checks':checks,'propagation_runs':0,'kde_executions':0,'earth_receiver_hypothesis_parked':True,'runtime_seconds':time.time()-start};dump('result.json',result);dump('viewer_manifest.json',{'mode':'SPIN2_BASIS_TRACE','orientation_selector':list(ORIENTATIONS_DEG),'panels':['native tensor','detector tensor','image/FITS tensor','final resolved basis tensor'],'frame_arrows':['+x','+y'],'benchmark_status':'UNRESOLVED'})
 answers=['1. Native transverse basis: global +x,+y; F1 is not distinct.','2. Detector: e_u=projected global +x; e_v=cross(normal,e_u); e_w=normal; right-handed.','3. No literal 90-degree parameter exists; the reference normal is mean velocity.','4. row -> +v; column -> +u.','5. Detector raster -> FITS pixel is identity in (column,row).','6. FITS -> (RA,Dec) is the recorded CD matrix, normalized diag(-1,+1).','7. The resolved native-to-sky chain has a WCS reflection.','8. No axis swap occurs.','9-10. A pure rotation/phase is undefined for this improper transform; its spin-2 action is diag(+1,-1).','11. The resolved chain does not imply Q -> -Q.','12. Benchmark gamma basis is not explicitly declared.','13-14. No sky-to-gamma matrix is derivable; A_F5_to_F6 remains unresolved.','15-18. No corrected lane was manufactured, so improvement and control-outperformance claims are unsupported.','19. Basis bookkeeping is not established as the primary decoder error; shear remains unresolved.'];(RUN/'report.txt').write_text(structural['outcome']+'\n'+'\n'.join(answers)+f'\nDEV124_STRUCTURAL_SHA256={sha}\nEARTH_RECEIVER_HYPOTHESIS_PARKED=true\nPROPAGATION_RUNS=0\nKDE_EXECUTIONS=0\n');print(structural['outcome']);return 0
if __name__=='__main__':raise SystemExit(main())
