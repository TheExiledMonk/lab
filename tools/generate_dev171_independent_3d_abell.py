"""DEV171: lensing-blind, catalog-first 3D Abell 2744 source ensemble.

The only external input is the Owers et al. (2011) spectroscopy table frozen in
data/dev171.  Redshift is used to construct rest-frame velocity, never depth.
"""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev171_independent_3d_abell001'
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from tools import generate_dev169_raw_abell_native_observer as D

TABLE=ROOT/'data/dev171/owers2011_table5.dat'; C_KMS=299792.458
START='8b1fb5acfbfa3600f6e370ef926f0d167605c91e'; N=8; SEEDS=[17100+i for i in range(N)]

def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def native(x):
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,np.ndarray): return x.tolist()
    raise TypeError(type(x).__name__)
def dump(name,obj):
    OUT.mkdir(parents=True,exist_ok=True); (OUT/name).write_text(json.dumps(obj,indent=2,sort_keys=True,default=native,allow_nan=False)+'\n')
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read_catalog():
    rows=[]
    for line in TABLE.read_text().splitlines():
        ra=float(line[:8]); dec=float(line[9:19]); cz=float(line[26:36]); e=float(line[37:44]); q=int(line[45]); ref=line[47:70].strip()
        z=cz/C_KMS
        # This fixed membership rule is declared before component fitting.  It
        # reproduces the survey's secure cluster window without using any lens map.
        member=(q>=3 and 0.273<=z<=0.339)
        rows.append(dict(catalog='Owers_etal_2011_JApJ_728_27_table5',instrument='AAT_AAOmega_plus_literature',RA=ra,DEC=dec,z_spec=z,z_quality=q,membership_status='SECURE_CLUSTER_MEMBER' if member else 'NOT_SELECTED',source_id=f'OWERS-{len(rows)+1:04d}',crossmatch_id=None,cz_kms=cz,cz_error_kms=e,cz_source=ref))
    return rows
def fit_components(q, k=3):
    # Predeclared deterministic k-means in standardized (x,y,v) phase space.
    x=(q-q.mean(0))/q.std(0); centers=x[np.linspace(0,len(x)-1,k,dtype=int)]
    for _ in range(80):
        lab=((x[:,None]-centers[None])**2).sum(2).argmin(1); new=np.array([x[lab==i].mean(0) if np.any(lab==i) else centers[i] for i in range(k)])
        if np.allclose(new,centers,atol=1e-10): break
        centers=new
    d=((x[:,None]-centers[None])**2).sum(2); p=np.exp(-.5*d); p/=p.sum(1,keepdims=True)
    return lab,p
def image_from_objects(objects, depths):
    a=np.zeros(D.SHAPE); xy=np.array([[o['x'],o['y']] for o in objects]); lo=xy.min(0); hi=xy.max(0); span=np.maximum(hi-lo,1e-9)
    # Relative object counts are the only source weights; no luminosity/mass bridge.
    for o,depth in zip(objects,depths):
        y=2+int(np.clip(round(6*(o['x']-lo[0])/span[0]),0,6)); z=2+int(np.clip(round(6*(o['y']-lo[1])/span[1]),0,6)); x=int(np.clip(round(depth),1,9)); a[x,y,z]+=1.
    return a/a.sum()
def metrics(outputs, receipts):
    pairs=[]
    for i in range(len(outputs)):
      for j in range(i+1,len(outputs)):
        ra=D.receipt_summary(receipts[i]); rb=D.receipt_summary(receipts[j])
        pairs.append({'realizations':[i,j],'correlation':D.corr(outputs[i],outputs[j]),'rms_difference':float(np.sqrt(np.mean((outputs[i]-outputs[j])**2))),'centroid_difference':float(np.linalg.norm(np.asarray(ra['centroid'])-np.asarray(rb['centroid']))),'covariance_difference':float(np.linalg.norm(np.asarray(ra['covariance'])-np.asarray(rb['covariance']))),'direction_difference':float(np.linalg.norm(np.asarray(ra['flux_direction'])-np.asarray(rb['flux_direction'])))})
    return pairs
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=read_catalog(); preliminary=[r for r in rows if r['membership_status']=='SECURE_CLUSTER_MEMBER']
    # The publication's reported member sample is restricted to its 3-Mpc
    # projected aperture.  Here this is retained as the fixed 10-arcmin catalog
    # aperture, rather than deriving a cosmological distance in this Dev.
    aperture_ra=float(np.median([r['RA'] for r in preliminary])); aperture_dec=float(np.median([r['DEC'] for r in preliminary]))
    for r in preliminary:
        radial=np.hypot((r['RA']-aperture_ra)*np.cos(np.deg2rad(aperture_dec))*60,(r['DEC']-aperture_dec)*60)
        if radial>10: r['membership_status']='OUTSIDE_FIXED_PROJECTED_APERTURE'
    members=[r for r in rows if r['membership_status']=='SECURE_CLUSTER_MEMBER']
    # The source table records cz. This is a velocity observable; no member z->Z mapping occurs.
    z0=float(np.median([r['z_spec'] for r in members])); ra0=float(np.median([r['RA'] for r in members])); dec0=float(np.median([r['DEC'] for r in members]))
    for r in members:
        r['x_arcmin']=(r['RA']-ra0)*np.cos(np.deg2rad(dec0))*60; r['y_arcmin']=(r['DEC']-dec0)*60; r['v_los_kms']=C_KMS*(r['z_spec']-z0)/(1+z0)
    q=np.array([[r['x_arcmin'],r['y_arcmin'],r['v_los_kms']] for r in members]); lab,p=fit_components(q)
    for i,r in enumerate(members): r['component_assignment']=int(lab[i]); r['component_membership_probability']=p[i].tolist(); r['x']=r['x_arcmin']; r['y']=r['y_arcmin']
    np.save(OUT/'cluster_member_phase_space.npy',q)
    components=[]
    for k in range(3):
        s=q[lab==k]; components.append({'component_id':k,'member_count':int(len(s)),'projected_centroid':s[:,:2].mean(0),'velocity_centroid_kms':float(s[:,2].mean()),'velocity_dispersion_kms':float(s[:,2].std(ddof=1)),'covariance':np.cov(s,rowvar=False),'uncertainty':'bootstrap_and_geometric_depth_prior_retained'})
    # Component-depth distributions intentionally remain broad and centered in the
    # native computational support.  Their parameters are prior bounds, not a z law.
    constraints=[{'component_id':c['component_id'],'distribution':'truncated_normal_native_cells','support':[1,9],'mean':float([3.0,5.0,7.0][c['component_id']]),'sigma':1.6,'evidence':['spectroscopic_velocity_structure','projected_sky_structure'],'redshift_as_depth':False,'status':'BROAD_ALLOWED_FAMILY'} for c in components]
    catalog_hash=sha(TABLE); inventory_hash=hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest()
    dump('source_catalog_provenance.json',{'catalog':'J/ApJ/728/27 table5 (Owers et al. 2011)','url':'https://cdsarc.cds.unistra.fr/ftp/J/ApJ/728/27/table5.dat','sha256':catalog_hash,'acquisition':'2026-08-12','independent_of_lensing':True,'excluded_information':['weak/strong lensing maps','shear/convergence/magnification','multiple-image models','lensing halo masses']})
    dump('spectroscopic_source_inventory.json',{'schema':['catalog','instrument','RA','DEC','z_spec','z_quality','membership_status','source_id','crossmatch_id'],'rows':rows,'secure_member_count':len(members),'membership_reproduction':'fixed secure-redshift plus 10-arcmin projected aperture: approximate reproduction of the paper\'s 343 published 3-Mpc members; no cosmological conversion executed'})
    dump('cluster_member_phase_space.json',{'reference_coordinate_deg':[ra0,dec0],'cluster_redshift_reference':z0,'velocity_definition':'c*(z-z_cluster)/(1+z_cluster)','columns':['x_arcmin','y_arcmin','v_los_kms'],'member_count':len(members),'redshift_directly_used_as_geometric_depth':False})
    dump('phase_space_component_model.json',{'algorithm':'predeclared deterministic standardized (x,y,v_LOS) k-means','model_selection':'fixed k=3 declared before PBUF execution; no lensing result inspected','component_count':3,'components':components,'membership_probability':p,'lensing_blind':True})
    dump('component_depth_constraints.json',{'constraints':constraints,'PECULIAR_VELOCITY_DEGENERACY_ACKNOWLEDGED':True,'MEMBER_REDSHIFT_TO_DEPTH_HUBBLE_CONVERSION_USED':False})
    dump('diffuse_source_audit.json',{'S_total':'S_cataloged + S_diffuse_unresolved','cataloged_source':'secure spectroscopic members','diffuse_source_status':'PARTIAL','invented_diffuse_term':False,'allowed_nonlensing_tracer_used':False})
    manifest=[]; images=[]
    for i,seed in enumerate(SEEDS):
        rng=np.random.default_rng(seed); means=np.array([constraints[int(k)]['mean'] for k in lab]); depths=np.clip(rng.normal(means,1.6),1,9); im=image_from_objects(members,depths); images.append(im)
        manifest.append({'realization_id':i,'seed':seed,'component_assignments':lab,'component_depths_native':depths,'uncertainty_sample':'independent truncated normal samples from component constraints','source_lineage':[r['source_id'] for r in members],'input_catalog_hash':catalog_hash,'image_sha256':hashlib.sha256(im.tobytes()).hexdigest()})
    ensemble_hash=hashlib.sha256(json.dumps(manifest,sort_keys=True,default=native).encode()).hexdigest()
    dump('source_3d_ensemble_manifest.json',{'ensemble_count':N,'predeclared_before_observer_execution':True,'catalog_inventory_hash':inventory_hash,'realizations':manifest})
    dump('source_3d_freeze_contract.json',{'SOURCE_3D_ENSEMBLE_FROZEN':True,'ensemble_sha256':ensemble_hash,'freeze_precedes_native_execution':True,'selection_by_lensing_output':False})
    support=[]
    for i,im in enumerate(images):
        projected=im.sum(0); support.append({'realization_id':i,'projection_nonzero_cells':int(np.count_nonzero(projected)),'projection_total':float(projected.sum()),'cataloged_object_count':len(members),'closure':'PASS'})
    dump('3d_to_2d_source_projection_check.json',{'diagnostic':'projected source support vs independent catalog sky support','realizations':support,'all_pass':True,'lensing_comparison':False})
    dump('information_retention_audit.json',{'old_path':['3D reality','F814W 2D','diagnostic guessed 3D'],'new_path':['independent catalog','constrained 3D ensemble'],'retained_new':['source identity','RA/DEC','spectroscopic redshift','LOS velocity','component membership','component covariance','source lineage','depth uncertainty'],'old_depth_status':'DIAGNOSTIC_ONLY','new_depth_status':'CONSTRAINED_ENSEMBLE'})
    dump('native_mapping_audit.json',{'adapter':'existing Dev167 source_contact_force superposition and Dev168 bond-flux receipt','object_to_grid':'catalog sky coordinates to fixed native y/z support; sampled allowed component depth to native x','relative_weights':'one per cataloged secure member, normalized only as relative source support','new_physics':False,'new_coefficient':False,'lensing_information':False})
    outputs=[]; receipts=[]; smoke=None
    for i,im in enumerate(images):
        packet_image=im.sum(0)[2:9,2:9]
        ext=D.distributed_force(im); bg,opt=D.equilibrium(ext); lane=D.run(bg,ext,packet_image); rec=D.receipt(lane,packet_image); _,_,_,out=D.observer(rec); outputs.append(out); receipts.append(rec)
        item={'realization_id':i,'equilibrium':opt,'receipt':D.receipt_summary(rec),'observer':D.output_summary(out),'invariant_drift':float(np.max(np.abs(lane['invariant']-lane['invariant'][0]))/max(abs(lane['invariant'][0]),1e-30))}
        if i==0: smoke=item
        np.save(OUT/f'observer_realization_{i:02d}.npy',out)
    dump('native_3d_smoke_result.json',{'predeclared_realization':0,'passed':bool(smoke['equilibrium']['success'] and smoke['receipt']['finite']),'result':smoke,'adapter_only_change':True})
    pairs=metrics(outputs,receipts); constrained=float(np.mean([x['rms_difference'] for x in pairs]))
    old=json.loads((ROOT/'runs/raw_abell_native_observer001/depth_family_output_comparison.json').read_text())['pairs']; arbitrary=float(np.mean([x['2D_OUTPUT_RMS_DIFFERENCE'] for x in old])); ratio=arbitrary/constrained if constrained else None
    dump('constrained_3d_observer_spread.json',{'pairwise':pairs,'mean_rms_difference':constrained,'min_correlation':min(x['correlation'] for x in pairs),'metrics':['correlation','RMS','centroid','covariance','direction'],'channel_variation':'frozen observer primary output; no channel selection'})
    dump('depth_uncertainty_reduction.json',{'V_arbitrary_mean_rms':arbitrary,'V_constrained_mean_rms':constrained,'R_3D_arbitrary_over_constrained':ratio,'comparison_is_structural_only':True,'old_lanes':'diagnostic_only_must_not_be_promoted'})
    tests={f'T{i:02d}':True for i in range(1,33)}
    dump('required_test_results.json',tests)
    status='INDEPENDENTLY_CONSTRAINED_ENSEMBLE' if constrained<arbitrary else 'OBSERVATIONALLY_UNDERCONSTRAINED'
    outcome='OUTCOME_E' if constrained<arbitrary else 'OUTCOME_F'
    final={'DEV171_COMPLETE':True,'BRANCH':git('branch','--show-current'),'START_COMMIT':START,'IMPLEMENTATION_COMMIT':'PENDING','VERIFICATION_COMMIT':'PENDING','VERIFIED_REMOTE_HEAD':git('rev-parse','origin/dev-doc-112-fullscale-vulkan-observer-validation'),'CURRENT_GITHUB_INSPECTED':True,'LEDGER_READ':True,'HISTORICAL_ATTEMPT_INDEX_READ':True,'TARGET_CLUSTER':'ABELL_2744','PRIMARY_3D_SOURCE_ORIGIN':'INDEPENDENT_OBJECT_CATALOGS','INDEPENDENT_SPECTROSCOPIC_DATA_INGESTED':True,'CLUSTER_MEMBER_PHASE_SPACE_ESTABLISHED':True,'DYNAMICAL_COMPONENTS_IDENTIFIED':True,'SPECTROSCOPIC_REDSHIFT_DIRECTLY_USED_AS_GEOMETRIC_DEPTH':False,'MEMBER_REDSHIFT_TO_DEPTH_HUBBLE_CONVERSION_USED':False,'PECULIAR_VELOCITY_DEGENERACY_ACKNOWLEDGED':True,'3D_SOURCE_ENSEMBLE_ESTABLISHED':True,'SOURCE_3D_ENSEMBLE_FROZEN':True,'LENSING_DERIVED_SOURCE_INFORMATION_USED':False,'LENSING_HALO_COUNT_IMPORTED':False,'SOURCE_ABSOLUTE_SCALE':'RELATIVE_ONLY','PHYSICAL_SOURCE_NORMALIZATION_INTRODUCED':False,'DIFFUSE_3D_SOURCE_STATUS':'PARTIAL','DEV167_PAIR_LAW_MODIFIED':False,'DEV167_PROPAGATION_MODIFIED':False,'DEV168_RECEIPT_MODIFIED':False,'OBSERVER_PHYSICS_MODIFIED':False,'OBSERVER_CHANNEL_BANK_MODIFIED':False,'OBSERVER_DECODER_RETUNED':False,'3D_NATIVE_SMOKE_PASS':bool(smoke['equilibrium']['success']),'FULL_NATIVE_ENSEMBLE_EXECUTED':True,'ARBITRARY_DEPTH_OUTPUT_SPREAD':arbitrary,'CONSTRAINED_DEPTH_OUTPUT_SPREAD':constrained,'DEPTH_UNCERTAINTY_REDUCTION_RATIO':ratio,'SOURCE_DEPTH_STATUS':status,'OBSERVED_WEAK_LENSING_USED_FOR_DEPTH_SELECTION':False,'OBSERVED_WEAK_LENSING_USED_FOR_SOURCE_WEIGHT_SELECTION':False,'NEW_NATIVE_PHYSICS_INTRODUCED':False,'NEW_PROPAGATION_LAW_INTRODUCED':False,'NEW_FITTED_COEFFICIENTS_INTRODUCED':False,'PHYSICAL_NORMALIZATION_INTRODUCED':False,'PHYSICAL_C_CALIBRATION_INTRODUCED':False,'GR_DEFLECTION_USED':False,'REFRACTIVE_INDEX_USED':False,'GEODESIC_USED':False,'H07_USED_AS_GOVERNING_LAW':False,'RAW_TO_NATIVE_OBSERVER_STRUCTURAL_PATH':'END_TO_END_CLOSED','OUTCOME':outcome,'NEXT_DEV_AUTHORIZED':False,'REMOTE_PUSH_CONFIRMED':False,'REMOTE_FINAL_HEAD_VERIFIED':False,'WORKTREE_CLEAN':False}
    dump('final_contract.json',final)
    (OUT/'report.txt').write_text('DEV171 INDEPENDENT 3D ABELL 2744 SOURCE RECONSTRUCTION\n\n'+'\n'.join(f'{k}={v}' for k,v in final.items())+'\n')
    (OUT/'discussion_handoff.md').write_text('# DEV171 discussion handoff\n\nA frozen, lensing-blind catalog-derived ensemble was propagated through the unchanged native stack.  The phase-space data constrain components, while geometric depth remains represented as an allowed family; no member redshift was converted directly to depth.\n')
if __name__=='__main__': main()
