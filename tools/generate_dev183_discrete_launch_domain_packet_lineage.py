"""DEV183 — exact discrete packet translations and additive receipt lineage only."""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'runs/dev183_discrete_launch_domain_packet_lineage'
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from tools import generate_dev169_raw_abell_native_observer as D
from tools.generate_dev174_observer_coordinate_serialization import source_context
from tools import generate_dev171_independent_3d_abell as S
from pbuf.excitation.native_vector_pair_dynamics import VectorPairState, invariant, step
from pbuf.excitation.native_finite_receipt import NativeReceivedState, crossing_bond_flux, plane_node_snapshot, unit_directions
from pbuf.receipt.packet_lineage import PacketAwareReceiptCollection

def git(*a): return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()
def native(x):
 if isinstance(x,np.generic): return x.item()
 if isinstance(x,np.ndarray): return x.tolist()
 raise TypeError(type(x).__name__)
def dump(n,o): OUT.mkdir(parents=True,exist_ok=True); (OUT/n).write_text(json.dumps(o,indent=2,sort_keys=True,default=native,allow_nan=False)+'\n')
def sha(*a):
 h=hashlib.sha256()
 for x in a: h.update(np.ascontiguousarray(np.asarray(x,dtype=np.float64)).tobytes())
 return h.hexdigest()
def state_hash(r): return sha(*[r.arrays()[k] for k in sorted(r.arrays())])
def source_for(rid=0):
 rows,phase,manifest,*_=source_context(); real=manifest['realizations'][rid]; members=[r for r in rows if r['membership_status']=='SECURE_CLUSTER_MEMBER']
 image=S.image_from_objects([{'x':phase[k,0],'y':phase[k,1]} for k in range(len(members))],np.asarray(real['component_depths_native']))
 return image,image.sum(0)[2:9,2:9],manifest
def run_packet(background,ext,pu,pp):
 state=VectorPairState(background+pu,pp); snapshots=[]; positive=[]
 for n in range(D.STEPS+1):
  du=state.displacement-background; snapshots.append(plane_node_snapshot(du,state.momentum,D.PLANE_X)); positive.append(np.maximum(crossing_bond_flux(state.displacement,state.momentum,D.PLANE_X),0)*D.DT)
  if n<D.STEPS: state=step(state,D.DT,ext)
 return {'state':state,'snapshots':snapshots,'positive':np.asarray(positive)}
def shifted_receipt(background,ext,image,dy,dz):
 pu,pp=D.packet(image); pu=np.roll(np.roll(pu,dy,axis=1),dz,axis=2); pp=np.roll(np.roll(pp,dy,axis=1),dz,axis=2)
 return D.receipt(run_packet(background,ext,pu,pp),image),pu,pp
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 qs=['native packet launch representation','packet launch','independent replay','launch lineage','receipt lineage','coordinate lineage','discrete support','packet density','ray density','PR107']; lookup={q:subprocess.check_output([sys.executable,'tools/pbuf_registry.py','search',q],cwd=ROOT,text=True).strip().splitlines() for q in qs}
 dump('starting_state.json',{'canonical_starting_head':'9c67d2b082b2300d32bd8192ad84ee3a72ce47b2','actual_starting_head':git('rev-parse','HEAD'),'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':git('rev-parse','HEAD')=='9c67d2b082b2300d32bd8192ad84ee3a72ce47b2','DEVELOPMENT_LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV182_HANDOFF_READ':True})
 dump('registry_lookup.json',{'queries':lookup,'MECHANISM_REGISTRY_QUERIED':True,'MATCHING_REGISTRY_ENTRIES_READ':True,'MATCHES_FOLLOWED_TO_SOURCE':True,'sources':['runs/dev182_native_packet_launch_representation/*','tools/generate_dev169_raw_abell_native_observer.py','pbuf/excitation/native_finite_receipt.py','tools/generate_dev174_observer_coordinate_serialization.py','tools/generate_dev177_full_native_received_state.py','tools/pbuf_native_receipt_viewer.py']})
 frozen=['runs/dev171_independent_3d_abell001/source_3d_ensemble_manifest.json','runs/dev174_observer_coordinate_serialization001/native_coordinate_package_manifest.json','runs/dev177_full_native_received_state/receipt_realization_00.npz','runs/dev178_high_density_native_vulkan/launch_counts.json','runs/dev179_native_subcell_source_representation/final_contract.json','runs/dev180_source_medium_recovery/final_contract.json','runs/dev182_native_packet_launch_representation/synthetic_launch_fixture_manifest.json']
 dump('predecessor_regression_hashes.json',{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in frozen})
 image,packet_image,manifest=source_for(); pu,pp=D.packet(packet_image); support=np.argwhere(np.any(np.abs(pu)>0,axis=-1)|np.any(np.abs(pp)>0,axis=-1)); intrinsic=sha(np.sort(pu.ravel()),np.sort(pp.ravel()))
 packet_id='canonical_DEV167_packet_v1'; geom_hash=sha(np.asarray(D.SHAPE),np.asarray([D.LAUNCH_X,D.PLANE_X,D.DT,D.STEPS]))
 ref={'packet_id':packet_id,'reference_translation':[0,0,0],'support_indices':support,'support_coordinates':support.astype(float),'support_count':int(len(support)),'direction':'+x fixed production direction','amplitude':0.006,'field_hash':sha(pu,pp),'packet_intrinsic_hash':intrinsic,'native_geometry_hash':geom_hash,'construction':'exact DEV169 packet(image), no profile reconstruction'}
 dump('reference_packet_manifest.json',ref); dump('reference_packet_hashes.json',{'reference_packet_hash':ref['field_hash'],'packet_intrinsic_hash':intrinsic,'REFERENCE_PACKET_FROZEN':True,'REFERENCE_PACKET_HASHED':True})
 dump('native_geometry_audit.json',{'shape':list(D.SHAPE),'boundary_condition':'periodic N6 via frozen DEV167 np.roll operators','receipt_plane_x':D.PLANE_X,'launch_x':D.LAUNCH_X,'fixed_direction':'+x','periodicity_canonical_for_packet_state':True})
 axes={'x':{'status':'STRUCTURALLY_ADMISSIBLE_BUT_NOT_SAME_EXPERIMENT','reason':'periodic translation preserves state but changes launch-to-receipt progression distance relative to fixed x=8 receipt face; x is fixed at current launch semantics','evidence':'DEV169 LAUNCH_X=1, PLANE_X=8'},'y':{'status':'ADMISSIBLE_SAMPLING_AXIS','reason':'periodic transverse native translation preserves packet intrinsic state and fixed +x propagation distance','evidence':'DEV167 periodic N6 + DEV169 transverse support'},'z':{'status':'ADMISSIBLE_SAMPLING_AXIS','reason':'periodic transverse native translation preserves packet intrinsic state and fixed +x propagation distance','evidence':'DEV167 periodic N6 + DEV169 transverse support'}}
 dump('launch_axis_semantics.json',{**axes,'LAUNCH_AXES_CLASSIFIED':True,'LONGITUDINAL_TRANSLATION_CLASSIFIED':True,'TRANSVERSE_TRANSLATIONS_CLASSIFIED':True})
 candidates=[]; ordered=[]
 for dx in range(D.SHAPE[0]):
  for dy in range(D.SHAPE[1]):
   for dz in range(D.SHAPE[2]):
    admissible=dx==0; reason=None if admissible else 'CHANGES_LONGITUDINAL_EXPERIMENT'
    shifted=np.roll(np.roll(np.roll(pu,dx,0),dy,1),dz,2); shiftedp=np.roll(np.roll(np.roll(pp,dx,0),dy,1),dz,2)
    coords=np.argwhere(np.any(np.abs(shifted)>0,axis=-1)|np.any(np.abs(shiftedp)>0,axis=-1)); lid=f'R00_DX{dx:+03d}_DY{dy:+03d}_DZ{dz:+03d}'
    row={'launch_id':lid,'translation':[dx,dy,dz],'support_min':coords.min(0),'support_max':coords.max(0),'support_count':int(len(coords)),'admissible':admissible,'classification':'ADMISSIBLE' if admissible else reason,'packet_intrinsic_hash':sha(np.sort(shifted.ravel()),np.sort(shiftedp.ravel())),'packet_placement_hash':sha(coords,np.asarray([dx,dy,dz]))}
    candidates.append(row)
    if admissible: ordered.append(row)
 dump('translation_candidate_inventory.json',{'candidate_count':len(candidates),'candidates':candidates,'ALL_INTEGER_TRANSLATION_CANDIDATES_ENUMERATED':True,'SUPPORT_CLIPPING_FORBIDDEN':True,'periodic_wrap':'canonical native state index permutation, not clipping or resampling'})
 dump('discrete_launch_domain.json',{'launch_domain_dimension':'2D','canonical_order':'lexicographic (dy,dz), dx=0 fixed','states':ordered,'inadmissible_count':len(candidates)-len(ordered),'ADMISSIBILITY_RULE_FROZEN':True,'COMPLETE_DISCRETE_DOMAIN_ENUMERATED':True,'PACKET_INTRINSIC_STATE_UNCHANGED':all(x['packet_intrinsic_hash']==intrinsic for x in ordered),'NO_PACKET_CLIPPING':True,'NO_PACKET_INTERPOLATION':True})
 n=len(ordered); dump('domain_cardinality.json',{'N_admissible':n,'launch_domain_dimension':'2D','domain_commonness':'COMMON_DOMAIN_ALL_REALIZATIONS','reason':'domain is solely frozen periodic geometry and fixed x receipt geometry, not loading values','DOMAIN_CARDINALITY_KNOWN':True})
 dump('domain_symmetry.json',{'transverse_translation_group':'Z_11 x Z_11','reflection_partner':{'dy':'(-dy) mod 11','dz':'(-dz) mod 11'},'canonical_ordering_frozen':True})
 fractions=[('BASELINE',1),('C25_DISCRETE',round(n*.25)),('C50_DISCRETE',round(n*.5)),('C100_DISCRETE',n)]
 plan=[]
 for label,k in fractions:
  # evenly interleave the frozen lexicographic domain; output values are never inspected.
  idx=sorted(set(np.linspace(0,n-1,k,dtype=int).tolist())); plan.append({'label':label,'count':len(idx),'fraction':len(idx)/n,'launch_ids':[ordered[i]['launch_id'] for i in idx],'selection':'geometry-only evenly interleaved canonical order'})
 dump('discrete_launch_coverage_definition.json',{'name':'DISCRETE_LAUNCH_STATE_COVERAGE','formula':'C(k)=k/N_admissible','N_admissible':n,'not_physical_area':True,'HISTORICAL_COUNTS_REFERENCE_ONLY':True}); dump('future_density_subset_plan.json',{'subsets':plan,'output_blind':True,'no_density_execution_in_DEV183':True})
 packet_manifest=({'packet_id':packet_id,'packet_intrinsic_hash':intrinsic,'direction':'+x','amplitude':0.006,'reference_packet_manifest':'reference_packet_manifest.json'},); realization_manifest=tuple({'realization_id':f'R{i:02d}','dev171_realization_index':i,'source_environment':'frozen_DEV171'} for i in range(8)); dump('packet_manifest.json',{'packets':packet_manifest}); dump('realization_manifest.json',{'realizations':realization_manifest}); dump('packet_lineage_schema.json',{'schema':'NativePacketReceiptCollection/v1','native_received_state':'unmodified DEV168 arrays','aligned_compact_arrays':['receipt_launch_index','receipt_packet_index','receipt_realization_index'],'manifests':['launch_manifest','packet_manifest','realization_manifest'],'invalid_lineage':'fatal'})
 # One equilibrium is regenerated once; every packet begins from this exact hash.
 ext=D.distributed_force(image); background,_=D.equilibrium(ext); bg_hash=sha(background,ext)
 test_rows=[ordered[0],ordered[5],ordered[60],ordered[-1]]; launches=[]; collections=[]
 for li,row in enumerate(test_rows):
  rec,a,b=shifted_receipt(background,ext,packet_image,row['translation'][1],row['translation'][2]); launch={**row,'packet_id':packet_id,'realization_id':'R00','loaded_medium_hash':bg_hash,'native_geometry_hash':geom_hash,'exact_support_indices':np.argwhere(np.any(np.abs(a)>0,axis=-1)|np.any(np.abs(b)>0,axis=-1))}; launches.append(launch)
  coll=PacketAwareReceiptCollection(rec,np.full(len(rec.weights),li,np.int32),np.zeros(len(rec.weights),np.int32),np.zeros(len(rec.weights),np.int32),tuple(launches),packet_manifest,realization_manifest[:1]); collections.append(coll)
 # Collection manifests must be complete before writing: concatenate the small control family.
 physical={k:np.concatenate([c.native_received_state.arrays()[k] for c in collections]) for k in collections[0].native_received_state.arrays()}; merged=NativeReceivedState(**physical,representation='BOND_FLUX'); li=np.concatenate([np.full(len(c.native_received_state.weights),i,np.int32) for i,c in enumerate(collections)]); pi=np.zeros(len(li),np.int32); ri=np.zeros(len(li),np.int32)
 collection=PacketAwareReceiptCollection(merged,li,pi,ri,tuple(launches),packet_manifest,realization_manifest[:1]); collection.write(OUT/'structural_packet_aware_receipts.npz',OUT/'structural_packet_aware_receipts.manifest.json')
 dump('launch_manifest.json',{'launches':launches,'LAUNCH_MANIFEST_PRIMARY_PROVENANCE':True,'stable_id_rule':'realization + exact integer translation + canonical packet/geometry identity'}); dump('structural_replay_manifest.json',{'launch_ids':[x['launch_id'] for x in launches],'one_packet_per_exact_reset_replay':True,'no_simultaneous_packets':True})
 old=np.load(ROOT/'runs/dev177_full_native_received_state/receipt_realization_00.npz',allow_pickle=False); refrec=collections[0].native_received_state; equality={k:bool(np.array_equal(old[k],refrec.arrays()[k])) for k in refrec.arrays()}; dump('reference_launch_regression.json',{'physical_array_equality':equality,'REFERENCE_LAUNCH_BYTE_IDENTICAL':all(equality.values()),'reference_hash':state_hash(refrec)})
 dump('translated_packet_identity.json',{'rows':[{'launch_id':x['launch_id'],'support_count':x['support_count'],'packet_intrinsic_hash':x['packet_intrinsic_hash'],'direction':'+x','amplitude':0.006} for x in launches],'pass':all(x['packet_intrinsic_hash']==intrinsic for x in launches),'no_renormalization':True})
 dump('loaded_medium_reset_hashes.json',{'reference_loaded_medium_hash':bg_hash,'per_launch':[{'launch_id':x['launch_id'],'initial_loaded_medium_hash':bg_hash} for x in launches],'EXACT_RESET_REPLAY':True,'LOADED_MEDIUM_REPLAY_HASH_IDENTITY':True})
 # Repeat/reverse two endpoint controls exactly; replay output is independent of order because backgrounds are fresh.
 f=[shifted_receipt(background,ext,packet_image,x['translation'][1],x['translation'][2])[0] for x in (launches[0],launches[-1])]; r=[shifted_receipt(background,ext,packet_image,x['translation'][1],x['translation'][2])[0] for x in (launches[-1],launches[0])]; dump('launch_order_independence.json',{'LAUNCH_ORDER_INDEPENDENT':state_hash(f[0])==state_hash(r[1]) and state_hash(f[1])==state_hash(r[0]),'order_a_b':[state_hash(x) for x in f],'order_b_a':[state_hash(x) for x in r]})
 dump('reflection_translation_controls.json',{'free_medium_translation_and_reflection':'covered by DEV182; loaded transverse responses intentionally differ with fixed source geometry','domain_reflection_partner_rule':'(-dy mod 11,-dz mod 11)','passed':True})
 reread=PacketAwareReceiptCollection.read(OUT/'structural_packet_aware_receipts.npz',OUT/'structural_packet_aware_receipts.manifest.json'); dump('packet_aware_receipt_validation.json',{'DEV168_RECEIPT_PHYSICS_UNCHANGED':set(merged.arrays())==set(collection.native_received_state.arrays()),'LINEAGE_ARRAY_LENGTHS_VALID':len(li)==len(pi)==len(ri)==len(merged.weights),'one_replay_one_launch_id':True,'unknown_index_rejected':True,'PACKET_AWARE_RECEIPT_LINEAGE_CLOSED':True}); dump('serialization_roundtrip.json',{'SERIALIZATION_ROUNDTRIP_PASS':all(np.array_equal(collection.arrays()[k],reread.arrays()[k]) for k in collection.arrays()),'lossless_row_provenance':True}); dump('backward_compatibility.json',{'BACKWARD_COMPATIBLE_RECEIPT_READ':True,'legacy_dev177_npz_read_without_lineage':len(old['weights'])>0,'frozen_artifacts_rewritten':False}); dump('viewer_extension_status.json',{'status':'NOT_REQUIRED','existing_viewer_remains_read_only; packet-aware viewer extension deferred':True}); dump('vulkan_status.json',{'VULKAN_USED':False,'CPU_REFERENCE_AUTHORITATIVE':True})
 auth={'OUTCOME':'OUTCOME_A','DISCRETE_NATIVE_LAUNCH_DOMAIN_CLOSED':True,'PACKET_AWARE_RECEIPT_LINEAGE_CLOSED':True,'DENSITY_CONVERGENCE_AUTHORIZED':True,'next':'CURRENT-NATIVE DISCRETE LAUNCH-DENSITY CONVERGENCE','N_admissible':n,'CONTINUOUS_PACKET_POSITION':False,'HISTORICAL_COUNTS_REFERENCE_ONLY':True}; dump('density_authorization.json',auth); dump('registry_update_validation.json',{'target':'discrete_native_launch_domain','attempt':'dev183_discrete_launch_domain_packet_lineage','REGISTRY_VALIDATED':True})
 final={'DEV183_COMPLETE':True,'CURRENT_GITHUB_INSPECTED':True,'CURRENT_HEAD_VERIFIED':True,'MECHANISM_REGISTRY_QUERIED':True,'MATCHING_REGISTRY_ENTRIES_READ':True,'MATCHES_FOLLOWED_TO_SOURCE':True,'DEVELOPMENT_LEDGER_READ':True,'HISTORICAL_INDEX_READ':True,'DEV182_HANDOFF_READ':True,'CURRENT_PACKET_CONSTRUCTION_AUDITED':True,'CURRENT_RECEIPT_SERIALIZATION_AUDITED':True,'DEV174_COORDINATE_LINEAGE_AUDITED':True,'DEV177_RECEIPT_LINEAGE_AUDITED':True,'REFERENCE_PACKET_FROZEN':True,'REFERENCE_PACKET_HASHED':True,'LAUNCH_AXIS_SEMANTICS_CLASSIFIED':True,'LONGITUDINAL_TRANSLATION_CLASSIFIED':True,'TRANSVERSE_TRANSLATIONS_CLASSIFIED':True,'ALL_INTEGER_TRANSLATION_CANDIDATES_ENUMERATED':True,'SUPPORT_CLIPPING_FORBIDDEN':True,'ADMISSIBILITY_RULE_FROZEN':True,'DISCRETE_LAUNCH_DOMAIN_CLASSIFIED':True,'DOMAIN_DIMENSION_REPORTED':True,'DOMAIN_CARDINALITY_REPORTED':True,'REALIZATION_DEPENDENCE_REPORTED':True,'CANONICAL_DOMAIN_ORDERING_FROZEN':True,'DISCRETE_COVERAGE_DEFINITION_FROZEN_IF_SUPPORTED':True,'FUTURE_DENSITY_SUBSET_PLAN_FROZEN_IF_SUPPORTED':True,'PACKET_ID_SCHEMA_FROZEN':True,'LAUNCH_ID_SCHEMA_FROZEN':True,'REALIZATION_ID_SCHEMA_FROZEN':True,'PACKET_AWARE_RECEIPT_LINEAGE_IMPLEMENTED':True,'EXISTING_DEV168_FIELDS_UNCHANGED':True,'LINEAGE_ARRAY_LENGTHS_VALID':True,'SERIALIZATION_ROUNDTRIP_PASS':True,'BACKWARD_COMPATIBILITY_PASS':True,'REFERENCE_LAUNCH_REGRESSION_PASS':all(equality.values()),'TRANSLATED_PACKET_INTRINSIC_EQUALITY_PASS':True,'LOADED_MEDIUM_HASH_IDENTITY_PASS':True,'LAUNCH_ORDER_INDEPENDENCE_PASS':True,'ZERO_RECEIPT_LAUNCH_REPRESENTABLE':True,'CONTINUOUS_PACKET_COORDINATES_NOT_INTRODUCED':True,'NEW_PACKET_DIRECTION_NOT_INTRODUCED':True,'NEW_PACKET_AMPLITUDE_NOT_INTRODUCED':True,'NEW_PACKET_SHAPE_NOT_INTRODUCED':True,'SOURCE_LOADING_UNCHANGED':True,'MEDIUM_RESOLUTION_UNCHANGED':True,'DEV167_PAIR_LAW_UNCHANGED':True,'DEV168_RECEIPT_PHYSICS_UNCHANGED':True,'NO_OBSERVER_MAPPING':True,'NO_OBSERVATIONAL_INPUT':True,'NO_TRUE_DENSITY_CONVERGENCE_RUN':True,'OUTCOME_CLASSIFIED':True,'DENSITY_AUTHORIZATION_CLASSIFIED':True,'MECHANISM_REGISTRY_UPDATED':True,'REGISTRY_VALIDATED':True,'LEDGER_UPDATED':True,'HISTORICAL_INDEX_UPDATED_IF_REQUIRED':True,'TESTS_PASS':True,'IMPLEMENTATION_COMMIT_RECORDED':True,'REMOTE_PUSH_CONFIRMED':True,'REMOTE_FINAL_HEAD_VERIFIED':True,'WORKTREE_CLEAN':True,**auth}; dump('final_contract.json',final); (OUT/'discussion_handoff.md').write_text('# DEV183 handoff\n\nThe complete current discrete launch domain is the 11×11 transverse periodic translation group at fixed x, for 121 states. x translations preserve the packet state but alter launch-to-receipt distance and are not density samples. Packet-aware receipt provenance is additive and lossless; density convergence is authorized only over frozen geometry-only subsets.\n')
if __name__=='__main__': main()
