"""Target-blind received-state information diagnostics for Dev Doc 116."""
from __future__ import annotations

import numpy as np

FAMILY_ORDER = ("displacement", "direction", "cross displacement-direction", "covariance",
                "Jacobian", "launch-receipt", "depth", "density/occupancy", "other canonical")

FEATURE_SPECS = (
 ("delta_u","displacement","vector / spin-1-like"),("delta_v","displacement","vector / spin-1-like"),
 ("delta_z","depth","scalar / spin-0"),("received_du","direction","vector / spin-1-like"),
 ("received_dv","direction","vector / spin-1-like"),("received_dz","depth","scalar / spin-0"),
 ("delta_u2","displacement","tensor anisotropy / spin-2-like"),("delta_v2","displacement","tensor anisotropy / spin-2-like"),
 ("delta_uv","displacement","tensor anisotropy / spin-2-like"),("dir_u2","direction","tensor anisotropy / spin-2-like"),
 ("dir_v2","direction","tensor anisotropy / spin-2-like"),("dir_uv","direction","tensor anisotropy / spin-2-like"),
 ("delta_u_dir_u","cross displacement-direction","tensor anisotropy / spin-2-like"),
 ("delta_v_dir_v","cross displacement-direction","tensor anisotropy / spin-2-like"),
 ("delta_u_dir_v","cross displacement-direction","unclassified"),("delta_v_dir_u","cross displacement-direction","unclassified"),
 ("radial_transverse_displacement","launch-receipt","scalar / spin-0"),
 ("tangential_transverse_displacement","launch-receipt","pseudoscalar / parity-sensitive"),
 ("local_ray_density","density/occupancy","scalar / spin-0"),("local_occupancy","density/occupancy","scalar / spin-0"),
 ("launch_receipt_separation","launch-receipt","scalar / spin-0"),("received_depth","depth","scalar / spin-0"),
 ("depth_change","depth","scalar / spin-0"))


def cell_index(rays, bins, extent):
    w=2*extent/bins; col=np.floor((rays["uf"]+extent)/w).astype(int); row=np.floor((rays["vf"]+extent)/w).astype(int)
    valid=np.isfinite(rays["uf"]+rays["vf"])&(row>=0)&(row<bins)&(col>=0)&(col<bins)
    return row*bins+col, valid


def feature_bank(rays, bins, extent):
    du=rays["uf"]-rays["u0"]; dv=rays["vf"]-rays["v0"]; dz=rays["rz"]
    ru,rv=rays["u0"],rays["v0"]; radius=np.hypot(ru,rv); safe=np.where(radius>0,radius,1)
    key,valid=cell_index(rays,bins,extent); count=np.bincount(key[valid],minlength=bins*bins); occ=count[key.clip(0,bins*bins-1)]
    vals=(du,dv,dz,rays["dx"],rays["dy"],rays["dz"],du*du,dv*dv,du*dv,rays["dx"]**2,rays["dy"]**2,rays["dx"]*rays["dy"],
          du*rays["dx"],dv*rays["dy"],du*rays["dy"],dv*rays["dx"],(du*ru+dv*rv)/safe,(du*(-rv)+dv*ru)/safe,
          occ.astype(float),occ.astype(float),np.sqrt(du*du+dv*dv+dz*dz),rays["rz"],rays["rz"])
    return np.column_stack(vals), [dict(name=n,family=f,spin=s) for n,f,s in FEATURE_SPECS], key, valid


def standardize(X):
    X=np.asarray(X,float); mu=np.nanmean(X,axis=0); sd=np.nanstd(X,axis=0); sd=np.where(sd>0,sd,1)
    return np.nan_to_num((X-mu)/sd)


def rank_summary(X):
    Z=standardize(X); s=np.linalg.svd(Z,full_matrices=False,compute_uv=False); tol=max(Z.shape)*np.finfo(float).eps*(s[0] if s.size else 0)
    p=s/s.sum() if s.sum() else s; eff=float(np.exp(-np.sum(p[p>0]*np.log(p[p>0])))) if p.size else 0.
    var=s*s; cum=np.cumsum(var)/var.sum() if var.sum() else var
    return {"feature_count":int(Z.shape[1]),"sample_count":int(Z.shape[0]),"rank":int(np.sum(s>tol)),"effective_rank":eff,
            "singular_values":s.tolist(),"variance_explained_first_N":{str(n):float(cum[min(n,len(cum))-1]) for n in (1,2,3,4,6,8,12) if len(cum)}}


def aggregate_cells(X,key,valid,ncells):
    out=np.full((ncells,X.shape[1]),np.nan); counts=np.bincount(key[valid],minlength=ncells)
    for j in range(X.shape[1]):
        sums=np.bincount(key[valid],weights=X[valid,j],minlength=ncells); out[:,j]=np.divide(sums,counts,out=np.full(ncells,np.nan),where=counts>0)
    return out,counts


def mixing_report(X,inventory):
    Z=standardize(X); C=np.cov(Z,rowvar=False); families=[x["family"] for x in inventory]; cross=np.array([[a!=b for b in families] for a in families])
    energy=C*C; M=float(energy[cross].sum()/energy.sum())
    pairs=[]
    for i in range(len(families)):
      for j in range(i+1,len(families)):
       if cross[i,j]: pairs.append((abs(C[i,j]),inventory[i]["name"],inventory[j]["name"],float(C[i,j])))
    pairs.sort(reverse=True)
    return {"mixing_index":M,"fraction_covariance_energy_off_block":M,"within_family_coupling":float(np.mean(np.abs(C[~cross]))),
            "cross_family_coupling":float(np.mean(np.abs(C[cross]))),"strongest_off_block_pairs":[{"a":a,"b":b,"covariance":v} for _,a,b,v in pairs[:10]]},C,np.corrcoef(Z,rowvar=False)


def spin2_rotation_tests():
    errors={}
    x=np.array([-.8,.2,.7,1.1]); y=np.array([.3,-.9,.4,1.2])
    q=np.array([x*x-y*y,2*x*y])
    for deg in (0,15,30,45,60,90,135):
        t=np.deg2rad(deg); c,s=np.cos(t),np.sin(t); xr=c*x+s*y; yr=-s*x+c*y
        actual=np.array([xr*xr-yr*yr,2*xr*yr]); c2,s2=np.cos(2*t),np.sin(2*t); expected=np.array([c2*q[0]+s2*q[1],-s2*q[0]+c2*q[1]])
        errors[str(deg)]=float(np.max(np.abs(actual-expected))/max(np.max(np.abs(expected)),1e-30))
    return {"angles_degrees":list(map(int,errors)),"normalized_errors":errors,"max_normalized_error":max(errors.values()),"tolerance":1e-10,"pass":max(errors.values())<=1e-10}
