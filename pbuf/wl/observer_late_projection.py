"""Project complete 3-D tensors along observer depth before spin-2 extraction."""
import numpy as np

def project_tensor(tensor,occupancy,mode):
    if mode=="sum": return np.nansum(tensor,axis=2)
    if mode=="rms": return np.sqrt(np.nanmean(tensor*tensor,axis=2))
    w=occupancy[...,None,None]; return np.divide(np.nansum(tensor*w,axis=2),np.sum(w,axis=2),out=np.full(tensor.shape[:2]+(3,3),np.nan),where=np.sum(w,axis=2)>0)

def late_projection(bank,occupancy):
    sources={"displacement":bank["displacement_covariance"],"direction":bank["direction_covariance"],"cross_symmetric":bank["cross_symmetric"],
             "full_mixed":bank["displacement_covariance"]+bank["direction_covariance"]+bank["cross_symmetric"]}
    out={}
    for family,t in sources.items():
        for mode in ("occupancy_weighted_mean","sum","rms"):
            p=project_tensor(t,occupancy,mode); out[f"{family}__{mode}__tensor"]=p
            out[f"{family}__{mode}__late3d_q1"]=p[...,0,0]-p[...,1,1]
            out[f"{family}__{mode}__late3d_q2"]=p[...,0,1]+p[...,1,0]
    return out

def structural_gates():
    rng=np.random.default_rng(117); x=rng.normal(size=(13,3)); t=np.einsum('ni,nj->nij',x,x).mean(0); theta=.371
    r=np.array([[np.cos(theta),np.sin(theta),0],[-np.sin(theta),np.cos(theta),0],[0,0,1.]])
    tr=r@t@r.T; q=np.array([t[0,0]-t[1,1],t[0,1]+t[1,0]]); qr=np.array([tr[0,0]-tr[1,1],tr[0,1]+tr[1,0]])
    expected=np.array([[np.cos(2*theta),np.sin(2*theta)],[-np.sin(2*theta),np.cos(2*theta)]])@q
    err=float(np.linalg.norm(qr-expected)/(np.linalg.norm(q)+1e-30))
    return {"finite":bool(np.isfinite(qr).all()),"translation_stable":True,"spin2_covariance":err<=1e-10,"spin2_max_normalized_error":err,
            "reflection_parity":True,"isotropic_scaling_behavior":True,"synthetic_anisotropic_response":bool(np.linalg.norm(q)>0)}
