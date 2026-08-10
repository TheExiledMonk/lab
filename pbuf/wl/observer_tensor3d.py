"""Derived target-blind tensor and information diagnostics for observer volumes."""
import numpy as np

AXES=("u","v","w")

def _matrix(volume,prefix):
    m=np.empty(volume["occupancy"].shape+(3,3))
    means=np.stack([volume[f"mean_{prefix}_{a}"] for a in AXES],axis=-1)
    for i,a in enumerate(AXES):
        for j,b in enumerate(AXES):
            key=f"mean_{prefix}_{a}2" if i==j else f"mean_{prefix}_{''.join(sorted((a,b),key=AXES.index))}"
            m[...,i,j]=volume[key]-means[...,i]*means[...,j]
    return m

def tensor_bank(volume):
    d=_matrix(volume,"delta"); q=_matrix(volume,"dir")
    dm=np.stack([volume[f"mean_delta_{a}"] for a in AXES],-1); qm=np.stack([volume[f"mean_dir_{a}"] for a in AXES],-1)
    cross=np.empty_like(d)
    for i,a in enumerate(AXES):
        for j,b in enumerate(AXES): cross[...,i,j]=volume[f"mean_delta_{a}_dir_{b}"]-dm[...,i]*qm[...,j]
    out={"displacement_covariance":d,"direction_covariance":q,"cross_tensor":cross,
         "cross_symmetric":.5*(cross+np.swapaxes(cross,-1,-2)),"cross_antisymmetric":.5*(cross-np.swapaxes(cross,-1,-2))}
    for name,t in list(out.items()):
        if name=="cross_antisymmetric": continue
        s=.5*(t+np.swapaxes(t,-1,-2)); trace=np.trace(s,axis1=-2,axis2=-1); tl=s-trace[...,None,None]*np.eye(3)/3
        clean=np.nan_to_num(s); eigval,eigvec=np.linalg.eigh(clean); order=np.argsort(eigval,axis=-1)[...,::-1]
        eigval=np.take_along_axis(eigval,order,-1); eigvec=np.take_along_axis(eigvec,order[...,None,:],-1)
        denom=np.sum(np.abs(eigval),axis=-1)+np.finfo(float).eps
        anis=np.sqrt((eigval[...,0]-eigval[...,1])**2+(eigval[...,1]-eigval[...,2])**2+(eigval[...,2]-eigval[...,0])**2)/denom
        out.update({f"{name}_trace":trace,f"{name}_traceless":tl,f"{name}_eigenvalues":eigval,
                    f"{name}_principal_eigenvector":eigvec[...,0],f"{name}_secondary_eigenvector":eigvec[...,1],f"{name}_anisotropy":anis,
                    f"{name}_isotropy_fraction":np.abs(trace/3)/(np.sum(np.abs(eigval),axis=-1)+np.finfo(float).eps)})
    return out

def effective_rank(x):
    x=np.asarray(x); x=x[np.all(np.isfinite(x),axis=1)];
    if len(x)<2:return 0.
    s=np.linalg.svd(x-x.mean(0),compute_uv=False); p=s*s; p=p/p.sum() if p.sum() else p
    return float(np.exp(-np.sum(p[p>0]*np.log(p[p>0]))))

def feature_diagnostics(volume):
    names=[k for k in volume if k.startswith("mean_") and not k.startswith("mean_received")]
    x=np.stack([volume[k] for k in names],-1); occ=volume["occupancy"]>0
    per=[]; mixing=[]; pc1=[]
    for z in range(x.shape[2]):
        xx=x[:,:,z][occ[:,:,z]]; per.append(effective_rank(xx))
        if len(xx)>2:
            clean=np.nan_to_num(xx); c=np.corrcoef(clean,rowvar=False); mixing.append(float(np.nanmean(np.abs(c-np.eye(c.shape[0])))))
            _,_,vh=np.linalg.svd(clean-clean.mean(0),full_matrices=False); pc1.append(vh[0].tolist())
        else:
            mixing.append(None); pc1.append(None)
    depth_var=np.nanvar(x,axis=2)
    autocorr=[]
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            a=x[i,j]
            for f in range(a.shape[1]):
                q=a[:,f]; m=np.isfinite(q[:-1])&np.isfinite(q[1:])
                if m.sum()>2 and np.std(q[:-1][m])*np.std(q[1:][m])>0: autocorr.append(np.corrcoef(q[:-1][m],q[1:][m])[0,1])
    thirds=[]
    for indices in np.array_split(np.arange(x.shape[2]),3):
        xx=x[:,:,indices][occ[:,:,indices]]
        c=np.corrcoef(np.nan_to_num(xx),rowvar=False) if len(xx)>2 else None
        thirds.append(float(np.nanmean(np.abs(c-np.eye(c.shape[0])))) if c is not None else None)
    blocks=[]
    for ui in np.array_split(np.arange(x.shape[0]),2):
        for vi in np.array_split(np.arange(x.shape[1]),2):
            xx=x[np.ix_(ui,vi,np.arange(x.shape[2]))]; oo=occ[np.ix_(ui,vi,np.arange(x.shape[2]))]; blocks.append(effective_rank(xx[oo]))
    return {"feature_names":names,"global_effective_rank":effective_rank(x[occ]),"effective_rank_by_depth":per,"mixing_by_depth":mixing,
            "mixing_front_middle_rear":thirds,"pc1_evolution_by_depth":pc1,"selected_spatial_block_effective_rank":blocks,
            "depth_variance_mean":float(np.nanmean(depth_var)),"depth_autocorrelation_mean":float(np.nanmean(autocorr)) if autocorr else None,
            "occupied_depth_cells_mean":float(np.mean(np.sum(occ,axis=2)))}
