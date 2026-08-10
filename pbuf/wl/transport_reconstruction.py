"""Target-blind reconstruction lanes for frozen launch/receipt transport."""
from __future__ import annotations
import numpy as np
from scipy import ndimage
from .transport_receiver_decode import rasterize
from .reconstructed_geometry import geometry_from_derivatives
from .transport_mesh import mapped_mesh,rasterize_mesh

PATTERNS=("uniform_lattice","checkerboard","isotropic_point_grid","horizontal_bars","vertical_bars","bars_45","bars_135","concentric_rings")

def diagnostic_patterns(shape):
    y,x=np.indices(shape);p=max(8,min(shape)//16);cx=(shape[1]-1)/2;cy=(shape[0]-1)/2
    return {"uniform_lattice":(((x%p)==0)|((y%p)==0)).astype(float),
            "checkerboard":(((x//p)+(y//p))%2).astype(float),
            "isotropic_point_grid":(((x%p)==p//2)&((y%p)==p//2)).astype(float),
            "horizontal_bars":((y%p)<max(1,p//4)).astype(float),
            "vertical_bars":((x%p)<max(1,p//4)).astype(float),
            "bars_45":((np.mod(x+y,p))<max(1,p//4)).astype(float),
            "bars_135":((np.mod(x-y,p))<max(1,p//4)).astype(float),
            "concentric_rings":((np.mod(np.hypot(x-cx,y-cy),p))<max(1,p//4)).astype(float)}

def _deposit_geometry(t,g,resolution):
    out={"transport_density":rasterize(t["uf"],t["vf"],np.ones_like(t["uf"]),resolution)}
    for k in ("local_area_change","transport_area_ratio","local_orientation","local_anisotropy","local_curvature","axis_ratio","principal_axis","spin2_shape_q1","spin2_shape_q2"):
        out[k]=rasterize(t["uf"],t["vf"],g[k],resolution)
    out["depth"]=rasterize(t["uf"],t["vf"],t["wf"],resolution)
    for k in ("dir_u","dir_v","dir_w"):out[k]=rasterize(t["uf"],t["vf"],t[k],resolution)
    return out

def reconstruct_endpoint_only(t,patterns=None,resolution=64):
    patterns=patterns or diagnostic_patterns(t["u0"].shape);out={}
    # Endpoint-only has deliberately discarded launch identity: evaluate pattern
    # on nearest launch-grid coordinate corresponding to the received endpoint.
    u=t["u0"][0];v=t["v0"][:,0];iu=np.clip(np.rint((t["uf"]-u[0])/(u[1]-u[0])).astype(int),0,len(u)-1);iv=np.clip(np.rint((t["vf"]-v[0])/(v[1]-v[0])).astype(int),0,len(v)-1)
    for n,z in patterns.items():out["reconstructed_intensity_"+n]=rasterize(t["uf"],t["vf"],z[iv,iu],resolution)
    out["transport_density"]=rasterize(t["uf"],t["vf"],np.ones_like(t["uf"]),resolution);return out

def reconstruct_raw(t,patterns=None,resolution=64):
    patterns=patterns or diagnostic_patterns(t["u0"].shape);out={}
    for n,z in patterns.items():out["reconstructed_intensity_"+n]=rasterize(t["uf"],t["vf"],z,resolution)
    out["transport_density"]=rasterize(t["uf"],t["vf"],np.ones_like(t["uf"]),resolution);return out

def reconstruct_order(t,first,second=None,patterns=None,resolution=64):
    patterns=patterns or diagnostic_patterns(t["u0"].shape);g=geometry_from_derivatives(first,second);out=_deposit_geometry(t,g,resolution)
    # Deterministic symmetric subcell quadrature. J/H alter sample receipt points;
    # source values remain launch-correspondent and no coefficients are fitted.
    du=float(t["launch_spacing_u"])/4;dv=float(t["launch_spacing_v"])/4
    acc={n:np.zeros((resolution,resolution)) for n in patterns};count=np.zeros((resolution,resolution))
    J=g["jacobian"]
    for a,b in ((-du,-dv),(-du,dv),(du,-dv),(du,dv)):
        ru=t["uf"]+J[...,0,0]*a+J[...,0,1]*b;rv=t["vf"]+J[...,1,0]*a+J[...,1,1]*b
        if second is not None:
            ru+=.5*(second["d_uu_delta_u"]*a*a+2*second["d_uv_delta_u"]*a*b+second["d_vv_delta_u"]*b*b)
            rv+=.5*(second["d_uu_delta_v"]*a*a+2*second["d_uv_delta_v"]*a*b+second["d_vv_delta_v"]*b*b)
        count+=rasterize(ru,rv,np.ones_like(ru),resolution)
        for n,z in patterns.items():acc[n]+=rasterize(ru,rv,z,resolution)
    for n in patterns:out["reconstructed_intensity_"+n]=np.divide(acc[n],count,out=np.zeros_like(count),where=count>0)
    return out

def reconstruct_mesh(t,patterns=None,resolution=64):
    patterns=patterns or diagnostic_patterns(t["u0"].shape);m=mapped_mesh(t["u0"],t["v0"],t["uf"],t["vf"]);out={}
    for n,z in patterns.items():out["reconstructed_intensity_"+n]=rasterize_mesh(m,z,resolution)
    for k in ("transport_area_ratio","orientation_change","orientation_spin2_q1","orientation_spin2_q2"):out[k]=rasterize(m["centroid_u"],m["centroid_v"],m[k],resolution)
    out["transport_density"]=rasterize(m["centroid_u"],m["centroid_v"],np.ones_like(m["centroid_u"]),resolution);return out,m

def reconstruct_patches(t,geometry,widths=(2,4,8,16),resolution=64):
    out={}
    for w in widths:
        for k in ("transport_area_ratio","local_orientation","local_anisotropy","local_curvature","spin2_shape_q1","spin2_shape_q2"):
            z=ndimage.uniform_filter(np.asarray(geometry[k],float),size=w,mode="nearest")
            out[f"patch{w}_{k}"]=rasterize(t["uf"],t["vf"],z,resolution)
    return out

def quadratic_taylor_error(alpha=.03,beta=.02,n=31):
    """Analytic key test: Hessian Taylor map is exact for the quadratic warp."""
    u=np.linspace(-1,1,n);v=np.linspace(-1,1,n);U,V=np.meshgrid(u,v);du=.37*(u[1]-u[0]);dv=-.31*(v[1]-v[0])
    uf=U+alpha*U**2;vf=V+beta*U*V
    exact_u=(U+du)+alpha*(U+du)**2;exact_v=(V+dv)+beta*(U+du)*(V+dv)
    first_u=uf+(1+2*alpha*U)*du;first_v=vf+beta*V*du+(1+beta*U)*dv
    second_u=first_u+alpha*du**2;second_v=first_v+beta*du*dv
    return float(np.sqrt(np.mean((first_u-exact_u)**2+(first_v-exact_v)**2))),float(np.sqrt(np.mean((second_u-exact_u)**2+(second_v-exact_v)**2)))
