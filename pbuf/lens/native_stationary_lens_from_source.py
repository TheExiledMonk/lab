"""Stationary N6 response for a distributed Dev159 source constraint."""
from __future__ import annotations
import numpy as np

from pbuf.excitation.native_bond_state import positive_gradient, relational_imbalance
from pbuf.source.native_distributed_source_constraint import distributed_source_imposed_excursion


def stationary_distributed_response(source: np.ndarray) -> np.ndarray:
    forcing = distributed_source_imposed_excursion(source)
    axes = [2*np.pi*np.fft.fftfreq(n) for n in forcing.shape]
    mesh = np.meshgrid(*axes, indexing="ij")
    stiffness = (2.0/3.0)*sum(np.sin(k/2.0)**2 for k in mesh)
    fhat = np.fft.fftn(forcing); qhat = np.zeros_like(fhat)
    mask = stiffness > 0
    qhat[mask] = fhat[mask]/stiffness[mask]
    q = np.fft.ifftn(qhat).real
    return q-q.mean()


def equilibrium_residual(q: np.ndarray, source: np.ndarray) -> np.ndarray:
    return relational_imbalance(q)/6.0 + distributed_source_imposed_excursion(source)


def response_inventory(q: np.ndarray) -> dict:
    bonds = positive_gradient(q); bounded = bonds/(1.0+np.abs(bonds))
    return {"node_excursion_peak":float(np.max(np.abs(q))),
            "bond_excursion_peak":float(np.max(np.abs(bonds))),
            "bounded_strain_peak":float(np.max(np.abs(bounded))),
            "bounded_stress_peak":float(np.max(np.abs(bounded))),
            "accumulated_deformation":float(np.sum(np.abs(bonds))),
            "medium_response_l1":float(np.sum(np.abs(relational_imbalance(q)/6.0)))}


def weighted_geometry(field: np.ndarray) -> dict:
    w = np.abs(np.asarray(field,float)); total=float(w.sum())
    if not total > 0: raise ValueError("field has no response")
    zz,yy,xx=np.indices(w.shape); coords=(zz,yy,xx)
    means=[float((w*a).sum()/total) for a in coords]
    cov=np.empty((3,3))
    for i,a in enumerate(coords):
        for j,b in enumerate(coords): cov[i,j]=float((w*(a-means[i])*(b-means[j])).sum()/total)
    transverse=cov[1:,1:]; eigvals,eigvecs=np.linalg.eigh(transverse); order=np.argsort(eigvals)[::-1]
    angle=float(np.degrees(np.arctan2(eigvecs[1,order[0]],eigvecs[0,order[0]])))
    return {"centroid_zyx":means,"r_z_rms":float(np.sqrt(max(cov[0,0],0))),
            "r_perp_rms":float(np.sqrt(max(np.trace(transverse),0))),
            "r_3d_rms":float(np.sqrt(max(np.trace(cov),0))),
            "transverse_principal_variances":eigvals[order].tolist(),"transverse_principal_angle_degrees":angle,
            "total_native_response":total,"peak_deformation":float(w.max())}
