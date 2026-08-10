"""Local derivative and conditioning diagnostics for Dev133 optics."""
from __future__ import annotations
import numpy as np
from pbuf.wl.reverse_transport import transport_diagnostics


def compose_transport(optical, pbuf, *, compatible=True):
    if not compatible:
        return None, "TRANSPORT_COMPOSITION_NOT_VALID"
    a, b = np.asarray(optical, float), np.asarray(pbuf, float)
    if a.shape[-2:] != (2, 2) or b.shape[-2:] != (2, 2):
        raise ValueError("transport matrices must end in (2,2)")
    return np.matmul(a, b), "TRANSPORT_COMPOSITION_VALID"


def derivative_invariants(jacobian):
    j = np.asarray(jacobian, float)
    diag = transport_diagnostics(j.reshape((-1, 2, 2)))
    sv = np.stack((diag["sigma_max"], diag["sigma_min"]), axis=-1)
    cond = np.divide(sv[:, 0], sv[:, 1], out=np.full(len(sv), np.inf), where=sv[:, 1] > 0)
    orientation = np.arctan2(j.reshape(-1, 2, 2)[:, 1, 0], j.reshape(-1, 2, 2)[:, 0, 0])
    anisotropy = np.divide(sv[:, 0]-sv[:, 1], sv[:, 0]+sv[:, 1], out=np.zeros(len(sv)), where=(sv[:, 0]+sv[:, 1])>0)
    return {"determinant": diag["detJ"], "singular_values": sv, "condition_number": cond,
            "orientation": orientation, "anisotropy": anisotropy, "classification": diag["classification"]}


def finite_difference_derivatives(function, points, *, step=1e-5):
    """Central-difference 2-D Jacobian and component Hessians."""
    x = np.asarray(points, float)
    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("points must have shape (N,2)")
    n = len(x); j = np.empty((n, 2, 2)); h = np.empty((n, 2, 2, 2)); eye = np.eye(2)*step
    f0 = np.asarray(function(x), float)
    for a in range(2):
        fp, fm = np.asarray(function(x+eye[a]), float), np.asarray(function(x-eye[a]), float)
        j[:, :, a] = (fp-fm)/(2*step); h[:, :, a, a] = (fp-2*f0+fm)/(step*step)
    fpp=np.asarray(function(x+eye[0]+eye[1]),float);fpm=np.asarray(function(x+eye[0]-eye[1]),float)
    fmp=np.asarray(function(x-eye[0]+eye[1]),float);fmm=np.asarray(function(x-eye[0]-eye[1]),float)
    h[:,:,0,1]=h[:,:,1,0]=(fpp-fpm-fmp+fmm)/(4*step*step)
    return j, h

