"""Exact regular launch-grid triangle topology and mapped geometry."""
from __future__ import annotations
import numpy as np
from .transport_receiver_decode import rasterize


def regular_triangles(shape):
    ny,nx=shape;ids=np.arange(ny*nx).reshape(shape)
    a=ids[:-1,:-1].ravel();b=ids[:-1,1:].ravel();c=ids[1:,:-1].ravel();d=ids[1:,1:].ravel()
    return np.vstack((np.column_stack((a,b,c)),np.column_stack((d,c,b))))


def _signed_area(x,y,tri):
    x0,x1,x2=(x[tri[:,i]] for i in range(3));y0,y1,y2=(y[tri[:,i]] for i in range(3))
    return .5*((x1-x0)*(y2-y0)-(y1-y0)*(x2-x0))


def mapped_mesh(u0,v0,uf,vf):
    shape=np.asarray(u0).shape;tri=regular_triangles(shape)
    a0=_signed_area(np.ravel(u0),np.ravel(v0),tri);af=_signed_area(np.ravel(uf),np.ravel(vf),tri)
    ratio=np.divide(af,a0,out=np.zeros_like(af),where=np.abs(a0)>0)
    cu=np.mean(np.ravel(uf)[tri],axis=1);cv=np.mean(np.ravel(vf)[tri],axis=1)
    e0=np.column_stack((np.ravel(u0)[tri[:,1]]-np.ravel(u0)[tri[:,0]],np.ravel(v0)[tri[:,1]]-np.ravel(v0)[tri[:,0]]))
    ef=np.column_stack((np.ravel(uf)[tri[:,1]]-np.ravel(uf)[tri[:,0]],np.ravel(vf)[tri[:,1]]-np.ravel(vf)[tri[:,0]]))
    t0=np.arctan2(e0[:,1],e0[:,0]);tf=np.arctan2(ef[:,1],ef[:,0]);dt=tf-t0
    return {"triangles":tri,"launch_area":a0,"received_area":af,"transport_area_ratio":ratio,
            "centroid_u":cu,"centroid_v":cv,"launch_orientation":t0,"received_orientation":tf,
            "orientation_change":dt,"orientation_spin2_q1":np.cos(2*dt),"orientation_spin2_q2":np.sin(2*dt)}


def rasterize_mesh(mesh,vertex_values,resolution=64,bounds=(-8.,8.)):
    """Conservative deterministic triangle-centroid rasterization.

    Each exact regular-grid triangle contributes its mean vertex value weighted
    by mapped triangle area.  This preserves topology and area without Delaunay.
    """
    tri=mesh["triangles"];z=np.mean(np.ravel(vertex_values)[tri],axis=1);w=np.abs(mesh["received_area"])
    num=rasterize(mesh["centroid_u"],mesh["centroid_v"],z*w,resolution,bounds)
    den=rasterize(mesh["centroid_u"],mesh["centroid_v"],w,resolution,bounds)
    return np.divide(num,den,out=np.zeros_like(num),where=den>0)
