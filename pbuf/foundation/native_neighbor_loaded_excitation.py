"""Loaded propagation experiments with no direct loading/excitation term."""
from __future__ import annotations
import numpy as np
from .native_neighbor_state import NativeNeighborState,local_link_frame
from .native_neighbor_dynamic_projection import progress,quadratic_norm

LOADS=("unloaded","weak compact","moderate compact","strong compact unsaturated","diffuse","asymmetric","two-center","shell-like","gradient")
EXCITATIONS=("single transverse A","single transverse B","45-degree combination","handedness +","handedness -","narrow-band packet","broad packet","two-packet interference")
def loading_profile(kind,n=64,strength=.5):
    x=np.arange(n); c=(n-1)/2; s=max(2,n/10)
    if kind==0:return np.zeros(n)
    if kind in (1,2,3): return strength*np.exp(-.5*((x-c)/s)**2)
    if kind==4:return strength*np.exp(-.5*((x-c)/(2*s))**2)
    if kind==5:return strength*np.exp(-.5*((x-.4*n)/s)**2)*(1+.2*(x-c)/n)
    if kind==6:return strength*(np.exp(-.5*((x-.35*n)/s)**2)+np.exp(-.5*((x-.65*n)/s)**2))/2
    if kind==7:return strength*np.exp(-.5*((abs(x-c)-2*s)/(s/2))**2)
    return strength*x/max(1,n-1)
def excitation(kind,n=64,amplitude=1.):
    x=np.arange(n); g=np.exp(-.5*((x-n/4)/(n/12 if kind!=6 else n/6))**2); p=[(1,0),(0,1),(2**-.5,2**-.5)][min(kind,2)]
    a=amplitude*g[:,None]*np.array(p)[None,:]
    if kind==7:a+=amplitude*np.exp(-.5*((x-n/3)/(n/12))**2)[:,None]*np.array([0.,1.])
    return a
def frames_from_loading(load):
    load=np.asarray(load,float); slope=np.gradient(load); return np.stack([local_link_frame([1.,s,0.]) for s in slope])
def run_case(load_kind,excitation_kind,n=64,strength=.5,amplitude=1.,steps=8):
    load=loading_profile(load_kind,n,strength); x=excitation(excitation_kind,n,amplitude); frames=frames_from_loading(load)
    st=NativeNeighborState(load,x,frames); n0=quadratic_norm(x); progress(st,steps); n1=quadratic_norm(st.transverse)
    weights=np.sum(st.transverse**2,axis=1); centroid=float(np.sum(np.arange(n)*weights)/max(np.sum(weights),1e-30))
    return {"load":f"LOAD{load_kind:02d}","excitation":f"EX{excitation_kind+1:02d}","norm_before":n0,"norm_after":n1,
            "centroid":centroid,"interaction":"GEOMETRIC_ONLY" if load_kind else "NO_INTERACTION","new_interaction_coefficients":0}
def run_matrix(n=64): return [run_case(i,j,n=n) for i in range(9) for j in range(8)]
