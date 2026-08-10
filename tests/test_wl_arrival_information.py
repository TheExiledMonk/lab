import numpy as np
from pbuf.wl.arrival_formation import form_arrival_events, native_receiver_plane
from pbuf.wl.arrival_information import (geometry_rank_ladder,reconstruction_audits,
    endpoint_comparison,information_preservation,survival_audit)


def events(side=9):
    y,x=np.mgrid[-1:1:complex(side),-1:1:complex(side)];uv=np.column_stack((x.ravel(),y.ravel()))
    pos=np.column_stack((uv,np.full(len(uv),4.7)));d=np.column_stack((.02*uv[:,0],-.03*uv[:,1],np.ones(len(uv))))
    return form_arrival_events(pos,d,native_receiver_plane(),launch_uv=uv,side=side,scales=(1,))


def test_rank_reconstruction_and_endpoint_comparison():
    e=events();ladder=geometry_rank_ladder(e);assert [q["stage"] for q in ladder]==["G0","G1","G2","G3","G4"]
    audits=reconstruction_audits(e);assert len(audits)==4
    p=native_receiver_plane();pos=np.column_stack((e.receiver_reference["launch_u"],e.receiver_reference["launch_v"],np.full(e.ray_count,4.7)))
    comp=endpoint_comparison(e,pos,p);assert comp["classification"]=="EXPLICIT_INTERSECTION_REQUIRED"


def test_foreign_key_preservation_and_depth_survival():
    e=events();primary={"ray_index":np.arange(e.ray_count),"receive_w":np.linspace(0,1,e.ray_count)}
    a=information_preservation(e,primary);assert a["DEV129_RECEIVER_FIELDS_LOST"]==0
    s=survival_audit(e,{"receive_w":primary["receive_w"]},"depth");assert s["available_channels"]==1
