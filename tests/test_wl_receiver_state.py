import numpy as np

from pbuf.wl.receiver_state import (build_receiver_state, channel_manifest,
    manifest_sha256, receiver_matrix)


def packet(side=9, transform=lambda u,v:(u,v), depth=None, direction=None):
    v,u=np.mgrid[:side,:side].astype(float);uf,vf=transform(u,v)
    w=np.zeros_like(u) if depth is None else depth(u,v)
    launch=np.column_stack((u.ravel(),v.ravel(),np.zeros(side*side)))
    receive=np.column_stack((np.asarray(uf).ravel(),np.asarray(vf).ravel(),np.asarray(w).ravel()))
    d=np.tile([0.,0.,1.],(side*side,1)) if direction is None else direction(u,v).reshape(-1,3)
    direct=np.linalg.norm(receive-launch,axis=1)
    return {"endpoint_launch_position":launch,"endpoint_receive_position":receive,
        "endpoint_initial_direction":np.tile([0.,0.,1.],(side*side,1)),"endpoint_final_direction":d,
        "path_path_length":direct,"path_straight_line_distance":direct,"path_path_excess":np.zeros(side*side),
        "path_total_direction_change":np.zeros(side*side),"path_net_direction_change":np.zeros(side*side),
        "path_path_curvature_integral":np.zeros(side*side)}


def center(state,name,side=9): return state.channel_bank["C7"]["s1_"+name].reshape(side,side)[side//2,side//2]


def test_identity_receiver():
    s=build_receiver_state(packet())
    assert np.all(s.channel_bank["C2"]["delta_u"]==0)
    assert np.all(s.channel_bank["C2"]["delta_v"]==0)
    np.testing.assert_allclose(center(s,"local_received_area_ratio"),1,atol=1e-14)
    np.testing.assert_allclose(center(s,"local_anisotropy"),0,atol=1e-14)


def test_translation_expansion_and_stretch():
    identity=build_receiver_state(packet())
    translated=build_receiver_state(packet(transform=lambda u,v:(u+3,v-2)))
    np.testing.assert_allclose(translated.channel_bank["C2"]["delta_u"],3)
    np.testing.assert_allclose(translated.channel_bank["C2"]["delta_v"],-2)
    np.testing.assert_allclose(center(translated,"position_cov_uu"),center(identity,"position_cov_uu"))
    expanded=build_receiver_state(packet(transform=lambda u,v:(2*u,2*v)))
    np.testing.assert_allclose(center(expanded,"local_received_area_ratio"),4)
    np.testing.assert_allclose(center(expanded,"local_anisotropy"),0,atol=1e-14)
    stretched=build_receiver_state(packet(transform=lambda u,v:(3*u,2*v)))
    assert center(stretched,"position_eigenvalue_1") > center(stretched,"position_eigenvalue_2")
    np.testing.assert_allclose(center(stretched,"local_received_area_ratio"),6)


def test_rotation_depth_and_direction_independence():
    theta=.37;c,s=np.cos(theta),np.sin(theta)
    rotated=build_receiver_state(packet(transform=lambda u,v:(c*u-s*v,s*u+c*v)))
    np.testing.assert_allclose(center(rotated,"local_received_area_ratio"),1,atol=1e-14)
    np.testing.assert_allclose(center(rotated,"local_anisotropy"),0,atol=1e-14)
    base=build_receiver_state(packet())
    depth=build_receiver_state(packet(depth=lambda u,v:u+2*v))
    np.testing.assert_array_equal(base.channel_bank["C0"]["receive_u"],depth.channel_bank["C0"]["receive_u"])
    assert center(depth,"var_receive_w")>center(base,"var_receive_w")
    varying=build_receiver_state(packet(direction=lambda u,v:np.stack((.01*u,.01*v,np.ones_like(u)),axis=-1)))
    np.testing.assert_array_equal(base.channel_bank["C0"]["receive_u"],varying.channel_bank["C0"]["receive_u"])
    assert center(varying,"direction_trace")>0


def test_validity_alias_constant_and_manifest_determinism():
    state=build_receiver_state(packet())
    assert not state.validity_masks["C7_s1"].reshape(9,9)[0,0]
    state.channel_bank["C9"]["duplicate_receive_u"]=state.channel_bank["C0"]["receive_u"].copy()
    state.channel_bank["C9"]["constant"]=np.ones(state.ray_count)
    _,_,_,audit=receiver_matrix(state,("C0","C9"))
    assert audit["constant_channels"]
    assert audit["alias_channels"]
    assert manifest_sha256(channel_manifest(state))==manifest_sha256(channel_manifest(state))
