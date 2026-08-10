import numpy as np

from pbuf.wl.bundle_transport_geometry import (
    deposit_bundle, fit_bundle_scale, reconstruct_launch_topology,
)


def regular_rays(side=6):
    y, x = np.mgrid[-1:1:complex(side), -1:1:complex(side)]
    n = side * side
    return {"u0": x.ravel(), "v0": y.ravel(), "launch_x": x.ravel(),
            "launch_y": y.ravel(), "uf": (2*x+.25*y).ravel(),
            "vf": (-.1*x+1.5*y).ravel(), "rx": (2*x+.25*y).ravel(),
            "ry": (-.1*x+1.5*y).ravel(), "rz": np.zeros(n),
            "dx": np.ones(n), "dy": np.ones(n), "dz": np.ones(n),
            "e1": np.array([1.,0.,0.]), "e2": np.array([0.,1.,0.])}


def test_exact_topology_and_affine_matrix_are_preserved():
    rays=regular_rays(); topology=reconstruct_launch_topology(rays)
    bundle=fit_bundle_scale(rays,topology,1,np.ones(36),np.zeros(36))
    assert topology["shape"] == (6,6)
    np.testing.assert_allclose(bundle["bundle_matrix_2d"][2,2],
                               [[2.,.25],[-.1,1.5]],atol=1e-12)
    assert bundle["valid"][2,2]
    assert bundle["parity_class"][2,2] == 1


def test_received_deposition_retains_component_and_spin2_fields():
    rays=regular_rays();topology=reconstruct_launch_topology(rays)
    bundle=fit_bundle_scale(rays,topology,1,np.ones(36),np.zeros(36))
    maps=deposit_bundle(bundle,8,bounds=(-3.,3.))
    assert maps["occupancy"].sum() == 36
    assert maps["launch_component_count"].max() >= 1
    assert np.all(np.isfinite(maps["bundle_q1"]))
