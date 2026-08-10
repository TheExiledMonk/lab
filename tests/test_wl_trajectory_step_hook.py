import hashlib

import numpy as np

from pbuf.labs.foundation import los_consistent_ray_geometry001 as GEO
from pbuf.wl.trajectory_state import TrajectoryStepAccumulator, summarize_trajectory


def _case():
    grid = np.linspace(-2, 2, 9)
    yy, xx = np.meshgrid(grid, grid, indexing="ij")
    field = {"xgrid": grid, "ygrid": grid,
             "rx": .01 + .002*xx, "ry": -.005 + .001*yy}
    return field, np.array([-.3, .2]), np.array([.1, -.25])


def _hash(result):
    h = hashlib.sha256()
    for name in ("x", "y", "z", "vx", "vy", "vz"):
        h.update(np.ascontiguousarray(result[1][name]).tobytes())
    return h.hexdigest()


def test_null_and_active_hook_endpoint_bitwise_parity_and_event_counts():
    field, x0, y0 = _case()
    null = GEO._propagate_g3d(field, .03, 160, x0, y0)
    observer = TrajectoryStepAccumulator(2, 160, diagnostic_indices=[0, 1])
    active = GEO._propagate_g3d(field, .03, 160, x0, y0, step_observer=observer)
    assert _hash(null) == _hash(active)
    receipt, paths = observer.finalize()
    assert np.all(receipt.path_summary["number_of_steps"] == 159)
    assert paths["positions"].shape == (160, 2, 3)


class _HostileObserver:
    def __init__(self): self.events = []
    def _attack(self, position, direction):
        self.events.append("state")
        for value in (position, direction):
            try: value[...] = 999
            except ValueError: pass
    def observe_launch(self, **kw): self.events.append("launch"); self._attack(kw["position"], kw["direction"])
    def observe_step(self, **kw): self.events.append(kw["step_index"]); self._attack(kw["position"], kw["direction"])
    def observe_termination(self, **kw): self.events.append("termination")


def test_mutation_isolation_and_exact_event_order():
    field, x0, y0 = _case(); null = GEO._propagate_g3d(field, .03, 160, x0, y0)
    hostile = _HostileObserver()
    active = GEO._propagate_g3d(field, .03, 160, x0, y0, step_observer=hostile)
    assert _hash(null) == _hash(active)
    assert hostile.events == (["launch", "state"] +
                              [item for k in range(1, 160) for item in (k, "state")] +
                              ["termination"])


def test_streaming_path_direction_curvature_and_native_parity():
    field, x0, y0 = _case(); observer = TrajectoryStepAccumulator(2, 160, [0, 1])
    GEO._propagate_g3d(field, .03, 160, x0, y0, step_observer=observer)
    receipt, paths = observer.finalize()
    for ray in range(2):
        native = {name.removeprefix("native_"): values[:, ray]
                  for name, values in paths.items() if name.startswith("native_")}
        expected, expected_native = summarize_trajectory(
            paths["positions"][:, ray], paths["directions"][:, ray], native)
        for key in ("path_length", "net_direction_change", "total_direction_change",
                    "path_curvature_integral", "path_curvature_squared_integral", "curvature_max"):
            np.testing.assert_allclose(receipt.path_summary[key][ray], expected[key], rtol=2e-14, atol=2e-14)
        for name in native:
            for key in ("integral", "absolute_integral", "square_integral", "minimum", "maximum", "rms"):
                np.testing.assert_allclose(receipt.native_path_summary[name][key][ray],
                                           expected_native[name][key], rtol=2e-14, atol=2e-14)
