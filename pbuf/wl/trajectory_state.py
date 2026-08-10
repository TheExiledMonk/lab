"""Target-blind, additive 3-D trajectory-state bookkeeping (Dev127).

Nothing in this module propagates a ray.  It consumes states produced by a
propagator and therefore cannot feed measurements back into the dynamics.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np

PATH_FRACTIONS = np.linspace(0.0, 1.0, 9)
BUNDLE_SCALES = (1, 2, 4, 8, 16, 32)
EPS = np.finfo(np.float64).eps


def diagnostic_ray_indices(side):
    """Deterministic target-blind center/cardinal/diagonal/radial/grid sample."""
    side=int(side)
    if side < 1: raise ValueError("side must be positive")
    c=(side-1)/2; chosen={(int(round(c)),int(round(c)))}
    edge=(0,side-1)
    chosen.update(((int(round(c)),j) for j in edge));chosen.update(((i,int(round(c))) for i in edge))
    chosen.update((i,j) for i in edge for j in edge)
    for radius in np.linspace(0,max(c,0),16):
        angle=2*np.pi*(len(chosen)%16)/16
        chosen.add((int(np.clip(round(c+radius*np.sin(angle)),0,side-1)),
                    int(np.clip(round(c+radius*np.cos(angle)),0,side-1))))
    flat=np.linspace(0,side*side-1,min(32,side*side),dtype=np.int64)
    chosen.update(divmod(int(q),side) for q in flat)
    return np.array(sorted(i*side+j for i,j in chosen),dtype=np.int64)


def native_field_manifest():
    """Fields naturally sampled by the frozen G3D propagator."""
    return [
        {"name":"rx_sample","source_module":"pbuf.labs.foundation.los_consistent_ray_geometry001._sample",
         "units":"native response per propagation length","where_sampled":"before direction update and at checkpoints",
         "currently_used_by_propagation":True},
        {"name":"ry_sample","source_module":"pbuf.labs.foundation.los_consistent_ray_geometry001._sample",
         "units":"native response per propagation length","where_sampled":"before direction update and at checkpoints",
         "currently_used_by_propagation":True},
    ]


def dev121_mapping():
    """Semantic mapping of the immutable Dev121 29-channel contract."""
    from .multiscale_transport_relations import canonical_manifest
    rows=[]
    for item in canonical_manifest():
        name=item["name"];order=item["derivative_order"]
        if name in ("u0","v0","uf","vf","wf","dir_u","dir_v","dir_w"):
            relation="exact equivalent"
        elif order == 1: relation="history-enhanced equivalent"
        elif order == 2: relation="history-enhanced equivalent"
        else: relation="endpoint-only"
        rows.append({"index":item["index"],"dev121_name":name,"mapping":relation,
                     "dev127_source":"endpoint" if order==0 else ("first_order_transport_history" if order==1 else "second_order_transport_history")})
    return rows


def _vectors(a, name):
    x = np.asarray(a, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3 or len(x) < 1:
        raise ValueError(f"{name} must have shape (samples, 3)")
    if not np.all(np.isfinite(x)):
        raise ValueError(f"{name} contains non-finite values")
    return x


def _unit(d):
    n = np.linalg.norm(d, axis=-1, keepdims=True)
    return np.divide(d, n, out=np.zeros_like(d), where=n > 0)


def _angles(a, b):
    return np.arccos(np.clip(np.sum(_unit(a) * _unit(b), axis=-1), -1.0, 1.0))


def sample_path_fractions(positions, directions, fractions=PATH_FRACTIONS):
    """Linearly sample state by normalized cumulative arc length.

    Interpolation is diagnostic-only and never modifies propagation state.
    """
    x, d = _vectors(positions, "positions"), _vectors(directions, "directions")
    if len(x) != len(d):
        raise ValueError("positions and directions must have equal sample counts")
    f = np.asarray(fractions, float)
    if np.any((f < 0) | (f > 1)) or np.any(np.diff(f) < 0):
        raise ValueError("fractions must be ordered in [0, 1]")
    seg = np.linalg.norm(np.diff(x, axis=0), axis=1)
    s = np.r_[0.0, np.cumsum(seg)]
    q = f * s[-1] if s[-1] > 0 else np.zeros_like(f)
    out_x = np.column_stack([np.interp(q, s, x[:, j]) for j in range(3)])
    out_d = _unit(np.column_stack([np.interp(q, s, d[:, j]) for j in range(3)]))
    return {"fractions": f, "positions": out_x, "directions": out_d,
            "interpolation": "linear_in_cumulative_arc_length"}


def summarize_trajectory(positions, directions, native_fields: Mapping[str, Sequence[float]] | None = None,
                         basis: np.ndarray | None = None):
    """Return the fixed-size endpoint/path/native summary for one ray."""
    x, d = _vectors(positions, "positions"), _vectors(directions, "directions")
    if len(x) != len(d):
        raise ValueError("positions and directions must have equal sample counts")
    basis = np.eye(3) if basis is None else np.asarray(basis, float)
    if basis.shape != (3, 3):
        raise ValueError("basis must contain e_u, e_v, e_w as rows")
    delta = np.diff(x, axis=0); ds = np.linalg.norm(delta, axis=1)
    theta = _angles(d[:-1], d[1:]) if len(d) > 1 else np.empty(0)
    length = float(ds.sum()); direct = float(np.linalg.norm(x[-1] - x[0]))
    cumulative = np.r_[0.0, np.cumsum(ds)]
    mid = cumulative[:-1] + .5 * ds
    curvature = np.divide(theta, ds, out=np.zeros_like(theta), where=ds > 0)
    abs_weight = np.abs(curvature) * ds
    k1 = float(abs_weight.sum()); k2 = float(np.sum(curvature * curvature * ds))
    kmean = float(np.sum(curvature * ds) / length) if length else 0.0
    kvar = float(np.sum((curvature-kmean)**2 * ds) / length) if length else 0.0
    kcentroid = float(np.sum(mid * abs_weight) / (length * k1)) if length and k1 else 0.0
    imax = int(np.argmax(theta)) if len(theta) else 0
    # Frozen straight reference uses the actual launch state and initial direction.
    ref = x[0] + cumulative[:, None] * _unit(d[:1])[0]
    transverse = (x-ref) @ basis[:2].T
    trans_norm = np.linalg.norm(transverse, axis=1)
    it = int(np.argmax(trans_norm))
    proj = _unit(d) @ basis.T
    # Trapezoids are the exact discrete equivalent for sampled path data.
    def trap(values):
        return float(np.sum(.5*(values[:-1]+values[1:])*ds)) if len(ds) else 0.0
    net = d[-1]-d[0]; net_angle = float(_angles(d[:1], d[-1:])[0])
    out = {
        "launch_x": float(x[0,0]), "launch_y": float(x[0,1]), "launch_z": float(x[0,2]),
        "receive_x": float(x[-1,0]), "receive_y": float(x[-1,1]), "receive_z": float(x[-1,2]),
        "initial_dir_x": float(d[0,0]), "initial_dir_y": float(d[0,1]), "initial_dir_z": float(d[0,2]),
        "final_dir_x": float(d[-1,0]), "final_dir_y": float(d[-1,1]), "final_dir_z": float(d[-1,2]),
        "net_direction_change_x": float(net[0]), "net_direction_change_y": float(net[1]),
        "net_direction_change_z": float(net[2]), "net_direction_change": net_angle,
        "path_length": length, "straight_line_distance": direct, "path_excess": max(0.0, length-direct),
        "total_direction_change": float(theta.sum()),
        "direction_history_ratio": float(theta.sum()/(net_angle+EPS)),
        "maximum_local_direction_change": float(theta[imax]) if len(theta) else 0.0,
        "step_of_maximum_direction_change": imax,
        "path_location_of_maximum_direction_change": float(mid[imax]/length) if length and len(mid) else 0.0,
        "path_curvature_integral": k1, "path_curvature_squared_integral": k2,
        "curvature_mean": kmean, "curvature_rms": float(np.sqrt(k2/length)) if length else 0.0,
        "curvature_max": float(np.max(np.abs(curvature))) if len(curvature) else 0.0,
        "curvature_variance": kvar, "curvature_path_centroid": kcentroid,
        "integral_du_ds": trap(proj[:,0]), "integral_dv_ds": trap(proj[:,1]),
        "integral_abs_du_ds": trap(np.abs(proj[:,0])), "integral_abs_dv_ds": trap(np.abs(proj[:,1])),
        "integral_du2_ds": trap(proj[:,0]**2), "integral_dv2_ds": trap(proj[:,1]**2),
        "max_abs_delta_u": float(np.max(np.abs(transverse[:,0]))),
        "max_abs_delta_v": float(np.max(np.abs(transverse[:,1]))),
        "max_transverse_norm": float(trans_norm[it]), "step_of_max_transverse_displacement": it,
        "path_location_of_max_transverse_displacement": float(cumulative[it]/length) if length else 0.0,
        "final_delta_u": float(transverse[-1,0]), "final_delta_v": float(transverse[-1,1]),
        "number_of_steps": max(0, len(x)-1), "termination_status": "provided_final_state",
    }
    native = {}
    for name, values in sorted((native_fields or {}).items()):
        v = np.asarray(values, float)
        if v.shape != (len(x),) or not np.all(np.isfinite(v)):
            raise ValueError(f"native field {name!r} must have one finite scalar per sample")
        wmean = trap(v)/length if length else float(v[0])
        ia = int(np.argmax(np.abs(v)))
        native[name] = {"integral": trap(v), "absolute_integral": trap(np.abs(v)),
                        "square_integral": trap(v*v), "minimum": float(v.min()),
                        "maximum": float(v.max()), "mean": wmean,
                        "rms": float(np.sqrt(trap(v*v)/length)) if length else abs(float(v[0])),
                        "path_location_of_max_abs": float(cumulative[ia]/length) if length else 0.0}
    return out, native


@dataclass
class TrajectoryReceipt3D:
    endpoint: dict
    path_summary: dict
    native_path_summary: dict = field(default_factory=dict)
    bundle_summary: dict = field(default_factory=dict)

    def as_dict(self):
        return {"endpoint": self.endpoint, "path_summary": self.path_summary,
                "native_path_summary": self.native_path_summary,
                "bundle_summary": self.bundle_summary}


class TrajectoryAccumulator:
    """Streaming facade; retains one ray only, never a production ray bank."""
    def __init__(self, basis=None):
        self.basis = basis; self.positions=[]; self.directions=[]; self.native={}

    def update(self, position, direction, native_fields=None):
        self.positions.append(np.asarray(position, float).copy())
        self.directions.append(np.asarray(direction, float).copy())
        supplied = native_fields or {}
        if self.native and set(supplied) != set(self.native):
            raise ValueError("native field set changed during accumulation")
        for key, value in supplied.items(): self.native.setdefault(key, []).append(float(value))

    def finalize(self):
        summary, native = summarize_trajectory(self.positions, self.directions, self.native, self.basis)
        endpoint_keys = {k for k in summary if k.startswith(("launch_", "receive_", "initial_dir_", "final_dir_"))}
        return TrajectoryReceipt3D({k: summary[k] for k in endpoint_keys},
                                   {k: v for k, v in summary.items() if k not in endpoint_keys}, native)


class TrajectoryStepAccumulator:
    """Passive vectorized observer for the frozen propagator.

    Fixed-size summaries are streamed for every ray.  Complete sample arrays
    are retained only for the deterministic diagnostic subset.  Callback
    return values are deliberately always ``None``.
    """

    def __init__(self, ray_count, expected_steps, diagnostic_indices=(), basis=None):
        self.ray_count = int(ray_count)
        self.expected_steps = int(expected_steps)
        self.diagnostic_indices = np.asarray(diagnostic_indices, dtype=np.int64)
        if self.ray_count < 1 or self.expected_steps < 1:
            raise ValueError("ray_count and expected_steps must be positive")
        if np.any((self.diagnostic_indices < 0) | (self.diagnostic_indices >= self.ray_count)):
            raise ValueError("invalid diagnostic ray index")
        self.basis = np.eye(3) if basis is None else np.asarray(basis, dtype=np.float64)
        if self.basis.shape != (3, 3):
            raise ValueError("basis must have shape (3, 3)")
        n = self.ray_count
        self._launched = False; self._terminated = False; self._last_step = 0
        self._path_length = np.zeros(n); self._total_turn = np.zeros(n)
        self._max_turn = np.zeros(n); self._max_turn_step = np.zeros(n, dtype=np.int64)
        self._k2 = np.zeros(n); self._kmax = np.zeros(n)
        self._k_weighted_mid = np.zeros(n)
        self._native = {}
        self._diag_position = []; self._diag_direction = []; self._diag_ds = []
        self._diag_native = {}
        self._fraction_steps = set()
        for value in PATH_FRACTIONS * max(0, self.expected_steps - 1):
            self._fraction_steps.update((int(np.floor(value)), int(np.ceil(value))))
        self._fraction_states = {}

    @staticmethod
    def _state(value, n, name):
        a = np.asarray(value, dtype=np.float64)
        if a.shape != (n, 3) or not np.all(np.isfinite(a)):
            raise ValueError(f"{name} must have shape ({n}, 3) and be finite")
        return a

    def observe_launch(self, *, ray_index, position, direction, launch_coordinates,
                       native_state=None):
        if self._launched:
            raise RuntimeError("duplicate launch event")
        np.testing.assert_array_equal(np.asarray(ray_index), np.arange(self.ray_count))
        self._previous_position = self._state(position, self.ray_count, "position").copy()
        self._initial_position = self._previous_position.copy()
        self._previous_direction = self._state(direction, self.ray_count, "direction").copy()
        self._initial_direction = self._previous_direction.copy()
        self._last_native = {k: np.asarray(v, dtype=np.float64).copy()
                             for k, v in (native_state or {}).items()}
        for key, value in self._last_native.items():
            if value.shape != (self.ray_count,):
                raise ValueError(f"native field {key!r} has wrong shape")
            self._native[key] = {q: np.zeros(self.ray_count) for q in
                                 ("integral", "absolute_integral", "square_integral")}
            self._native[key].update(minimum=value.copy(), maximum=value.copy(),
                                     max_abs=np.abs(value), argmax_path=np.zeros(self.ray_count))
        q = self.diagnostic_indices
        self._diag_position.append(self._previous_position[q].copy())
        self._diag_direction.append(self._previous_direction[q].copy())
        for key, value in self._last_native.items(): self._diag_native.setdefault(key, []).append(value[q].copy())
        if 0 in self._fraction_steps:
            self._fraction_states[0] = (self._previous_position.copy(), self._previous_direction.copy())
        self._launched = True

    def observe_step(self, *, ray_index, step_index, position, direction,
                     native_state=None, ds=None):
        if not self._launched or self._terminated or int(step_index) != self._last_step + 1:
            raise RuntimeError("invalid step event order")
        pos = self._state(position, self.ray_count, "position")
        direction = self._state(direction, self.ray_count, "direction")
        segment = np.linalg.norm(pos - self._previous_position, axis=1)
        turn = _angles(self._previous_direction, direction)
        midpoint = self._path_length + .5 * segment
        curvature = np.divide(turn, segment, out=np.zeros_like(turn), where=segment > 0)
        self._path_length += segment; self._total_turn += turn
        replace = turn > self._max_turn
        self._max_turn[replace] = turn[replace]; self._max_turn_step[replace] = int(step_index)
        self._k2 += curvature * curvature * segment
        self._kmax = np.maximum(self._kmax, np.abs(curvature))
        self._k_weighted_mid += midpoint * np.abs(curvature) * segment
        supplied = {k: np.asarray(v, dtype=np.float64) for k, v in (native_state or {}).items()}
        if set(supplied) != set(self._last_native):
            raise ValueError("native field set changed during accumulation")
        for key, value in supplied.items():
            if value.shape != (self.ray_count,): raise ValueError(f"native field {key!r} has wrong shape")
            prior = self._last_native[key]; bank = self._native[key]
            bank["integral"] += .5 * (prior + value) * segment
            bank["absolute_integral"] += .5 * (np.abs(prior) + np.abs(value)) * segment
            bank["square_integral"] += .5 * (prior*prior + value*value) * segment
            bank["minimum"] = np.minimum(bank["minimum"], value)
            bank["maximum"] = np.maximum(bank["maximum"], value)
            replace_native = np.abs(value) > bank["max_abs"]
            bank["max_abs"][replace_native] = np.abs(value[replace_native])
            bank["argmax_path"][replace_native] = self._path_length[replace_native]
            self._last_native[key] = value.copy()
        q = self.diagnostic_indices
        self._diag_position.append(pos[q].copy()); self._diag_direction.append(direction[q].copy())
        self._diag_ds.append(segment[q].copy())
        for key, value in supplied.items(): self._diag_native[key].append(value[q].copy())
        self._previous_position = pos.copy(); self._previous_direction = direction.copy()
        self._last_step = int(step_index)
        if self._last_step in self._fraction_steps:
            self._fraction_states[self._last_step] = (pos.copy(), direction.copy())

    def observe_termination(self, *, ray_index, termination_status, final_step_index):
        if not self._launched or self._terminated or int(final_step_index) != self._last_step:
            raise RuntimeError("invalid termination event")
        self.termination_status = str(termination_status); self._terminated = True

    def finalize(self):
        if not self._terminated: raise RuntimeError("trajectory has not terminated")
        direct = np.linalg.norm(self._previous_position - self._initial_position, axis=1)
        net_angle = _angles(self._initial_direction, self._previous_direction)
        k1 = self._total_turn
        length = self._path_length
        path = {
            "path_length": length.copy(), "straight_line_distance": direct,
            "path_excess": np.maximum(0., length-direct), "net_direction_change": net_angle,
            "total_direction_change": self._total_turn.copy(),
            "maximum_local_direction_change": self._max_turn.copy(),
            "path_location_of_maximum_direction_change": np.divide(
                self._max_turn_step, max(1, self._last_step)),
            "path_curvature_integral": k1.copy(),
            "path_curvature_squared_integral": self._k2.copy(),
            "curvature_mean": np.divide(k1, length, out=np.zeros_like(k1), where=length>0),
            "curvature_rms": np.sqrt(np.divide(self._k2, length, out=np.zeros_like(length), where=length>0)),
            "curvature_max": self._kmax.copy(),
            "curvature_path_centroid": np.divide(self._k_weighted_mid, length*k1,
                                                   out=np.zeros_like(length), where=(length*k1)>0),
            "number_of_steps": np.full(self.ray_count, self._last_step, dtype=np.int64),
        }
        endpoint = {"launch_position": self._initial_position.copy(),
                    "receive_position": self._previous_position.copy(),
                    "initial_direction": self._initial_direction.copy(),
                    "final_direction": self._previous_direction.copy()}
        native = {}
        for key, bank in self._native.items():
            item = {k: v.copy() for k, v in bank.items() if k != "max_abs"}
            item["mean"] = np.divide(bank["integral"], length, out=np.zeros_like(length), where=length>0)
            item["rms"] = np.sqrt(np.divide(bank["square_integral"], length,
                                             out=np.zeros_like(length), where=length>0))
            item["path_location_of_max_abs"] = np.divide(bank["argmax_path"], length,
                                                           out=np.zeros_like(length), where=length>0)
            native[key] = item
        return TrajectoryReceipt3D(endpoint, path, native), self.diagnostic_full_paths()

    def fixed_fraction_states(self):
        """Return Dev127's nine states using its frozen linear interpolator."""
        positions=[]; directions=[]; last=max(0, self.expected_steps-1)
        for fraction in PATH_FRACTIONS:
            value=float(fraction*last); lo=int(np.floor(value)); hi=int(np.ceil(value)); weight=value-lo
            p0,d0=self._fraction_states[lo]; p1,d1=self._fraction_states[hi]
            positions.append((1-weight)*p0+weight*p1)
            directions.append(_unit((1-weight)*d0+weight*d1))
        return {"fractions": PATH_FRACTIONS.copy(), "positions": np.asarray(positions),
                "directions": np.asarray(directions),
                "interpolation": "DEV127 linear_in_cumulative_arc_length; constant-step frozen path"}

    def diagnostic_full_paths(self):
        positions = np.asarray(self._diag_position)
        directions = np.asarray(self._diag_direction)
        return {"ray_indices": self.diagnostic_indices.copy(), "positions": positions,
                "directions": directions, "ds": np.asarray(self._diag_ds),
                **{f"native_{k}": np.asarray(v) for k, v in self._diag_native.items()}}


def bundle_history(positions, launch_uv, fractions=PATH_FRACTIONS):
    """Local-cell area and first/second-order transport at fixed depths.

    ``positions`` is (depth, row, column, 3); transverse coordinates are the
    first two components.  Derivatives use numpy's deterministic centered-grid
    convention, including exact mixed Hessian components.
    """
    p=np.asarray(positions,float); uv=np.asarray(launch_uv,float)
    if p.ndim != 4 or p.shape[-1] != 3 or uv.shape != p.shape[1:3]+(2,):
        raise ValueError("incompatible bundle position/launch grid")
    f=np.asarray(fractions,float)
    if len(p) != len(f): raise ValueError("one bundle state is required per fraction")
    du=float(np.median(np.diff(uv[0,:,0]))); dv=float(np.median(np.diff(uv[:,0,1])))
    if du == 0 or dv == 0: raise ValueError("degenerate launch grid")
    J=np.empty(p.shape[:3]+(2,2)); H=np.empty(p.shape[:3]+(2,3))
    for k in range(len(f)):
        for c in range(2):
            gy,gx=np.gradient(p[k,...,c],dv,du,edge_order=2)
            J[k,...,c,0]=gx;J[k,...,c,1]=gy
            _,gxx=np.gradient(gx,dv,du,edge_order=2); gyy,_=np.gradient(gy,dv,du,edge_order=2)
            gyx,_=np.gradient(gx,dv,du,edge_order=2)
            H[k,...,c,:]=np.stack((gxx,gyx,gyy),axis=-1)
    sv=np.linalg.svd(J,compute_uv=False);det=np.linalg.det(J);area=det
    initial=area[0]; ratio=np.divide(area,initial,out=np.full_like(area,np.nan),where=np.abs(initial)>EPS)
    return {"fractions":f,"jacobian":J,"hessian":H,"trace":np.trace(J,axis1=-2,axis2=-1),
            "determinant":det,"singular_value_1":sv[...,0],"singular_value_2":sv[...,1],
            "anisotropy":(sv[...,0]-sv[...,1])/(sv.sum(-1)+EPS),
            "orientation":np.arctan2(J[...,1,0]-J[...,0,1],J[...,0,0]+J[...,1,1]),
            "second_order_norm":np.linalg.norm(H,axis=(-2,-1)),"signed_area":area,
            "area_ratio":ratio,"area_min":np.nanmin(area,axis=0),"area_max":np.nanmax(area,axis=0),
            "area_final":area[-1],"minimum_abs_area":np.nanmin(np.abs(area),axis=0),
            "orientation_change":np.any(np.signbit(area)!=np.signbit(initial),axis=0),
            "signed_area_crossing":np.any(area[:-1]*area[1:]<=0,axis=0)}


def bundle_separation_history(positions, pairs, fractions=PATH_FRACTIONS):
    """Summarize separation for explicitly launch-defined ray pairs."""
    p=np.asarray(positions,float);pairs=np.asarray(pairs,np.int64);f=np.asarray(fractions,float)
    if p.ndim != 3 or p.shape[-1] != 3 or len(p) != len(f):
        raise ValueError("positions must have shape (depth, rays, 3)")
    if pairs.ndim != 2 or pairs.shape[1] != 2 or np.any(pairs<0) or np.any(pairs>=p.shape[1]):
        raise ValueError("pairs must be valid ray-index pairs")
    sep=np.linalg.norm(p[:,pairs[:,0]]-p[:,pairs[:,1]],axis=-1)
    imin=np.argmin(sep,axis=0);imax=np.argmax(sep,axis=0);q=np.arange(len(pairs))
    return {"fractions":f,"pairs":pairs,"separation":sep,"initial_separation":sep[0],
            "minimum_separation":sep[imin,q],"maximum_separation":sep[imax,q],
            "final_separation":sep[-1],"separation_rms":np.sqrt(np.mean(sep*sep,axis=0)),
            "path_location_of_min_separation":f[imin],"path_location_of_max_separation":f[imax]}


def effective_rank(bank, variance=0.95):
    x=np.nan_to_num(np.asarray(bank,float));x=x-x.mean(0,keepdims=True)
    s=np.linalg.svd(x,compute_uv=False); power=s*s
    if not power.sum(): return {"effective_rank":0.0,"variance_count":0}
    p=power/power.sum(); return {"effective_rank":float(np.exp(-np.sum(p[p>0]*np.log(p[p>0])))),
                                 "variance_count":int(np.searchsorted(np.cumsum(p),variance)+1)}


def reconstruction_r2(source, predictors):
    y=np.asarray(source,float);x=np.asarray(predictors,float)
    if y.ndim==1:y=y[:,None]
    if x.ndim==1:x=x[:,None]
    X=np.column_stack((np.ones(len(x)),x));coef=np.linalg.lstsq(X,y,rcond=None)[0];pred=X@coef
    den=np.sum((y-y.mean(0))**2,axis=0); residual=np.sum((y-pred)**2,axis=0)
    ratio=np.divide(residual,den,out=np.ones(y.shape[1]),where=den>0)
    return np.where(den>0,1-ratio,0.)
