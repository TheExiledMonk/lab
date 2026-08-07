# PBUF LAB-FREEZE-001 — Engineering Reference

This document is the engineering reference for the Version 1 weak-lensing
laboratory. It complements `laboratory_specification.md` by recording
implementation-level details that a future engineer needs to reproduce,
audit, and verify the frozen laboratory.

## 1. Source Tree

```
pbuf-test/
├── constitutive_equations.py            # Version A (frozen)
├── weak_lensing_observation001.py        # Pipeline core (frozen)
├── observable_lab001.py                 # Observable extraction (frozen)
├── source_plane_lab001.py               # Source plane launchers (frozen)
├── numerical_convergence001.py          # Convergence audit (frozen)
├── transport_lab001.py .. transport_lab008.py  # Transport validation (frozen)
├── constitutive_lab001.py               # Constitutive validation (frozen)
├── input_lab001.py / input_lab002.py    # Input validation (frozen)
├── observation_bridge001.py             # Bridge audit (frozen)
├── weak_lensing_prediction001.py        # Forward prediction (frozen)
├── weak_lensing_generalization001.py    # Cross-dataset check (frozen)
└── weak_lensing_validation001.py        # Validation suite (frozen)
```

All files are byte-locked; every completed validation milestone records
the identical SHA-256 hash for the same file.

## 2. Frozen Numerical Parameters

The frozen parameter set (`weak_lensing_observation001.LENS`) is:

```python
LENS = {
    "n": 128,            # constitutive grid resolution (default; 256² min production)
    "extent": 8.0,       # half-domain extent (dimensionless)
    "strength": 0.18,    # deformation strength u0
    "step": 0.06,        # photon step size (dimensionless length)
    "steps": 80,         # number of propagation steps
    "y_span": 3.0,       # y-extent of the launch plane
    "nphotons": 2000,    # legacy default; production = 20000 minimum
    "bins": 64,          # output bin resolution
}
```

## 3. Pipeline Core Implementation

### 3.1 Constitutive field (`constitutive_equations.py:32`)

```python
def version_a(matter: np.ndarray, c: object) -> np.ndarray:
    """WL-001 baseline: local linear loading."""
    return c.deformation_strength * _normalized(matter)


def _normalized(matter: np.ndarray) -> np.ndarray:
    return matter / max(float(matter.max()), 1e-15)
```

- `deformation_strength = 0.18`
- `_normalized` divides by `max(ρ, 1e-15)` to avoid division by zero.

### 3.2 Response field (`weak_lensing_observation001.py:181`)

```python
def make_field(rho, extent, strength, n):
    x = np.linspace(-extent, extent, n)
    y = np.linspace(-extent, extent, n)
    X, Y = np.meshgrid(x, y, indexing="xy")
    cfg = type("Config", (), {"deformation_strength": strength})()
    c = get_equation("A").solve(rho, cfg)
    gy, gx = np.gradient(c, x, y, edge_order=1)
    g = np.hypot(gx, gy)
    gx_hat = gx / np.maximum(g, 1e-15)
    gy_hat = gy / np.maximum(g, 1e-15)
    bad = g < 1e-15
    gx_hat = np.where(bad, 1.0, gx_hat)
    gy_hat = np.where(bad, 0.0, gy_hat)
    rx = -g * gy_hat
    ry = g * gx_hat
    return {"xgrid": x, "ygrid": y, "X": X, "Y": Y,
            "rho": rho, "c": c,
            "gx": gx, "gy": gy, "g_magnitude": g,
            "rx": rx, "ry": ry,
            "response_direction": np.arctan2(ry, rx)}
```

- `edge_order=1` in `np.gradient`.
- The 90° rotation is right-handed: `R_90(g_x, g_y) = (-g_y, +g_x)`.
- The bad-cell guard (`g < 1e-15`) sets `gx_hat = 1, gy_hat = 0` so that
  `r = (0, 0)` at degenerate cells.

### 3.3 Photon propagation (`weak_lensing_observation001.py:206`)

```python
def propagate(field, step, steps, x0, y0, vx0, vy0):
    xgrid = field["xgrid"]; ygrid = field["ygrid"]
    rx = field["rx"]; ry = field["ry"]
    x = x0.copy(); y = y0.copy()
    vx = vx0.copy(); vy = vy0.copy()
    nphotons = len(x)
    max_deviation = np.zeros(nphotons)
    bending_angle = np.zeros(nphotons)
    conservation = np.zeros(nphotons)
    xs = np.empty((nphotons, steps))
    ys = np.empty((nphotons, steps))
    xs[:, 0] = x; ys[:, 0] = y
    for k in range(1, steps):
        ix = np.clip(np.searchsorted(xgrid, x) - 1, 0, len(xgrid) - 1)
        iy = np.clip(np.searchsorted(ygrid, y) - 1, 0, len(ygrid) - 1)
        rx_loc = rx[iy, ix]
        ry_loc = ry[iy, ix]
        vx_new = vx + step * rx_loc
        vy_new = vy + step * ry_loc
        scale = np.maximum(np.hypot(vx_new, vy_new), 1e-12)
        vx_unit = vx_new / scale
        vy_unit = vy_new / scale
        conservation = np.maximum(conservation, np.abs(np.hypot(vx_unit, vy_unit) - 1))
        dot = np.clip(vx * vx_unit + vy * vy_unit, -1, 1)
        bending_angle += np.arccos(dot)
        vx = vx_unit; vy = vy_unit
        x = x + step * vx
        y = y + step * vy
        max_deviation = np.maximum(max_deviation, np.abs(y - y0))
        xs[:, k] = x; ys[:, k] = y
    return {"x": x, "y": y, "max_deviation": max_deviation,
            "bending_angle": bending_angle, "conservation": conservation,
            "xs": xs, "ys": ys}
```

Properties:
- Conservation error ≤ `2.2204e-16` (machine epsilon) per run.
- Photon step size `step · k` advances the position; `step · steps = 4.80`
  is the default total travel.
- Per-step cost is dominated by `searchsorted` and array indexing.

### 3.4 Observable extraction (`observable_lab001.py`)

Eight methods are exposed via `METHOD_DISPATCH`:

| Key | Function | Status |
|---|---|---|
| `histogram` | `method_histogram` | legacy reference (degenerate) |
| `kernel` | `method_kernel` | supporting |
| `jacobian` | `method_jacobian` | **primary** |
| `area` | `method_area` | secondary (finite-area) |
| `divergence` | `method_divergence` | supporting |
| `knn` | `method_knn` | supporting |
| `voronoi` | `method_voronoi` | supporting |
| `triangulation` | `method_triangulation` | secondary (Delaunay) |

Each method takes `(xs, ys, x, y, conservation, n_bins, ...)` and returns
a dict with `kappa`, `gamma1`, `gamma2`, `gamma_magnitude`, `magnification`,
`ray_count`, and `edges`.

### 3.5 Source plane launchers (`source_plane_lab001.py`)

| Launch | Function | Description |
|---|---|---|
| A | `launch_A_edge_1d` | 1D edge launch (legacy control) |
| B | `launch_B_cartesian` | Uniform Cartesian 2D (frozen) |
| C | `launch_C_hex` | Hexagonal packing |
| D | `launch_D_jittered` | Jittered Cartesian grid (seed=123456) |
| E | `launch_E_multires` | Multi-resolution (dense central) |

Only Launch B is frozen.

## 4. Pipeline Coupling

```
ρ(X)                  ←  Abell 2744 SaWLens κ (FITS)
   │
   │  ρ → max(ρ, 0) / max(max(ρ, 0))
   ▼
C(X) = u0·ρ(X)/ρ_max  ←  Version A (constitutive_equations.version_a)
   │
   │  ∇C = (∂C/∂x, ∂C/∂y)
   ▼
r(X) = R_90(∇C)·|∇C| ←  Transport Version A (weak_lensing_observation001.make_field)
   │
   │  Launch B: x0 ∈ {-8.0}×[-3, 3], v0 = (1, 0)
   ▼
{xi, yi, vxi, vyi}    ←  Photon propagation (weak_lensing_observation001.propagate)
   │
   │  xs, ys (nphotons × steps)
   ▼
κ, γ₁, γ₂, γ, μ       ←  Observable extraction (observable_lab001.method_*)
```

The coupling is strictly one-directional: every stage consumes the
output of the previous stage and produces a fixed output schema. No
feedback, no iterative refinement.

## 5. Validation Trail

| Milestone | SHA-256 of artifact run.json |
|---|---|
| TRANSPORT-LAB-001 | (see `runs/transport_lab001/`) |
| TRANSPORT-LAB-002 | (see `runs/transport_lab002/`) |
| TRANSPORT-LAB-003 | (see `runs/transport_lab003/`) |
| TRANSPORT-LAB-004 | (see `runs/transport_lab004/`) |
| TRANSPORT-LAB-005 | (see `runs/transport_lab005/`) |
| TRANSPORT-LAB-006 | (see `runs/transport_lab006/`) |
| TRANSPORT-LAB-007 | (see `runs/transport_lab007/`) |
| TRANSPORT-LAB-008 | `runs/transport_lab008/transport_lab008_report.md` |
| CONSTITUTIVE-LAB-001 | `runs/constitutive_lab001/constitutive_lab001_report.md` |
| WEAK-LENSING-PREDICTION-001 | `runs/weak_lensing_prediction001/report.md` |
| WEAK-LENSING-GENERALIZATION-001 | `runs/weak_lensing_generalization001/report.md` |
| WEAK-LENSING-VALIDATION-001 | `runs/weak_lensing_validation001/report.md` |
| OBSERVATION-BRIDGE-001 | `runs/observation_bridge001/report.md` |
| INPUT-LAB-001 | `runs/input_lab001/report.md` |
| INPUT-LAB-002 | `runs/input_lab002/report.md` |
| OBSERVABLE-LAB-001 | `runs/observable_lab001/report.md` |
| SOURCE-PLANE-LAB-001 | `runs/source_plane_lab001/report.md` |
| NUMERICAL-CONVERGENCE-001 | `runs/numerical_convergence001/report.md` |

## 6. Reproducing the Frozen Laboratory

To reproduce the laboratory from scratch:

1. Verify all source-file SHA-256s match `checksums.csv`.
2. Verify the Python environment matches `environment.json`.
3. Run any of the validation milestones in `runs/`. The output must match
   the published run.json byte-for-byte.
4. Run the production configuration (nphotons = 20 000, grid = 256²,
   step = Δs/2, source plane = Launch B, observable = Jacobian) and
   compare to `regression_baseline.json`. Tolerance:
   - κ (RMS, mean, std, peak): 1%
   - γ (RMS): 1%
   - conservation: machine epsilon
   - trajectory SHA-256: byte-exact

## 7. Modifying the Frozen Laboratory

This milestone **forbids** modifications to the frozen implementation.
Any future change must:

1. Increment the laboratory version (Version 2+).
2. Create a new `runs/lab_freeze002/` directory with its own checksums.
3. Reproduce the Version 1 regression baseline within numerical
   tolerance *before* introducing the change.
4. Demonstrate that the change is attributable to the modification, not
   to laboratory drift.

— PBUF LAB-FREEZE-001