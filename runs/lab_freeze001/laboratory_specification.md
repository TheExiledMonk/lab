# PBUF LAB-FREEZE-001 — Version 1 Weak-Lensing Laboratory Specification

## Status

Engineering freeze milestone. Documentation only. No code modifications.

This milestone establishes the **canonical reference implementation** of the
Version 1 weak-lensing laboratory. Every subsequent PBUF physics investigation
SHALL use this laboratory as the reference measurement instrument unless a
verified implementation defect is discovered.

---

## 1. Frozen Physics

The laboratory implements the following exact, validated physical model.
Every component below is a result of the completed validation programme.

### 1.1 Constitutive (Version A)

The constitutive field `C(X)` is defined by:

```
C(X) = u0 · ρ(X) / ρ_max
```

with

| Symbol | Meaning | Frozen value |
|---|---|---|
| `u0` | Deformation strength | `0.18` |
| `ρ(X)` | Matter input field | per cluster (Section 2) |
| `ρ_max` | Maximum of the input field | computed per input |

Reference implementation: `constitutive_equations.py`, function `version_a`
(SHA-256: `e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f`).

### 1.2 Transport (Version A)

The frozen transport law is the **neighbour-to-neighbour** propagation rule
with the following components:

- **Response direction:** 90° transverse (right-handed rotation of `∇C`).
- **Response magnitude:** linear amplitude `A = |∇C|`.
- **Update rule:** direct vector addition `v_new = v + step · r`.
- **Normalisation:** per-step unit-speed renormalisation `|v| := 1`.
- **Propagation:** velocity stored on a uniform Cartesian grid;
  photon position updated by `x := x + step · v` per step.
- **Conservation:** maximum drift of `|v|` from 1 is `2.2204e-16`
  (machine epsilon) for every frozen run.

These five components are the *exact* law that was frozen by TRANSPORT-LAB-001
through TRANSPORT-LAB-008 and CONSTITUTIVE-LAB-001.

### 1.3 Source Plane (Configuration B)

The frozen photon launch is a **uniform two-dimensional Cartesian grid**:

| Property | Value |
|---|---|
| Plane | `x ∈ [-extent, -extent + y_span]`, `y ∈ [-y_span, y_span]` |
| Density | uniform |
| Launch velocity | `(1, 0)` for every photon |
| Geometric seed | none (deterministic; no jitter) |
| Random seed | not applicable (Launch B is deterministic) |

Launch B is selected over the other four source-plane candidates by the
SOURCE-PLANE-LAB-001 ranking (most stable, fewest degeneracies).

### 1.4 Observable Extraction

The frozen observable extraction scheme is a hierarchy of methods implemented
in `observable_lab001.py` (SHA-256: `2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132`).

| Rank | Method | Use |
|---|---|---|
| **Primary** | `jacobian` | the frozen reference observable |
| Secondary | `area` (finite-area) | cross-check |
| Secondary | `triangulation` (Delaunay) | cross-check |
| Supporting | `kernel` | KDE density |
| Supporting | `knn` | k-nearest-neighbour density |
| Supporting | `voronoi` | Voronoi-area method |
| Supporting | `divergence` | displacement divergence |
| Legacy | `histogram` | retained as historical reference; **suppressed** (degenerate) |

The primary observable is the **Jacobian** method (linear fit per bin).

### 1.5 Frozen Numerical Configuration

#### Minimum validated production configuration

| Parameter | Value |
|---|---|
| Photon count | 20 000 |
| Constitutive grid | 256² |
| Step size | Δs / 2 |
| Source plane | Cartesian 2D |
| Observable | Jacobian |

Validated by NUMERICAL-CONVERGENCE-001 (Outcome A): at these settings the
Jacobian κ field changes by less than 1% relative to the next finer
configuration.

#### Higher-accuracy configuration

| Parameter | Value |
|---|---|
| Photon count | 50 000 |
| Constitutive grid | 512² |
| Step size | Δs / 4 |
| Source plane | Cartesian 2D |
| Observable | Jacobian |

Validated by NUMERICAL-CONVERGENCE-001 (Outcome A): at these settings the
Jacobian κ field changes by less than 0.1% relative to the next finer
configuration.

---

## 2. Inputs

### 2.1 Matter field `ρ(X)`

The laboratory accepts a dimensionless matter field on the pipeline grid
`X ∈ [-extent, extent]²` with `extent = 8.0` and grid size `n = 128`
(= 16 384 cells).

For validation against published weak-lensing maps, the canonical matter
input is:

```
ρ(X) = max(κ_observed(X), 0) / max(max(κ_observed(X), 0))
```

from the Frontiers-Fields SaWLens reconstruction of the Abell 2744 cluster
(file: `PBUF_benchmark/WL-001_Abell2744/hlsp_frontier_model_abell2744_merten_v1_kappa.fits`).

Per OBSERVATION-BRIDGE-001, this substitution is an *approximation* (the
published κ is a `Σ/Σ_crit` map, not a mass density) but it is the frozen
control and is reproduced byte-for-byte by every subsequent run.

### 2.2 Constitutive field `C(X)`

Generated from `ρ(X)` via Version A (Section 1.1).

### 2.3 Response field `r(X)`

```
g_x = ∂C/∂x ,    g_y = ∂C/∂y
|r| = sqrt(g_x² + g_y²)
|r_hat_x| = g_x / max(|r|, 1e-15)
|r_hat_y| = g_y / max(|r|, 1e-15)
r_x = -|r| · r_hat_y ,    r_y = +|r| · r_hat_x
```

i.e. `r = R_90(∇C) · |∇C|` where `R_90` is the right-handed 90° rotation.

### 2.4 Source plane

Launch B (Cartesian 2D; Section 1.3).

### 2.5 Pipeline parameters (frozen)

| Parameter | Symbol | Frozen value |
|---|---|---|
| Grid resolution | `n` | 128 |
| Domain half-extent | `extent` | 8.0 |
| Deformation strength | `u0` | 0.18 |
| Launch y-span | `y_span` | 3.0 |
| Step size (Δs) | `step` | 0.06 |
| Number of steps | `steps` | 80 |
| Total photon travel | `step · steps` | 4.80 |
| Output bin resolution | `bins` | 64 |
| Photon count (default) | `nphotons` | 2000 (legacy) / 20 000 (min) / 50 000 (high) |

---

## 3. Processing

### 3.1 Constitutive generation

```
def C_field(rho):
    return u0 * rho / max(rho.max(), 1e-15)
```

Function: `constitutive_equations.version_a(matter, cfg)`.

### 3.2 Neighbour transport

```
def response(C_field, xgrid, ygrid):
    gy, gx = np.gradient(C_field, xgrid, ygrid, edge_order=1)
    g = np.hypot(gx, gy)
    gx_hat = np.where(g < 1e-15, 1.0, gx / np.maximum(g, 1e-15))
    gy_hat = np.where(g < 1e-15, 0.0, gy / np.maximum(g, 1e-15))
    rx = -g * gy_hat
    ry = +g * gx_hat
    return rx, ry
```

Function: `weak_lensing_observation001.make_field`, lines 181-203.

### 3.3 Photon propagation

```
def propagate(rx, ry, xgrid, ygrid, x0, y0, vx0, vy0, step, steps):
    for k in 1..steps-1:
        ix = clip(searchsorted(xgrid, x) - 1, 0, N-1)
        iy = clip(searchsorted(ygrid, y) - 1, 0, N-1)
        rx_loc = rx[iy, ix]; ry_loc = ry[iy, ix]
        vx_new = vx + step * rx_loc
        vy_new = vy + step * ry_loc
        scale = max(hypot(vx_new, vy_new), 1e-12)
        vx_unit = vx_new / scale; vy_unit = vy_new / scale
        # accumulate bending angle
        dot = clip(vx*vx_unit + vy*vy_unit, -1, 1)
        bending_angle += arccos(dot)
        vx, vy = vx_unit, vy_unit
        x += step * vx; y += step * vy
        store (xs[:,k], ys[:,k])
    return xs, ys, max_deviation, bending_angle, conservation
```

Function: `weak_lensing_observation001.propagate`, lines 206-245.

Properties:
- Conservation error `||v| - 1| ≤ 2.2204e-16` for every step (machine epsilon).
- Photon total travel = `step · steps = 4.80` (default) or `2.40`
  (step/2) or `1.20` (step/4).
- Runtime scales linearly in `nphotons × steps`.

### 3.4 Observable extraction

#### 3.4.1 Jacobian (primary)

A ray-bundle Jacobian is constructed per output bin from the photon
displacement field `(xf - x0, yf - y0)`. The Jacobian matrix is computed by
a linear fit, and the convergence and shear are recovered from its trace and
off-diagonal components.

Function: `observable_lab001.method_jacobian` (SHA-256: `2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132`).

#### 3.4.2 Secondary methods

- `area` — finite-area distortion (initial vs final spread).
- `triangulation` — Delaunay triangulation area method.
- `kernel` — Gaussian KDE (Scott bandwidth).
- `knn` — adaptive k-nearest-neighbour density.
- `voronoi` — Voronoi area method.
- `divergence` — displacement divergence (`∇·d`).
- `histogram` — histogram occupancy (legacy; degenerate by construction).

---

## 4. Outputs

The laboratory produces the following observables:

| Observable | Symbol | Frozen definition |
|---|---|---|
| Convergence | κ | `0.5 * (N_final/N_initial - 1)` for `N_initial > 0`; or Jacobian-trace equivalent for the `jacobian` method |
| Shear component 1 | γ₁ | `0.5 * (∂d_x/∂x - ∂d_y/∂y)` |
| Shear component 2 | γ₂ | `0.5 * (∂d_x/∂y + ∂d_y/∂x)` |
| Shear magnitude | γ | `sqrt(γ₁² + γ₂²)` |
| Magnification | μ | `1 / ((1 - κ)² - |γ|²)` |
| Deflection | (d_x, d_y) | `(xf - x0, yf - y0)` |
| Photon trajectories | (xs, ys) | full `(nphotons × steps)` array |

All observables are evaluated on the output bin grid `bins × bins = 64 × 64`
mapped onto the same `[-8, 8]²` domain as the constitutive field.

---

## 5. Frozen Algorithms

The canonical algorithm registry is `algorithm_registry.csv`. Every stage
records:

- algorithm name
- implementation file
- function
- SHA-256 checksum
- dependencies

Summary:

| Stage | Algorithm | File | Function | SHA-256 |
|---|---|---|---|---|
| Constitutive | Version A | `constitutive_equations.py` | `version_a` | `e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f` |
| Response | 90° transverse | `weak_lensing_observation001.py` | `make_field` | `a5c3632fec9adfc2659d5c283d07c599db6db7edb2e34e4aebd84e35434642bc` |
| Transport | direct addition + renormalisation | `weak_lensing_observation001.py` | `propagate` | `a5c3632fec9adfc2659d5c283d07c599db6db7edb2e34e4aebd84e35434642bc` |
| Source plane | Launch B (Cartesian 2D) | `source_plane_lab001.py` | `launch_B_cartesian` | `efa9d74924cb61a3b48a69fa075055512d86391d03194be342597420bc353de4` |
| Observable (primary) | Jacobian | `observable_lab001.py` | `method_jacobian` | `2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132` |
| Observable (secondary) | Area | `observable_lab001.py` | `method_area` | `2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132` |
| Observable (secondary) | Triangulation (Delaunay) | `observable_lab001.py` | `method_triangulation` | `2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132` |
| Observable (legacy) | Histogram | `observable_lab001.py` | `method_histogram` | `2867c0bf94fabe3fbba0264d5d272ecadbdf45fabb341a6a8fe77972fbaec132` |

This set of algorithms is the **canonical reference implementation**.
Reproducibility is verified by identical SHA-256 across every validation
milestone (see `checksums.csv`).

---

## 6. Validation Summary

### 6.1 Transport (TRANSPORT-LAB-001 .. TRANSPORT-LAB-008)

The transport law is uniquely selected: **neighbour-to-neighbour, 90°
transverse response, linear magnitude, direct addition, per-step
renormalisation**. Multiple equivalent variants were eliminated:

| Lab | Variable | Outcome |
|---|---|---|
| LAB-001 | propagation kernel | equivalent (5 candidates within 5%) |
| LAB-002 | interpolation scheme | equivalent (bilinear, nearest) |
| LAB-003 | boundary condition | equivalent (clip, wrap) |
| LAB-004 | gradient scheme | equivalent (gradient, finite-diff) |
| LAB-005 | launch position | equivalent (left, right) |
| LAB-006 | normalisation | equivalent (rescale, replace) |
| LAB-007 | update rule | equivalent (direct, projected) |
| LAB-008 | response amplitude law | **linear** selected (Cand 1 / Cand 4 within 7%) |

The amplitude-law ablation (LAB-008) is the decisive selection: only the
linear law `A = |∇C|` (Cand 1) and the logarithmic law `A = log(1 + |∇C|)`
(Cand 4) are statistically equivalent. All other candidates (sqrt, quadratic,
saturating, piecewise, threshold) produce bending that deviates from the
control by −100% to +3 161%. The frozen choice is the linear law.

### 6.2 Observable (OBSERVABLE-LAB-001)

| Finding | Status |
|---|---|
| `histogram` method suppresses information | **YES** (Q5) — degenerate by construction |
| `jacobian` method is physically meaningful | **VALIDATED** (std(κ) = 0.15, finite pixel count) |
| `area` (finite-area) method is meaningful in 2D | **VALIDATED** (std(κ) = 0.14 in Launch B) |
| `triangulation` (Delaunay) method is meaningful in 2D | **VALIDATED** (std(κ) = 0.12 in Launch B) |

The histogram method is retained only as a legacy reference; it produces a
constant κ = -0.5 in every frozen trajectory set because photons leave the
launch column entirely within `step · steps = 4.80` units.

### 6.3 Source plane (SOURCE-PLANE-LAB-001)

| Finding | Status |
|---|---|
| 1D edge launch (Launch A) is degenerate | **CONFIRMED** — 5 of 8 methods degenerate |
| 2D Cartesian launch (Launch B) removes degeneracy | **CONFIRMED** — 3 of 3 previously-degenerate methods become valid |

Launch ranking (most stable first): **Launch B**, Launch D, Launch C, Launch E, Launch A.

### 6.4 Numerical convergence (NUMERICAL-CONVERGENCE-001)

| Group | Converging? | Convergence order (Jacobian, field) |
|---|---|---|
| A (photon count) | YES | p_obs = +0.98 (R²=0.99) |
| B (grid) | YES | p_obs = +0.59 (R²=0.94) |
| C (step) | YES (already converged at Δs = 0.06) | p_obs = +0.85 (R²=0.89) |
| D (domain) | N/A (consistency, not refinement) | p_obs = +2.24 |
| E (kNN Jacobian neighbourhood) | YES | p_obs = +1.14 |

Converged solution (Jacobian method, RMS):

- κ = 0.134 ± 0.001 at nphotons = 100 000
- |γ| = 0.084 ± 0.001 at nphotons = 100 000

Uncertainty ranking (largest first): **E, D, A, B, C**.

### 6.5 Machine-precision conservation

Maximum conservation error across **every** frozen run is
`2.220446049250313e-16` (= 2²⁻⁵², IEEE-754 double-precision machine epsilon).
No run exceeds this.

---

## 7. Laboratory Limits

The frozen laboratory has the following experimentally established
limitations. (No speculation.)

1. **No cosmological bridge.** κ is a local photon-density distortion, not
   `Σ/Σ_crit`. There is no dependence on lens redshift, source redshift, or
   the distance ratio `D_ls/D_s`. (OBSERVATION-BRIDGE-001)
2. **Dimensionless constitutive field.** `C(X)` is dimensionless; no
   physical length, mass, or energy is attached.
3. **No Σ_crit.** The laboratory never computes a critical surface-mass
   density.
4. **No source-redshift dependence.** Photons are launched with velocity
   `(1, 0)`; there is no source-plane redshift.
5. **No physical angular scale.** The pipeline grid is dimensionless
   Cartesian on `[-8, 8]²`; no WCS, no RA/Dec, no arcsec-per-pixel.
6. **Benchmark comparison pending.** Across all 5 clusters and 14 input
   candidates (INPUT-LAB-001), Pearson(κ) is `nan` because the histogram
   method produces a constant κ = -0.5. γ Pearson is below 0.06 in every
   case. No quantitative agreement with published observations has been
   demonstrated.
7. **Version A physics only.** Versions B, C, D in `constitutive_equations.py`
   are NOT frozen. The frozen laboratory uses Version A exclusively.
8. **Mirror asymmetry is a property of the transport.** The 90° transverse
   response has a definite handedness, so the mirror-symmetry test fails
   by construction. (WEAK-LENSING-VALIDATION-001)
9. **Histogram method is degenerate.** The legacy `histogram` extraction
   produces a constant κ = -0.5 by construction. It is retained as
   historical reference only.
10. **No domain rescaling.** Domain-size variations (Group D in
    NUMERICAL-CONVERGENCE-001) rescale the apparent cluster size, so RMS κ
    changes monotonically with domain extent. This is a consistency check,
    not a refinement study.

See `known_limitations.md` for the full discussion.

---

## 8. Regression Baseline

The regression baseline (see `regression_baseline.json`) records the exact
predicted observables for the **20 000-photon production configuration**
(Launch B, grid 256², Δs/2, Jacobian method). Future versions MUST reproduce:

| Quantity | Frozen value (20k) | Frozen value (50k) |
|---|---|---|
| RMS κ (Jacobian) | `1.3522e-01` | `1.3430e-01` |
| RMS γ (Jacobian) | `8.6084e-02` | `8.4649e-02` |
| Mean κ (Jacobian) | `-4.1174e-03` | `-4.8573e-03` |
| std κ (Jacobian) | `1.3515e-01` | `1.3421e-01` |
| Peak \|κ\| | `4.5195e-01` | `4.7014e-01` |
| Conservation | `2.2204e-16` | `2.2204e-16` |
| Trajectory SHA-256 | `80d8fe47bd0d45672e16c659fc1b095114106eacd6a52634fccde1fc81b41b3f` | n/a (different seed) |

Numerical tolerance: **1%** on RMS κ and RMS γ; **machine epsilon** on
conservation; **byte-exact** on the trajectory checksum.

---

## 9. Reproducibility

### 9.1 Python version and packages

See `environment.json`. The frozen pipeline was validated on:

- Python 3.14.x
- numpy ≥ 2.0
- matplotlib ≥ 3.10
- astropy ≥ 6.0 (FITS reading only)
- scipy ≥ 1.13 (`map_coordinates`)

### 9.2 Random seeds

| Component | Seed |
|---|---|
| Launch B (Cartesian 2D) | none (deterministic) |
| Launch D (jittered) | `123456` (recorded; not used by the frozen launch) |
| Photon initial position | deterministic from launch geometry |
| Photon initial velocity | `(1, 0)` (deterministic) |

### 9.3 Numerical precision

All laboratory computations use IEEE-754 double precision (`float64`).
Conservation error = `2.220446049250313e-16` (= 2⁻⁵²) for every frozen run.

### 9.4 SHA-256 of frozen sources

See `checksums.csv`. All frozen source files reproduce the same hash
across every completed validation milestone.

---

## 10. Performance

Measured by NUMERICAL-CONVERGENCE-001 on the Abell 2744 cluster at the
default grid (n=128), source-plane Launch B, Jacobian method.

### 10.1 20 000-photon production configuration

(Full 20 000-photon run; propagated + extracted by all 8 methods.)

| Quantity | Value |
|---|---|
| Propagation runtime | `7.10e-02 s` |
| Jacobian extraction runtime | `2.06e-01 s` |
| Total runtime | `5.80e-01 s` |
| Peak `tracemalloc` | `3.75e+07 B` ≈ 35.8 MB |
| Peak RSS | `2.08e+05 kB` ≈ 203 MB |

### 10.2 50 000-photon production configuration

| Quantity | Value |
|---|---|
| Propagation runtime | `1.78e-01 s` |
| Jacobian extraction runtime | `3.13e-01 s` |
| Total runtime | `1.10 s` |
| Peak `tracemalloc` | `9.34e+07 B` ≈ 89.1 MB |
| Peak RSS | `3.30e+05 kB` ≈ 323 MB |

### 10.3 Scaling

Runtime scales linearly with `nphotons` (propagation) and with `nphotons²`
for the n_neighbours-dependent kNN Jacobian (Group E audit only; not the
primary observable).

Memory scales linearly with `nphotons × steps` (trajectory array) and with
`n²` (constitutive grid; production uses n=256, ≈ 65 536 cells = 0.5 MB).

---

## 11. Closing Statement

The Version 1 weak-lensing laboratory has completed the validation
programme:

- TRANSPORT-LAB-001 through TRANSPORT-LAB-008
- CONSTITUTIVE-LAB-001
- WEAK-LENSING-PREDICTION-001
- WEAK-LENSING-GENERALIZATION-001
- WEAK-LENSING-VALIDATION-001
- OBSERVATION-BRIDGE-001
- INPUT-LAB-001
- INPUT-LAB-002
- OBSERVABLE-LAB-001
- SOURCE-PLANE-LAB-001
- NUMERICAL-CONVERGENCE-001

It is hereby **frozen** as the canonical reference implementation for all
subsequent PBUF development. Any future change to the constitutive law,
transport, source-plane geometry, observable extraction, or numerical
algorithm MUST be introduced in a new (Version 2) laboratory; it SHALL NOT
modify this frozen laboratory.

— PBUF LAB-FREEZE-001