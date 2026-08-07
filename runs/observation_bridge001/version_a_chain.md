# Version A Pipeline - Stage-by-Stage Audit

All numerical parameters are frozen (identical to weak_lensing_prediction001.py and weak_lensing_generalization001.py):

- `n_grid = 128`
- `extent = 8.0`
- `strength = 0.18`
- `step = 0.06`
- `steps = 80`
- `nphotons = 2000`
- `bins = 64`

## Matter density ρ(X)

- Mathematical quantity: `ρ : R^2 -> R`
- Physical meaning: Input scalar field declared to be a non-negative dimensionless matter density on the pipeline grid.
- Units: dimensionless (ρ_max = 1 after normalisation)
- Assumptions: ρ >= 0; field supplied externally; no cosmological, redshift, or physical-unit content; no baryonic/dark partition.

## Constitutive field C(X) = 0.18 · ρ(X) / ρ_max

- Mathematical quantity: `C : R^2 -> R,    C = 0.18 * ρ / ρ_max`
- Physical meaning: Local linear scalar proxy of the medium's deformation by matter. C is the output of Version A's constitutive equation. No physical dimension is attached; C exists only as a scalar intermediate field.
- Units: dimensionless (bounded by 0.18)
- Assumptions: Local, linear, isotropic response; no propagation; the coefficient 0.18 is the frozen 'deformation strength' of Version A.

## Gradient field ∇C

- Mathematical quantity: `∇C = (gx, gy) = (∂C/∂x, ∂C/∂y)`
- Physical meaning: Gradient of the constitutive scalar field. Has the magnitude of a response per dimensionless length.
- Units: dimensionless per dimensionless length
- Assumptions: Finite differences with edge_order=1; the grid [-8, 8] x [-8, 8] carries no physical length scale.

## Response field r = (rx, ry)

- Mathematical quantity: `A = |∇C|;    r_x = -A * (∂C/∂y)/A;    r_y = A * (∂C/∂x)/A  (90-degree transverse rotation of the unit gradient)`
- Physical meaning: Vector field representing the local transverse response of the medium. The 90-degree rotation is the frozen 'transport' choice; the response amplitude equals the gradient magnitude.
- Units: dimensionless per dimensionless length (input to velocity update)
- Assumptions: Neighbour-to-neighbour coupling, direct addition, instantaneous renormalisation; no retardation, no falloff, no medium rigidity.

## Photon propagation

- Mathematical quantity: `v_{k+1} = (v_k + step * r) / |v_k + step * r|;    x_{k+1} = x_k + step * v_{k+1}`
- Physical meaning: Iterative ray-tracing through the response field. Photons start at x = -8 with v = (1, 0); the pipeline runs for steps = 80 iterations with step = 0.06, so the maximum propagation distance is 4.8 dimensionless units (≪ 2*extent = 16).
- Units: step in dimensionless length units; velocity renormalised to unit speed per step
- Assumptions: Frozen parameters n_grid = 128, extent = 8, strength = 0.18, step = 0.06, steps = 80, nphotons = 2000, bins = 64; identical to weak_lensing_prediction001 and weak_lensing_generalization001.

## Predicted convergence κ

- Mathematical quantity: `κ(x,y) = 0.5 * (N_final(x,y) / N_initial(x,y) - 1)`
- Physical meaning: Local photon-count ratio in a 64 x 64 histogram of (x_f, y_f) after propagation versus the initial (x_0, y_0) histogram. The pipeline only produces finite κ on bins where N_initial > 0.
- Units: dimensionless
- Assumptions: Bins with N_initial = 0 are filled with NaN; bins where photons left the initial x = -8 column but never returned show the constant value -0.5. No physical unit is attached; κ here is not the surface-mass-density-to-critical-density ratio.

## Predicted shear γ_1, γ_2

- Mathematical quantity: `α = mean photon displacement in each bin;    γ_1 = 0.5 * (∂α_x/∂x - ∂α_y/∂y);    γ_2 = 0.5 * (∂α_x/∂y + ∂α_y/∂x)`
- Physical meaning: Components of the 2 x 2 Jacobian of the photon-displacement field, evaluated on the same 64 x 64 grid. No convergence-to-shear correction is applied; no reduced-shear division by (1 - κ) is performed.
- Units: dimensionless
- Assumptions: NaN-filled bins propagate into α and hence γ; the reported fields contain a mixture of finite values (where photons landed) and NaN/0 (elsewhere).

## Predicted magnification μ

- Mathematical quantity: `μ = 1 / ((1 - κ)^2 - |γ|^2)`
- Physical meaning: Magnification derived from the standard lensing identity, using the predicted κ and γ. NaN wherever the denominator is non-positive or one of the inputs is undefined.
- Units: dimensionless
- Assumptions: Same as κ and γ above; no cosmological distance factor.

