# PBUF CORE-001 microscopic state and coarse graining

## Result

CORE-001 is complete as a **conditional microscopic model**. It supplies explicit inputs with which MB-001 can be revisited. The construction formalizes the requested premise; it does not establish that spacetime actually has this microstructure, derive the number three, or explain the numerical value `1/137`.

## Microscopic state and dynamics

Put a dimensionless state `q_i=(q_i^1,q_i^2,q_i^3) in R^3` at each site of an isotropic lattice of spacing `a`. The three components are the three stipulated fundamental degrees of freedom. Let `e` be a unit vector in this internal state space selected by the matter coupling, `eta_i=rho_i/rho_*`, and `g_dev=1/137`. Define

`F = epsilon_* sum_i [kappa_0 |q_i|^2/2 + kappa_1 sum_<ij> |q_j-q_i|^2/2 - g_dev eta_i e.q_i]`.

Here `epsilon_*` is an energy, while `q`, `kappa_0`, `kappa_1`, `eta`, and `g_dev` are dimensionless. The matter vertex is normalized directly by the PBUF coupling `g_dev`; there is no auxiliary coupling or source multiplier. The local evolution is overdamped relaxation

`tau dq_i/dt = -d(F/epsilon_*)/dq_i + xi_i`,

where `tau` is a time and optional noise has zero mean. Static CORE-001 uses `xi=0`. Positive `kappa_0` and `kappa_1` make the unloaded quadratic energy bounded below. Matter perturbs the state through the explicitly assumed linear interaction; no quadratic density response is derived here.

## Coarse-graining map

Choose a nonnegative radial kernel `W_L` of width `L`, normalized so `sum_i a^d W_L(x-x_i)=1`, and define

`u(x) = C_L[q](x) = e . sum_i a^d W_L(x-x_i) q_i`.

Thus `u` is dimensionless. Normalization preserves uniform states and, with periodic or vanishing-flux boundaries, preserves the projected spatial mean. Radial `W_L` preserves spatial rotations; simultaneous rotation of `q` and `e` preserves internal-basis covariance. The source selects `e`, so full internal `O(3)` symmetry is explicitly broken by matter, not accidentally by coarse graining.

## Continuum limit and connection to MB-001

Assume `a << L << L_macro`, fields vary little between sites, transverse components have relaxed, fluctuations have finite short-range correlations, and boundary terms vanish. Taylor expansion of neighbor differences gives

`F_cont = integral [K u^2/2 + G |grad u|^2/2 - s(rho)u] d^d x`,

with

- `K = epsilon_* kappa_0/a^d` (energy per volume),
- `G = epsilon_* kappa_1 a^(2-d)` (energy per volume times length squared),
- `s(rho) = epsilon_* g_dev (rho/rho_*)/a^d` (energy per volume),
- `ell=sqrt(G/K)=a sqrt(kappa_1/kappa_0)` (length).

Stationarity gives `K u-div(G grad u)=s(rho)`, exactly the conditional WL-003/MB-001 continuum form. These equations show how microscopic parameters *could* supply `s`, `K`, `G`, and `ell`; they do not determine their values. A nonzero fixed `ell` as `a->0` requires the renormalized ratio `kappa_1/kappa_0 ~ (ell/a)^2`.

Limits are explicit: `kappa_1->0` gives independent local response, `kappa_0->infinity` suppresses deformation, `rho->0` gives the unloaded state, and `L/a->infinity` with `L/L_macro->0` suppresses microscopic fluctuations without erasing macroscopic variation.

## Traceability matrix

| ID | Definition/equation | Status | Assumption or derivation boundary |
|---|---|---|---|
| CORE-001-A01 | `three real components q_i in R^3` | working premise | The number three is stipulated by CORE-001; no microscopic derivation is claimed. |
| CORE-001-A02 | `g_dev = 1/137` | working premise | Dimensionless matter--state coupling scale; its value and interpretation are assumed. |
| CORE-001-A03 | `periodic isotropic lattice with spacing a` | modeling assumption | Provides a regulator and a controlled long-wavelength expansion. |
| CORE-001-E01 | `F=epsilon_* sum_i [kappa_0|q_i|^2/2 + kappa_1 sum_<ij>|q_j-q_i|^2/2 - g_dev eta_i e.q_i]` | corrected microscopic free energy | Defines recovery, nearest-neighbor transmission, and direct PBUF matter loading. |
| CORE-001-E02 | `tau dq_i/dt = -d(F/epsilon_*)/dq_i + xi_i` | defined local evolution | Overdamped relaxation; zero-mean noise xi may be omitted in the static limit. |
| CORE-001-E03 | `C_L[q](x)=e . sum_i a^d W_L(x-x_i) q_i` | defined coarse graining | W_L is nonnegative, rotationally symmetric, normalized, and a << L. |
| CORE-001-E04 | `u(x)=C_L[q](x)` | definition | Dimensionless scalar deformation accepted by the existing continuum interface. |
| CORE-001-E05 | `K=epsilon_* kappa_0/a^d` | conditional derivation | Local stiffness in the aligned, long-wavelength sector. |
| CORE-001-E06 | `G=epsilon_* kappa_1 a^(2-d)` | conditional derivation | Gradient stiffness from the nearest-neighbor term; convention-dependent O(1) factors are absorbed in kappa_1. |
| CORE-001-E07 | `s(rho)=epsilon_* g_dev eta/a^d` | corrected conditional derivation | eta=rho/rho_* is dimensionless; the direct linear source follows from the assumed interaction. |
| CORE-001-E08 | `ell=sqrt(G/K)=a sqrt(kappa_1/kappa_0)` | conditional derivation | Propagation length; a finite continuum ell requires parameter scaling under a->0. |
| CORE-001-E09 | `K u-div(G grad u)=s(rho)` | macroscopic limit | WL-003 form after alignment, scale separation, isotropy, and static relaxation. |

## Consistency checks

All executable checks pass: **True**. The normalized kernel preserves constants and the periodic mean; spatial and internal covariance errors are below floating-point tolerances. The discrete Laplacian errors at `n=32,64,128` are `[0.027174656444203145, 0.006818738140463345, 0.001706255665739933]`, demonstrating the expected long-wavelength convergence. These are mathematical consistency checks, not observational validation.

## Remaining theoretical gaps

1. PBUF has not derived why the microscopic state has exactly three components or why `g_dev=1/137` governs this coupling.
2. The lattice/regulator, energy scale `epsilon_*`, spacing `a`, relaxation time `tau`, couplings, noise statistics, and internal direction `e` are not predicted.
3. The assumed linear matter interaction is not derived from a covariant action; it does not derive Version D's quadratic normalized source.
4. A relativistic evolution law, tensorial metric map, causality, gauge behavior, and coupling to stress-energy rather than a static density are unresolved.
5. Renormalization and universality beyond the quadratic, isotropic, short-correlation regime remain to be shown.
6. No equality between `ell` and the observed baryonic width is implied, and no weak-lensing parameter has been fitted or changed.
