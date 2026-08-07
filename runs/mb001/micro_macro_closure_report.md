# PBUF MB-001 micro--macro closure report

## Result: Outcome C

Existing supplied PBUF theory is insufficient to close the frozen macroscopic law. No supplied document defines a microscopic deformation variable, a quantitative matter--medium interaction, a free energy, a correlation function, or a dispersion relation from which `s(rho)`, `K`, or `G` can be calculated. The named V11 quantity `alpha_T(a)`, together with `epsilon_0,T(a)`, `k_max(a)`, and thermal rigidity, is explicitly only a candidate input in the Physics Starting Pack and is supplied there without equations, definitions, numerical values, or dimensions. This reference does not identify `alpha_T(a)` with any microscopic matter--state vertex.

Accordingly, no constitutive code was changed. Version D remains empirical, and its unchanged rerun is a reproducibility control—not validation of a new closure.

## What can be established conditionally

For the scalar laboratory variable `u`, static local balance may be written as MB-001-E01. Defining the flux in MB-001-E02 gives conservation form. If `K` and `G` are positive constants, division by `K` gives `(1 - ell^2 Laplacian)u = s/K` and MB-001-E03 follows. These statements determine the form and dimensions of a possible closure but do not derive its coefficients from PBUF.

## Research-task disposition

1. **Microscopic quantity:** absent. The laboratory's dimensionless scalar `u` is only a proxy; no supplied PBUF equation maps a microstate, strain, occupancy, or metric perturbation to it.
2. **Matter loading:** absent. “Mass loads spacetime” is qualitative; no action, conjugate force, susceptibility, or response function derives `s(rho)`. MB-001-E04 remains empirical.
3. **Effective stiffness:** absent. Thermal rigidity is named but never quantitatively related to `K`, and its dimensions are not supplied.
4. **Propagation coefficient:** absent. No gradient-energy coefficient, microscopic coupling, correlation function, or dispersion relation derives `G`.
5. **Propagation length:** only the conditional ratio `ell=sqrt(G/K)` emerges. Its value does not; MB-001-E05 remains empirical.
6. **Conservation form:** MB-001-E02 supplies the conditional continuum form.

## Equation-to-PBUF traceability matrix

| ID | Equation | Status | Origin | Unclosed relationship |
|---|---|---|---|---|
| MB-001-E01 | `K u - div(G grad u) = s(rho)` | conditional continuum balance; not a closed PBUF law | none supplied; WL-003 conditional balance | PBUF supplies no map from microscopic variables to s, K, or G |
| MB-001-E02 | `J = -G grad u; K u + div J = s(rho)` | conservation form of MB-001-E01 | algebraic rewriting of MB-001-E01 | the microscopic transported quantity and coefficient G are undefined |
| MB-001-E03 | `ell^2 = G/K` | derived algebraically, conditional on MB-001-E01 | MB-001-E01 with constant positive K and G | neither K nor G is supplied by PBUF, so ell has no predicted value |
| MB-001-E04 | `s(rho)/K = u0 (rho/rho_max)^2` | empirical Version-D identification; not derived | WL-002 Version D only | microscopic matter-medium interaction/action or response law is absent |
| MB-001-E05 | `ell = sigma_rho` | empirical Version-D identification; not derived | WL-002 Version D only | PBUF gives no correlation-length or dispersion relation |

Full physical interpretations, dimensional justifications, assumptions, and limiting cases are in `equation_traceability.csv` and `closure_equations.json`.

## Frozen-laboratory validation against archived Version D

The rerun used the archived configuration: `True`. Equation ID unchanged: `True`. All validation gates pass: `True`.

| Artifact | Maximum absolute change | Delta RMSE | Topology correlation |
|---|---:|---:|---:|
| deformation | 0 | 0 | 0.9999999999999999 |
| gradient_x | 0 | 0 | 1.0 |
| gradient_y | 0 | 0 | 1.0 |
| pbuf | 0 | 0 | 1.0 |
| residual | 0 | 0 | 1.0 |
| observation_minus_pbuf | 0 | 0 | 1.0 |
| pbuf_minus_lcdm | 0 | 0 | 1.0 |
| photon trajectories | 0 | 0 | 1.0 |

PBUF RMSE was `0.00019210569324` archived and `0.00019210569324` on rerun (change `0`). Detailed deformation, gradient, trajectory, topology, residual, RMSE, stability, and gate comparisons are in `validation.json`.

## Exact missing physical law and next milestone

PBUF must supply a dimensionally explicit microscopic energy/action or response functional `F[microstate, rho]` together with a coarse-graining definition `u=C[microstate]`. Its long-wavelength expansion must independently predict the matter-conjugate source `s(rho)`, local curvature `K`, and gradient coefficient `G`. This would make `ell=sqrt(G/K)` a prediction and test—rather than assume—the Version-D relations `s/K=u0(rho/rho_max)^2` and `ell=sigma_rho`.

Until that law is supplied, retain Version D only as the leading empirical candidate and do not update the constitutive law.
