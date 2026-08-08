# PBUF Foundation — Deformed Medium Geometry Propagation 001

## Purpose

Test whether propagation can arise from the **form of the already-supported accumulated medium deformation itself**, without coupling light directly to stiffness and without introducing GR/LCDM potential machinery or an observationally fitted coefficient.

Frozen chain:

```text
rho
 -> existing A8 transport
 -> raw c_state
 -> six-neighbor bounded-strain equilibrium
 -> accumulated medium state u
 -> deformed-medium geometry
 -> geodesic propagation
```

## Geometric candidate

Take the y=0 cross-section of the scalar accumulated state and interpret it as a deformed sheet:

```text
X(x,z) = (x,z,u(x,z))
```

The induced metric follows directly:

```text
g_ij = delta_ij + partial_i(u) partial_j(u)
```

Propagation is then computed from the geodesic equation of that induced geometry. There is **no free propagation coefficient**.

The bounded-strain stiffness remains upstream only: it determines the equilibrium shape `u`; it is not present in the propagation equation.

## Frozen inputs

- existing raw `c_state` source path;
- existing bounded-strain accumulation implementation;
- `K0=1`;
- `epsilon_max=1`;
- source radius `3.5`;
- reference mass `2.0`;
- mass ladder `0.5,1,2,4,8`;
- impact parameters `6,7,8,9,10`;
- ray interval `z=-20..20`;
- integration step `0.10`.

## Predeclared structural checks

1. zero source gives zero geometric deflection;
2. centered path gives zero transverse deflection;
3. reflection antisymmetry;
4. deflection points toward the source;
5. weak response scales approximately as `mass^1`;
6. impact response scales approximately as `b^-1`.

## Interpretation rule

This lab tests one explicit candidate: **scalar-height induced geometry**.

If it fails, do not modify or tune the supported accumulation bridge. The correct conclusion is that scalar `u` alone is insufficient to specify the propagation geometry and that a vector/tensor deformation description must be derived.

## Guardrails

No G. No GR potential decomposition. No Weyl machinery. No LCDM. No free propagation coefficient. No observational amplitude calibration. No native rescaling. No fitted/tuned K. No inserted `1/r` response. No inserted `1/b` response. No Rmax. No cosmology. No observed lensing target. No kappa/shear observation. No Quantum Engine. No Planck input.

## Valid outcomes

- `DEFORMED_MEDIUM_GEOMETRY_PROPAGATION_SUPPORTED`
- `DEFORMED_MEDIUM_GEOMETRY_PROPAGATION_PARTIAL_SUPPORT`
- `DEFORMED_MEDIUM_GEOMETRY_PROPAGATION_NOT_SUPPORTED`

Partial/null outcomes are scientifically valid and must not be repaired or tuned.

## Runner

```bash
PYTHONPATH=. python pbuf/labs/foundation/deformed_medium_geometry_propagation001.py
```

Return the current branch/HEAD, exit code, complete raw stdout/stderr, `git status --short`, and `git stash list`. Do not modify, repair, tune, rescale, reinterpret, or merge anything.
