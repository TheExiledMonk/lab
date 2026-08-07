# PBUF PHOTON-001 — Microscopic-to-photon coupling derivation

## Result: exact coupling remains underdetermined

The established PBUF chain defines the coarse scalar deformation `u=C_L[q]`, but it contains no electromagnetic action, effective metric map, refractive response, or symmetry principle that fixes how photons read `u`. Therefore no unique photon coupling—and no numerical coupling coefficient—can be derived without adding a theoretical postulate.

The strongest explicit conditional result is the isotropic geometrical-optics family `S_ray=E0 integral n(u) ds`. Its Euler-Lagrange equation is `dt/ds=(I-tt^T) grad ln n`. In the small-deformation limit, `n=1+beta u+O(u^2)` and ray curvature is `beta(I-tt^T)grad u`. PBUF does not determine `beta=(dn/du)|_0`; in particular, the matter-to-microstate premise `g_dev=1/137` does not establish a photon coupling.

No lensing run was executed, no parameter was fitted, and no propagator or constitutive ranking was changed.

## Conditional derivation

Locality, stationarity, spatial isotropy, parity, losslessness, and absence of dispersion/birefringence reduce a scalar optical response to a positive index `n(u)`. Varying the Fermat functional with fixed endpoints gives `d(n t)/ds=grad n`. Taking the component perpendicular to the unit tangent yields the curvature equation. Consequently deformation itself controls phase, while its transverse gradient controls bending; an accumulated phase is an output of the same response rather than an independent fundamental coupling.

An effective optical metric `ds_eff^2=-c^2dt^2/n(u)^2+dx^2`, up to conformal freedom for null paths, reproduces this ray law. It is an equivalent parametrization under the stated static isotropic assumptions, not a microscopic PBUF derivation.

## Coupling assumptions catalogue

| ID | Assumption | Status | Needed for | If false |
|---|---|---|---|---|
| P001-A01 | The CORE-001 coarse field u is the photon-accessible scalar. | conditional | scalar optical reduction | Photon coupling may depend on other projections, tensor modes, or derivatives. |
| P001-A02 | Geometrical optics applies: wavelength is short compared with the variation scale of u. | standard approximation | ray action and path equation | A wave equation and diffraction must replace rays. |
| P001-A03 | The medium is local, static, isotropic, parity even, nondispersive, and nonbirefringent in its rest frame. | conditional symmetry choice | one scalar refractive index n(u) | Direction, frequency, polarization, or history dependent optical tensors are allowed. |
| P001-A04 | Photon number/frequency is conserved in the static background and propagation is lossless. | conditional conservation requirement | real Fermat functional | Complex response and absorption/emission terms are required. |
| P001-A05 | n(u) is differentiable near the unloaded state, with n(0)=1. | normalization plus regularity | small-deformation expansion | No linear weak-field limit follows. |
| P001-A06 | A covariant completion, if used, provides an effective metric whose null rays reproduce the same optical index in the static isotropic limit. | conditional | effective-metric interpretation | Fermat optics remains phenomenological rather than covariantly derived. |

## Equation traceability matrix

| ID | Equation | Status | Origin | Assumptions | Meaning/boundary |
|---|---|---|---|---|---|
| P001-E01 | u(x)=C_L[q](x) | established definition | CORE-001-E03/E04 | CORE-001 scale separation and scalar projection | Supplies a possible photon input, not a coupling law. |
| P001-E02 | S_ray[x]=E0 integral n(u(x)) \|dx/dlambda\| dlambda | conditional general form | spatial isotropy, locality, stationarity, losslessness | P001-A01--A05 | The scalar field can influence rays only through an undetermined optical response n(u). |
| P001-E03 | d(n t)/ds=grad n | derived from P001-E02 | Euler-Lagrange variation | arc length s; unit tangent t | Only transverse index gradients bend a ray. |
| P001-E04 | dt/ds=(I-t t^T) grad ln n(u) | derived from P001-E03 | projection perpendicular to t | n>0 | A uniform deformation changes optical phase but not the ray path. |
| P001-E05 | n(u)=1+beta u+O(u^2), beta=(dn/du)\|_0 | conditional expansion; coefficient missing | P001-A05 | \|u\| small and differentiability | beta is a dimensionless photon-coupling response not fixed by g_dev=1/137. |
| P001-E06 | dt/ds=beta (I-t t^T) grad u+O(u grad u) | derived from P001-E04/E05 | first-order expansion | weak deformation | The deformation gradient drives bending only after beta is supplied. |
| P001-E07 | Delta Phi=(E0/hbar c) integral [n(u)-1] ds | derived conditionally from P001-E02 | eikonal phase | coherent monochromatic wave and geometrical optics | Uniform u can be phase-visible even when it causes no deflection. |
| P001-E08 | ds_eff^2=-c^2 dt^2/n(u)^2+dx^2 (up to conformal factor) | equivalent representation, not PBUF-derived | null condition gives \|dx/dt\|=c/n | P001-A03 and P001-A06 | A metric interpretation adds no prediction until n(u) is known; null paths do not fix the conformal factor. |
| P001-E09 | current WL: update v proportional to -grad u, then normalize | empirical compatibility target | pbuf_experiment.py propagate | implicit beta/sign/normalization; x component is additionally weighted by 0.15 | Matches the structure of P001-E06 only at leading order after projection, but is not uniquely derived and is not rotationally isotropic as coded. |

## Compatibility and alternatives

The frozen WL routine samples `grad u`, updates a direction approximately proportional to `-grad u`, and renormalizes it. Renormalization supplies the transverse projection to first order, so its broad structure can approximate the linearized conditional ray law with an implicit negative `beta`. This is only structural compatibility: the code gives the x-gradient an extra factor `0.15`, so it is not the rotationally invariant scalar law above, and neither its sign nor normalization follows from PBUF. PHOTON-001 therefore does not authorize changing it.

| Hypothesis | Allowed | Distinctive effect | PBUF status |
|---|---|---|---|
| direct deformation n(u) | yes, conditionally | uniform u changes phase; gradients bend | not uniquely implied |
| gradient-only n(\|grad u\|^2) or higher-derivative coupling | yes | uniform u invisible; typically nonlinear and introduces extra boundary sensitivity | not selected by current axioms |
| nonlocal accumulated response | yes | history/path dependence | requires a kernel absent from PBUF |
| effective metric g_eff[u] | yes | covariant null geodesics and possible time delay | tensor map, causal dynamics, and normalization absent |
| direct vector/tensor microstate coupling | yes | polarization or direction dependence | hidden by CORE-001 scalar projection; not excluded by three-component ontology |

The three-component ontology alone does not choose among these hypotheses. A scalar coarse projection can make three-component, scalar, and generic-N microscopic models photon-equivalent.

## Predicted observables (conditional)

| ID | Observable | Conditional prediction | Dependencies | Future pass/fail criterion |
|---|---|---|---|---|
| P001-O01 | ray deflection/curvature | proportional to the transverse gradient of u in the weak scalar limit | A01--A05 and unknown beta | PHOTON-002 numerical curvature agrees with beta(I-tt)grad u to preregistered truncation tolerance |
| P001-O02 | relative phase/time delay | line integral of n(u)-1; can be nonzero without bending | coherent timing/phase readout and n(u) | phase scales linearly with path length in a uniform-u slab and shows zero transverse deflection |
| P001-O03 | frequency dependence | none under A03; chromaticity signals dispersion or failure of scalar nondispersive closure | multi-frequency propagation | paths agree across frequency within tolerance, or A03 is rejected |
| P001-O04 | polarization dependence | none under A03; splitting signals birefringent/tensor coupling | polarization-resolved propagation | polarization paths agree within tolerance, or scalar coupling is rejected |
| P001-O05 | rotational covariance | rotating u and initial ray rotates the output path identically | isotropic scalar hypothesis | rotated/unrotated solutions agree after inverse rotation; current 0.15 x weighting is expected not to satisfy this gate |

These are discrimination tests for a chosen coupling family, not observationally validated PBUF predictions.

## PHOTON-002 implementation specification

Implement a new, isolated candidate interface `optical_response(u)->n and grad_log_n(u,grad_u); integrate dt/ds=(I-tt^T)grad_log_n` while preserving the frozen propagator. Record `u`, coordinates, independently justified `n(u)`/`beta`, initial ray data, numerical controls, and optional frequency/polarization labels. Emit paths, tangents, curvature, phase/time delay, convergence tests, and complete coupling provenance.

Required gates:

1. u=0 and spatially uniform u produce straight rays
2. constant transverse grad u produces the analytic small-beta curvature
3. longitudinal gradients do not create first-order transverse curvature
4. step refinement converges at the integrator's declared order
5. rotations and translations commute with propagation
6. polarization and frequency null tests hold for scalar nondispersive n
7. phase accumulates through a uniform-u slab while deflection remains zero
8. comparison to frozen WL is diagnostic only and changes no existing artifact or propagator

Observational comparison remains forbidden until a photon/electromagnetic action independently fixes `n(u)` or at least `beta`. A future implementation may compare candidate paths with the frozen interface as a compatibility diagnostic, but must expose the sign, units, projection, and anisotropic x weighting rather than absorb them into a fitted constant.

## Completion assessment

PHOTON-001 satisfies the document's alternative completion route: the conditional interaction is mathematically explicit and traceable, and the exact missing theoretical step is identified as a microscopic/covariant map from `q` or `u` to the electromagnetic action (equivalently `n(u)` or `g_eff[u]`). Automated completeness checks pass: **True**.
