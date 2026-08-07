# PBUF LENS-LOAD-001 — Lens 001 native-load reconstruction

## Decision

The requested observed-lens loading field is **not identifiable from the
supplied inputs under the frozen framework**. This is a structural result, not
a numerical failure. The frozen elastic map determines deformation from a
given generalized load, but no frozen member of the medium-to-metric family
maps that deformation to weak-lensing observables. Its normalization and
spatial support are also unselected. Observational uncertainties were not
provided.

Consequently no reconstructed-load array, load map, mass correlation, or
forward optical validation is asserted. Writing one would introduce exactly
the metric/coupling assumption forbidden by this milestone.

## Local domain and boundary

The executable gate uses the square data footprint `[-8,8] x [-8,8]`, which
is the smallest region evidenced by the supplied archived baryonic and optical
arrays and contains the archived ray window. It declares zero placement on the
four sides. Full Dirichlet data remove rigid modes and give a well-posed local
elliptic restriction. This is explicitly a finite isolated-boundary
approximation: the frozen framework provides no falloff law, so real-data use
must enlarge the footprint until the reconstruction is insensitive to the
boundary. `computational_domain.json` records the exact grid and caveat.

## Frozen forward and conditional inverse

For the authorized placement realization,

\[
 F=\operatorname{Grad}_0y,\quad C=F^\sharp F,\quad
 E=(C-\mathbf1)/2,
\]

\[
 P_C={K_0\over2}\operatorname{tr}(E)\mathbf1+\mu_0E_{TF},
 \qquad P_F=2FP_C,
\]

and the static weak problem is

\[
 \int_\Omega P_F[y]:\operatorname{Grad}_0\eta\,dV_0
 =\langle\Pi,\eta\rangle+\int_{\Gamma_N}\bar t\cdot\eta\,dA_0.
\]

No constitutive equation is changed. If an admissible placement and boundary
traction are independently known, the exact local inverse is

\[
 \langle\Pi_{\rm req},\eta\rangle
 =\int_\Omega P_F[y]:\operatorname{Grad}_0\eta\,dV_0,
 \qquad b_{\rm req}=-\operatorname{Div}_0P_F[y]
\]

where the pointwise form requires regularity. The program can evaluate and
forward-check a scalar anti-plane restriction when passed `--placement`; it
labels every such artifact **conditional** and **not inferred from lensing**.

## Why the optical inverse cannot run

The complete observation operator would be

\[
 \mathcal H_{G,\beta,z}(\Pi)=M R_z G[q(S_\beta(\Pi)),C(S_\beta(\Pi))].
\]

The frozen corpus fixes the elastic part `S_beta`, but does not select `G`.
Therefore the Jacobian factor

\[
 D\mathcal H=DM\,DR_z\,DG\,\mathcal L_0^{-1}
\]

is unknown before ordinary finite-field and line-of-sight null directions are
even considered. Baryonic matter cannot close this gap without a forbidden
matter–load hypothesis. The archived `observation.csv` is an image-intensity
proxy, not a calibrated weak-shear catalogue with covariance.

## Deliverable disposition

| Requested item | Result |
|---|---|
| local computational domain | delivered, with boundary limitation |
| local inverse formulation | delivered as a set-valued/conditional inverse |
| reconstructed loading field and load maps | not identifiable; not fabricated |
| mass correlation | not defined without a load estimate |
| forward validation against lensing uncertainty | impossible without `G`, measurement operator, and covariance |
| exportable load dataset | schema represented by the conditional output; no observed Lens 001 dataset claimed |

`reconstruction_status.json` is the machine-readable readiness decision and
`load_characterization.csv` records every unavailable characterization rather
than encoding missing values as physical zeros.

## Completion status

**BLOCKED_BY_FROZEN_IDENTIFIABILITY.** To complete the mission without changing
the constitutive law or V11, supply an independently authorized, frozen-compatible
metric map and optical measurement operator, calibrated single-lens weak-lensing
data with covariance, and boundary/background data. These are closure and
problem-data requirements, not requests for cosmology or a universal
matter–medium interaction law.
