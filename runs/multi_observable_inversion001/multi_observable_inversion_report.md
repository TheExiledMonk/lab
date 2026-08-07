# PBUF MULTI-OBSERVABLE-INVERSION-001 — Reconstruction of the Universal Matter-to-Load Operator

## 0. Decision

Combining independent matter, geometry, kinematic, thermodynamic, and optical
measurements makes the inverse problem **strictly better constrained whenever
their linearized observation operators remove different null directions**.
It does not, under the frozen framework alone, make the unrestricted
matter-to-load operator unique.

The decisive separation is

\[
 \boxed{
 \text{primitive measurements}
 \longrightarrow (T,z,\beta,d_{\rm dyn},d_{\rm gas},d_{\rm opt})
 \longrightarrow \mathfrak P_{\rm joint}
 \longrightarrow \Pi_{\rm req}
 \longrightarrow \mathfrak L? }
 \tag{MOI-001}
\]

Matter observations constrain the input descriptor `T`; geometry constrains
the reference region, optical paths, and boundary/background bookkeeping;
non-optical response measurements can constrain deformation or the effective
field through distinct forward maps; and lensing constrains projected optical
action.  These are not interchangeable roles.

For a fixed system, fixed forward closures, and known boundary work, the
native load remains uniquely determined by an admissible placement:

\[
 \langle\Pi_{\rm req}[y],\eta\rangle
 =\int_\Omega P_F[y]:\operatorname{Grad}_0\eta\,dV_0,
 \qquad \eta\in V_0.                                      \tag{MOI-002}
\]

Additional observations improve recovery of `y` only through independently
specified observation maps.  All lensing channels still contain the same
unselected medium-to-metric map `G`; stellar/gas dynamics require an
independently specified tracer-response law; and bulk/boundary separation
requires independent boundary data.  Data abundance cannot select any of
these missing maps by itself.

The smallest operator family authorized after the joint audit is therefore
not a formula.  It is the intersection class

\[
 \boxed{\mathfrak C_*=
 \{\mathfrak L:\ T\mapsto V_0^*\mid
 \mathfrak L\text{ is gauge-basic, natural/covariant, dimensionally and
 mechanically admissible, and intersects every independent joint-data
 fiber}\}.}                                               \tag{MOI-003}
\]

The frozen corpus does not select within `C_*` the input component, index
bridge, normalization, differential order, locality/support, additivity,
species rule, or configuration dependence.  Finite observations determine an
operator only on the observed input set (and only modulo output null spaces);
infinitely many extensions remain unless an operator class and all forward
maps are fixed independently.

Thus the milestone outcome is:

1. multi-observable reconstruction is mathematically and observationally
   preferable to weak-lensing-only reconstruction;
2. it can remove specific optical, support, boundary, and model ambiguities
   when an independent channel is sensitive to them;
3. it can falsify and strongly reduce a predeclared operator class across
   heterogeneous systems; but
4. it cannot uniquely recover an unrestricted universal interaction law from
   the frozen framework without new physical assumptions.

No observable is fitted here.  No metric map is derived or selected.  No
ontology, state variable, constitutive law, V11 statement, or empirical
coupling constant is added.

## 1. Observable catalogue

### 1.1 Classification rules

An **observable** below means a recorded instrument-level quantity or a
calibrated quantity obtained without using the elastic reconstruction.  A
**descriptor** is a model-mediated representation constructed from such
measurements.  A **derived proxy** is not an independent datum merely because
it is published as a catalogue column.

Independence from elastic reconstruction is weaker than statistical
independence.  Two channels may both be non-circular yet share calibration,
selection, distance, membership, or line-of-sight errors.  The archive must
retain primitive measurements and their joint covariance.

### 1.2 Complete information catalogue

| Category | Primitive measured information | Common derived information | Role in the joint inverse | Independence qualification |
|---|---|---|---|---|
| sky position and morphology | angular coordinates, image pixels, isophotes, sizes, orientations | centroids, ellipticity, position angle, light profile, projected axis ratio | tracer support, projected geometry, sampling and masks | independent of elastic reconstruction if not lensing-corrected with the tested model; deprojection is model-dependent |
| distances and redshifts | spectroscopic/photometric redshift indicators, standard-ruler/candle measurements where available | angular-diameter/luminosity distance, lens/source ordering, physical scale | ray geometry `z`, volume conversion, reference-domain scale | distance products share background calibration; a lensing-inferred distance is optical, not independent geometry |
| stellar light | multi-band fluxes, spectra, surface-brightness maps | luminosity, colour, stellar population descriptors | spatial tracer for baryonic matter | direct light is independent; dust, population and distance conversions are correlated systematics |
| stellar mass | stellar spectra/photometry plus population model | mass-to-light ratio and stellar-mass map | matter descriptor `T_*` | model-mediated, not primitive; independent of lensing only if no lensing/dynamical mass calibration enters |
| atomic/molecular gas | line intensity, frequency, width, absorption, resolved emission | column density, gas mass, velocity field | matter support, baryonic descriptor, tracer kinematics | conversion/excitation/opacity and distance dependent; independent if conversion is not calibrated by the tested elastic/optical model |
| ionized/hot gas | X-ray photon positions, energies, counts, spectra; thermal/radio observables | emissivity, temperature, abundance, electron density, pressure, gas mass | thermodynamic matter descriptor and geometry | surface brightness and spectrum are independent; deprojection and plasma calibration are model-mediated |
| total baryonic matter | independently constructed stellar, cold-gas, hot-gas, dust and compact-object inventories | baryonic density/stress-energy descriptor with covariance | principal candidate input `T_b` | derived sum with shared distance, aperture and membership covariance; omitted sectors remain support uncertainty |
| galaxy/tracer positions | angular positions, redshifts, membership indicators | number-density field, substructure, projected/3-D geometry | environment, support and independent tracer distribution | membership and projection are correlated; positions do not directly measure native load |
| stellar/galaxy kinematics | spectral-line centroids and shapes, proper motions where available | line-of-sight velocities, velocity-dispersion profile, higher moments, rotation | independent response channel `d_dyn` | primitive velocities are independent; dynamical mass/potential is circular unless a declared independent dynamics law is used only in the forward model |
| gas kinematics | Doppler shifts, line profiles, resolved velocity maps | rotation curve, dispersion, inflow/outflow descriptors | independent response and equilibrium diagnostic | conversion to enclosed mass assumes dynamics and often equilibrium |
| gas thermodynamics | X-ray/SZ/radio spectral and intensity information | temperature, pressure, entropy and density profiles | independent matter/response constraints | hydrostatic mass is not primitive and must not be treated as independent matter |
| time variability | repeated flux, position, spectrum and delay measurements | variability curves, pattern speeds, transient delays | history/time-dependent response when frozen duration/kinetic closure exists | static milestone can retain raw data but may not interpret it without authorized dynamics |
| weak lensing | source image pixels, shapes, sizes, fluxes, redshifts | reduced shear, magnification, flexion, convergence or aperture maps | projected optical constraint `d_W` | raw optical channel; convergence/mass maps are model inversions and are not matter inputs |
| strong lensing | multiple-image positions, resolved arcs/rings, relative shapes and fluxes, source redshifts | critical curves, lens potential or enclosed lensing mass | high-resolution optical constraint `d_S` | independent measurements but correlated with weak lensing through `G`, line-of-sight structure and lens model |
| Einstein radius | angular radius or critical-curve geometry | enclosed projected lensing mass | compressed strong-lensing constraint | derived from strong-lensing images; never an additional independent datum if those pixels/positions are already used |
| lensing time delays | measured light curves and relative arrival delays | time-delay distance or Fermat-potential differences | optical-plus-clock constraint `d_\Delta` | measured delays add information, but derived distances/potentials assume `G`, a lens model and background calibration |
| non-lensing time delays | independently timed propagation signals, if applicable | path-integrated propagation response | possible distinct propagation constraint | only independent if its forward law and calibration do not reuse the tested optical reconstruction |
| external environment | positions, redshifts, light/gas measurements of surrounding structures | external convergence/shear, tidal/environment model | boundary/background descriptor `\beta` | raw environment data are independent; lensing-derived external convergence is part of the optical inverse |
| domain and boundary data | observed truncation geometry, asymptotic/background conditions, independently justified displacement/traction information | boundary package `\beta=(\bar y,\bar t,\beta_{bg})` | separates bulk load from boundary work | geometry alone does not measure traction; most exact elastic boundary values remain latent |
| instrument/selection metadata | point-spread response, exposure, noise, masks, selection, calibration standards | measurement operators and covariance | defines every `M_a` and uncertainty | not physical evidence for an operator, but indispensable to identifiability |

Strong lensing, weak lensing, Einstein radius, magnification, and time delay are
distinct summaries or samples of the optical field, not automatically
independent physical channels.  Stellar mass, gas mass, baryonic mass,
hydrostatic mass, and dynamical mass likewise name different inference levels;
only their non-reused primitive measurements can establish independence.

## 2. Information independence and non-reuse

### 2.1 Independence labels

Use the following labels:

- `I`: physically distinct primitive information, conditional on recorded
  cross-covariance;
- `C`: non-circular but statistically or systematically correlated;
- `D`: deterministically/model-derived from another entry and not additional
  information when its parent data are used;
- `F`: shares an essential forward factor, so it may shrink a fiber but cannot
  independently identify that factor;
- `X`: circular for the stated use and forbidden as simultaneous input and
  validation.

### 2.2 Independence matrix

Rows are proposed inputs; columns are validation or reconstruction channels.

| Input \ validation | stellar/gas matter inventory | geometry/environment | kinematics | X-ray/SZ thermodynamics | weak lensing | strong lensing | Einstein radius | lensing delays |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| primitive stellar/gas measurements | C | C | I/C | I/C | I/C | I/C | I/C | I/C |
| population/conversion-derived baryonic map | — | C | I/C | C | I/C | I/C | I/C | I/C |
| sky/redshift/environment geometry | C | — | C | C | C | C | C | C |
| primitive kinematic spectra/velocities | I/C | C | — | I/C | I/C | I/C | I/C | I/C |
| dynamical mass/potential product | C | C | X | C | X if lensing-calibrated | X if lensing-calibrated | X if lensing-calibrated | X if lens-model calibrated |
| primitive X-ray/SZ photons/spectra | I/C | C | I/C | — | I/C | I/C | I/C | I/C |
| hydrostatic mass product | C | C | C | X | X if lensing-calibrated | X if lensing-calibrated | X if lensing-calibrated | X if lens-model calibrated |
| weak-lensing image catalogue | I/C | C | I/C | I/C | X | F | D/F | F |
| weak-lensing convergence/mass map | X as matter input | C | C | C | D/X | F/X | D/F/X | F/X |
| strong-lensing images/positions | I/C | C | I/C | I/C | F | X | D | F |
| Einstein radius | I/C | C | I/C | I/C | F | D | X | F |
| measured lensing light curves/delays | I/C | C | I/C | I/C | F | F | F | X |
| lens-model potential/time-delay distance | X as matter input | C | C | C | F/X | D/X | D/X | D/X |

`I/C` means the physical channel is distinct but known shared nuisance
parameters must appear in the joint covariance.  All optical columns share
`G`, ray geometry, and line-of-sight structure; their errors need not be
identical, but they cannot select `G` merely by being combined.

### 2.3 Ledger rule

For each system define disjoint raw-data index sets

\[
 R_s^{\rm input},\quad R_s^{\rm reconstruct},\quad
 R_s^{\rm validate},\qquad
 R_s^{\rm input}\cap R_s^{\rm validate}=\varnothing,
 \quad R_s^{\rm reconstruct}\cap R_s^{\rm validate}=\varnothing. \tag{MOI-004}
\]

Disjoint catalogue labels are insufficient: two products derived from the
same pixels or spectra have the same raw-data lineage.  Every derived product
must carry a provenance map to its primitive measurements, assumptions, and
nuisance parameters.  A datum used to choose geometry, regularization,
operator support, normalization, or candidate class is training information
and cannot be called validation.

## 3. Joint inverse formulation

### 3.1 Latent variables and forward maps

For system `s`, let

\[
 \theta_s=(T_s,z_s,\beta_s,y_s,G,\kappa_s,\nu_s),            \tag{MOI-005}
\]

where `T_s` is an authorized matter descriptor, `z_s` is observed geometry,
`beta_s` is boundary/background data, `y_s` is placement, `G` is an admissible
metric map, `kappa_s` collects independently specified tracer-dynamics or gas
response closures, and `nu_s` contains nuisance/calibration variables.  These
symbols organize existing information; they add no PBUF state variables.

For measurement families `a`, write

\[
 d_{a,s}=\mathcal H_{a,s}(\theta_s)+\epsilon_{a,s}.          \tag{MOI-006}
\]

Representative factors are

\[
\begin{aligned}
 \mathcal H_{m,s}&=M_{m,s}\,Q_m(T_s,z_s,\nu_s),\\
 \mathcal H_{g,s}&=M_{g,s}\,Q_g(z_s,\beta_s,\nu_s),\\
 \mathcal H_{o,s}&=M_{o,s}R_{o,s}(z_s)G[q(y_s),C(y_s)],\\
 \mathcal H_{k,s}&=M_{k,s}K_s[y_s,G,T_{\rm tracer,s};\kappa_s],\\
 \mathcal H_{x,s}&=M_{x,s}X_s[T_{\rm gas,s},y_s,G;\kappa_s].
                                                               \tag{MOI-007}
\end{aligned}
\]

`K_s` and `X_s` are placeholders for independently authorized forward laws,
not laws supplied by this milestone.  If no such law is frozen or declared,
the corresponding raw data remain descriptive and cannot be converted into a
constraint on `y` or `Pi`.

Mechanical feasibility is

\[
 \mathcal A_{\Omega_s}(y_s)=\Pi_s+\mathcal T_{\bar t_s},
 \qquad \Pi_s\in V_{0,s}^*.                                \tag{MOI-008}
\]

For acceptance regions `A_{a,s}` determined from measurement uncertainty, the
joint inverse is the set

\[
 \boxed{
 \mathfrak P_s^{\rm joint}=
 \{(\theta_s,\Pi_s):
 \mathcal H_{a,s}(\theta_s)\in A_{a,s}\ \forall a,
 \ \mathcal A_{\Omega_s}(y_s)=\Pi_s+\mathcal T_{\bar t_s}\}.}
                                                               \tag{MOI-009}
\]

The requested generalized inverse is therefore

\[
 \boxed{
 \mathcal I:(O_1,\ldots,O_n)\longmapsto
 \operatorname{proj}_{\Pi}\mathfrak P_s^{\rm joint}
 \subseteq V_{0,s}^*.}                                      \tag{MOI-010}
\]

It is set-valued unless the joint forward map is injective on the admissible
quotient and stable against noise.

### 3.2 Better-posedness criterion

Let `J_a=D\mathcal H_a` at an admissible solution after nuisance, gauge, and
boundary variables are included, and define the stacked derivative

\[
 \mathcal J_{\rm joint}\delta\theta
 =(J_1\delta\theta,\ldots,J_n\delta\theta).                 \tag{MOI-011}
\]

Then

\[
 \ker\mathcal J_{\rm joint}=\bigcap_a\ker J_a.              \tag{MOI-012}
\]

The joint inverse is locally more identifiable than weak lensing alone iff

\[
 \ker\mathcal J_{\rm joint}
 \subsetneq \ker J_W                                      \tag{MOI-013}
\]

on the physical quotient, and it is locally stable only if the stacked
operator has an appropriate positive lower bound on the remaining observable
subspace.  Correlated data change the weighting and effective rank, not the
kernel-intersection identity.  Repeated measurements proportional to an
existing row improve precision but do not add identifiability.

Since matter observations primarily constrain `T`, not `y`, their addition
does not by itself reconstruct `Pi`.  Their essential contribution is to make
the operator comparison `T -> Pi` testable once independent response channels
constrain the load fiber.

## 4. Identifiability analysis

### 4.1 Null-direction audit

| Ambiguity | Can joint observables reduce it? | Condition for reduction | Irreducible remainder under frozen framework |
|---|---|---|---|
| optical projection/line-of-sight nulls | yes | distinct source planes, strong-lensing localization, magnification/delay information, or non-optical response has a derivative not vanishing on the optical null | unsampled regions and perturbations invisible to every retained channel |
| reduced-shear/mass-sheet/source-position families | conditionally | absolute magnification, independent geometry/distance, time delays, boundary information, or non-optical response breaks the same transformation | transformations that can be absorbed jointly into `G`, source, distance, or boundary variables |
| coordinate and rigid-motion gauge | no physical removal is needed | impose quotient/gauge convention consistently | gauge-equivalent representatives remain intentionally identical |
| unselected metric map `G` | optical multiplicity alone: no; heterogeneous channels: only conditionally | `G` belongs to a predeclared identifiable class and a non-optical map constrains `y` independently | unrestricted factorization ambiguity between `G`, deformation, and loading |
| boundary displacement/traction | yes, conditionally | independent environment/asymptotic data or direct boundary response | unknown boundary work can always trade against bulk load |
| exterior/background influence | yes, conditionally | complete observed environment and selected support rule | unobserved exterior and unrestricted nonlocal `G`/operator support |
| matter support | yes | resolved stellar/gas/tracer observations at sufficient depth and independent membership | unseen matter-bearing distinctions, deprojection, finite resolution, and below-threshold support |
| line-of-sight deprojection | yes, but not generally uniquely | multiple viewing constraints, symmetry independently justified, kinematics/thermodynamics with fixed forward laws | single-view 3-D arrangements producing the same projections |
| body load versus singular/interface/boundary load | conditionally | resolution and test functions sensitive to each support class plus known boundary work | sub-resolution distributions with identical generalized work |
| load representative | no unique pointwise representative is required | compare in `V_0^*` | distributional representatives with identical action on all admissible variations are mechanically equivalent |
| operator input choice | yes across heterogeneous systems | independently vary density, stress, momentum, species and geometry | collinear descriptors on the observed sample cannot be separated |
| differential order/locality | conditionally | spatially resolved inputs and outputs over multiple scales/geometries | operators agreeing on measured bandwidth and support remain equivalent |
| normalization | conditionally | absolute input, response, distance and boundary calibration with fixed `G` | composite rescalings among unselected forward factors |
| species/additivity | conditionally | systems with independently varying matter composition and common calibration | sectors that never vary independently are not identifiable |
| off-sample operator extension | no from finite samples | requires independent regularity/complexity assumptions or sufficiently rich continuum excitation | infinitely many operators coincide on all observed `T_s` and differ elsewhere |

### 4.2 Structural rank statement

Let `N_gauge` denote the authorized gauge/rigid subspace.  Conditional local
identifiability of `Pi` requires

\[
 \ker\mathcal J_{\rm joint}
 \subseteq N_{\rm gauge}+\ker D\mathcal A_\Omega,           \tag{MOI-014}
\]

together with independent boundary work and local injectivity of the
constitutive branch.  Conditional identifiability of an operator perturbation
`delta L` over systems `s=1,...,S` further requires

\[
 \delta\mathfrak L[T_s]
 \notin N_{\Pi,s}\ \text{for some }s
 \quad\text{for every nonzero admissible }\delta\mathfrak L, \tag{MOI-015}
\]

where `N_{Pi,s}` is the residual load-observation null space.  If a nonzero
operator perturbation maps every observed input into the corresponding null
space, it is unidentifiable regardless of measurement precision.

The mathematical gain can be reported without data as kernel codimension,
rank of the stacked derivative on a finite candidate tangent space, singular
value bounds, and intersections of admissible fibers.  No numerical value is
claimed here.

## 5. Reduction of the universal operator class

### 5.1 Frozen admissibility envelope

Every surviving universal loading operator must satisfy all of the following:

1. **Codomain:** `L[T;q_0,beta]` is one generalized placement-load covector in
   `V_0*`; it is not `W`, stress, deformation, or a second field.
2. **Gauge and covariance:** it descends from the complete physical-state
   quotient and transforms naturally under authorized changes of reference
   description.
3. **Work and units:** its pairing with a placement variation has energy
   dimension; a regular body representative has dimension `E L^-4`.
4. **Mechanical admissibility:** it is continuous/distributionally admitted,
   satisfies force/moment compatibility where required, and permits a solution
   in the frozen constitutive domain.
5. **Universality:** the same rule, normalization, index bridge, support rule,
   differential order, species rule, and configuration-dependence rule applies
   to every system.
6. **Source separation:** it enters only on the load side of frozen balance
   and does not alter `W`, `P_C`, `P_F`, `C`, or V11.
7. **Isotropic obstruction:** a nonzero parity-even point-local algebraic map
   from a symmetric rank-two matter tensor to a spatial load covector is
   excluded on the isotropic unloaded branch unless an already-authorized
   direction or derivative is present.
8. **Boundary/support declaration:** bulk, interface, singular, boundary, and
   nonlocal contributions are explicitly distinguished.

### 5.2 Data-fiber reduction

For each system define the admissible input set `T_s` and load fiber

\[
 \mathcal F_s=\operatorname{proj}_{\Pi}
 \mathfrak P_s^{\rm joint}.                                \tag{MOI-016}
\]

For any predeclared structural family `C_0`, observations reduce it to

\[
 \boxed{
 \mathfrak C_S
 =\{\mathfrak L\in\mathfrak C_0:
 \mathfrak L[T_s;q_{0,s},\beta_s]\in\mathcal F_s
 \text{ for every }s\}.}                                  \tag{MOI-017}
\]

This is the smallest empirically admissible class without fitting.  Empty
`C_S` falsifies `C_0`; a strict subset is genuine reduction; a singleton is
identification only relative to `C_0` and the fixed forward closures.

With no independently predeclared `C_0`, the smallest frozen family remains
`C_*` in MOI-003.  It includes, without selecting among them, finite-jet local
operators using authorized gradients/directions, boundary-supported
functionals, singular/interface loads, configuration-dependent work forms,
and support-declared nonlocal functionals.  Multiple observations may reject
members whose predictions miss a joint fiber, but cannot make this unrestricted
family unique because an operator may be changed away from the observed input
set or within every `N_{Pi,s}`.

## 6. Cross-system universality protocol

For galaxies A and B and cluster C, or any larger heterogeneous collection:

### Phase 1 — preregistration

1. freeze `C_0`, including locality/support, derivative order, input tensor,
   normalization convention, species/additivity rule, and covariance rule;
2. freeze all observation maps, `G` class, dynamics/gas closures, boundary
   conventions, gauge quotient, and acceptance criteria;
3. record raw-data provenance and partition systems and channels into discovery
   and validation sets before inspecting validation residuals.

### Phase 2 — independent system reconstruction

For each `s`, construct matter-input uncertainty `T_s`, geometry/boundary
uncertainty, and the joint load fiber `F_s` independently.  Do not force a
common operator during this phase.  Compare generalized virtual work rather
than coordinate components, and retain multimodality and null directions.

### Phase 3 — intersection test

Compute the logical intersection MOI-017.  A universal candidate exists iff

\[
 \exists\mathfrak L\in\mathfrak C_0\quad
 \mathfrak L[T_A]\in\mathcal F_A,\quad
 \mathfrak L[T_B]\in\mathcal F_B,\quad
 \mathfrak L[T_C]\in\mathcal F_C.                          \tag{MOI-018}
\]

System geometry and boundary data may enter only through the same declared
natural rule; system-specific retuning is failure of the universal candidate.

### Phase 4 — held-out forward validation

Use the unchanged chain

\[
 T_s\xrightarrow{\mathfrak L}\Pi_s
 \xrightarrow{\text{frozen balance}}y_s
 \xrightarrow{F,C}C_s
 \xrightarrow{G}g_s^{\rm eff}
 \xrightarrow{R,M}d_s^{\rm pred}.                          \tag{MOI-019}
\]

Test withheld systems and withheld channels without changing `L`, `G`,
boundary rules, normalization, or nuisance priors.  Validate in primitive
observable space.  Failure in a measured direction falsifies the candidate;
residuals lying wholly in an admitted null direction are inconclusive.

### Phase 5 — universality stress tests

Require the collection to vary independently in scale, morphology,
composition, environment, and matter-state gradients.  Leave-one-system-class
out tests distinguish interpolation from universality.  A class surviving
galaxies but failing clusters is not universal; separate per-class rules are
not permitted to masquerade as one operator.

Agreement supports only `C_S`.  It does not prove uniqueness against operator
families not placed in `C_0`.

## 7. Uncertainty propagation and robustness

### 7.1 Block propagation

Let the primitive data vector be `d`, with full covariance `Sigma_d`, and let
`r(theta,L)` stack measurement and mechanical residuals.  Linearization gives

\[
 D_\theta r\,\delta\theta
 +D_{\mathfrak L}r\,\delta\mathfrak L
 +D_d r\,\delta d=0.                                      \tag{MOI-020}
\]

On an identified quotient and for a declared generalized inverse `+`,

\[
 \delta\theta
 =-(D_\theta r)^+
   (D_d r\,\delta d+D_{\mathfrak L}r\,\delta\mathfrak L)
   +n,\qquad n\in\ker D_\theta r.                          \tag{MOI-021}
\]

The load variation is

\[
 \delta\Pi
 =D_y\mathcal A_\Omega\,\delta y
 -\delta\mathcal T_{\bar t},                              \tag{MOI-022}
\]

and, where covariance propagation is meaningful,

\[
 \Sigma_\Pi\approx J_{\Pi d}\Sigma_dJ_{\Pi d}^*
 +J_{\Pi\beta}\Sigma_\beta J_{\Pi\beta}^*
 +\Sigma_{\rm structural}.                                \tag{MOI-023}
\]

`Sigma_structural` must not hide discrete uncertainty over `G`, dynamics/gas
closures, support class, deprojection, or regularization; those are propagated
as separate fibers/model branches, not converted automatically into Gaussian
noise.

### 7.2 Error budget by stage

| Stage | Principal uncertainties | Effect on reconstructed loading |
|---|---|---|
| matter measurement | flux/spectral calibration, distance, conversion factors, population/plasma models, incompleteness, membership | moves and blurs operator input `T`; chiefly limits discrimination among input components and species rules |
| geometry | distance/background calibration, centering, inclination, deprojection, line-of-sight structure, domain choice | changes physical scale, support, ray map and derivative estimates; couples to both input and response |
| boundary/background | exterior tides, truncation, unknown traction/displacement, asymptotic prescription | trades directly with inferred bulk load and can dominate large-scale/low-order modes |
| deformation/response reconstruction | optical projection, tracer anisotropy/equilibrium, gas equilibrium, sampling and resolution | controls which components of `y` and hence `Pi` are observed; null directions produce unbounded rather than merely large variance |
| metric map | admissible form, support, normalization and derivatives of `G` | structural factorization ambiguity; typically dominant until independently fixed |
| mechanical inversion | branch selection, ellipticity margin, gauge/rigid quotient, differentiation of noisy `y` | `Pi=A(y)-T` amplifies poorly resolved spatial modes; conditioning worsens near loss of ellipticity |
| operator comparison | finite system diversity, correlated descriptors, bandwidth, off-sample extension | limits rank of operator identification even when each `Pi_s` is well reconstructed |
| optical prediction | PSF, shape/noise calibration, source redshifts, ray/background model | affects validation power but cannot repair upstream structural non-identifiability |

There is no universal ordering of numerical dominance without data.  The
structurally dominant uncertainties are those generating exact or near-null
directions: unrestricted `G`, unknown boundary work, incomplete 3-D support,
and unselected response laws.  Once those are fixed, measurement-specific
dominance is assessed through singular spectra or generalized Fisher/normal
operators, not asserted in advance.

### 7.3 Required robustness tests

Vary, without fitting an interaction law: data partitions; spatial resolution
and masks; line-of-sight/environment models; boundary/domain placement;
deprojection branches; admissible `G` branches; tracer anisotropy/equilibrium
branches; load support class; gauge conventions; and regularization used only
to display representatives.  A claimed operator restriction is robust only if
it persists at the fiber level across these variations.

## 8. Physical interpretation of a surviving class

If a universal class survives MOI-017 and held-out validation, the following
interpretations remain compatible with the frozen one-medium ontology and
local balance.  They are classifications, not selections.

| Mathematical signature | Compatible one-medium interpretation |
|---|---|
| regular bulk covector | internal transfer of generalized placement work from a matter-bearing aspect of the same medium |
| divergence/eigenstress form | stress polarization or equilibrium-mismatch representation, subject to its compatibility and gauge freedom |
| boundary-supported functional | exchange across a regional accounting boundary or representation of an excluded exterior |
| singular/interface support | localized defect-like or interface organization within the one medium |
| finite-jet dependence on matter descriptors | local interaction using authorized spatial variation and index bridges |
| configuration-dependent work form | follower interaction or state-dependent internal exchange without changing `W` |
| support-declared nonlocal functional | extended bookkeeping of one-medium interaction, only if observationally established and explicitly bounded |
| one additive rule across sectors | universal resultant work channel with sector contributions |
| one non-additive rule | collective interaction among matter-bearing distinctions, still producing one resultant load |

Elastic and optical equivalence cannot determine microscopic mechanism.
Different internal accounts that induce the same generalized work are
physically indistinguishable by this reconstruction.  None introduces a
second medium, force carrier, new constitutive response, or independent
preferred-state field.

## 9. Readiness assessment

### 9.1 What is ready

The frozen framework is ready for a **multi-observable, set-valued
reconstruction program** with:

- a complete primitive-observable and provenance catalogue;
- explicit joint likelihood/acceptance interfaces and cross-covariance;
- canonical conditional reconstruction `y -> Pi_req`;
- quotient-aware kernel and stability tests;
- operator-class intersection across systems; and
- held-out forward validation with the unchanged mechanical pipeline.

This program can determine whether a specified local/nonlocal, tensorial,
support, species, or additivity class is incompatible with all observations.

### 9.2 What is not ready

An observationally unique unrestricted operator is not ready because the
frozen theory does not supply:

1. a selected and identifiable medium-to-metric map `G`;
2. complete independent boundary work;
3. selected tracer-dynamical and gas-response maps for converting those raw
   channels into deformation constraints;
4. a unique 3-D matter descriptor and pullback/index bridge; or
5. an independently bounded operator class whose off-sample extension can be
   identified.

### 9.3 Completion theorem

**Theorem — status of MULTI-OBSERVABLE-INVERSION-001.**  Let a collection of
independent observables have forward derivatives `J_a`.  Their joint inverse
has null space equal to `intersection_a ker J_a`; hence it is strictly better
posed than the weak-lensing-only inverse precisely when that intersection is
a strict subset of the weak-lensing kernel on the authorized quotient.  With
known boundary work and sufficient independent constraints on placement, the
frozen constitutive first variation then yields a unique generalized load.

For multiple systems, a declared universal family reduces exactly to the
intersection class MOI-017.  It can be strongly constrained or falsified when
systems excite independent input directions and their residual load null
spaces do not contain the corresponding operator differences.  However,
because the metric, response, boundary, support, and unrestricted
off-sample-extension ambiguities remain unfrozen, the complete corpus does not
select a unique universal matter-to-load operator.  The smallest remaining
admissible family is MOI-003, further restricted by the eight frozen
conditions in Section 5.1 and by all joint-data fibers.

Therefore multiple independent observations provide a **significantly better
constrained reconstruction methodology**, but not a derivation or unique
identification of the missing interaction law without additional physical
closure.  This conclusion introduces no empirical constant and makes no
constitutive, ontological, metric-map, or V11 modification.

## 10. Traceability

| Result used | Frozen source |
|---|---|
| one-medium ontology, complete gauge-basic state, no new ontology | FOUNDATION-001; STATE-003; LOCAL-STATE-001 |
| placement, `F`, objective `C`, and admissible configuration | DEFORMATION-001 |
| frozen stored energy, stress evaluation, and tangent | HYPER-001; CONSTITUTIVE-CONSTRUCTION-001 |
| weak balance, generalized work, boundary separation, and regional locality | BALANCE-001; LOCALITY-001; WEAK-LENSING-LOCALITY-001 |
| load codomain, regular representative, dimensions, and compatibility | NATIVE-SOURCE-001 |
| absent matter projection, isotropic tensor obstruction, and residual operator family | SOURCE-PROJECTION-001 |
| internal one-medium work-transfer interpretation | MATTER-MEDIUM-INTERACTION-001 |
| conditional mismatch/eigenstress representation | EQUILIBRIUM-MISMATCH-001 |
| conditional `y -> Pi`, optical non-identifiability, and cross-system fiber test | INVERSE-SOURCE-001 |
