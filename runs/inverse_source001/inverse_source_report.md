# PBUF INVERSE-SOURCE-001 — Reconstruction of Native Elastic Loading from Weak-Lensing Systems

## 0. Decision

The frozen framework admits a **well-posed mechanical inverse only after the
placement/deformation and boundary work are known**.  It does not admit a
well-posed end-to-end inverse from weak-lensing observations alone.

For a known admissible placement `y` on a bounded reference region `Omega`,
the frozen constitutive law fixes `P_F[y]`.  The unique total generalized load
on the zero-Dirichlet test space `V_0` is

\[
 \boxed{\langle\Pi_{\rm req}[y],\eta\rangle
 =\int_\Omega P_F[y]:\operatorname{Grad}_0\eta\,dV_0,
 \qquad \eta\in V_0.}                                      \tag{IS-001}
\]

If the Neumann traction is independently known, the unique bulk generalized
load is

\[
 \boxed{\langle b_{\rm req},\eta\rangle
 =\int_\Omega P_F[y]:\operatorname{Grad}_0\eta\,dV_0
 -\int_{\Gamma_N}\bar t\cdot\eta\,dA_0.}                  \tag{IS-002}
\]

When regular enough, `b_req=-Div_0 P_F[y]`; its natural traction is
`P_F[y]N`.  Equations IS-001--IS-002 are reconstructions, not physical
interpretations or new interaction laws.

Weak lensing does not determine `y`.  It observes a projected optical
functional of an effective metric, while the frozen corpus has not selected
the map `G` from medium deformation/state to that metric.  Even for a selected
`G`, ordinary weak-lensing data retain gauge, line-of-sight, boundary, and
lensing degeneracies.  Hence the actual frozen inverse is

\[
 \boxed{\text{set-valued and underdetermined; regularization or additional
 data can select a representative but cannot create uniqueness.}}          \tag{IS-003}
\]

Cross-system reconstruction can nevertheless test a *declared candidate
class* of universal matter-to-load operators.  It cannot nonparametrically
prove a unique universal law unless matter descriptors, the metric map,
boundary/background data, identifiability conditions, and the candidate
operator class are fixed independently.

## 1. Spaces and the complete inverse problem

### 1.1 Frozen mechanical maps

Let `Y_bar_y` be the admissible placement space, `V_0` its zero-Dirichlet
tangent/test space, and

\[
 F[y]=\operatorname{Grad}_0y,\qquad C[y]=F[y]^\sharp F[y],
 \qquad P_F[y]=2F[y]D_CW(C[y]).                              \tag{IS-004}
\]

Here `W` is the frozen interior law and the selected solution must remain in
the authorized domain.  Define the internal first-variation operator

\[
 \langle\mathcal A_\Omega(y),\eta\rangle
 :=\int_\Omega P_F[y]:\operatorname{Grad}_0\eta\,dV_0.
                                                                    \tag{IS-005}
\]

For boundary package `beta=(bar y,bar t,beta_bg)` define

\[
 \Pi_{\rm bulk}=\mathcal A_\Omega(y)-\mathcal T_{\bar t},
 \qquad
 \langle\mathcal T_{\bar t},\eta\rangle
 =\int_{\Gamma_N}\bar t\cdot\eta\,dA_0.                    \tag{IS-006}
\]

This is exactly the frozen balance rearranged.  No inverse constitutive law is
needed: stress is evaluated forward from the reconstructed deformation.

Let `S_beta` denote the (possibly local-branch) elastic solution operator,

\[
 y=S_\beta(\Pi),\qquad
 \mathcal A_\Omega(y)=\Pi+\mathcal T_{\bar t}.                \tag{IS-007}
\]

On a uniformly strongly elliptic branch with gauge/rigid modes removed,
`S_beta` is locally single-valued under the same assumptions used by the
frozen regional-sufficiency theorem.

### 1.2 Optical maps

Let `G` be an admissible effective-metric map, `R_z` the ray/bundle operator
for known lens-source-observer geometry `z`, and `M` the measurement operator
(sampling, masks, point-spread response, and the actually reported reduced
shear/magnification/position information).  Their composition is

\[
 \mathcal H_{G,\beta,z}(\Pi)
 :=M\,R_z\,G\bigl[q(S_\beta(\Pi)),C(S_\beta(\Pi));\mathcal D\bigr].
                                                                    \tag{IS-008}
\]

For observed lensing datum `d`, covariance or weighting operator `N`, and
admissible set `K`, the inverse problem is the feasibility problem

\[
 \boxed{\mathfrak P(d)=
 \{(G,\beta,y,\Pi):G\in\mathfrak G_{\rm frozen},\ y\in K,
 \ \mathcal A_\Omega(y)=\Pi+\mathcal T_{\bar t},\
 \ M R_zG[q(y),C(y)]=d\}.}                                  \tag{IS-009}
\]

With noisy data, equality is replaced by a declared acceptance region, for
example `||H(Π)-d||_{N^{-1}} <= delta`.  This does not constitute fitting in
this milestone; it states the future data interface.

If `G` and `beta` are independently fixed, one may write the representative
regularized problem

\[
 \Pi_\alpha\in\arg\min_{\Pi\in\mathcal P_{\rm adm}}
 \left\{\|\mathcal H_{G,\beta,z}(\Pi)-d\|_{N^{-1}}^2
 +\alpha\mathcal R(\Pi)\right\}.                            \tag{IS-010}
\]

`R` is a reconstruction convention (smoothness, sparsity, minimum norm, or
another declared prior), not a new constitutive term, state variable, or
physical coupling.  No particular `R`, `alpha`, or empirical coefficient is
selected here.

### 1.3 Weak-field linearization

At an unloaded background, write

\[
 \delta d=\mathcal J\,\delta\Pi,
 \qquad
 \mathcal J:=DM\,DR_z\,DG_{q_0}\,\mathcal L_0^{-1},          \tag{IS-011}
\]

where `L_0` is fixed by `(K_0,mu_0)`.  This factorization isolates all inverse
failures.  The frozen framework fixes `L_0`, but does not select `DG`; finite
and noisy optical sampling also makes `DM DR_z` non-injective.  The retained
V11 condition constrains the composite response but does not identify its
individual factors.

## 2. Identifiability theorem

### Theorem 1 — conditional load reconstruction and end-to-end non-identifiability

Under the frozen milestones:

1. Given an admissible placement `y` and complete boundary work, IS-001 gives
   one and only one total generalized load in `V_0*`.  Given `bar t`, IS-002
   gives one and only one bulk generalized load in `V_0*`.
2. A pointwise body-load representative is unique only up to equality as a
   distribution and exists only with the required regularity.  Without known
   boundary traction, bulk and boundary loading are not separately
   identifiable.
3. Given only weak-lensing observations and geometry, the native load is not
   identifiable.  The solution set contains freedom from the unselected
   metric map and, even after selecting it, from the null space of the optical
   and projection operators and from incomplete boundary/background data.
4. Therefore the frozen end-to-end inverse is underdetermined, not
   overdetermined.  Multiple measurements may make a finite candidate model
   statistically overconstrained, but they do not remove structural null
   directions.  Regularization is normally required for a stable
   representative and cannot turn a structural null direction into observed
   information.

**Proof.**  For fixed `y`, the right side of IS-001 is a specified continuous
linear functional on `V_0`, hence defines exactly one element of `V_0*`.
Subtracting specified traction proves the second statement.  If traction is
not specified, any admissible boundary functional `Delta T` produces the same
total work after the compensating change `Pi_bulk -> Pi_bulk-Delta T`.

For lensing, IS-011 shows that every element of `ker J` is invisible to first
order.  More fundamentally, `DG` is not selected: admissible pairs `(DG,
delta Pi)` can have the same composite optical response.  Optical projection
and finite sampling add further kernels.  Thus distinct loads occur in the
same data fiber.  More data can reduce a kernel only after all forward factors
and the tested model class are fixed.  A penalty chooses among elements of a
fiber but supplies no observation that distinguishes them.  QED.

### 2.1 Observable and hidden quantities

| Quantity | Status from weak-lensing data |
|---|---|
| image ellipticities / reduced shear, positions, flux or time-delay data actually measured | observable, subject to calibration and sampling |
| lens, source, observer geometry and declared background distances | supplied problem data, not inferred native state |
| optical tidal action integrated along sampled null bundles | constrained through the chosen lensing model |
| full four-dimensional effective metric | hidden; only restricted projected combinations are probed, modulo coordinate/gauge freedom |
| metric perturbations outside the sampled optical tube | hidden |
| medium placement `y`, rotation-bearing `F`, and objective deformation `C` | hidden until a particular injective `G` and sufficient boundary/gauge data are supplied |
| stresses `P_C,P_F` | derived from `C`, hence hidden at the lensing-only stage |
| generalized load `Pi`, regular body load `b`, and boundary traction | hidden; only conditionally reconstructible from `y` and boundary data |
| matter descriptor `T` | must be obtained independently for a loading-law test; lensing must not be used as both matter input and load output |
| physical origin/species partition of a reconstructed resultant load | not observable from elastic response alone |

Familiar lensing degeneracies—mass-sheet/source-position transformations,
unmeasured line-of-sight structure, finite-field boundaries, and reduced-shear
rather than absolute-shear information—are instances of non-injectivity at the
optical stage.  The theorem does not depend on any one named degeneracy.

## 3. Reconstruction target

The primary target shall be

\[
 \boxed{\Pi_{\rm req}\in V_0^*,\quad
 \text{the equivalent generalized placement-load functional}.}   \tag{IS-012}
\]

This target is already authorized, is invariant under changes of a regular or
singular representative, pairs directly with the frozen variation, and
includes loads that are not smooth body densities.  It also avoids assigning
physical origin.

Derived representations may be reported with explicit qualifications:

1. `b_req=-Div_0P_F` when it is a regular bulk covector density and traction is
   separately fixed;
2. `t_req=P_FN` for the boundary component;
3. the stress polarization `P*` only when the load is known to lie in the
   already-authorized eigenstress form `Pi=B*P*`.  Since `P*` has a
   divergence-free/gauge freedom and not every generalized load has this
   form, it is not the universal primary target;
4. a preferred strain only under an independently declared equilibrium-
   mismatch hypothesis.  It is not reconstructed as a new state variable.

Two loads are observationally equivalent for a data set when they lie in the
same fiber of `H`.  They are mechanically equivalent when they have identical
action on every admissible test variation.  These equivalence relations must
not be conflated.

## 4. Multi-system universality methodology

For independent systems `s=1,...,S`, use a common frozen convention for the
reference state, constitutive branch, gauge quotient, and definition of the
generalized-load target.

### Stage A — independent reconstruction

For every system archive separately:

\[
 D_s=(d_s,z_s,\beta_s,\text{selection/calibration metadata}),
 \qquad
 \mathfrak P_s=\mathfrak P(D_s).                              \tag{IS-013}
\]

Propagate observational, geometric, boundary, metric-map, and regularization
uncertainty into a set or posterior over `Pi_s`; never compare only one
regularized point estimate.  Use held-out lensing channels or source redshift
planes for forward validation.

### Stage B — independent matter descriptors

Construct only already-authorized matter descriptors `T_s` from observations
independent of the lensing inversion where possible.  Record frame, pullback,
support, resolution, and uncertainty.  Do not call the reconstructed load a
function of mass merely because mass was used in geometry or priors.

### Stage C — state a candidate operator class before comparison

A universal hypothesis is a **single natural operator**, not equality of raw
source maps across differently shaped systems:

\[
 H_{\rm univ}:\quad \Pi_s=\mathfrak L[T_s;q_{0,s},\beta_s]
 \quad\text{for all }s,                                      \tag{IS-014}
\]

where the same functional form, normalization, differential/support class,
covariance rule, and species rule apply to every system.  `q_0` and `beta`
are displayed so environmental or reference dependence cannot masquerade as
non-universality.  No new coefficient is introduced by this definition.

### Stage D — quotient-aware comparison

1. Pull every input and output into the corresponding intrinsic reference
   geometry; compare invariant virtual work on matched dimensionless test
   families, not coordinate components.
2. Test whether one `L` intersects every reconstruction fiber:

\[
 \boxed{\exists\mathfrak L\in\mathfrak C\quad\text{such that}\quad
 \mathfrak L[T_s;q_{0,s},\beta_s]\in
 \operatorname{proj}_{\Pi}\mathfrak P_s\quad\forall s.}      \tag{IS-015}
\]

3. Require cross-validation: determine the law/class on a training subset,
   predict the load equivalence classes and lensing of held-out systems, and
   run the frozen forward pipeline without refitting.
4. Stratify by matter type and environment to test the claimed universality,
   while retaining one resultant placement-load slot.
5. Diagnose residuals in observable space and generalized-work space.  A
   discrepancy aligned with an inverse null direction neither confirms nor
   refutes the law.

IS-015 can falsify a specified class `C`.  Agreement can support that class,
but cannot prove uniqueness against unrestricted operators: finitely many
input-output pairs admit infinitely many extensions unless locality,
regularity, symmetry, and complexity are fixed independently.

## 5. Dimensional analysis and cross-system invariants

For length-valued placement and reference volume `L^3`,

\[
 [C]=[F]=1,\qquad [W]=[P_F]=[K_0]=[\mu_0]=E L^{-3},
 \qquad [b]=E L^{-4}.                                      \tag{IS-016}
\]

The generalized-load functional has

\[
 [\langle\Pi,\eta\rangle]=E,
 \qquad [\Pi]\ \text{intrinsically means energy per placement variation}.
                                                                    \tag{IS-017}
\]

Boundary traction has `[t]=E L^-3`.  An eigenstress polarization has
`[P*]=E L^-3`, while an infinitesimal preferred strain is dimensionless.

For a system with a declared geometric length `L_s` (for example a domain or
lens scale fixed by the geometry, not a fitted coupling), natural normalized
fields are

\[
 \widehat X=X/L_s,\qquad
 \widehat P=P_F/K_0,\qquad
 \widehat b=L_s b/K_0,\qquad
 \widehat t=\bar t/K_0,\qquad
 \chi=\mu_0/K_0.                                           \tag{IS-018}
\]

The load may be compared invariantly through normalized work

\[
 \widehat{\mathcal W}_s[\widehat\eta]
 ={\langle\Pi_s,L_s\widehat\eta\rangle\over K_0L_s^3}.
                                                                    \tag{IS-019}
\]

Other dimensionless quantities already available include principal strains,
`t=tr E`, `E_TF:E_TF`, ratios of geometric distances, and normalized optical
observables such as shear.  If a matter descriptor has stress-energy dimension
`E L^-3`, then `T/K_0` is dimensionless and `L_s b/K_0` can be compared with
dimensionless spatial derivatives of `T/K_0`.  This is a dimensional test,
not authorization of a gradient loading law.

The choice of `L_s` must be declared and varied in robustness checks.  The
frozen framework contains no intrinsic material length and therefore does not
select one universal system scale.  Normalization by an observed lens scale
introduces no coupling constant, but scale-dependent collapse can be evidence
only for a specified covariant law.

## 6. Physical-interpretation framework

Only after a universal operator survives IS-015 may its mathematical
properties be classified.  Consistent interpretations within the frozen
one-medium ontology include:

| Reconstructed mathematical signature | Permitted interpretation (not a selection) |
|---|---|
| regular resultant placement load | internal matter-to-elastic generalized-work transfer |
| boundary-supported load | exchange represented through regional boundary work or excluded exterior influence |
| distributional localized load | singular/defect-like organization of the same medium, without a second substance |
| self-equilibrated `B*P*` form | stress polarization or weak-field equilibrium mismatch |
| configuration-dependent generalized load | follower interaction or finite prestrain representation |
| finite-jet dependence on `T` | local matter-medium interaction using authorized derivatives |
| support-controlled functional dependence | non-point-local interaction representation, only if its support is independently established |
| one rule for all matter sectors | universal coupling into the single resultant load channel |
| partition-dependent contributions summing to one load | species contributions within one-medium accounting |

None of these interpretations turns matter into a second medium, changes `W`,
or makes a preferred-state field independent of the complete state.  Elastic
response alone cannot distinguish physically different mechanisms that induce
the same generalized work.

## 7. Forward consistency and validation

For any candidate reconstructed operator `L`, validation uses the existing
pipeline without alteration:

\[
 T_s\xrightarrow{\mathfrak L}\Pi_s
 \xrightarrow[\beta_s]{\mathcal A_\Omega(y_s)=
 \Pi_s+\mathcal T_{\bar t_s}}y_s
 \xrightarrow{F,C}C_s
 \xrightarrow{G}g_s^{\rm eff}
 \xrightarrow{R_{z_s},M}d_s^{\rm pred}.                      \tag{IS-020}
\]

The validation protocol is:

1. verify covariance, gauge-null work, dimensions, regularity, admissible
   support, and pure-Neumann force/moment compatibility of `L[T_s]`;
2. solve the unchanged frozen weak balance with the unchanged `W`, `K_0`,
   `mu_0`, state domain, reference, and boundary data;
3. verify equilibrium residuals and that `C_s` stays in the authorized
   strongly elliptic branch;
4. apply one independently selected, frozen-compatible `G` with the retained
   V11 metric/cone gate;
5. propagate the declared null bundles and compare only to held-out optical
   observations;
6. repeat across systems with the same `L` and no system-specific
   normalization.

This reinsertion neither modifies V11 nor derives a new metric or constitutive
law.  The exact weak-field consistency gate remains

\[
 [DG_{q_0}\mathcal L_0^{-1}D\mathfrak L_{T_0}[\delta T]]_{\rm gauge}
 =[h^{\rm V11}[\delta T]]_{\rm gauge}.                       \tag{IS-021}
\]

## 8. Readiness assessment

Three logically different outcomes must be separated.

1. **Calibration only.**  If reconstructions depend materially on an
   unselected `G`, regularizer, boundary model, or lensing null direction,
   they are equivalent-source maps useful for forward calibration, not an
   identified interaction law.
2. **Constraint on the interaction postulate.**  If multiple independent
   systems, independently measured matter descriptors, and a fixed forward
   closure restrict the same covariant operator class, reconstruction can
   reject support, tensor, species, sign, or functional alternatives.  This
   is the realistic success level of the present methodology.
3. **Replacement of an independent postulate.**  This would require the data
   and frozen structures to identify one operator uniquely (including
   normalization, index bridge, locality/support, boundary separation, and
   species universality) and to validate it out of sample.  Even then the law
   is empirically reconstructed, not derived from the prior ontology.  Its
   adoption remains an empirical interaction principle, although no separate
   guessed formula would be needed.

Consequently successful multi-system reconstruction would ordinarily
**constrain and empirically calibrate** the missing interaction postulate.  It
does not logically eliminate the need to state an interaction law.  Only a
uniqueness result over a predeclared admissible operator class could replace
an independently conjectured postulate by an empirically identified one.

## 9. Completion theorem

### Theorem 2 — status of INVERSE-SOURCE-001

The frozen PBUF framework supplies a canonical conditional inverse
`y -> Pi_req` through its constitutive first variation and balance.  Its most
natural target is the generalized placement-load covector.  The framework
does not supply an injective map from observed weak lensing to `y`, because
the effective metric map is unselected and the optical observation operator
has nontrivial equivalence classes.  Thus it does **not presently admit a
mathematically well-posed reconstruction of native loading from weak-lensing
systems alone**.

With an independently fixed metric map, adequate boundary/background data,
additional observations sufficient for the declared candidate class, and
explicit regularization, it admits a conditional set-valued reconstruction.
Cross-system comparison can then falsify or constrain universal operator
classes via IS-015 and validate survivors through IS-020.  It cannot infer a
unique unrestricted loading law from finitely many lensing systems.

This milestone introduces no state variable, empirical coupling, fitted
value, constitutive law, ontology, or V11 modification.

## 10. Traceability

| Result used | Frozen source |
|---|---|
| one-medium ontology and gauge-basic physical state | FOUNDATION-001; LOCAL-STATE-001 |
| objective dimensionless `C=F^sharp F` and placement realization | DEFORMATION-001 |
| frozen hyperelastic variation and tangent | HYPER-001; CONSTITUTIVE-CONSTRUCTION-001 |
| generalized balance and boundary work | BALANCE-001; LOCALITY-001 |
| local elliptic regional solve and missing metric locality/selection | WEAK-LENSING-LOCALITY-001; LOCAL-STATE-001 |
| load codomain, work pairing, dimensions, and regular representative | NATIVE-SOURCE-001 |
| missing matter-to-load operator and composite V11 gate | SOURCE-PROJECTION-001; MATTER-MEDIUM-INTERACTION-001 |
| conditional eigenstress/mismatch representation | EQUILIBRIUM-MISMATCH-001 |
