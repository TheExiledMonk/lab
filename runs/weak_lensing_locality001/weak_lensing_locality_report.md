# PBUF WEAK-LENSING-LOCALITY-001 — Minimal Local Information Required for Weak Lensing

## 0. Decision and scope

The frozen architecture gives a **local static elastic boundary-value problem**, but it does **not yet give a closed local weak-lensing problem**.

More precisely:

1. In the authorized placement realization, the internal response at (X) depends on the local first jet of (y), and static balance is a second-order elliptic system on any bounded operating region on which strong ellipticity holds.
2. The influence of the excluded exterior on that regional elastic solution can be represented by admissible boundary displacement or traction data (and by prescribed body loading in the region). The complete instantaneous state of the universe is therefore not required to solve the **native static deformation subproblem**.
3. End-to-end weak lensing is not presently computable because two maps remain open: the projection of a specified mass-energy distribution into the native source (b), and selection of the medium-to-metric map (G). In particular, METRIC-001 admits both finite-jet local maps and functional kernels of unrestricted spatial support. The frozen results therefore do not prove that (g^{\rm eff}) in a bounded region is determined by regional state data and boundary traces.
4. Once a local source projection and a local finite-jet (G) satisfying the frozen V11 gates are supplied, a bounded weak-lensing calculation needs only regional source data, frozen constitutive data, admissible boundary data, the local reference/background prescription, and the metric along the relevant null bundle. It does not need the complete universal state (q).

Thus the answer to the central question is conditional:

\[
\boxed{\text{bounded local sufficiency is proved for static native deformation, but not yet for complete weak lensing.}}
\tag{WL-001}
\]

This report adds no ontology, state variable, microscopic constituent, coefficient, fit, evolution law, or modification of V11.

## 1. Complete information-sufficiency audit

The classifications below concern a static or quasistatic weak-lensing calculation around one isolated or weakly coupled object. “Required” means required in every such calculation. “Conditional” means required only for a selected realization, boundary treatment, or time-dependent calculation. “Unnecessary” means the quantity need not be supplied as independent computational input, even when it exists globally or is derived during the solve.

| Frozen quantity or structure | Status | Justification for one-object weak lensing |
|---|---|---|
| One continuous spacetime medium (mathcal M) and FP-1--FP-6 | required as constraints, not field data | They constrain admissible models and preserve the one-metric/V11 limit; a numerical solver does not require a separate value for the global ontology. |
| Complete universal state (q=[\widehat q]_{\mathcal G}\in\mathcal Q_{\rm phys}) | unnecessary as a complete global input | Ontically every observable is a functional of (q), but functional dependence does not imply that the whole argument must be known. The local placement realization and boundary trace can determine the regional elastic solution. End-to-end local sufficiency remains conditional on the support of (G). |
| Restriction (q|_\Omega), represented locally by (y|_\Omega) | required | This is the unknown regional medium state in the authorized bounded-domain realization. It is solved for, not prescribed throughout. |
| Gauge group (mathcal G=\operatorname{Diff}(\mathcal M)\ltimes\mathcal G_{\rm int}) | required as quotient/gauge handling | A representative or gauge condition is needed to solve component equations and remove rigid/gauge null modes; gauge-related data are not additional physical information. |
| Unloaded comparison state (q_0), (C[q_0,q_0]={\bf1}) | required | Deformation is relational. A regional reference/background prescription is needed to form (C). The complete exterior realization of (q_0) is not needed when its restriction and boundary compatibility are specified. |
| Choice among fixed, instantaneous-natural, or cosmological reference families | conditionally required | It matters when cosmological expansion could be counted as strain. For an isolated calculation, one must declare which already-admissible reference/background is used; no new family is selected here. |
| Reference carrier (mathcal B_0), region (Omega\Subset\mathcal B_0), measure (dV_0), operators (operatorname{Grad}_0,operatorname{Div}_0) | required | They define the local placement BVP and its weak form. |
| Placement (y(X)) and (F=\operatorname{Grad}_0y) | required in the selected local realization | They provide the explicit local unknown and first derivative. Intrinsically one may instead use (q) and ((D_qC)^*), but some concrete realization is mandatory. |
| Objective deformation (C=F^\sharp F\in\overline{\mathcal D_C}\subset\operatorname{Sym}^+(3)) | required | It is the argument of the frozen stored energy and the input to the metric-map family. |
| (E=(C-{\bf1})/2), (t=\operatorname{tr}E), (E_{\rm TF}) | required only as derived coordinates | They are convenient for the selected constitutive formula, but are not independent input or state variables. |
| (I_1,I_2,I_3) or principal stretches | unnecessary as separate input | They are algebraic functions of (C). They are alternative invariant coordinates, not extra information. |
| Admissible spectral domain (mathcal D_C) and its hard endpoint | conditionally required | The domain is required for admissibility in every solve; its endpoint reaction matters only if the solution reaches the boundary. Weak lensing in the weak-deformation interior does not activate it. |
| Frozen moduli (K_0>0,mu_0>0) | required | They fix the selected minimal interior response and reference strong ellipticity. No numerical values are introduced here. |
| Selected Candidate-A/B interior energy (Q(C)=K_0t^2/2+\mu_0E_{\rm TF}:E_{\rm TF}) | required | It closes the internal elastic operator. A and B have identical interior equations; only endpoint bookkeeping differs. |
| (P_C=DW), (P_F=2FP_C), and tangent/acoustic form | required as derived response | These give stress transmission, the regional PDE, ellipticity checks, and natural traction. They are computed from (C,K_0,mu_0), not supplied independently. |
| Normal cone (N_{\mathcal K}) / tangent-cone inequality | conditionally required | It is active only at the hard admissibility boundary and is not a new multiplier field. |
| Static balance (-\operatorname{Div}_0P_F=b) and weak form | required | This is the closed local constitutive-balance skeleton for a prescribed native source and boundary data. |
| Native body source (b) or source covector (mathcal S) | required | A source must load the medium. Its value in (Omega) is problem data, but the map from physical mass-energy to (b) is not frozen. |
| Boundary partition `Gamma_D union Gamma_N`, displacement `bar y`, traction `bar t` | required up to an admissible well-posed choice | A bounded elliptic problem needs sufficient essential/natural data and removal of rigid/gauge modes. Both full Dirichlet and admissible mixed data are allowed; pure Neumann requires compatibility and quotienting rigid modes. |
| Kinetic/inertial operator (mathcal K_\tau), momentum, and initial data | unnecessary for a strictly static calculation; conditionally required otherwise | They remain open and are needed for retardation, transients, waves, or time-dependent lenses, but not for the frozen static equilibrium equation. |
| Ordered history `q(s)`, duration functional, and calibrated duration `tau` | unnecessary for static deformation; conditionally required for propagation/time dependence | The order label is never physical time. A stationary metric and a V11 null trajectory can be treated without solving universal evolution; time-dependent lensing needs duration and kinetic closure. |
| Balance currents, total energy, momentum, Noether charges | unnecessary for static lensing unless used as diagnostics | BALANCE-001 does not freeze nontrivial universal charges. They are not additional lensing inputs. |
| Effective metric family (g^{\rm eff}=G[q,C;\mathcal D]) | required | Weak lensing is operationally computed from the Lorentzian metric and its null cone. The type and constraints of (G) are frozen, but no member is selected. |
| Locality class, derivative order, and normalization of (G) | required closure, presently missing | Without these, regional deformation does not determine regional metric perturbation. This is the decisive locality gap. |
| V11 weak matching ([D_qG|_{q_0}\delta q]_{\rm gauge}=[h^{\rm V11}]_{\rm gauge}), Lorentzian signature, and cone match | required as gates | They constrain any selected metric map and photon propagation without selecting that map. |
| Clock/ruler calibration and synchronization convention | conditionally required | Proper null curves are invariant, but angular positions, distances, and time-dependent observables require an observer/chart calibration. This is operational data, not microscopic state. |
| Photon trajectory/null propagation rule | required but not a new PBUF state | In the retained V11 regime, photons follow the universal null propagation structure of the one effective metric. The relevant ray/bundle initial or endpoint conditions must be specified. |
| Source, lens, and observer event/position data | required problem data | They choose which null bundle and observable are computed; they are not universal state variables. |
| Cosmological homogeneous variables and complete saturation history | unnecessary for an isolated local calculation | They are required only if the desired observable includes cosmological distance/background evolution. Frozen milestones prohibit identifying them directly with local strain. |
| External masses and tidal environment outside (Omega) | conditionally required only through boundary/background data | Weak coupling can be summarized by admissible boundary displacement/traction or a declared background metric. Their full internal states are unnecessary. |
| Full past history of the universe | unnecessary | STATE-003 distinguishes a present state from a history, and the frozen elastic energy is state-local and history-free. |
| Rank-four clock/ruler realization | unnecessary for the selected native elastic solve; conditionally required by a future concrete (G) | GOVERNING-EQUATION-001 uses the authorized rank-three placement realization. A rank-four realization may not be silently added. |
| Barrier profile, gradient length, nonlocal constitutive kernel, dissipation, internal variables | unnecessary and not frozen in the canonical baseline | LOCALITY-001 and CONSTITUTIVE-CONSTRUCTION-001 eliminate them from the minimal internal operator. |

## 2. Local continuum audit

Let (Omega\subset\mathcal B_0) be bounded and let

\[
\mathcal E_\Omega[y]=\int_\Omega W_A(F^\sharp F)\,dV_0.
\tag{WL-002}
\]

For every admissible variation (eta),

\[
\delta\mathcal E_\Omega[y;\eta]
=\int_\Omega P_F:\operatorname{Grad}_0\eta\,dV_0
=-\int_\Omega\operatorname{Div}_0P_F\cdot\eta\,dV_0
+\int_{\partial\Omega}(P_FN)\cdot\eta\,dA_0.
\tag{WL-003}
\]

This identity establishes the continuum structure:

| Sector | Bounded-domain status |
|---|---|
| deformation | local: (C(X)) is obtained from the first jet (j_X^1y) |
| stress | local: (P_C(X)=DW(C(X))), (P_F(X)=2F(X)P_C(X)) |
| constitutive response | closed in the weak interior by (K_0,mu_0,Q) and the hard domain |
| balance | local second-order divergence equation with body loading and boundary traction |
| effective metric | not closed: (G) may be a local finite-jet operator or a nonlocal functional kernel |

Near the unloaded state, the static operator is

\[
-\mu_0\Delta u-\left(K_0+{\mu_0\over3}\right)
\nabla(\nabla\!\cdot u)=f,
\tag{WL-004}
\]

which is strongly elliptic because (K_0,mu_0>0). At finite deformation the same conclusion holds only on the declared connected region where the lifted Legendre--Hadamard form is positive. Consequently the architecture naturally admits local continuum mechanics, but a theorem of global existence or uniqueness on the entire finite-deformation domain is not frozen.

## 3. Locality theorem

**Theorem 1 (regional elastic sufficiency).** Let (Omega\Subset\mathcal B_0) be bounded. Assume an admissible placement realization, prescribed (b|_\Omega), admissible Dirichlet/Neumann data on (partial\Omega), and a solution branch wholly contained in a uniformly strongly elliptic subset of (operatorname{int}\mathcal D_C). After gauge/rigid modes are fixed or quotiented, any locally unique static solution in (Omega) is determined by

\[
\mathfrak I_\Omega=
\{\Omega,q_0|_\Omega,K_0,\mu_0,\mathcal D_C,
b|_\Omega,\bar y|_{\Gamma_D},\bar t|_{\Gamma_N}\},
\tag{WL-005}
\]

and not by an independently supplied complete exterior or universal state.

**Proof.** The weak equation on (Omega) is

\[
\int_\Omega P_F(y):\operatorname{Grad}_0\eta\,dV_0
=\langle b,\eta\rangle+\int_{\Gamma_N}\bar t\cdot\eta\,dA_0
\tag{WL-006}
\]

for every zero-Dirichlet admissible variation. Every coefficient and term in this equation belongs to the regional data set (WL-005) or is an algebraic/local differential function of the unknown regional placement. No value of the placement, deformation, or stress outside the closure of `Omega` occurs. If two global states induce the same data in (WL-005), their regional restrictions solve the same weak problem; local uniqueness on the quotient makes those restrictions identical. Thus exterior distinctions not represented in the boundary data cannot alter the regional solution. ∎

This theorem is an information statement, not a claim that elliptic influence is point-local. A change of source anywhere in (Omega), or a change of boundary data anywhere on (partial\Omega), can change the solution throughout (Omega).

**Corollary.** One global continuous medium does not imply that a regional computation requires the whole medium's instantaneous state. The medium's exterior influence enters through boundary traces exactly as it does in classical continua.

## 4. Boundary-condition analysis

### 4.1 Infinite isolated idealization

An isolated solution may be posed by requiring approach to the unloaded/background class:

\[
C\to{\bf1},\qquad P_F\to0
\quad\text{toward the asymptotic exterior},
\tag{WL-007}
\]

with a representative/gauge normalization sufficient to remove rigid modes. This is a symbolic asymptotic condition, not a numerical falloff law. A falloff rate is not frozen and must not be invented.

### 4.2 Finite computational domain

The weakest frozen boundary vocabulary is:

* essential data `y = bar y` on `Gamma_D`;
* natural data `P_F N = bar t` on `Gamma_N`;
* an admissible mixed partition with enough constraint to remove gauge/rigid null modes; or
* pure traction data only when force/moment compatibility holds and solutions are taken modulo the corresponding rigid modes.

For an isolated truncation, one may prescribe boundary values induced by the undeformed/background exterior or use traction-free data only where that approximation is physically declared. The frozen theory does not establish that a finite traction-free boundary exactly represents infinity. Robin or absorbing conditions are not native consequences: they require a boundary constitutive law, boundary energy, or dynamic closure.

### 4.3 External loading and weak coupling

The exterior need not be reconstructed. Its effect can be supplied as:

\[
\bar y_{\rm ext}\text{ on }\Gamma_D,
\qquad
\bar t_{\rm ext}\text{ on }\Gamma_N,
\qquad\text{and/or}\qquad
b_{\rm ext}|_\Omega.
\tag{WL-008}
\]

These data may represent an external tidal environment, but the frozen framework supplies neither their values nor a multipole formula. “Weakly coupled” justifies treating them as prescribed problem data; it does not authorize discarding them.

### 4.4 Cosmological background

There are two mathematically distinct permitted treatments:

1. **Background as reference/boundary data:** choose an already-admissible homogeneous/isotropic (q_0) prescription and solve for a regional perturbation relative to it.
2. **Cosmological observable:** if angular-diameter distances or evolution along a long ray are part of the requested observable, the background effective metric along that ray is additional problem data and generally lies outside the isolated local solve.

The first does not require the complete cosmological state. The second requires the relevant background solution along the optical path, not the complete instantaneous state of the universe. No frozen result permits identifying V11 homogeneous quantities with local (C) or supplying their evolution from this milestone.

## 5. Domain-of-dependence results

Two meanings must be separated.

**Theorem 2 (constitutive support).** For the canonical baseline,

\[
\operatorname{supp}\delta\mathcal A[y]
\subseteq\operatorname{supp}\eta
\tag{WL-009}
\]

up to the boundary trace generated by integration by parts. Equivalently, (mathcal A[y](X)) is a differential operator of spatial order two whose coefficients depend on (j_X^1y). It has no integral kernel or dependence on distant constitutive state.

**Proof.** `W_A` is pointwise in `C`, `C` is pointwise in the first jet of `y`, and its Euler operator is `-Div_0 P_F`. Differential operators do not enlarge support except through the closure of the support and boundary traces. ∎

This theorem proves **local constitutive dependence**, not a finite-radius static domain of dependence. For an elliptic equation the solution at a point generally has a Green representation involving all regional sources and the whole boundary:

\[
u(X)=\int_\Omega \mathbb G(X,Y)f(Y)\,dY
+\int_{\partial\Omega}\mathbb B(X,Y)[\text{boundary data}]\,dA_Y.
\tag{WL-010}
\]

Thus “nearby” means data in the selected bounded problem, not merely an arbitrarily small neighborhood of (X).

**Theorem 3 (no frozen end-to-end domain-of-dependence theorem).** The frozen architecture does not imply that (g^{\rm eff}(X)), and hence lensing at (X), depends only on (q|_\Omega) and boundary traces.

**Proof.** METRIC-001 permits

\[
\delta g^{\rm eff}_{\mu\nu}(X)
=\int R_{\mu\nu A}(X,Y)\,\delta q^A(Y)\,d\mu(Y).
\tag{WL-011}
\]

No support restriction on `R(X,·)` is frozen. Choose two admissible global states that agree on the closure of `Omega` and induce identical elastic boundary data but differ on a set `U` outside `Omega`. An admissible kernel with `R(X,Y) != 0` for `X` in `Omega` and `Y` in `U` can assign different effective metrics at `X` to the two states. Therefore regional metric equality does not follow from the frozen assumptions. ∎

If a later authorized closure selects a finite-jet (G) of order (r), then (g^{\rm eff}(X)) depends only on (j_X^rq,j_X^rC). Boundary data must then provide the regularity and derivative traces needed to evaluate that jet. If it selects a causal nonlocal (G), its kernel support—not global ontology—determines the required domain.

For dynamics, finite propagation and a causal domain of dependence cannot be proved until (mathcal K_\tau) is selected as a positive hyperbolic closure and its characteristic cone is matched to (g^{\rm eff}). The frozen static result does not require such an evolution law.

## 6. Comparison with classical field theories

| Theory | Local law | Why distant matter can matter | Regional data needed | Comparison with PBUF |
|---|---|---|---|---|
| elasticity | local stress-strain law and elliptic/hyperbolic balance | solution propagates boundary/source influence through the PDE | material law, regional loads, boundary/initial data | PBUF's frozen internal operator has exactly this local first-gradient structure |
| electromagnetism | local Maxwell equations | constraints and retarded fields carry source/boundary influence | regional sources plus boundary/initial/radiation data | one global electromagnetic field does not require instantaneous knowledge everywhere; PBUF is analogous at the field-equation level once its missing maps are selected |
| fluid mechanics | local balance and constitutive fluxes | pressure constraints, inflow/outflow, and initial data couple a domain | equation of state/transport laws, regional and boundary/initial data | PBUF likewise separates constitutive response, balance, and problem data, but has no frozen kinetic closure |
| General Relativity | local covariant field equations; constraints on hypersurfaces | boundary/asymptotic data and causal past determine a solution branch | stress-energy, constraint-satisfying initial/boundary data, gauge | PBUF retains the one-metric V11 limit but has not frozen the native source-to-metric field equation |

The presence of one global continuous medium is no more an obstruction to local field solutions than the existence of one spacetime metric, one electromagnetic field, or one connected elastic body. Locality is a property of operators and their support, not of whether the ontology contains one connected entity. PBUF matches classical continuum locality in its elastic sector; it has not yet matched GR's closed source-to-metric sector.

## 7. Minimal weak-lensing information flow

The smallest honest chain is

```text
prescribed mass-energy distribution in/near Omega
        |
        v   [MISSING: native source projection]
native body/source covector b|Omega + external boundary loading
        |
        v   [FROZEN: static local balance + W_A]
regional placement y -> C -> P_F
        |
        v   [MISSING: selected G and its support/normalization]
effective metric g_eff on a tube containing the relevant rays
        |
        v   [RETAINED V11: universal null propagation]
null ray / optical bundle with source-observer endpoint data
        |
        v
weak-lensing observables (angular distortion, convergence/shear,
and any requested timing/magnification quantity)
```

Information requirements by step are:

1. **Mass-energy to native load:** regional mass-energy is sufficient only after the source projection is supplied. The complete universal distribution is unnecessary if excluded matter is represented by boundary/background data.
2. **Load to deformation:** bounded regional data suffice under Theorem 1.
3. **Deformation to metric:** bounded regional data suffice only if the selected (G) has local/controlled support. This is not frozen.
4. **Metric to photon trajectory:** only the metric on a neighborhood of the relevant null ray or bundle, plus endpoint/initial ray data, is needed. A complete metric on all spacetime is unnecessary.
5. **Trajectory to observable:** observer/source calibration and, for cosmological observables, the relevant background distance information are needed; the full universe is not.

No step logically requires “the complete instantaneous state of the universe” as such. The current framework nevertheless cannot prove that the metric step has bounded support, so it cannot yet replace that global argument by a finite data set.

## 8. Computational feasibility assessment

### What can be computed now

Given symbolic (K_0,mu_0,mathcal D_C), a prescribed native source (b), and admissible boundary data, the following are numerically formulable without new physics:

* the static weak problem (WL-006);
* its reference-linearized strongly elliptic form (WL-004);
* the finite-deformation variational inequality when the hard domain is active;
* (C,P_C,P_F) and admissibility/ellipticity diagnostics.

Standard finite-element or other elliptic discretizations are structurally applicable. Actual well-posedness still requires a nonempty weakly closed admissible class, lower semicontinuity/coercivity in the chosen realization, compatible loads, and operation inside a uniformly elliptic branch.

### What cannot be computed from the frozen data

A weak-lensing prediction cannot yet be produced, even symbolically as a unique functional of a given mass-energy distribution, because the composite response

\[
\boxed{
\mathscr H: T^{\rm matter}
\longmapsto b
\longmapsto q|_\Omega
\longmapsto G[q,C]
\longmapsto g^{\rm eff}
}
\tag{WL-012}
\]

is undefined in its first and third arrows. The decisive missing mathematical object for **local lensing** is therefore a support-controlled, normalized source-to-metric response operator (mathscr H). In the frozen factorization it consists of:

1. the native source projection (T^{\rm matter}\mapsto b); and
2. one selected (G\in\mathfrak G), including its locality/kernel support, derivative order, normalization, and V11 cone match.

Calling only (G) missing would ignore that a physical mass distribution cannot yet be inserted as (b). Calling only the source map missing would ignore that solved deformation cannot yet be converted into a unique lensing metric. Kinetic closure is an additional obstruction only for time-dependent lensing, not for a static lens.

## 9. Global ontology versus local prediction

The categories are mathematically distinct:

| Global ontology / full-universe description | Local prediction / one lens |
|---|---|
| (q\in\mathcal Q_{\rm phys}), a global gauge class on (mathcal M) | equivalence class of restrictions (q|_\Omega) that produce the same regional source, boundary traces, and relevant metric tube |
| complete global history ([\gamma]) | no history for a static solve; local initial history data only if dynamics is selected |
| global reference/background realization | (q_0|_\Omega) plus compatible boundary/background prescription |
| all external sources and fields | their induced regional body load, boundary traction/displacement, and relevant background metric |
| global effective metric, if a closure produces one | (g^{\rm eff}) only on a neighborhood of the source-observer null bundle |
| homogeneous cosmological evolution and global saturation history | relevant background metric/distances only when the observable spans cosmological scales |

Formally, let (pi_\Omega) map a global state to the regional data tuple in (WL-005), augmented by whatever metric data a selected (G) requires. A local observable has the form

\[
\mathcal O_{\rm WL}=\widetilde{\mathcal O}_{\rm WL}\circ\pi_\Omega
\tag{WL-013}
\]

exactly when it is constant on each fiber of (pi_\Omega). Theorem 1 proves this factorization for regional static deformation. Theorem 3 proves that the frozen admissible family (mathfrak G) does not guarantee it for (g^{\rm eff}) or lensing. This is the precise mathematical distinction between ontic global completeness and predictive local sufficiency.

## 10. Quantities explicitly unnecessary for local static weak-lensing prediction

Subject to the missing local source and metric closures identified above, the following are not independent inputs to a bounded one-object calculation:

* the complete instantaneous state (q) outside the computational region and relevant optical tube;
* the complete past or future history of the universe;
* a global evolution law for the universe;
* global values of (C,W,P_C,P_F) away from the region;
* all three invariants in addition to (C), or strain coordinates in addition to (C);
* kinetic momentum, inertia, initial velocity, and duration calibration for a strictly static lens;
* a barrier profile when the weak solution remains in the admissible interior;
* a gradient constitutive term, kernel communication law, intrinsic length, dissipation, memory, or internal variables;
* a rank-four primitive elastic state;
* universal Noether charges or a total universal energy;
* the internal state of distant matter whose influence is already encoded in boundary/background data;
* the global metric away from the relevant null bundle;
* V11 saturation history or a reconstruction of homogeneous cosmological state, unless the selected observable explicitly requires cosmological background propagation.

## 11. Completion statement

The frozen PBUF framework possesses the same bounded-domain locality as classical first-gradient elasticity for its native static deformation sector. A single global medium does not force use of the complete instantaneous universal state: regional loads plus admissible boundary data are sufficient for that sector.

The framework does not yet possess a complete local weak-lensing boundary-value formulation. The source projection and support-controlled effective metric map are absent, and the admissible metric family still includes nonlocal kernels. Therefore the strongest authorized conclusion is neither “the whole universe is required” nor “local weak lensing is already closed.” It is:

\[
\boxed{
\begin{array}{c}
\text{local native deformation: yes;}\\
\text{complete local weak-lensing prediction: conditional and presently unclosed;}\\
\text{complete-universe knowledge: not intrinsically required once the missing maps localize.}
\end{array}}
\tag{WL-014}
\]
