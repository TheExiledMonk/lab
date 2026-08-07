# PBUF FND-003 — Three-dimensional microscopic state justification

## Result: Outcome C

The existing ontology is compatible with a three-component spatial-vector microscopic state, but it does not uniquely require one. Three spatial dimensions determine the dimension of tangent vectors, not the representation carried by an otherwise unspecified microscopic state. Exactly three follows only after adding that the state transforms linearly and faithfully as a spatial vector and independently resolves all directions. Those are new postulates, so Outcome A and Outcome B are not justified.

## Ontology-to-mathematics mapping

| ID | Claim | Classification | Source | Reason/boundary |
|---|---|---|---|---|
| FND-003-T01 | Physical space has three spatial dimensions | working premise | FND-003 mission; stated PBUF ontology | Taken as the starting ontology, not derived by this milestone. |
| FND-003-T02 | The microscopic state exists | explicitly stated in PBUF | CORE-001 and FND-002 | A PBUF premise; existence alone fixes neither representation nor component count. |
| FND-003-T03 | The microscopic state has exactly three real components | new assumption introduced after V11 | CORE-001 working premise | Three-dimensional space permits scalar, vector, tensor, and larger microscopic states; dimension matching is not a theorem. |
| FND-003-T04 | Each component is associated with one spatial direction | new assumption introduced after V11 | FND-003 candidate interpretation | Requires q to be a spatial vector/covector and a choice of basis; components are basis-dependent, while q is geometric. |
| FND-003-T05 | A faithful linear SO(3) vector realization needs at least three real components | mathematically derived | conditional representation argument | SO(3) has no faithful real representation in dimensions 1 or 2; its defining three-dimensional representation is faithful. The premise that q must be faithful and linear is additional. |
| FND-003-T06 | The three-component realization is unique | not established | model comparison | Even after selecting the vector representation, dynamics, parity, locality, couplings, and field content remain open; non-vector realizations also satisfy the base ontology. |
| FND-003-T07 | g_dev=1/137 is an intrinsic microscopic coupling | working premise | CORE-001 and FND-002 | No supplied symmetry, normalization, action, or renormalization prescription derives or operationally isolates it. |
| FND-003-T08 | g_dev directly normalizes the corrected CORE-001 linear source | mathematically derived | corrected CORE-001 energy | The source term is -epsilon_* g_dev eta e.q; no independent coupling multiplier or inverse-rescaling freedom remains. |
| FND-003-T09 | The current theory identifies 1/137 as a separately measurable parameter | not established | corrected identifiability audit | The normalized microscopic equation is sensitive directly to g_dev, but no supplied principle derives its numerical value or links it to an independently completed downstream observable. |
| FND-003-T10 | u=C_L[e.q] is a scalar continuum deformation | mathematically derived | CORE-001 conditional map | It is rotationally scalar only when q and the matter-selected e transform together and the kernel is isotropic and normalized. |
| FND-003-T11 | A direction-free nonzero linear map from a spatial vector q to a scalar u exists | not established | SO(3) invariance | No nonzero SO(3)-invariant linear functional on the vector representation exists; one must add e, use derivatives such as div(q), or use a nonlinear invariant such as \|q\|. |
| FND-003-T12 | The long-wave scalar equation K u-div(G grad u)=s(rho) follows | mathematically derived | CORE-001/FND-002 | Conditional on the projection, stable analytic energy, locality, isotropy, scale separation, and decoupling of other modes; coefficients and source law are not predicted. |

## Three-dimensional state definition

The minimal conditional realization is `q(x) in T_x Sigma (or T_x^* Sigma), with dim(Sigma)=3`. In an orthonormal frame its coordinates are `q^a`, `a=1,2,3`, with `q'^a=R^a_b q^b for R in SO(3)`. The directions label geometric basis components, not three invariant substances. Changing frame mixes them. No extra internal degree of freedom is introduced in this realization, but excluding scalar, tensor, spinorial, or additional sectors is itself a model choice.

The conditional representation argument in `mathematical_derivation.md` proves that three is minimal under the faithful-linear-vector premises. It also records why those premises do not follow from dimensionality alone.

## Mapping to continuum deformation

| Map | Definition | Advantage | Boundary |
|---|---|---|---|
| Matter-selected projection | u(x)=integral W_L(x-y) e(y).q(y) dy | Linear and matches CORE-001 | Requires an additional transforming vector e; a fixed background e breaks isotropy |
| Longitudinal scalar | u(x)=L integral W_L(x-y) div(q(y)) dy | Rotation scalar without an internal direction | Introduces a derivative and length L; does not equal the CORE-001 map |
| Magnitude scalar | u(x)=integral W_L(x-y) \|q(y)\| dy | Direction-free rotation scalar | Nonlinear, nonnegative, and changes the weak-field/source expansion |

The existing CORE-001 interface selects the first map. With normalized isotropic `W_L`, a transforming unit vector `e`, scale separation, stable local response, and decoupled transverse modes, it yields the scalar effective energy and stationary equation `K u-div(G grad u)=s(rho)`. This is a conditional coarse-graining derivation, not a derivation of the vector ontology or the coefficients. A spatially fixed `e` would add preferred-direction structure absent from the base ontology.

## Treatment of 1/137

In the corrected microscopic energy, `g_dev=1/137` appears directly in the matter vertex. The previous inverse-rescaling argument was an artefact of an auxiliary modelling choice and is withdrawn. The current record therefore classifies `1/137` as a direct working premise, not a mathematically derived value. Deriving it requires a microscopic principle that fixes the coupling and its applicable scale behavior, plus a completed response/access chain for empirical identifiability.

## Remaining irreducible postulates

P1. The microscopic state carries the defining spatial-vector (or covector) representation of SO(3), rather than a scalar, tensor, spinorial, or unrelated internal representation. Needed for: identifying three components with three spatial directions.
P2. The rotation action is linear and faithful and all three directional responses are independent; no additional microscopic sectors are required. Needed for: conditional minimality of three.
P3. A covariant scalarization mechanism is selected: a matter-provided vector e, a longitudinal derivative, or a nonlinear invariant. Needed for: mapping q to the scalar continuum field u.
P4. The stable, local, isotropic long-wave expansion and decoupling assumptions of CORE-001/FND-002 hold. Needed for: recovering K u-div(G grad u)=s(rho).
P5. If g_dev=1/137 is to be derived rather than postulated, a microscopic principle must fix its value and specify any applicable scale behavior. Needed for: turning 1/137 from a direct PBUF premise into a prediction.

## Recommendation for FND-004

FND-004 should formulate the covariant microscopic action and representation content. It should choose and justify the SO(3) representation, derive rather than select the scalarization/source coupling, test whether a matter-provided `e` creates forbidden preferred-direction effects, and state any applicable scale behavior for `g_dev`. Its gates should require (1) a unique representation or an explicit postulate, (2) a rotation-covariant map to `u`, (3) a mode-decoupling calculation, and (4) a direct, operationally testable coupling. No empirical fitting or weak-lensing changes should occur until those gates pass.

## Completion checks

All checks pass: **True**. The work is theory-only and imports or modifies no frozen-laboratory module.
