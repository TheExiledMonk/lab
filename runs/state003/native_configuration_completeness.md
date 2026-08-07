# PBUF STATE-003 — Completeness of the Native Configuration Variable

## 0. Decision

The frozen milestones select **Outcome A in the ontic sense**, with an essential
qualification:

\[
\boxed{
q\text{ represents the entire instantaneous physical state of the one medium,}
\quad q\not\equiv\text{placement data alone}.}
\tag{S3-001}
\]

This follows directly from STATE-002 S2-002--S2-003, where it is an accepted
definition rather than a derived mechanical theorem. There is no additional
instantaneous physical state information outside \(q\). Any later instantaneous
physical observable must be a gauge-invariant functional of \(q\), or a
relational functional of \(q\) and an explicitly declared comparison state.

However, the frozen milestones do **not** establish either of the stronger
claims

\[
 q\text{ has a known component realization},\qquad
 q\text{ is sufficient Cauchy data for unique continuation}.
\tag{S3-002}
\]

Ontic completeness, coordinatization completeness, and dynamical/Markov
completeness are distinct. Only the first is frozen. Accordingly, Outcome B is
false as a characterization of the accepted ontology, while Outcome C applies
only to the unproved stronger questions in S3-002—not to whether \(q\) denotes
the complete physical state.

The opposite-wave argument in EVOLUTION-PRINCIPLE-001 does not follow from the
frozen ontology. If opposite propagation orientations are physically distinct
at one order state, completeness requires them to correspond to distinct
\(q\)'s, even though the frozen realization does not show how. Two classical
solutions with the same placement but opposite velocity show only that
**placement** is incomplete; they do not show that the PBUF \(q\) is incomplete.
The core successor underdetermination survives independently through STATE-002
S2-014 and the absence of a history-selection law.

## 1. Formal meaning of complete

STATE-002 fixes

\[
 \mathcal Q_{\rm phys}=\mathcal A/\mathcal G,
 \qquad q=[\widehat q]_{\mathcal G},
 \tag{S3-003}
\]

where \(\mathcal A\) consists of admissible **complete-configuration
representatives** and \(\mathcal G\) contains accepted descriptive
redundancies. Its explicit completeness clause is

\[
 \mathcal O=\mathcal O[q],\qquad C=C[q,q_0].
 \tag{S3-004}
\]

The strongest justified definition is therefore:

> A complete physical configuration is the gauge-equivalence class of an
> admissible representative of the one medium that separates all physical
> distinctions belonging to one instantaneous order state; every instantaneous
> physical object is determined by that class, directly or relative to an
> authorized comparison state.

“Instantaneous” here means “belonging to one order state,” not simultaneous in
a fundamental or coordinate time. The definition makes no claim that the
components of \(\widehat q\) are known.

### Three non-equivalent completeness notions

| Notion | Mathematical meaning | Frozen status |
|---|---|---|
| ontic/state completeness | no additional instantaneous physical distinction lies outside \(q\) | accepted by FP-4 and S2-003 |
| realization completeness | an explicit list of fields/components gives an injective representation of every \(q\) | not established; STATE-002 expressly leaves the realization undiscovered |
| dynamical completeness | one \(q\), plus prescribed boundary/source data, determines a unique future history | expressly not established; S2-014 denies the implication |

Calling \(q\) a “configuration” therefore must not import the narrower
classical meaning “placement coordinate.” Nor may “complete” be expanded into a
Markov or determinism axiom.

## 2. Information inventory

The inventory must distinguish abstract containment from an available formula.

| Candidate information | Status in \(q\) | Frozen justification |
|---|---|---|
| identity of the one medium at the order state | necessarily represented | FP-1, FP-4; S2-002 |
| admissibility, orientation/signature branch | necessarily represented/testable | S2-004, S2-007--S2-009 |
| geometric placement | permitted representative, not the universal definition of \(q\) | STATE-002 calls \(q\) a gauge class of sections, not an embedding; GOVERNING-EQUATION-001 uses placement only for its authorized local realization |
| relative deformation \(C\) | necessarily derived from \(q,q_0\), but does not exhaust \(q\) | S2-003, S2-008, S2-018; the map need not be injective |
| strain coordinates | derived reparametrizations of \(C\), not extra information | S2-008--S2-010 |
| stored energy | necessarily derived once the frozen constitutive choice is applied | HYPER-001 and S2-018: \(W=W(C[q,q_0])\) |
| constitutive stress/elastic response | necessarily derived after the constitutive derivative and realization are specified | \(P_C=DW(C)\); placement stress requires the pullback through \(D_qC\) |
| total stress-energy | not yet defined | BALANCE-001 and INERTIA-001 require kinetic/source and effective mapping structure |
| wave content | not presently decidable as a mode decomposition; any instantaneous physical wave distinction must nevertheless be represented by \(q\) | no closed wave operator or normal-mode decomposition exists; S2-003 governs any later physical observable |
| propagating-mode amplitudes | not presently defined; conditional functionals of \(q\) if they are instantaneous observables | LOCALITY-001 supplies only conditional waves after kinetic closure |
| spatial phase/profile | represented insofar it is a physical distinction in the global field configuration | STATE-002 makes \(q\) a global infinite-dimensional field configuration, not merely three numbers |
| temporal phase or propagation orientation | no explicit extractor is frozen; if physically distinct at one state it must distinguish \(q\)'s | S2-003 plus absence of a kinetic realization |
| momentum or momentum-equivalent quantity | neither established nor excluded as information encoded by \(q\); no momentum functional is derived | INERTIA-001 says a kinetic momentum principle is missing; S2-003 requires any later instantaneous physical momentum to be functional of the complete state |
| arbitrary past history | not necessarily encoded | STATE-002 separately defines histories as curves in \(\mathcal Q_{\rm adm}\); no injective history-to-state encoding is stated |
| duration, coordinate time, or order-label rate | not contained as primitive state information | FP-3; S2-012--S2-013; DURATION-001 |
| gauge representative data | deliberately not physical information | S2-001--S2-002, S2-006 |

Two cautions follow. First, “not yet defined” is not “physically outside \(q\).”
Second, a quantity can depend on a history segment without being additional
instantaneous state information. The frozen texts do not decide whether a
future momentum description is an instantaneous functional of a richer
realization of \(q\), a relational functional of a history germ, or unavailable
until kinetic closure. They only prohibit declaring a new independent state
sector here.

## 3. Classical comparison

Let a standard material continuum use placement \(y(X)\), deformation
\(F=\operatorname{Grad}y\), strain \(C=F^\sharp F\), velocity \(v\), momentum
\(p\), and internal variables \(\alpha\). The frozen correspondence is:

| Classical object | PBUF status |
|---|---|
| placement \(y\) | one permitted local representative of \(q\), not the intrinsic definition |
| deformation \(F,C\) | \(C[q,q_0]\) is derived; \(F\) is representative-dependent |
| strain | derived from \(C\) |
| elastic stress | derived from \(W(C)\) and the kinematic pullback |
| velocity | no primitive analogue; calibrated tangent only after a history and duration are available |
| momentum | no frozen map; cannot be presumed independent of \(q\) and cannot be presumed already extractable |
| kinetic energy | not frozen |
| thermodynamic/internal variables | none authorized |
| phase-space point \((y,p)\) | no established PBUF correspondence; \(q\neq y\) does not imply \(q=(y,p)\) |

Thus classical placement space is demonstrably narrower than the declared PBUF
state concept. Classical phase space is not demonstrably identical to it. The
only safe relation is

\[
q\longmapsto C[q,q_0]\longmapsto(W,P_C),
\tag{S3-005}
\]

with possible placement representatives in the selected local realization.
No reverse reconstruction from these classical quantities to \(q\) is frozen.

## 4. State-completeness theorem

### Theorem — Ontic completeness without dynamical sufficiency

Under FP-4 and STATE-002 S2-002--S2-003:

1. no two distinct instantaneous physical states are represented by the same
   \(q\);
2. no additional independent instantaneous physical state datum may be needed
   to specify which physical state is occupied; but
3. neither uniqueness of future continuation nor an explicit complete
   coordinate realization follows.

**Proof.** By S2-002, points of \(\mathcal Q_{\rm phys}\) are equivalence classes
only under declared descriptive redundancies. By S2-003, completeness assigns
every later physical object to \(q\), directly or relationally with a declared
comparison state. If two physically distinct instantaneous states had the same
\(q\), some physical distinction between them would not be a functional of
\(q\), contradicting the accepted completeness definition. This proves (1) and
(2) as consequences of the frozen premise. STATE-002 simultaneously says that
the component realization has not been discovered, that \(C\) need not be
injective, and in S2-014 that one state does not imply a unique later state.
Therefore neither realization completeness nor deterministic continuation
follows, proving (3). \(\square\)

The theorem establishes Outcome A exactly as defined by the frozen ontology.
It does not derive what the hidden components of a representative are; “hidden
components” would itself be misleading because no additional field or sector
is authorized. It establishes only the abstract separation property of the
already accepted state variable.

## 5. Opposite-wave counterexample audit

EVOLUTION-PRINCIPLE-001 stated conditionally that a reversible wave may pass
through the same configuration with opposite calibrated tangents. That is true
for a classical configuration interpreted as placement or strain alone. It is
not a theorem about the frozen PBUF \(q\).

Assume two universes at one order state have identical \(q\) but physically
different propagation orientation. If propagation orientation is an
instantaneous physical distinction, this assumption contradicts the theorem in
Section 4. At least one of the following must instead hold in any future
realization:

1. the two universes have different \(q\)'s although their placement,
   deformation, and stored energy agree;
2. propagation orientation is relational history information rather than a
   distinction intrinsic to that single order state; or
3. the proposed wave distinction is not physically meaningful before a kinetic
   law defines it.

The frozen milestones do not select among these cases. Therefore the specific
same-\(q\), opposite-wave construction is **not an authorized counterexample**.
It imported the classical identification \(q=y\) or \(q=C\), which STATE-002
expressly does not make.

This correction does not establish deterministic evolution. S2-014 directly
states that \(q(s_1)\) does not imply a unique \(q(s_2)\). Multiple future
histories from one complete state are logically compatible with a
nondeterministic or non-Markovian history law; completeness alone is not a
uniqueness axiom.

## 6. Gauge audit

The frozen physical state space is the quotient by

\[
\mathcal G=\operatorname{Diff}(\mathcal M)\ltimes\mathcal G_{\rm int}.
\tag{S3-006}
\]

Within the accepted definition, \(\mathcal G\) contains the descriptive
redundancies already authorized, so gauge-related representatives describe one
physical configuration. Conversely, distinct points of
\(\mathcal Q_{\rm phys}\) are intended to represent distinct physical states.
This is the quotient's defining role.

What is not established is that an explicit future representative or gauge
group realization has been completely constructed. STATE-002 gives the
abstract quotient and objective tests, not a theorem that a concrete gauge
fixing is globally unique, free of stabilizers, or separates every orbit.
Accordingly:

* there are no authorized physical distinctions *within* one abstract gauge
  class;
* physical distinctions not visible in \(C\) may remain between different
  \(q\)'s, since \(q\mapsto C\) need not be injective; and
* failure of a concrete coordinate model to distinguish states would be a
  failure of that realization, not evidence for hidden physical variables.

## 7. History audit

STATE-002 places histories in a separate mathematical space:

\[
 [\gamma]\in C^0_{\rm ord}(S,\mathcal Q_{\rm adm})/
 \operatorname{Homeo}_+(S).
\tag{S3-007}
\]

A point \(q\) and a curve \([\gamma]\) are therefore different mathematical
types. No frozen injective map from full past histories into a current \(q\),
and no reconstruction of the past from \(q\), is supplied. Arbitrary history is
not necessarily encoded in one configuration.

This does not make history an additional instantaneous state variable. A
history is an ordered relation among complete states. The current \(q\) may
encode physical traces of the past—deformation patterns, wave profiles, or
other presently meaningful features—to the extent they are features of the
current complete state. But two distinct past curves may reach the same \(q\)
unless a future injectivity theorem says otherwise. The frozen rate-free,
history-free stored energy positively supplies no memory functional.

Hence the permitted conclusion is:

\[
\boxed{\text{present physical information is in }q;
\quad\text{the entire history is not proven to be.}}
\tag{S3-008}
\]

No memory variable or hidden degree of freedom is introduced by this
distinction.

## 8. Consequences for evolution theory

EVOLUTION-PRINCIPLE-001 requires a targeted reformulation, not reversal.

### Conclusions that remain valid

1. The continuous order has no intrinsic immediate neighbor, so a discrete
   \(q_n\to q_{n+1}\) remains a sampling unless transition structure is added.
2. Admissibility, \(W\), balance templates, and boundary data do not select a
   unique history.
3. Duration calibration does not itself construct a momentum covector or
   kinetic operator.
4. Stable reversible waves and V11 causal compatibility still require a closed
   dynamical structure.
5. S2-014 independently establishes that unique continuation is not frozen.

### Conclusion that must be narrowed

The claim that one PBUF \(q\) admits two physically distinct opposite-wave
states is unsupported. It must be restated conditionally:

> If a chosen realization identifies \(q\) only with placement/deformation,
> opposite wave directions demonstrate that realization's incompleteness. They
> do not demonstrate incompleteness of the abstract frozen state.

### Revised location of the missing information

Because \(q\) is ontically complete, this milestone may not identify momentum,
velocity, a prior configuration, or a hidden sector as missing **state
information**. What is missing is mathematical structure, not an additional
state datum:

* an explicit realization showing how all instantaneous distinctions are
  encoded in \(q\); and
* an evolution/kinetic rule acting on the already complete state space.

Whether a future realized \(q\) makes evolution first-order, whether evolution
is relational in a history germ, and whether momentum is an instantaneous
derived functional remain open. No kinetic law is derived here.

## 9. Deliverable summary

| Deliverable | Result |
|---|---|
| formal definition | \(q\) is the gauge class separating all instantaneous physical distinctions of the one medium |
| information inventory | \(C,W,P_C\) are derived; placement is only a realization; wave/momentum extractors are unconstructed but any instantaneous physical distinctions cannot lie outside \(q\) |
| classical comparison | PBUF \(q\) is broader than classical placement but is not proven identical to classical phase space |
| completeness theorem | Outcome A ontically; realization and dynamical completeness remain unproved |
| counterexample analysis | same-\(q\), opposite-wave universes do not follow and would contradict ontic completeness if orientation is instantaneous |
| gauge audit | abstract gauge classes identify physical identity; \(C\) alone does not separate all \(q\)'s |
| history audit | present physical traces belong to \(q\), but an entire past curve is not proven encoded |
| evolution consequence | successor non-closure remains; its classical opposite-wave rationale must be conditionalized |

## 10. Traceability and status

| Result | Frozen source |
|---|---|
| entire physical state, functional completeness | FOUNDATION-001 FP-4; STATE-002 S2-002--S2-003 |
| no known component realization | STATE-002 canonical definition and completion boundary |
| deformation derived but non-injective | STATE-002 S2-008--S2-011, S2-018 |
| global minimum degrees of freedom only | STATE-002 S2-016--S2-017 |
| histories separate from states; no deterministic implication | STATE-002 S2-012--S2-015 |
| duration/rates not primitive | DURATION-001 |
| kinetic and momentum maps open | GOVERNING-EQUATION-001; INERTIA-001 |
| conditional wave structure | LOCALITY-001; INERTIA-001 |
| successor underdetermination | EVOLUTION-PRINCIPLE-001 Theorems 1--2, excluding its conditional classical wave example as a PBUF proof |

**Status: complete.** The frozen \(q\) is the entire instantaneous physical
state by accepted definition, not merely a geometric placement. Its explicit
realization and its sufficiency for deterministic continuation are not frozen.
The EVOLUTION-PRINCIPLE-001 non-closure conclusion remains, but no missing
instantaneous momentum-like datum may be inferred from the same-configuration,
opposite-wave analogy.
