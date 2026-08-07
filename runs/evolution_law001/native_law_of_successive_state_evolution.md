# PBUF EVOLUTION-LAW-001 — Native Law of Successive State Evolution

## 0. Result

With every listed milestone held fixed, the strongest native evolution object is
not an immediate-successor map. It is a selection of oriented, unparameterized
admissible histories of the already complete state:

\[
 \mathcal H_{\rm adm}:=
 C^0_{\rm ord}(S,\mathcal Q_{\rm adm})/\operatorname{Homeo}_+(S),
 \qquad
 \boxed{\mathcal H_{\rm phys}\subseteq\mathcal H_{\rm adm}}.                 \tag{EL-001}
\]

Here \(\mathcal H_{\rm phys}\) is the class selected by the still-absent
evolution law. It is mathematical structure on the state space, not a new state
variable, field, constituent, momentum sector, or microscopic ontology.

For \(q\in\mathcal Q_{\rm adm}\), let \(\gamma\sim_0q\) mean that the oriented
history \(\gamma\) contains a chosen event whose state is \(q\), and let
\(\gamma_+^0\) be its future restriction from that event. The most general
state-indexed continuation correspondence induced by EL-001 is

\[
 \boxed{
 \mathcal E(q):={[\gamma_+^0]:[\gamma]\in\mathcal H_{\rm phys},\
                                  \gamma\sim_0q\}
 \subseteq \mathcal H^+_{\rm adm}(q).}                                      \tag{EL-002}
\]

Its future-state reachability projection is

\[
 \mathcal R(q):=\{q'\in\mathcal Q_{\rm adm}:q'\text{ occurs after }q
                   \text{ on some }[\gamma]\in\mathcal H_{\rm phys}\}.      \tag{EL-003}
\]

Thus the requested notation \(\mathcal E:q\leadsto\) future admissible states
must be read as a set-valued history/continuation correspondence. EL-003 alone
is generally too lossy to be an evolution law: it forgets ordering,
intermediate states, multiplicity, and correlations along a continuation.

No frozen result selects \(\mathcal H_{\rm phys}\), makes EL-002 single-valued,
assigns probabilities, or makes it depend only on its endpoint state. This is
the exact underdetermination result of this milestone.

## 1. Structural theorem

### Theorem 1 — Maximal frozen-compatible evolution framework

Every evolution law compatible with the frozen milestones can be represented
by a subclass \(\mathcal H_{\rm phys}\subseteq\mathcal H_{\rm adm}\), or by an
equivalent continuation rule whose concatenated solutions define such a
subclass, subject to all of the following and no stronger unconditional
requirements:

1. **Admissibility:** every state on a selected history belongs to
   \(\mathcal Q_{\rm adm}\), including the declared endpoint/tangent-cone
   restrictions.
2. **Gauge basicness:** selection is defined on
   \(\mathcal Q_{\rm phys}=\mathcal A/\mathcal G\); changing representatives
   cannot change whether a physical history is selected.
3. **Positive-order reparameterization invariance:** if a representative is
   selected, every increasing relabeling represents the same selected physical
   history. No step size or immediate successor is physical.
4. **Constitutive compatibility:** along every sufficiently regular selected
   history, \(C=C[q,q_0]\), \(W=W(C)\), \(P_C=DW\), and
   \(dW=P_C:dC\). Selection may not alter the frozen constitutive response.
5. **Locality compatibility:** the elastic part of any localized governing
   representative is the frozen local variational pullback
   \((D_qC)^*P_C\) (or \(-\operatorname{Div}_0P_F\) in the authorized placement
   realization). No gradient or nonlocal constitutive enrichment follows.
6. **Balance compatibility:** whenever a balance channel has separately been
   closed, every selected history satisfies its balance equation, sources, and
   boundary fluxes. The open balance templates do not themselves select a
   history.
7. **Energy consistency:** the stored-energy chain rule is mandatory. Total
   energy conservation is mandatory only for a separately closed conservative
   channel with the required full-law symmetry and closed/no-work boundary
   conditions.
8. **Duration compatibility:** any duration-parametrized representative is a
   gauge fixing of an already selected oriented history using the DURATION-001
   positive additive functional. The selection cannot depend on the arbitrary
   order label.
9. **V11/metric gate:** if the selected law has propagating characteristics in
   the V11 regime, their causal structure must pass the single-effective-metric
   compatibility gate. The metric map does not create the evolution law.

**Proof.** Items 1--3 are the definitions of the frozen state and history
spaces. Item 4 follows by composition with the frozen statewise hyperelastic
functional and its chain rule. Item 5 is the frozen LOCALITY-001 and
GOVERNING-EQUATION-001 elastic operator, which supplies restoring response but
not selection. Item 6 is exactly BALANCE-001's distinction between a balance
template and a closed channel. Item 7 separates the unconditional chain rule
from BALANCE-001's conditional Noether conservation. Item 8 follows from
FP-3 and DURATION-001. Item 9 is FP-5 and METRIC-001 applied conditionally to a
law possessing characteristics. None of the cited results adds uniqueness,
probability, memorylessness, reversal closure, or a kinetic principal symbol.
Hence those properties cannot be appended to the unconditional list. \(\square\)

The theorem states compatibility conditions, not a proposed physical law.

### Theorem 2 — No immediate-successor operator

On an interval-ordered continuous history, there is no invariant “next” state.

**Proof.** Between any two distinct labels of an interval lies another. An
increasing relabeling and insertion of further samples leave the physical
history unchanged. Therefore a map \(q_n\mapsto q_{n+1}\) requires a preferred
sampling or step structure absent from the quotient in EL-001. \(\square\)

## 2. Determinism audit

Use the following precise definitions:

* **unique continuation:** for each admissible conditioning datum, exactly one
  future history germ is selected, modulo gauge and increasing relabeling;
* **multiple continuation:** more than one such germ is selected;
* **probabilistic continuation:** a normalized probability measure or kernel is
  supplied on the continuation space;
* **history-dependent continuation:** continuations depend on an admitted past
  history segment/germ, not merely its endpoint;
* **relational continuation:** selection is stated through relations among two
  or more complete states or history segments rather than a unary flow.

### Theorem 3 — Determinism is undecidable from the frozen milestones

The frozen milestones require none of unique, multiple, probabilistic,
history-dependent, or relational continuation. They permit deterministic and
nondeterministic selections, but supply neither.

**Proof.** STATE-003 establishes separation of instantaneous physical states by
\(q\), not a selection on paths. The frozen restrictions define
\(\mathcal H_{\rm adm}\) and filters on its members. Provided at least two
compatible histories exist, choosing one, choosing both, or conditioning their
selection on a past segment leaves the same state space, constitutive law, and
compatibility filters. A probability law additionally requires a sigma-algebra
and normalized measure/kernel; none is frozen. Relational continuation is a
permitted representation of a selection but is not forced by the existence of
relations such as \(C[q,q_0]\). Therefore no listed continuation type is a
consequence. \(\square\)

“Nondeterministic” must not be silently replaced by “probabilistic.” A bare
set-valued correspondence supplies possibilities, not chances.

## 3. Markov audit

Ontic completeness and the Markov property have different mathematical types.
Ontic completeness says that \(q\) separates all physical distinctions at one
order state. A Markov evolution law additionally says that the conditional law
of the future factors through the current state.

For a set-valued history law this would require a continuation map \(K\) such
that, for every admitted past \(h_-\) ending at \(q\),

\[
 \mathcal E_+(h_-)=K(q),                                                     \tag{EL-004}
\]

together with a compatible composition/concatenation property. In a
probabilistic formulation it would require a measurable state space and a
transition kernel \(K(q,\cdot)\) whose conditional future distribution is
independent of the past given \(q\); a duration-homogeneous semigroup would be
an additional, still stronger structure.

### Theorem 4 — Ontic completeness does not imply Markov evolution

**Proof.** Completeness constrains functions of an instantaneous physical
state. It says nothing about whether a rule on curves factors through the
endpoint evaluation map \(h_-\mapsto q\). A law can condition admissible future
curves on relational properties of the past without adding an instantaneous
state variable. Conversely, a state-only kernel can be Markov. Both act on the
same complete state space. The precise missing structure is therefore a
history-selection/conditional-continuation rule plus, for Markovity, the
factorization and composition property EL-004. For stochastic Markovity a
sigma-algebra and transition probability kernel are additionally required.
\(\square\)

The state-indexed EL-002 is the union of continuations over all pasts reaching
\(q\). For a history-dependent law this union is only an envelope; it does not
retain which future is admissible after which past.

## 4. Reversibility

Constitutive reversibility and evolution reversibility are independent.

The frozen constitutive branch is state-local, rate-free, single-valued, and
hyperelastic. Therefore for a smooth path contained in that branch,

\[
 \int_{\gamma_{AB}}P_C:dC=W(C_B)-W(C_A),\qquad
 \int_{\gamma_{BA}}P_C:dC=W(C_A)-W(C_B).                                   \tag{EL-005}
\]

This is exact recovery of stored work. It excludes constitutive hysteresis,
plastic memory, and dissipation within the frozen branch.

Evolution reversibility would instead require a reversal operation \(\Theta\)
on physical histories, including any transformation of instantaneous features
already encoded by \(q\), such that

\[
 [\gamma]\in\mathcal H_{\rm phys}
 \Longrightarrow [\Theta\gamma]\in\mathcal H_{\rm phys}.                  \tag{EL-006}
\]

### Theorem 5 — Evolution reversibility is not fixed

The frozen constitutive reversibility EL-005 does not imply EL-006. Reversible,
irreversible, and branchwise/partially reversible history selections can all
use the same \(W\). No frozen reversal involution or reversal-closed selection
is supplied. Strict stable lossless wave propagation would conditionally
require a reversible conservative propagating sector, but it does not prove
global reversal invariance. \(\square\)

## 5. Conservation theorems

### Theorem 6 — Admissibility preservation

Every admissible PBUF evolution law preserves admissibility in the definitional
sense that its selected histories take values in \(\mathcal Q_{\rm adm}\).
This does not prove that an arbitrary candidate differential vector field is
tangent to the admissible domain; tangency or a viable endpoint rule is a test
that the candidate must pass.

### Theorem 7 — Gauge-class well-definedness

Every admissible evolution law preserves physical equivalence: representative
changes cannot change the selected quotient history. This is gauge basicness,
not conservation of a numerical “gauge charge.” It does not require a fixed
representative along a history.

### Theorem 8 — Balance compatibility is conditional

Every selected history must satisfy each separately closed applicable balance
law. No nontrivial conserved scalar, momentum, propagation charge, or total
energy follows from an unclosed balance template. Integrated conservation
additionally needs zero production and zero net boundary flux.

### Theorem 9 — Stored and total elastic energy

Every smooth selected history satisfies the identity
\(dW=P_C:dC\). It need not preserve \(W\). In a reversible wave, stored energy
normally exchanges with the still-unclosed kinetic channel. Conservation of
total elastic-plus-kinetic energy is required only after such a total energy,
its exchange/flux law, suitable boundary conditions, and the relevant symmetry
have been supplied. Consequently imposing \(W(q')=W(q)\) is neither a frozen
law nor a substitute for total-energy conservation.

These four theorems exhaust the unconditional conservation content relevant to
an evolution selector.

## 6. Continuous-wave compatibility

The frozen local elastic operator supplies restoring communication and an
acoustic stiffness. It does not supply a frequency, characteristic cone,
domain of dependence, or continuation selector.

### Theorem 10 — Conditional wave restriction

Continuous stable lossless wave propagation narrows \(\mathcal H_{\rm phys}\)
only conditionally: a wave-supporting member must possess a well-posed
propagating history law, a positive nondegenerate kinetic response on physical
wave modes, conservative exchange with \(W\), and characteristics compatible
with the V11 effective cone. These are tests on the absent law. The frozen
results neither construct such a member nor require all admissible evolution
laws to be wave supporting.

**Proof.** LOCALITY-001 fixes only the restoring spatial symbol. INERTIA-001
proves that stable reversible waves require a positive momentum-carrying
kinetic structure but that \(W\) does not determine it. METRIC-001 supplies a
causal compatibility gate, not a principal symbol. Hence wave propagation
constrains a selected law's outputs and symbol but cannot select the law from
the frozen data. \(\square\)

No independent momentum, kinetic energy, or kinetic operator is introduced by
this statement; the missing role is only identified.

## 7. Emergent-duration parameterization

Let \([\gamma]\in\mathcal H_{\rm phys}\) first be selected without a
fundamental time coordinate. DURATION-001 assigns, on a propagation-bearing
segment, the invariant accumulation

\[
 \tau(s)=\tau_0+\int_{s_0}^{s}F(q,\dot q,\dot z)\,ds,
 \qquad F(q,a\dot q,a\dot z)=aF(q,\dot q,\dot z),\quad a>0.                  \tag{EL-007}
\]

Where \(F>0\), \(\tau\) is monotone and may parameterize the same curve:

\[
 q_\tau:=\gamma\circ\tau^{-1}.                                              \tag{EL-008}
\]

EL-008 is a calibrated representative of an already selected physical
history. It does not turn \(\tau\) into a fundamental coordinate and does not
allow the evolution selector to depend on an arbitrary \(s\). In a constant
complete history, \(F=0\), so no physical-duration flow parameter exists; the
history is simply constant. Since the frozen theory selects neither \(F\) nor
its normalization uniquely, it also does not select a unique duration-rate
form of the evolution law.

## 8. Classical comparison

Classical Hamiltonian notation uses a phase-space state and a flow,

\[
 (q_{\rm cl},p_{\rm cl})\longmapsto(q'_{\rm cl},p'_{\rm cl}).               \tag{EL-009}
\]

PBUF instead declares \(q\) to be the complete instantaneous physical state
and writes the abstract history selection

\[
 q\leadsto\mathcal E(q).                                                     \tag{EL-010}
\]

The notation does not identify PBUF \(q\) with classical placement. If a future
realization derives gauge-invariant functionals

\[
 Q_{\rm cl}=Q[q],\qquad P_{\rm cl}=P[q]                                     \tag{EL-011}
\]

on an invariant sector and the selected PBUF histories push forward to a
closed deterministic flow of \((Q_{\rm cl},P_{\rm cl})\), then EL-010 can
reproduce EL-009 on that sector. The derived pair need not reconstruct all of
\(q\); exact equivalence would additionally require the pushforward to be
injective on the relevant physical states and to intertwine the two evolution
laws.

No such functionals, closure, injectivity, symplectic form, Legendre map, or
Hamiltonian flow are frozen. Therefore classical dynamics is a possible
derived representation, not a consequence and not an input. In particular,
EL-011 introduces no independent momentum state: any eventual momentum-like
quantity must be a functional of the already complete \(q\), or a relational
functional of a selected history as allowed by STATE-003.

## 9. Minimal missing mathematical structure

### Theorem 11 — Exact closure boundary

The minimal structure still absent is:

> a gauge-basic, positive-order-reparameterization-invariant selection rule on
> admissible oriented histories, equivalently a conditional continuation
> correspondence on admissible past germs, with enough composition structure
> to determine which continuations form one physical history.

This one role must decide which elements of \(\mathcal H_{\rm adm}\) are
physical. If unique evolution is desired, it must additionally make the
conditional continuation single-valued and well posed. If Markov evolution is
desired, it must factor through the current-state endpoint and obey the
appropriate composition law. If probabilistic evolution is desired, it must
also provide measurable continuation spaces and normalized kernels. If stable
reversible waves are required, it must further possess the conditional
kinetic, conservative, and causal properties in Theorem 10.

None of these optional strengthenings is part of the minimal role, and no
formula for the rule is proposed here.

**Proof.** Statewise structures determine the feasible points and elastic
response but not a subset of curves. Any actual evolution law must distinguish
physical curves from merely admissible ones, so a selection rule is necessary.
Once such a subset is supplied, EL-002 reconstructs every permitted future and
all frozen filters can be tested, so it is sufficient at the abstract operator
level. STATE-003 forbids interpreting this missing selector as additional
instantaneous information. \(\square\)

## 10. Deliverable summary

| Deliverable | Strongest justified result |
|---|---|
| general framework | selected unparameterized history class EL-001; set-valued continuation EL-002; reachability EL-003 |
| structural theorem | Theorem 1 plus the dense-order no-successor Theorem 2 |
| determinism | none of unique, multiple, probabilistic, hereditary, or relational continuation is forced |
| Markov audit | ontic completeness does not imply endpoint factorization or composition |
| reversibility | constitutive work recovery is frozen; evolution reversal is undecided |
| conservation | admissibility and quotient well-definedness are mandatory; balance and total-energy conservation are conditional |
| waves | only conditional kinetic, stability, well-posedness, and V11-cone gates |
| duration | duration may gauge-parameterize a selected history; it is not a fundamental coordinate |
| classical comparison | classical phase flow may be a derived pushforward only after unprovided closure and intertwining |
| missing structure | a gauge-basic reparameterization-invariant history-selection/conditional-continuation rule |

## 11. Scope and status

No ontology was reviewed. No state variable, hidden momentum sector,
microscopic constituent, auxiliary field, kinetic operator, kinetic energy,
fundamental coordinate time, constant, observational fit, V11 change, or
weak-lensing modification has been introduced.

**Status: complete negative-and-structural derivation.** The maximal compatible
operator class and every unconditional property have been identified; unique,
Markov, probabilistic, reversible, conservative, and wave-supporting evolution
remain additional mathematical properties of the absent history-selection
rule, not consequences of instantaneous state completeness.
