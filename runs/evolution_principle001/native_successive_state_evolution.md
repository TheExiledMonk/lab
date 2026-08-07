# PBUF EVOLUTION-PRINCIPLE-001 — Native Successive-State Evolution

## 0. Verdict

The frozen ontology defines the **space and order type of admissible histories**,
but it does not define a native immediate-successor map or relation:

\[
 \boxed{
 (\mathcal Q_{\rm adm},W,\text{balance templates},\text{boundary data})
 \not\Rightarrow q_n\longmapsto q_{n+1}.}
 \tag{SE-001}
\]

This is stronger than the earlier statement that elasticity does not determine
an inertial operator. Before duration calibration, STATE-002 represents a
history by a continuous order-preserving curve and quotients it by every
increasing relabeling,

\[
 [\gamma]\in\mathcal H_{\rm adm}
 =C^0_{\rm ord}(S,\mathcal Q_{\rm adm})/\operatorname{Homeo}_+(S).
 \tag{SE-002}
\]

When (S) is an interval, its order is dense: between two distinct ordered
labels lies another. It therefore has no immediate successor. A sequence
\((q_n)\) is only a sampling of an unparameterized history, and which sampled
state is called “next” changes when samples are inserted or removed. The frozen
relabeling quotient contains no preferred sampling or step size.

Even if a discrete sampling is imposed representationally, admissibility,
stored energy, balance templates, and boundary conditions filter possible
states but do not select one continuation. Momentum and inertia consequently
do not emerge from those data. The exact missing mathematical role is a
**gauge-basic kinetic/history-selection structure**: structure that compares
physical change between configurations, selects admissible unparameterized
histories (or their local directions), and becomes a nondegenerate
tangent-to-cotangent/Legendre map after duration calibration. This identifies
the role only; no new principle is formulated here.

## 1. Evolution without fundamental time

The native time-free object already authorized is not a discrete recurrence but
the oriented geometric history \([\gamma]\) in SE-002. It contains:

1. the complete configurations traversed;
2. their orientation/order; and
3. continuity in \(\mathcal Q_{\rm adm}\).

It does not contain an origin, numerical interval, rate, duration, or immediate
adjacency. Thus a continuous medium with emergent time **admits successive
configuration histories**, but “successive” means order along an
unparameterized curve, not a native map between neighboring elements.

For any \(q_-\prec q_+\) on a nonconstant continuous interval history, every
intermediate label may be inserted without changing the physical history
class. Consequently a purported local rule

\[
 q_{n+1}=T(q_n)
 \tag{SE-003}
\]

is not invariant under the frozen history equivalence unless it merely
describes a chosen section/sampling or carries extra adjacency structure. A
parameter-free law may select whole unparameterized curves, but the frozen
inputs do not supply such a law either.

## 2. Admissible-successor impossibility theorem

### Theorem 1 — No native immediate successor

Let the frozen state and history spaces be those of STATE-002. Then no
nontrivial, sampling-independent immediate-successor operation on a continuous
history follows from the ordered-history structure.

**Proof.** An interval is densely ordered. If \(s_0<s_1\), there exists
\(s_*\) with \(s_0<s_*<s_1\). Continuity supplies the intervening state
\(q(s_*)\), possibly equal to an endpoint on a constant segment but not
removable as an ordered label by a physical step-size criterion. Moreover,
\(\operatorname{Homeo}_+(S)\) preserves order while changing every numerical
spacing. Hence no pair of distinct labels is intrinsically adjacent, and no
immediate-successor operation descends from label space to the quotient
history. \(\square\)

### Theorem 2 — Frozen statewise data do not select continuation

Let \(\mathcal A_B\subseteq\mathcal Q_{\rm adm}\) denote configurations
satisfying the prescribed boundary conditions. Let \(W:\mathcal A_B\to
\mathbb R\cup\{+\infty\}\) be the frozen rate-free stored energy. Admissibility,
\(W\), BALANCE-001, and boundary conditions do not determine a unique
continuation from a general \(q\in\mathcal A_B\).

**Proof.**

* Admissibility and boundary conditions define a feasible set, not a direction
  or map on it.
* \(W\) assigns a scalar to each feasible configuration and its variation gives
  the elastic restoring covector. A scalar state function does not determine a
  tangent, orientation, transition cost, or map.
* BALANCE-001 provides identities and conditional balance templates. It
  explicitly leaves densities, fluxes, sources, and any momentum current open.
* At a stable unloaded state, minimizing \(W\) returns equilibrium; away from
  it, following \(-DW\) would require an unprovided mobility/metric and would
  select relaxation rather than conservative propagation. Requiring equal
  stored energy instead leaves, generically, an entire level set and therefore
  does not select a successor.

These structures are shared by mutually inequivalent continuations: constant
histories, quasistatic equilibrium families when boundary data vary, arbitrary
continuous feasible curves, and—after adding different kinetic closures—many
different propagating histories. Since identical frozen inputs admit different
continuations, no unique successor follows. \(\square\)

The theorem also excludes a nonunique “allowed-successor relation” as a
dynamical closure. The maximal relation inferred from statewise data is merely
feasibility, possibly supplemented by specified endpoint or boundary work. It
contains no criterion for immediacy and no physical selection among its pairs.

## 3. Conservation audit

| Quantity | What successive admissible states preserve | What does not follow |
|---|---|---|
| admissibility | Preserved if every member is required to lie in \(\mathcal Q_{\rm adm}\); continuity also preserves path-component membership | No theorem says an unconstrained successor map is tangent to the admissible domain or respects a hard endpoint |
| gauge equivalence | Preserved only if a proposed relation is well-defined on \(\mathcal Q_{\rm phys}=\mathcal A/\mathcal G\) | A representative-level update need not descend to the quotient |
| balance laws | Each selected history must satisfy any separately closed balance law | Balance templates alone supply neither a history nor the density, flux, and source needed to test it |
| stored energy | The chain rule \(dW=P_C:dC\) holds along a smooth sampled history | State succession alone does not imply \(W(q_{n+1})=W(q_n)\) |
| total energy | Nothing unconditional | Kinetic storage, work/flux pairing, a full evolution law, and calibrated-duration translation symmetry are missing |

Demanding preservation of \(W\) cannot replace total-energy conservation. In a
reversible elastic wave, stored energy generally changes and is exchanged with
kinetic energy. A configuration-only rule that fixes \(W\) forbids that exchange;
a rule that allows \(W\) to change has no frozen complementary energy account.

Gauge preservation is a **well-definedness condition**, not a conservation law:
if \([q]=[q']\), a physical relation must give the same successor classes from
both representatives. The frozen quotient requires this of any future rule but
does not construct the rule.

## 4. Can momentum emerge after duration calibration?

Suppose an unparameterized physical history has been selected by some means and
DURATION-001 then supplies a monotone physical duration \(\tau\). The history
acquires a calibrated tangent

\[
 v_\tau={Dq\over D\tau}\in T_q\mathcal Q_{\rm phys}.
 \tag{SE-004}
\]

SE-004 is a rate, not momentum. Momentum is a covector and requires a map

\[
 T_q\mathcal Q_{\rm phys}\longrightarrow T_q^*\mathcal Q_{\rm phys}
 \tag{SE-005}
\]

or the corresponding variational/history operator. Neither \(W(C)\),
\(D_qC\), duration calibration, nor a pair \((q_n,q_{n+1})\) supplies this
typed identification. In particular:

* a finite difference of configurations is generally coordinate-dependent and
  is not a cotangent vector;
* its normalization depends on calibrated duration and its comparison across
  different tangent spaces requires geometric structure;
* multiplying or otherwise remapping any candidate momentum by different
  positive gauge-basic operators leaves the same configuration curve while
  changing momentum and inertia.

Therefore momentum does **not** emerge from succession plus duration alone.
The precise missing ingredient is a nondegenerate, gauge-basic kinetic duality
or equivalent Legendre structure on physical configuration change. Duration
calibrates “how fast”; it does not determine “how much covector momentum.”

## 5. Can inertia emerge?

Resistance to changing a proposed successor sequence has no mathematical
meaning until neighboring histories can be compared and a cost, action, or
momentum response is assigned to their change. The frozen elastic Hessian
measures resistance to changing **configuration**; it does not measure resistance
to changing the history tangent or momentum.

Accordingly, successive configurations alone do not reproduce INERTIA-001's
structural theorem. They can impose only the kinematic gates already present:
admissibility, quotient covariance, order invariance, and static consistency.
They do not imply:

* a cotangent-valued inertial response;
* positive nondegenerate kinetic storage on propagating modes;
* a symmetric conservative kinetic principal part;
* a momentum/total-energy balance; or
* causal, well-posed propagation.

If the missing kinetic/history-selection structure is later supplied, its
duration-calibrated Legendre derivative can reproduce all ten structural gates
of INERTIA-001. In that case inertia is **derived from that structure**, not
from bare succession or elasticity. Relabeling the same missing content as
“resistance between states” does not remove the underdetermination.

## 6. Wave-propagation compatibility

A successor formulation is not intrinsically incompatible with waves, but the
frozen configuration-only data do not make it wave supporting.

1. **Reversibility.** A relation must admit the reversed oriented history (with
   any required state involution). Feasibility and hyperelastic path recovery do
   not imply this. A dissipative selection and a reversible selection can use
   the same \(W\).
2. **Finite propagation speed.** This is a statement about domains of dependence
   or a characteristic/causal relation. No principal symbol or transition cone
   is defined by a statewise energy and balance template.
3. **Stable waves.** The frozen acoustic stiffness gives the restoring symbol,
   but real stable frequencies require a positive kinetic symbol. With no such
   symbol, neither frequency nor propagation speed is defined.
4. **V11 compatibility.** METRIC-001 supplies the effective cone as a gate. It
   does not prove that a proposed successor relation has characteristics, nor
   that they lie on or within that cone.

There is also a state-completeness obstruction to a deterministic one-step map
on configuration alone. A reversible wave at the same configuration may have
opposite calibrated tangents. If both physical possibilities exist, \(q\) alone
cannot choose between them. A two-configuration recurrence can encode the
missing direction through \((q_{n-1},q_n)\), but that is second-order history
data plus a transition rule, not the requested map from one complete
configuration. In the continuum limit it is precisely tangent/kinetic
structure in another representation.

## 7. Comparison with classical dynamics

The two diagrams are not equivalent under the frozen inputs.

### A. Configuration \(\rightarrow\) momentum \(\rightarrow\) evolution

In a regular classical kinetic theory, configuration and calibrated tangent
determine momentum through a Legendre map; Hamiltonian or Euler--Lagrange
structure then determines evolution. Strictly, configuration alone is not
complete phase-space data, so the diagram suppresses tangent/momentum initial
data and constraints.

### B. Configuration \(\rightarrow\) successor \(\rightarrow\) derived momentum

A deterministic successor map presupposes enough information to select a
direction and amount of change. Deriving momentum afterward remains impossible
without a cotangent pairing/Legendre map. For reversible second-order systems,
one prior configuration or an equivalent tangent is also needed.

The formulations become **partially equivalent only after** a regular discrete
or continuous variational kinetic structure is supplied. Then a two-point
transition description and a momentum description can be related by Legendre
transforms, subject to constraints and gauge degeneracy. They can fail to be
equivalent when that transform is singular, the relation is multivalued,
dissipative, non-Markovian, or loses history information. None of the
equivalence hypotheses is frozen here.

Thus formulation B does not eliminate the kinetic principle. At best it moves
the same mathematical content into a two-point/history selector.

## 8. Minimal additional principle: exact mathematical role

The minimal remaining principle must perform all and only the following role:

1. select physical oriented unparameterized histories, or equivalently define
   permissible local physical directions, rather than merely feasible states;
2. compare change between configurations in a way that descends to the gauge
   quotient and is invariant under monotone order relabeling;
3. supply, after duration calibration, a nondegenerate physical
   tangent-to-cotangent relation on propagating modes, while retaining genuine
   gauge/constraint null directions;
4. provide the kinetic storage and exchange channel needed for reversible total
   energy accounting; and
5. determine a causal principal/transition structure testable against the
   single V11 effective cone.

This is one closure role, not five independent new fields or constants. It may
have continuous variational, constrained, or two-point discrete
representations, but choosing among those representations is not authorized by
this milestone. Boundary/source data remain problem data and cannot substitute
for it.

## 9. Deliverable summary

| Deliverable | Result |
|---|---|
| successive-state evolution analysis | Native object is an oriented unparameterized continuous history, not an immediate-successor sequence |
| admissible-successor theorem | No immediate successor exists on the dense frozen order; statewise constraints do not select continuation |
| conservation audit | Only imposed admissibility, quotient well-definedness, chain rule, and path-component persistence are automatic; no total-energy conservation follows |
| momentum emergence | Duration gives a tangent rate, not a momentum covector; a kinetic duality/Legendre map is missing |
| inertia emergence | Bare succession supplies none of the positive kinetic response required by INERTIA-001 |
| wave compatibility | Possible only conditionally after kinetic/history closure; reversibility, stability, finite speed, and V11 cone compliance are not derived |
| classical comparison | Partially equivalent only with a regular kinetic variational structure; otherwise fundamentally underdetermined |
| minimal remaining principle | A gauge-basic, reparametrization-invariant kinetic/history-selection role |

## 10. Traceability and milestone status

| Result | Frozen source |
|---|---|
| SE-002, dense-order obstruction, no deterministic continuation | STATE-002 S2-012--S2-015 |
| feasible domain and gauge quotient | STATE-002; DEFORMATION-001; HYPER-001 |
| statewise rate-free energy and restoring covector | HYPER-001; ENERGY-PRINCIPLE-001; CONSTITUTIVE-CONSTRUCTION-001 |
| conservation limitations | BALANCE-001 B-001--B-008 |
| duration supplies calibration but not momentum | DURATION-001 |
| elastic operator complete, kinetic slot open | GOVERNING-EQUATION-001 GE-001--GE-009 |
| locality supplies restoring communication but not kinetic propagation | LOCALITY-001 L-001--L-003 |
| kinetic structural gates and independence result | INERTIA-001 Theorems I--III |
| V11 causal gate | FOUNDATION-001; METRIC-001; INERTIA-001 |

**Status: complete negative derivation.** The frozen ontology naturally admits
ordered continuous histories but does not naturally define immediate successive
states or select an evolution. Momentum and inertia do not emerge from sampled
configuration change after duration calibration. Their common missing
mathematical source is the kinetic/history-selection role identified in Section
8. No ontology, constitutive law, state variable, constant, V11 result,
weak-lensing calculation, or observational fit has been changed.
