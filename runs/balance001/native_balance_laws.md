# PBUF BALANCE-001 — Native Balance Laws of the Spacetime Medium

## 0. Decision and scope

FOUNDATION-001, STATE-002, DEFORMATION-001, HYPER-001,
ENERGY-PRINCIPLE-001, DYNAMICS-001, DURATION-001, and METRIC-001 are fixed.
This milestone derives balance structure, not constitutive or evolution equations.

The central result is deliberately restrictive:

\[
\boxed{\text{The accepted architecture fixes balance templates and identities,
but no nontrivial universal scalar, energy, momentum, or propagation charge.}}
\tag{B-001}
\]

Ordinary conservation follows only after the full chosen action has the relevant
global symmetry, or after a closure explicitly declares a source-free current.
Treating locality, objectivity, covariance, or the existence of one medium as a
conservation theorem would add an assumption.

## 1. Native conservation hierarchy

### Tier 0 — exact kinematic and gauge facts

1. **One-state persistence.** Evolution is a continuous curve
   \(q:S\to\mathcal Q_{\rm adm}\) in the state space of the same one medium.
   This is persistence of the system represented, not conservation of an
   extensive “amount of medium”; no such density or measure-valued content was
   defined.
2. **Gauge/objectivity identities.** Gauge-related representatives are one
   state. Their variations are null directions of the physical first variation.
   These are identities/constraints, not charges.
3. **Order-relabeling identity.** For every selected differentiable native
   degree-one action,
   \[
      p=D_v\mathscr L,\qquad \langle p,v\rangle-\mathscr L\equiv0 .
      \tag{B-002}
   \]
   This is D-006--D-007. It is the parameter-Hamiltonian constraint, not energy
   conservation.
4. **Stored-energy chain rule.** On a smooth elastic branch,
   \[
      dW=P_C:dC,\qquad P_C=DW(C).
      \tag{B-003}
   \]
   Along any clock gauge, \(\dot W=P_C:\dot C\). This is exact reversible
   storage accounting, not a closed energy balance.
5. **Topological continuity.** A continuous history remains in its path
   component of \(\mathcal Q_{\rm adm}\). Hence any already-defined discrete
   invariant constant on path components is unchanged. No particular winding,
   defect, or homology charge has been defined, and therefore none is invented.
   The fixed carrier \(\mathcal M\) is not dynamically replaced. Jumps and
   topology-changing transition rules are outside the minimal history space.

### Tier 1 — universal local balance form, conditional on localization

Let \(a\) denote any scalar, vector, covector, or tensor quantity *already
derived from the selected realization*. A local reference-carrier description
has the bookkeeping form

\[
 \boxed{\partial_\tau \rho_a+\operatorname{Div}_0 J_a=\sigma_a,}
 \tag{B-004}
\]

after the emergent duration \(\tau\) has been calibrated. Here \(\rho_a\) is
the amount per reference volume \(dV_0\), \(J_a\) is outward transport per
reference area per duration, and \(\sigma_a\) is production or external supply
per reference volume per duration. Equivalently, for every material region
\(U\subset\mathcal B_0\),

\[
 {d\over d\tau}\int_U\rho_a\,dV_0
 =-\int_{\partial U}J_a\!\cdot N\,dA_0+\int_U\sigma_a\,dV_0 .
 \tag{B-005}
\]

The sign convention makes \(J_a\cdot N>0\) outward. B-004 follows from B-005
only with the regularity needed for localization. Before a clock gauge is
chosen, the invariant density form is

\[
 {d\over ds}\rho_a+\operatorname{Div}_0\mathcal J_a=\mathcal S_a,
 \quad
 (\rho_a,\mathcal J_a,\mathcal S_a)'
 =(\rho_a,\mathcal J_a/f',\mathcal S_a/f') ,
 \tag{B-006}
\]

under \(s'=f(s)\). Thus the spatial flux and source are order-density objects;
their integrals are invariant, while a numerical rate with respect to \(s\) is
not physical.

After METRIC-001 supplies an effective spacetime representation, the same
template is

\[
 \boxed{\nabla^{\rm eff}_\mu J_a^{\mu}=\Sigma_a,}
 \qquad
 \int_{\partial\Omega}J_a^\mu d\Sigma_\mu
 =\int_\Omega\Sigma_a\,dV_{g^{\rm eff}} .
 \tag{B-007}
\]

This is a covariant representation of a balance, not an independent law and
not a selection of \(G\), \(J_a\), or \(\Sigma_a\). A conservation law is the
special case \(\sigma_a=0\) (or \(\Sigma_a=0\)), plus vanishing boundary flux
for conservation of the integrated charge.

### Tier 2 — conditional Noether balances

For a chosen local action and an actual continuous variational symmetry with
generator \(\xi\), the associated current has the on-shell form

\[
 \nabla^{\rm eff}_\mu J_\xi^\mu=\Sigma_\xi .
 \tag{B-008}
\]

Exact invariance and no external/boundary breaking give \(\Sigma_\xi=0\).
Explicit inhomogeneity, prescribed loading, boundary work, or interaction with
an excluded subsystem appears in \(\Sigma_\xi\). B-008 states the Noether
structure only; no Euler--Lagrange equation is derived here.

## 2. Audit of requested candidate quantities

| Candidate | Native status | Density / flux / source status |
|---|---|---|
| medium content | The same one medium persists, but no extensive content density is defined; no continuity equation follows | \(\rho_m,J_m,\sigma_m\) require a realization and a declared additive content measure |
| elastic stored energy | \(W(C)\) and \(P_C=DW\) exist; B-003 is exact | total energy density and energy flux require kinetic/inertial structure, clock gauge, work pairing, and boundary loading; conservation additionally requires clock-shift symmetry |
| propagation | propagation generates duration and obeys the effective causal cone | no propagation-number/amplitude density was accepted, so no propagation conservation law follows; ray/phase transport needs process closure |
| generalized momentum | \(p=D_v\mathscr L\) exists after an action is chosen | conserved momentum/current requires translation symmetry of the full action; objectivity alone is insufficient |
| angular momentum | only an objective momentum map is structurally available | a conserved angular-momentum-type current requires a continuous global rotation action and invariant full action |
| action | \(\mathfrak S\) is reparameterization invariant and stationary on admissible histories | action is a functional, not a transported charge; B-002, not “action conservation,” is the native consequence |
| topology | path-component membership is preserved by continuous admissible histories | no universal local density/flux exists; a local topological current requires a separately defined invariant from existing \(q\) data |

## 3. Flux and source classification

A flux is defined operationally by B-005: it is the boundary term whose outward
integral changes the regional amount. The foundations do not fix component
formulae. Admissible source classes are:

| Class | Meaning | Integrated status |
|---|---|---|
| \(\sigma_a=0\) | locally conservative channel | global charge constant only for zero net boundary flux |
| \(\sigma_a^{\rm ext}\) | prescribed exchange with what the regional accounting excludes, including boundary/body loading | changes the accounted regional total |
| \(\sigma_a^{\rm int}\) | exchange between already-defined internal partitions | must sum to zero over a complete closed partition: \(\sum_A\sigma_{a,A}^{\rm int}=0\) |
| \(\sigma_a^{\rm diss}\) | irreversible conversion compatible with a declared dissipation inequality | requires constitutive/dynamical closure; not supplied by the accepted \(W\) |
| \(\sigma_a^{\rm geom}\) | connection terms seen after splitting a covariant tensor balance into chart components | not physical production; retained inside the covariant divergence in B-007 |
| \(\sigma_a^{\rm break}\) | explicit breaking of a conditional action symmetry | prevents the corresponding Noether conservation law |

No internal partition is present in the minimal state definition, so
\(\sigma^{\rm int}\) is only a rule for a later derived decomposition, not a
new sector. For the complete one-medium accounting, exchanges among such
partitions cancel; exchange with an external physical substrate is forbidden by
FP-1, though boundary flux through a selected subregion is allowed.

## 4. Symmetry--Noether audit

| Structure | Unconditional consequence | Extra premise needed for a conserved current |
|---|---|---|
| objectivity / internal gauge | physical quantities descend to the quotient; symmetry directions annihilate the first variation | continuous global group action on the chosen full action, invariant measure and boundary data; then momentum-map/angular type current |
| material relabeling / covariance | Noether identities and constraints among representative equations | a non-gauge global subgroup with a well-defined charge; gauge covariance alone gives no new observable charge |
| order reparameterization | degree-one integrand, B-002, radial Legendre degeneracy | none for B-002; it never implies physical energy conservation |
| emergent clock translation | nothing before a clock and clock-gauge action exist | full action independent of calibrated \(\tau\), with invariant boundary conditions; then energy current |
| spatial translation | nothing from isotropy alone | homogeneous realization/reference and invariant full action; then momentum current |
| spatial rotation | objective constitutive response | actual global rotational symmetry of full action and boundaries; then angular current |
| effective diffeomorphism covariance | tensorial B-007 and differential/gauge identities | covariance alone does not make a stress tensor divergence-free; a selected coupled action and source accounting are required |
| local Lorentz invariance | V11 matching constraint on \(G\) and observables | selected effective action with Lorentz symmetry; only then corresponding currents/identities |

## 5. Interaction, reversibility, and dissipation

The accepted state-local hyperelastic law is reversible on one elastic branch:
for any closed smooth cycle, \(\oint P_C:dC=\oint dW=0\). It contains no memory
or hysteresis. A stationary degree-one native action supplies a reversible
variational family only when the chosen integrand and admissible boundary data
also admit the reversed history of S2-015; stationarity by itself is not a
proof of time reversal.

Irreversible evolution is not ontologically forbidden, but it is not selected.
The mathematically admissible closure families are:

1. **reversible variational:** choose an admissible member of
   \(\mathcal L_{\rm nat}\), or after duration calibration the conditional
   clock-gauge family, and impose reversal-compatible action/boundaries;
2. **constrained reversible:** the same with gauge, hard-boundary, or other
   ideal constraints that do no dissipative work;
3. **dissipative rate closure:** after duration calibration, choose an
   objective, gauge-basic dissipation functional of existing \((q,v,C)\) and
   require a nonnegative dissipation rate \(\mathcal D\ge0\);
4. **rate-independent irreversible closure:** an objective positively
   homogeneous dissipation potential may generate path dependence, but genuine
   memory/hysteresis beyond dependence on the complete \(q\) would require
   additional state structure and is not presently authorized.

Families 3--4 require constitutive selection, normalization and a balance of
power/energy. They are admissible future closures, not derived laws. No entropy
field or temperature field is introduced, so no entropy balance is available.

## 6. Weak-field and V11 compatibility

Near \(q_0\),

\[
 C=\mathbf1+2\varepsilon+O(\varepsilon^2),\quad
 W={\lambda\over2}(\operatorname{tr}\varepsilon)^2
 +\mu\operatorname{tr}(\varepsilon^2)+O(|\varepsilon|^3),
 \tag{B-009}
\]

and METRIC-001 gives

\[
 g^{\rm eff}=\eta+R_0[\delta q]+O(\|\delta q\|^2),\qquad
 [R_0\delta q]_{\rm gauge}=[h^{\rm V11}]_{\rm gauge}.
 \tag{B-010}
\]

B-003 then reduces to the standard quadratic elastic power identity, while
B-007 reduces in a local inertial frame at a point to
\(\partial_\mu J_a^\mu=\Sigma_a\). DURATION-001 supplies the invariant clock
parameter and \(d\tau^2=-c^{-2}g^{\rm eff}_{\mu\nu}dx^\mu dx^\nu\); it never
identifies \(s\) with time. Thus the balance architecture admits the V11 local
Lorentz/weak-field limit without fixing \(\lambda,\mu,R_0\), a source map, or a
field equation. V11's homogeneous quantities are not reinterpreted as local
densities or strains, and V11 is unchanged.

## 7. Prerequisites for governing evolution equations

The exact dependency graph is

```text
q in Q_adm
  -> C[q,q0], v-ray, admissible history
  -> balance identities/templates (B-002--B-007)
  -> constitutive closure:
       W=Phi(I1,I2,I3), P_C=DW,
       selected kinetic/inertial and optional dissipative response,
       selected flux/source and boundary-work maps
  -> metric/clock closure G[q,C;D] and tau calibration
  -> source projection into the effective representation
  -> governing evolution equations + constraints + initial/boundary data
```

Balance restricts the form of any future equation: internal exchanges must
cancel, conservative channels must be divergences, gauge identities must be
respected, dissipation must have the declared sign, and the effective form must
be covariant. Constitutive laws determine the fluxes and productions; the
metric map determines operational spacetime representation; stored energy
supplies only \(P_C=DW\). None of these ingredients alone is an evolution
equation.

Mandatory remaining closure data are: a local realization and function spaces
for \(q\); one action or nonvariational kinetic law; emergent-duration gauge;
inertial normalization; constitutive \(\Phi\); any dissipative potential;
boundary/loading data; one \(G\in\mathfrak G\); native-to-effective source map;
and well-posedness/causality and V11 matching proofs.

## 8. Equation traceability

| Equation | Content | Fixed premises | Status |
|---|---|---|---|
| B-001 | no universal nontrivial charge follows | D-001--D-012; absence of content/propagation densities and clock translation | derived negative result |
| B-002 | parameter-Hamiltonian identity | D-001, D-006--D-007 | unconditional for a selected differentiable native action |
| B-003 | elastic storage identity | H-001, H-007 | unconditional on smooth elastic branch |
| B-004--B-006 | local/integral material balance template | locality, \(dV_0\), DU-002--DU-003, sufficient regularity | structural family; quantities require closure |
| B-007 | effective covariant balance template | M-001--M-002 plus B-004 | representation/matching family |
| B-008 | conditional Noether current | chosen symmetric action and admissible boundaries | conditional; no field equation derived |
| B-009 | weak elastic limit | DEFORMATION-001; H-011--H-014; EP-005--EP-006 | compatibility expansion |
| B-010 | weak metric limit | M-009--M-010 | V11 matching condition |

## 9. Completion boundary

BALANCE-001 supplies the native hierarchy, the most general local and covariant
balance forms, flux definition, complete source classification at the available
level, Noether audit, interaction/dissipation families, weak-field audit, and
closure dependency graph. It derives no Einstein equation, Euler--Lagrange
equation, constitutive model, new field, coupling, constant, observation fit, or
V11 modification.
