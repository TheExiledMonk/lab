# PBUF Native Microphysics Working System

## What exists today

The implemented system has two adjacent but uncoupled lanes. Static source arrays initialize frozen fast/slow scalar fields, whose mean `c_state` drives a nonlinear six-neighbor bounded-strain equilibrium solve for accumulated deformation `u`. Dynamic excitation is a signed two-component array on a one-dimensional periodic support. A progression step moves every two-vector to exactly one neighbor; the frame-aware form orthogonally changes its two coordinates. No neighbor branching, ray solver, renormalization, or loaded-state transfer factor participates.

Packet density is the squared two-component magnitude. Its density-weighted index is the centroid, and the centroid history is the diagnostic native packet trajectory. The current excitation observer does not compute direction or curvature. Dev151 unifies an audit longitudinal array `L` and the two excitation components as `(L,X1,X2)`, but `L` is not derived from production `c_state` in that path. Dev153 establishes no loaded/unloaded cross-effect.

## One-page equations

- **E01 A8 N6 mean:** $\bar u_a=\frac16\sum_{b\in N_6(a)}u_b$ — `pbuf/models/a8_state.py:40-77`
- **E02 fast update:** $u_f'=clip[u_f+0.03((\bar u_f-u_f)+0.5(u_s-u_f)),-5,5]$ — `pbuf/models/a8_state.py:101-112`
- **E03 slow update:** $u_s'=clip[u_s+0.003((\bar u_s-u_s)+0.5(u_f-u_s)),-5,5]$ — `pbuf/models/a8_state.py:101-112`
- **E04 combined loading state:** $c_{state}=(u_s+u_f)/2$ — `pbuf/models/a8_state.py:89-114`
- **E05 bond strain:** $\epsilon_{ab}=(u_b-u_a)/\Delta x$ — `pbuf/labs/foundation/c_state_bounded_strain_bridge001.py:128-133`
- **E06 bounded energy:** $W=-\frac{K\epsilon_{max}^2}{2}\ln(1-\epsilon^2/\epsilon_{max}^2)$ — `pbuf/wl/native_incremental_elastic_energy.py:18-26`
- **E07 bounded stress:** $\sigma=K\epsilon/(1-\epsilon^2/\epsilon_{max}^2)$ — `pbuf/wl/native_incremental_elastic_energy.py:29-33`
- **E08 tangent stiffness:** $K_{tan}=K(1+q)/(1-q)^2,\ q=(\epsilon/\epsilon_{max})^2$ — `pbuf/wl/native_incremental_elastic_energy.py:36-41`
- **E09 nonlinear equilibrium:** $A(u)_a=\Delta x^{-2}\sum_{b\in N_6(a)}K\frac{u_a-u_b}{1-\epsilon_{ab}^2/\epsilon_{max}^2}=c_{state,a}$ — `pbuf/labs/foundation/c_state_bounded_strain_bridge001.py:136-220`
- **E10 packet initialization:** $X_a=A\exp[-\tfrac12((a-a_0)/w)^2]\,p/\|p\|$ — `pbuf/excitation/native_excitation_state.py:24-32`
- **E11 free progression:** $X_a^{n+1}=X_{(a-d)\bmod N}^{n},\ d\in\{-1,+1\}$ — `pbuf/excitation/native_excitation_transfer.py:24-37`
- **E12 frame-aware progression:** $X_j^{n+1}=R_{ij}X_i^n,\ j=(i+1)\bmod N$ — `pbuf/foundation/native_neighbor_mixed_state.py:17-33`
- **E13 orthogonal frame map:** $R_{ij}=UV^T,\quad B_{ij}=E_{\perp,j}E_{\perp,i}^T=U\Sigma V^T$ — `pbuf/foundation/native_neighbor_state.py:48-52`
- **E14 excitation norm:** $N_X=\sum_a(X_{1,a}^2+X_{2,a}^2)$ — `pbuf/excitation/native_excitation_invariants.py:7-15`
- **E15 excitation density:** $\rho_{X,a}=X_{1,a}^2+X_{2,a}^2$ — `pbuf/foundation/native_neighbor_mixed_observer.py:6-9`
- **E16 packet centroid:** $\bar r_n=\sum_a a\rho_{X,a}^{(n)}/\sum_a\rho_{X,a}^{(n)}$ — `pbuf/foundation/native_neighbor_mixed_observer.py:6-9`
- **E17 FFT wavelength:** $m_n=argmax_{m\ge1}|FFT(X_1)_m|^2,\ \lambda_n=N/max(m_n,1)$ — `pbuf/foundation/native_neighbor_mixed_observer.py:10-12`
- **E18 unified state:** $\mathcal S=(L,X_1,X_2)$ — `pbuf/foundation/native_neighbor_state.py:54-70`

## The system without equations

Imagine two ledgers for the same conceptual medium. The static ledger spreads and balances persistent loading across six-faced grid neighbors, while enforcing a strain barrier. The dynamic ledger stores two signed disturbance numbers and shifts each pair one site per progression step. When local link frames differ, those two numbers are rotated without changing their combined squared size. Afterward, an observer squares them, locates their center, and records how that center changes. The code does not currently let the static ledger modify the dynamic update.

CARPET_FIBER_ANALOGY: A link's longitudinal entry can be pictured as fiber condition and its two transverse entries as sideways disturbance. The implementation supports the bookkeeping analogy; it does not prove a literal carpet-fiber ontology.

## Definitive flowcharts

```text
STATIC: rho -> u_slow/u_fast N6 updates -> c_state -> embed -> bounded N6 equilibrium -> accumulated u -> bond epsilon -> sigma/W/Ktan
DYNAMIC: packet controls -> X(N,2) -> one-neighbor periodic permutation -> optional orthogonal frame map -> X history -> rho_X -> centroid history (trajectory)
COMBINED: (static L, dynamic X1,X2) share rank-3 records; [NO ESTABLISHED L -> X UPDATE EDGE]
```

## ONE EXCITATION PROGRESSION STEP

1. Read `old[i]`, a two-number signed vector, and source frame `E_i`.
2. Select exactly `j=(i+1) mod N` (or the opposite direction in the basic transfer). There is no candidate list, branching, or weight.
3. Form raw transverse overlap `B=E_perp,j E_perp,i^T`; take its SVD `B=U Sigma V^T`; set `R=UV^T`.
4. Write `next[j]=R old[i]`. Each destination has exactly one writer.
5. Replace the whole transverse array with `next`; leave `L` unchanged; append a copy to history.
6. No renormalization, clipping, averaging, timestep, or ray update occurs.

## ONE STATIC EQUILIBRIUM UPDATE

For current `u`, calculate positive-axis bond differences divided by `DX`. Form secant bond weights `K/(1-(|epsilon|/epsilon_max)^2)`. Assemble the six-neighbor divergence operator `A(u)`, solve `A(candidate)=c_state` by conjugate gradient with zero boundary values, then damp: `u_new=0.65 candidate+0.35 u`. Stop when relative RMS change is at most `2e-7`, or after 30 Picard iterations.

## FROM EXCITATION STATE TO TRAJECTORY POINT

At every site, calculate `rho=X1^2+X2^2`. Calculate `rbar=sum(a*rho)/sum(rho)`. Repeating that observer operation for each stored state produces the centroid sequence. There is no current excitation function implementing `rbar_n,rbar_(n+1)->direction` or `direction->curvature`.

## D01–D30 domain reconstruction

### D01 — native spatial topology

#### Name

native spatial topology

#### Status

ESTABLISHED

#### Purpose

Provide Cartesian nearest-neighbor support

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

N_6=\{\pm\hat x,\pm\hat y,\pm\hat z\}

#### Update/evolution equation

N_6=\{\pm\hat x,\pm\hat y,\pm\hat z\}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Provide Cartesian nearest-neighbor support

#### Code implementation

pbuf/core/conventions.py:177-191; pbuf/models/a8_state.py:40-77

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 1

#### What it DOES establish

Provide Cartesian nearest-neighbor support

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Provide Cartesian nearest-neighbor support using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D02 — node/link representation

#### Name

node/link representation

#### Status

REPRESENTATIONAL

#### Purpose

Distinguish node arrays from audit link state

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

u,c: grid scalars; \mathcal S_i=(L_i,X_{1i},X_{2i})

#### Update/evolution equation

u,c: grid scalars; \mathcal S_i=(L_i,X_{1i},X_{2i})

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Distinguish node arrays from audit link state

#### Code implementation

pbuf/excitation/native_excitation_state.py:34-47; pbuf/foundation/native_neighbor_state.py:54-70

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 2

#### What it DOES establish

Distinguish node arrays from audit link state

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Distinguish node arrays from audit link state using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D03 — static loading input

#### Name

static loading input

#### Status

ESTABLISHED

#### Purpose

Turn rho into fast/slow initialization

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

u_s^0=u_f^0=\rho;\ c=(u_s+u_f)/2

#### Update/evolution equation

u_s^0=u_f^0=\rho;\ c=(u_s+u_f)/2

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Turn rho into fast/slow initialization

#### Code implementation

pbuf/models/a8_state.py:117-136

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 3

#### What it DOES establish

Turn rho into fast/slow initialization

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Turn rho into fast/slow initialization using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D04 — local medium loading/state

#### Name

local medium loading/state

#### Status

ESTABLISHED

#### Purpose

Store fast, slow, and combined scalar response

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

c_{state}=(u_s+u_f)/2

#### Update/evolution equation

c_{state}=(u_s+u_f)/2

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Store fast, slow, and combined scalar response

#### Code implementation

pbuf/models/a8_state.py:80-114

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 4

#### What it DOES establish

Store fast, slow, and combined scalar response

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Store fast, slow, and combined scalar response using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D05 — neighbor equilibrium

#### Name

neighbor equilibrium

#### Status

ESTABLISHED

#### Purpose

Solve discrete nonlinear force balance

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

A(u)=c_{state}

#### Update/evolution equation

A(u)=c_{state}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Solve discrete nonlinear force balance

#### Code implementation

pbuf/labs/foundation/c_state_bounded_strain_bridge001.py:128-220

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 5

#### What it DOES establish

Solve discrete nonlinear force balance

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Solve discrete nonlinear force balance using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D06 — bounded strain

#### Name

bounded strain

#### Status

ESTABLISHED

#### Purpose

Bound bond deformation constitutively

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\epsilon=(u_b-u_a)/\Delta x;\ W,\sigma,K_{tan}\text{ as E06--E08}

#### Update/evolution equation

\epsilon=(u_b-u_a)/\Delta x;\ W,\sigma,K_{tan}\text{ as E06--E08}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Bound bond deformation constitutively

#### Code implementation

pbuf/wl/native_incremental_elastic_energy.py:18-41

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 6

#### What it DOES establish

Bound bond deformation constitutively

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Bound bond deformation constitutively using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D07 — accumulated deformation response

#### Name

accumulated deformation response

#### Status

ESTABLISHED

#### Purpose

Produce long-range u from local c_state

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

u=A^{-1}_{bounded}(c_{state})

#### Update/evolution equation

u=A^{-1}_{bounded}(c_{state})

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Produce long-range u from local c_state

#### Code implementation

pbuf/labs/foundation/c_state_bounded_strain_bridge001.py:223-230

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 7

#### What it DOES establish

Produce long-range u from local c_state

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Produce long-range u from local c_state using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D08 — dynamic excitation state

#### Name

dynamic excitation state

#### Status

ESTABLISHED

#### Purpose

Store two signed components per site

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

X\in\mathbb R^{N\times2}

#### Update/evolution equation

X\in\mathbb R^{N\times2}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Store two signed components per site

#### Code implementation

pbuf/excitation/native_excitation_state.py:34-47

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 8

#### What it DOES establish

Store two signed components per site

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Store two signed components per site using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D09 — excitation initialization

#### Name

excitation initialization

#### Status

ESTABLISHED

#### Purpose

Create a localized Gaussian vector packet

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

X_a=Ae^{-((a-a_0)/w)^2/2}\hat p

#### Update/evolution equation

X_a=Ae^{-((a-a_0)/w)^2/2}\hat p

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Create a localized Gaussian vector packet

#### Code implementation

pbuf/excitation/native_excitation_state.py:24-32

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 9

#### What it DOES establish

Create a localized Gaussian vector packet

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Create a localized Gaussian vector packet using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D10 — excitation neighbor transfer

#### Name

excitation neighbor transfer

#### Status

ESTABLISHED

#### Purpose

Move all state to one periodic neighbor

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

X_a^{n+1}=X_{(a-d)\bmod N}^n\text{ or }R X

#### Update/evolution equation

X_a^{n+1}=X_{(a-d)\bmod N}^n\text{ or }R X

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Move all state to one periodic neighbor

#### Code implementation

pbuf/excitation/native_excitation_transfer.py:24-37; pbuf/foundation/native_neighbor_mixed_state.py:17-33

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 10

#### What it DOES establish

Move all state to one periodic neighbor

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Move all state to one periodic neighbor using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D11 — two transverse modes

#### Name

two transverse modes

#### Status

ESTABLISHED

#### Purpose

Represent independent transverse coordinates

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

X=(X_1,X_2)

#### Update/evolution equation

X=(X_1,X_2)

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Represent independent transverse coordinates

#### Code implementation

pbuf/excitation/native_excitation_invariants.py:24-26

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 11

#### What it DOES establish

Represent independent transverse coordinates

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Represent independent transverse coordinates using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D12 — excitation norm

#### Name

excitation norm

#### Status

ESTABLISHED

#### Purpose

Audit the quadratic invariant

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

N_X=\sum_a\|X_a\|^2

#### Update/evolution equation

N_X=\sum_a\|X_a\|^2

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Audit the quadratic invariant

#### Code implementation

pbuf/excitation/native_excitation_invariants.py:7-15

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 12

#### What it DOES establish

Audit the quadratic invariant

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Audit the quadratic invariant using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D13 — spatial wavelength

#### Name

spatial wavelength

#### Status

DIAGNOSTIC

#### Purpose

Estimate spatial periodicity from X1

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\lambda=N/argmax_{m\ge1}|FFT(X_1)_m|^2

#### Update/evolution equation

\lambda=N/argmax_{m\ge1}|FFT(X_1)_m|^2

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Estimate spatial periodicity from X1

#### Code implementation

pbuf/foundation/native_neighbor_mixed_observer.py:10-12

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 13

#### What it DOES establish

Estimate spatial periodicity from X1

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Estimate spatial periodicity from X1 using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D14 — native k

#### Name

native k

#### Status

UNRESOLVED

#### Purpose

Record a wave-number quantity

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

none

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

NO\ IMPLEMENTED\ k_n\ IN\ CANONICAL\ OBSERVER

#### Update/evolution equation

None implemented

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Record a wave-number quantity

#### Code implementation

NO_CANONICAL_SYMBOL / no current implementation

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 14

#### What it DOES establish

Record a wave-number quantity

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

The code does not currently calculate this quantity for excitation packets.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D15 — interference

#### Name

interference

#### Status

ESTABLISHED

#### Purpose

Demonstrate signed linear superposition

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

P(a+b)=P(a)+P(b)

#### Update/evolution equation

P(a+b)=P(a)+P(b)

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Demonstrate signed linear superposition

#### Code implementation

pbuf/excitation/native_excitation_invariants.py:16-19

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 15

#### What it DOES establish

Demonstrate signed linear superposition

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Demonstrate signed linear superposition using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D16 — polarization

#### Name

polarization

#### Status

ESTABLISHED

#### Purpose

Select a normalized two-vector orientation

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\hat p=p/\|p\|;\ X\mapsto QX

#### Update/evolution equation

\hat p=p/\|p\|;\ X\mapsto QX

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Select a normalized two-vector orientation

#### Code implementation

pbuf/excitation/native_excitation_state.py:24-32

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 16

#### What it DOES establish

Select a normalized two-vector orientation

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Select a normalized two-vector orientation using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D17 — handedness

#### Name

handedness

#### Status

ESTABLISHED

#### Purpose

Represent phase-quadrature sign classes in audit packets

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

X_2=\pm\,env\sin(ka)\text{ in Dev152 fixtures}

#### Update/evolution equation

X_2=\pm\,env\sin(ka)\text{ in Dev152 fixtures}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Represent phase-quadrature sign classes in audit packets

#### Code implementation

pbuf/foundation/native_neighbor_loaded_excitation.py:10-21

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 17

#### What it DOES establish

Represent phase-quadrature sign classes in audit packets

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Represent phase-quadrature sign classes in audit packets using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D18 — rank-3 unified neighbor representation

#### Name

rank-3 unified neighbor representation

#### Status

REPRESENTATIONAL

#### Purpose

Place L and X in one rank-3 record

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\mathcal S=(L,X_1,X_2)

#### Update/evolution equation

\mathcal S=(L,X_1,X_2)

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Place L and X in one rank-3 record

#### Code implementation

pbuf/foundation/native_neighbor_state.py:54-70

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 18

#### What it DOES establish

Place L and X in one rank-3 record

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Place L and X in one rank-3 record using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D19 — local frames

#### Name

local frames

#### Status

REPRESENTATIONAL

#### Purpose

Construct deterministic link triads

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

E=(e_\parallel,e_1,e_2)

#### Update/evolution equation

E=(e_\parallel,e_1,e_2)

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Construct deterministic link triads

#### Code implementation

pbuf/foundation/native_neighbor_state.py:26-32

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 19

#### What it DOES establish

Construct deterministic link triads

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Construct deterministic link triads using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D20 — frame transport

#### Name

frame transport

#### Status

ESTABLISHED

#### Purpose

Map transverse coordinates between frames

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

R=UV^T,\ R^TR=I

#### Update/evolution equation

R=UV^T,\ R^TR=I

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Map transverse coordinates between frames

#### Code implementation

pbuf/foundation/native_neighbor_state.py:48-52

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 20

#### What it DOES establish

Map transverse coordinates between frames

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Map transverse coordinates between frames using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D21 — mixed loaded+excited state

#### Name

mixed loaded+excited state

#### Status

REPRESENTATIONAL

#### Purpose

Run L and X together without backreaction

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

L^{n+1}=L^n;\ X_j^{n+1}=R_{ij}X_i^n

#### Update/evolution equation

L^{n+1}=L^n;\ X_j^{n+1}=R_{ij}X_i^n

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Run L and X together without backreaction

#### Code implementation

pbuf/foundation/native_neighbor_mixed_state.py:17-33

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 21

#### What it DOES establish

Run L and X together without backreaction

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Run L and X together without backreaction using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D22 — currently absent cross-effect

#### Name

currently absent cross-effect

#### Status

UNRESOLVED

#### Purpose

Record failure to derive C(L)

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\mathcal C_{ij}=I\text{ is the established null}

#### Update/evolution equation

\mathcal C_{ij}=I\text{ is the established null}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Record failure to derive C(L)

#### Code implementation

pbuf/foundation/native_loaded_link_response.py:52-64

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 22

#### What it DOES establish

Record failure to derive C(L)

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Record failure to derive C(L) using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D23 — excitation density

#### Name

excitation density

#### Status

DERIVED

#### Purpose

Measure packet morphology

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\rho_X=X_1^2+X_2^2

#### Update/evolution equation

\rho_X=X_1^2+X_2^2

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Measure packet morphology

#### Code implementation

pbuf/foundation/native_neighbor_mixed_observer.py:6-9

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 23

#### What it DOES establish

Measure packet morphology

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Measure packet morphology using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D24 — packet evolution

#### Name

packet evolution

#### Status

ESTABLISHED

#### Purpose

Apply permutation and optional frame map

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

X^{n+1}=P_R X^n

#### Update/evolution equation

X^{n+1}=P_R X^n

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Apply permutation and optional frame map

#### Code implementation

pbuf/foundation/native_neighbor_mixed_state.py:17-33

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 24

#### What it DOES establish

Apply permutation and optional frame map

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Apply permutation and optional frame map using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D25 — packet centroid

#### Name

packet centroid

#### Status

DERIVED

#### Purpose

Locate density-weighted packet center

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\bar r=\sum a\rho_a/\sum\rho_a

#### Update/evolution equation

\bar r=\sum a\rho_a/\sum\rho_a

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Locate density-weighted packet center

#### Code implementation

pbuf/foundation/native_neighbor_mixed_observer.py:6-9

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 25

#### What it DOES establish

Locate density-weighted packet center

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Locate density-weighted packet center using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D26 — native trajectory

#### Name

native trajectory

#### Status

DIAGNOSTIC

#### Purpose

Store centroid sequence as diagnostic path

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

X arrays/history

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\Gamma=(\bar r_0,\ldots,\bar r_N)

#### Update/evolution equation

\Gamma=(\bar r_0,\ldots,\bar r_N)

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Store centroid sequence as diagnostic path

#### Code implementation

pbuf/foundation/native_neighbor_mixed_observer.py:6-9

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 26

#### What it DOES establish

Store centroid sequence as diagnostic path

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Store centroid sequence as diagnostic path using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D27 — native direction

#### Name

native direction

#### Status

UNRESOLVED

#### Purpose

Would derive motion direction

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

none

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

NO\ EXCITATION\ DIRECTION\ ESTIMATOR

#### Update/evolution equation

None implemented

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Would derive motion direction

#### Code implementation

NO CURRENT EXCITATION IMPLEMENTATION

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 27

#### What it DOES establish

Would derive motion direction

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

The code does not currently calculate this quantity for excitation packets.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D28 — native curvature

#### Name

native curvature

#### Status

UNRESOLVED

#### Purpose

Would derive turning of direction

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

none

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

NO\ EXCITATION\ CURVATURE\ ESTIMATOR

#### Update/evolution equation

None implemented

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Would derive turning of direction

#### Code implementation

NO CURRENT EXCITATION IMPLEMENTATION

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 28

#### What it DOES establish

Would derive turning of direction

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

The code does not currently calculate this quantity for excitation packets.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D29 — historical geometric trajectory machinery

#### Name

historical geometric trajectory machinery

#### Status

HISTORICAL

#### Purpose

Analyze pre-excitation geometric path state

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

\kappa=\Delta\theta/\Delta s\text{ for geometric paths}

#### Update/evolution equation

\kappa=\Delta\theta/\Delta s\text{ for geometric paths}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Analyze pre-excitation geometric path state

#### Code implementation

pbuf/wl/trajectory_state.py:100-148

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 29

#### What it DOES establish

Analyze pre-excitation geometric path state

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Analyze pre-excitation geometric path state using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

### D30 — weak-lensing propagation machinery

#### Name

weak-lensing propagation machinery

#### Status

ESTABLISHED

#### Purpose

Propagate geometric optical rays/bundles

#### Physical meaning

Implementation-level native quantity or diagnostic; not an identification with established external physics.

#### Stored variables

scalar grid fields or geometric records

#### Units/dimensional status

native dimensionless/index units; absolute physical closure unresolved

#### Mathematical definition

r_{n+1}=G(r_n,d_n,medium)\text{ in separate WL stack}

#### Update/evolution equation

r_{n+1}=G(r_n,d_n,medium)\text{ in separate WL stack}

#### Inputs

upstream quantities listed in dependency graph

#### Outputs

Propagate geometric optical rays/bundles

#### Code implementation

pbuf/wl/propagation.py; pbuf/wl/trajectory_state.py

#### Upstream dependencies

See native_dependency_graph.json

#### Downstream consumers

See native_dependency_graph.json

#### Conservation/invariants

Quadratic X norm for dynamic permutation/orthogonal transport; static solver residual for equilibrium

#### Known controls/tests

Dev145-153 contracts and Dev154 source reconstruction; domain 30

#### What it DOES establish

Propagate geometric optical rays/bundles

#### What it DOES NOT establish

New cross-coupling, absolute units, or identification with gravity/EM/QM

PLAIN_LANGUAGE:

Propagate geometric optical rays/bundles using the concrete code path cited above.

CARPET_FIBER_ANALOGY:

A grid bond may be pictured as a fiber connection; this is analogy only, not an ontological claim.

## Historical trajectory distinction

| System | What evolves | Driver | Physical excitation? | Path imposed/extracted | Consumer |
|---|---|---|---|---|---|
| Historical trajectory state | positions/directions | geometric propagation | No | imposed/evolved geometry | receiver diagnostics |
| Dev148/149 excitation | signed `X(N,2)` | periodic permutation | Yes, native candidate | centroid extracted | excitation audits |
| Dev153 packet path | centroid history | same excitation history | diagnostic of X | extracted | loaded/unloaded comparison |
| Weak-lensing 3D rays | ray position/direction/bundle | frozen WL field | optical tracer, separate | solver-evolved | shear/receiver stack |

## User-supplied Dev153 statement verification

- `X_j = R_ij X_i` — **EXACT**: Exact for frame-aware one-neighbor progression.
- `R_ij^T R_ij = I` — **EXACT**: SVD polar factor enforces this numerically.
- `|X_j|^2 = |X_i|^2` — **EXACT**: Follows per transported vector.
- `C_ij = I in established Dev153 result` — **APPROXIMATE_DESCRIPTION**: Dev153 establishes absence of a derived nontrivial C; identity is the null control, not a separately discovered constitutive law.
- `rho_a = X1_a^2 + X2_a^2` — **EXACT**: Named amp in observer.
- `centroid is density-weighted spatial mean` — **EQUIVALENT**: It is a 1D index-weighted linear mean, not a general spatial vector or circular mean.
- `trajectory is centroid sequence` — **EQUIVALENT**: The code returns centroid history; packet_path artifacts interpret it as trajectory.
- `direction comes from centroid differences` — **INCORRECT**: No excitation direction estimator exists.
- `curvature comes from tangent differences` — **INCORRECT**: No excitation curvature estimator exists; separate WL geometry code has curvature.
- `loaded and unloaded packet paths are identical in Dev153` — **EXACT**: Dev153 execute applies the same fixed rotation irrespective of load.

## Parameters and hidden numerical operations

The machine-readable register contains 11 active/historical constants. Important hidden operations are the N6 factor `1/6`, fast/slow mean `1/2`, polarization normalization, SVD orthogonalization, periodic modulo indexing, centroid denominator guard, FFT exclusion of DC, static clipping, Picard damping, strain barrier guard, and solver tolerances.

## What we must stop re-deriving

- Native signed two-component excitation state and its two transverse modes.
- Source-free periodic permutation and exact quadratic-norm conservation.
- FFT-bin wavelength diagnostic, polarization covariance, and handedness fixtures.
- Rank-3 `(L,X1,X2)` representation and orthogonal F02–F06 frame transport.
- Excitation-density centroid trajectory extraction.
- Static N6 fast/slow response and bounded-strain accumulation.

## Real remaining gaps

- **ACTUAL_MISSING_PHYSICS** — loaded longitudinal state -> transverse evolution operator: No mixed Hessian, metric, allocation, or conservative exchange map is implemented.
- **ACTUAL_MISSING_PHYSICS** — excitation centroid history -> direction: No excitation observer implements the finite difference or zero-step rule.
- **IMPLEMENTED_UNDER_DIFFERENT_NAME** — excitation direction -> curvature: Curvature exists only in separate geometric WL trajectory code.
- **ACTUAL_MISSING_PHYSICS** — FFT wavelength -> native k: The canonical observer emits wavelength but not k.

## Retired mechanisms

Do not revive `Rmax`, historical `strength=0.18`, trajectory-as-physical-excitation, fundamental native time/`T0`, reflective boundaries as conservative zero flux, or fitted loading-to-excitation transfer factors.
