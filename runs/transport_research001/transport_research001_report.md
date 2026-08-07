# PBUF TRANSPORT-RESEARCH-001 -- Comparative Analysis of Native Wave Transport Mechanisms

## 0. Decision

**Outcome B.** No existing physical transport mechanism maps naturally
onto the existing PBUF ontology with minimal additional assumptions.

The six required systems fall into two structural families:

* **Family A (mechanical):** water surface waves (S1), elastic solids
  (S2), acoustic waves (S3), spin waves (S5), plasma / MHD (S6). All
  five share the architecture `local state + neighbour coupling +
  restoring mechanism + inertial resistance`, with the propagation
  speed set by the ratio of restoring stiffness to resistance.
* **Family B (electromagnetic):** Maxwell vacuum (S4). Unique in
  having neither a restoring mechanism nor an inertial resistance;
  propagation arises from the mutual curl coupling of `E` and `B`.

PBUF matches **Family A on the spatial half only**:

| Architecture slot        | PBUF status                |
|--------------------------|----------------------------|
| local state              | present                    |
| neighbour coupling       | present (static)           |
| local transfer           | present (static)           |
| restoring mechanism      | present                    |
| **resistance / inertia** | **MISSING**                |
| propagation speed        | not defined                |
| local steering           | available in principle     |
| wave equation            | not present                |

Family B is structurally incompatible: CORE-001's microscopic energy
has a mass-like onsite term `kappa_0|q|^2`, a scalar gradient term
`kappa_1|q_j - q_i|^2` (not a curl), and a single scalar triplet `q`
(not an `(E,B)` pair). V11's `alpha_resolved ~ 3 alpha_EM` is a
numerical identity and a dimensional-counting argument; it does not
identify `q` with an EM vector potential.

The closest structural match is therefore the elastic-solid pattern
(S2). Its completion requires adding the kinetic sector - exactly
the closure gap already identified by INERTIA-001 and EM-TRANSPORT-001.
Adding this sector is not a parameter renaming or a coordinate choice;
it is a new local physical principle, which the brief's "minimal
additional assumptions" criterion does not permit.

No ontology, field, coupling, constant, transport equation, V11
change, CORE-001 change, FOUNDATION-001 change, constitutive-law
change, cosmological result, metric construction, or weak-lensing
fit is introduced.

## 1. Method

For each of the six required systems the milestone records, in
physical language first and equations last:

1. The physical setting of the medium.
2. The disturbed quantity.
3. The local neighbour-to-neighbour interaction.
4. The quantity that is locally transferred.
5. The restoring mechanism.
6. The resistance to the response.
7. How the propagation speed is determined.
8. The local steering mechanism.
9. The minimum governing equations (no derivation).

The full per-system records are in `system_mechanism_audit.csv`.

The comparison table is built directly from these records
(`comparison_table.csv`) and is the only input to the pattern
identification step. The common architecture is then derived, not
assumed (`common_architecture.csv` and `common_architecture.json`).
Only after this comparison is complete is the common architecture
compared with the existing PBUF framework
(`pbuf_architecture_comparison.csv` and `pbuf_comparison.json`).

## 2. Per-system physical mechanism

### 2.1 Water surface waves (S1)

* **Setting.** A continuous fluid with a free surface.
* **Disturbed quantity.** Vertical displacement `eta(x,t)`.
* **Local interaction.** A displaced column exerts horizontal pressure
  on its neighbour; the neighbour accelerates and transfers momentum
  further along the surface.
* **Local transfer.** Hydrostatic pressure excess proportional to the
  height difference between adjacent columns.
* **Restoring.** Gravity (`rho g eta`).
* **Resistance.** Inertia of the water mass.
* **Speed.** `c^2 = (g/k) tanh(kh)` (deep-water limit `c^2 = g/k`).
* **Steering.** Bathymetry variations change the local `k tanh(kh)/h`
  and refract the wavefront.
* **Equations (minimum).** `d^2 eta/dt^2 = (g tanh(kh)/k) nabla^2 eta`.
  Exact: Euler + incompressibility + free-surface BC.

### 2.2 Elastic solids (S2)

* **Setting.** Continuous solid with placement `y(X,t)` and
  deformation `C = Grad y^sharp Grad y`.
* **Disturbed quantity.** Displacement `u(x,t) = y - x` (vector).
* **Local interaction.** Each material volume element pulls on its
  neighbours through the traction `P_F N` across the shared
  interface; the internal force is `Div P_F`.
* **Local transfer.** Traction across each material interface.
* **Restoring.** Gradient of the hyperelastic stored energy `W(C)`:
  `-Div P_F = -(D_y C)^* D_C W`.
* **Resistance.** Mass density `rho` of the medium.
* **Speed.** `c_alpha^2(n) = A_{iJkL} n_J n_L / rho` (acoustic tensor).
* **Steering.** Spatial variation of elastic moduli `C_ijkl(x)`.
* **Equations (minimum).** `rho u_tt = Div(A : sym grad u)`;
  nonlinear: `rho u_tt = Div P_F`.

### 2.3 Acoustic waves (S3)

* **Setting.** Continuous fluid with bulk modulus `B` and density
  `rho`.
* **Disturbed quantity.** Pressure perturbation `delta p(x,t)`.
* **Local interaction.** Compression of one volume pushes its
  neighbour through `B`; the neighbour compresses in turn.
* **Local transfer.** Pressure perturbation between adjacent volumes.
* **Restoring.** Compressibility (bulk modulus `B`).
* **Resistance.** Mass density `rho`.
* **Speed.** `c^2 = B/rho`.
* **Steering.** Spatial variation of `B(x)` or `rho(x)`.
* **Equations (minimum).** `d^2 p / dt^2 = c^2 nabla^2 p`.

### 2.4 Electromagnetic waves (S4)

* **Setting.** Vacuum with `E(x,t)` and `B(x,t)`.
* **Disturbed quantity.** `E` and `B`; equivalently the
  antisymmetric `F_{mu nu}`.
* **Local interaction.** `dE/dt` locally generates `curl B`
  (Ampere-Maxwell with `j = 0`); `dB/dt` locally generates
  `curl E` (Faraday).
* **Local transfer.** The fields themselves.
* **Restoring.** **NONE in vacuum.** No scalar potential energy
  storage.
* **Resistance.** **NONE in vacuum.** No inertia, no mass.
* **Speed.** `c^2 = 1/(epsilon_0 mu_0)` from the curl coupling.
* **Steering.** Spatial variation of `epsilon(x)`, `mu(x)`, or (in
  GR) the metric.
* **Equations (minimum).** `curl E = -dB/dt`, `curl B =
  mu_0 epsilon_0 dE/dt`, `div E = 0`, `div B = 0`. Each is
  first-order in time; the wave equation emerges only after taking
  the curl of one equation and substituting the other.

### 2.5 Spin waves / magnons (S5)

* **Setting.** Ordered magnetic lattice with local magnetization
  `m(x,t)`.
* **Disturbed quantity.** Direction of the local magnetization
  `m(x,t)`.
* **Local interaction.** Heisenberg exchange `-J sum_<ij> S_i . S_j`
  couples neighbouring spins; misalignment costs energy.
* **Local transfer.** Spin angular momentum: a precessing spin
  transfers angular momentum to its neighbour.
* **Restoring.** Anisotropy field `H_anis` and any external field
  `H_ext`.
* **Resistance.** Gyromagnetic precession: a spin cannot change
  direction instantaneously.
* **Speed.** Dispersion `omega(k) = gamma sqrt(H_eff (H_eff + D k^2))`;
  group velocity `v_g = d omega/dk`.
* **Steering.** Spatial gradients of `H_anis(x)` or `H_ext(x)`,
  spatial variation of `J(x)` or `D(x)`.
* **Equations (minimum).** `dm/dt = -gamma m x H_eff` (Landau-Lifshitz
  without damping).

### 2.6 Plasma waves / MHD (S6)

* **Setting.** Conducting fluid with density `n`, velocity `v`,
  pressure `p`, magnetic field `B`.
* **Disturbed quantity.** A combined state `{n, v, B, p}`; different
  MHD modes disturb different subsets.
* **Local interaction.** Charge separation creates `E`; currents
  modify `B`; the `J x B` force couples back to momentum; magnetic
  tension `B . nabla B / mu_0` transmits force along `B`.
* **Local transfer.** EM force (`J x B`), pressure gradients,
  magnetic tension.
* **Restoring.** Three mechanisms: electrostatics (charge
  separation), pressure gradients, magnetic tension.
* **Resistance.** Ion (and electron) inertia; for slow MHD modes
  the dominant mass is `rho_i`.
* **Speed.** Alfven `v_A = B/sqrt(mu_0 rho)`; sound
  `c_s = sqrt(gamma p/rho)`; magnetosonic `c = sqrt(c_s^2 + v_A^2)`.
* **Steering.** Spatial gradients of `rho`, `B`, `p`.
* **Equations (minimum).** Continuity, momentum (`rho dv/dt =
  -nabla p + J x B`), induction (`dB/dt = curl(v x B)`), energy,
  with `J = curl B / mu_0`, `div B = 0`.

## 3. Comparison table

| System        | Disturbed quantity     | Local interaction                | Restoring            | Resistance / inertia | Speed mechanism                         | Steering                                  | Equations (min)                                  |
|---------------|------------------------|----------------------------------|----------------------|----------------------|------------------------------------------|-------------------------------------------|--------------------------------------------------|
| Water surface | `eta(x,t)` scalar      | pressure between adjacent columns | gravity `rho g eta`  | fluid mass density   | `c^2 = g tanh(kh)/k`                     | depth variation `h(x)`                    | `eta_tt = c^2 nabla^2 eta`                       |
| Elastic solid | `u(x,t)` vector        | traction `P_F N`                 | `-Div P_F`           | mass density `rho`   | `c^2 = A/rho` (acoustic tensor)          | spatial variation of `C_ijkl(x)`          | `rho u_tt = Div(A : sym grad u)`                 |
| Acoustic      | `delta p(x,t)` scalar  | pressure between adjacent volumes| bulk modulus `B`     | mass density `rho`   | `c^2 = B/rho`                            | spatial variation of `B, rho`             | `p_tt = c^2 nabla^2 p`                           |
| Maxwell EM    | `E, B` vector pair     | curl coupling between `E, B`     | NONE                 | NONE                 | `c^2 = 1/(eps mu)` from curl constants   | spatial variation of `eps, mu` / metric   | `curl E = -B_t`; `curl B = eps mu E_t`           |
| Spin wave     | `m(x,t)` unit vector   | Heisenberg exchange `-J S_i.S_j` | anisotropy `H_anis`  | gyromagnetic precession | `omega = gamma sqrt(H_eff(H_eff+D k^2))` | gradients of `H_anis, H_ext, J, D`        | `dm/dt = -gamma m x H_eff`                       |
| Plasma / MHD  | `{n,v,B,p}` combined   | `J x B`, pressure, magnetic tension | electrostatics / pressure / magnetic tension | ion mass `rho_i` | `v_A`, `c_s`, magnetosonic combinations | gradients of `rho, B, p`                  | continuity + momentum + induction + energy       |

## 4. Common abstract transport architecture

After laying the systems side by side, a common architecture is
present in **five of the six** systems (S1, S2, S3, S5, S6). The
slots are:

1. A local state defined at every point.
2. A neighbour-coupling mechanism between adjacent infinitesimal
   regions.
3. A locally transferred quantity between neighbours.
4. A restoring mechanism that pulls the state back.
5. A resistance to the response (inertia).
6. A propagation speed set by the ratio of restoring strength to
   resistance.
7. A local steering mechanism through spatial variation of medium
   parameters.

The wave equation emerges as a second-order-in-time PDE from the
balance of restoring and resistance. This is Family A.

**Family B (S4 Maxwell) is structurally different.** It shares slots
1, 2, 3, 6, 7 with Family A but **lacks slots 4 and 5 entirely**.
Its wave equation emerges from the mutual first-order curl coupling
of `E` and `B`, with no scalar potential energy and no mass term.
This makes S4 a different architectural family, not a member of
Family A.

The full common-architecture record is in
`common_architecture.csv` and `common_architecture.json`. The key
finding, stated without PBUF-specific interpretation, is recorded
under `key_structural_finding` in `common_architecture.json`.

## 5. PBUF comparison

The common architecture is now compared, slot by slot, with the
existing PBUF framework. Every conclusion cites the relevant
frozen artifact (FOUNDATION-001, V11, CORE-001) explicitly. The
detailed evidence is in `pbuf_architecture_comparison.csv` and
`pbuf_comparison.json`. The structural conclusion:

| Architecture slot        | PBUF status                                  |
|--------------------------|----------------------------------------------|
| local state              | present (CORE-001-E03 / E04: `q in R^3`, `u`) |
| neighbour coupling       | present (CORE-001-E01: `kappa_1 |q_j-q_i|^2`; LOCALITY-001: `Div P_F`) |
| local transfer           | present (continuum: `P_F N` traction; `-Div P_F`) |
| restoring mechanism      | present (CORE-001-E01: `kappa_0|q|^2`; continuum `-Div P_F`) |
| **resistance / inertia** | **MISSING** (INERTIA-001 closure gap)       |
| propagation speed        | not defined (no time structure)              |
| local steering           | available in principle (moduli `C_ijkl(x)`)  |
| wave equation            | not present (CORE-001-E09 is elliptic)       |

* **Family A match.** PBUF matches the spatial half of Family A.
  It does not match Family A on the temporal half because it has no
  kinetic sector. This is the same closure gap already identified
  by INERTIA-001 and re-derived from the EM side in EM-TRANSPORT-001.
* **Family B match.** Family B requires (a) two coupled vector
  fields, (b) a curl operator in the kinetic structure, (c) gauge
  invariance, and (d) no mass-like onsite term. CORE-001-E01 has
  none of (a), (b), (c) and has a positive `kappa_0|q|^2` mass-like
  term that contradicts (d). V11's `alpha_resolved ~ 3 alpha_EM` is
  a numerical identity and a dimensional-counting argument; it does
  not identify `q` with an EM vector potential. Family B is therefore
  structurally incompatible with PBUF without rewriting the
  microscopic energy.

## 6. Closest structural match

Of the six systems, the elastic-solid pattern (S2) is the closest
match:

* **Shared slots.** local state (`q` vs `u`), neighbour coupling
  (`kappa_1 |q_j - q_i|^2` vs stress `Div P_F`), restoring
  (`-Div P_F` from stored energy), steering (spatial variation of
  moduli).
* **Missing slot.** resistance / inertia.
* **What completion requires.** A kinetic sector supplying positive
  momentum density or an equivalent symplectic structure. With this
  slot filled, the speed would be `c^2 = G/K` from the frozen elastic
  2-jet.

**Why this is NOT a clean mapping under "minimal additional
assumptions".** The kinetic sector is precisely the closure gap
INERTIA-001 left open: it cannot be derived from `F`, from
`alpha_EM`, from `alpha_resolved`, or from `g_dev`. Adding it is a
new local physical principle, not a parameter renaming or a
coordinate choice. Per the brief, this falls outside the
"minimal additional assumptions" criterion, and the mapping is
therefore not clean.

## 7. Compliance with the milestone brief

| Constraint                                              | Status |
|---------------------------------------------------------|--------|
| No new physics                                          | yes    |
| No new constants                                        | yes    |
| No new transport equations                              | yes    |
| No assumption that Maxwell applies directly             | yes    |
| No assumption that elasticity applies directly          | yes    |
| No assumption that magnetism is the answer              | yes    |
| No modification of V11                                  | yes    |
| No modification of FOUNDATION-001                       | yes    |
| No modification of CORE-001                             | yes    |
| No constitutive-law change                              | yes    |
| No cosmological result introduced                       | yes    |
| No weak-lensing fitting                                 | yes    |
| No metric construction                                  | yes    |
| No quantum-gravity or dark-sector / dark-energy content  | yes    |
| Local-mechanism focus only                              | yes    |
| Every PBUF claim cited to a frozen artifact             | yes    |

## 8. Closure

**Outcome B.** The six required systems fall into two structural
families. PBUF matches Family A (mechanical) on the spatial half
only; the kinetic / inertial resistance slot is missing. Family B
(Maxwell) is structurally incompatible with CORE-001. The closest
match is the elastic-solid pattern (S2), whose completion requires
the kinetic sector already flagged by INERTIA-001. The decision is
recorded in `decision.json`; the completion record is in
`validation.json`.
