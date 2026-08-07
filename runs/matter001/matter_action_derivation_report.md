# PBUF MATTER-001 — Derivation of the matter action for the spacetime medium

## Decision

**Outcome B, with an embedded Outcome C at the unresolved projection step.** V11 determines the operator carried by matter at the authoritative relativistic level: the full stress-energy tensor `T_m^{mu nu}`, obtained by varying the standard matter action with respect to the metric. It does not determine which scalar projection of that tensor loads the PBUF medium, because it supplies no local medium field, no map from that field to the metric, and no normalized coarse-graining map to MB-001's `u`.

Thus the matter action is partially derived without a new coupling, while the continuum source remains a family of admissible projections. The exact missing principle is a **covariant constitutive identification of the medium with geometry**, including `chi`, `g[chi]`, and `u=C[chi,g]`. Once that map is given, the source normalization follows by functional differentiation; it must not be separately assigned.

## 1. Authoritative boundary

V11 section 1.1 says that special relativity, general relativity, quantum mechanics, Einstein's equations, and standard quantum dynamics remain intact. The elastic medium is a constitutive interpretation. V11 then supplies a homogeneous background elastic density `Omega_sigma(a)` through its alpha/thermal pipeline, but it never displays a covariant local elastic action, a local deformation field, or a matter-medium vertex. Equations (16)--(17) normalize cosmological density parameters; they are not local source laws.

ARCH-001 removes the invented `g_dev` vertex. MB-001 retains only the unclosed balance `K u-div(G grad u)=s[rho]`. PHOTON-001 independently leaves the optical readout unclosed. None of the V11 alphas has the role of a matter charge.

## 2. Most general V11-consistent matter action

At sector level the admissible action is

`S = S_EH[g] + S_sigma[g,chi; V11 microphysics] + S_m[g,Psi]`.

`S_EH` and `S_m` are the standard generally covariant gravity and matter sectors; the cosmological constant is absent. `S_sigma` denotes the as-yet unknown local completion whose homogeneous stress-energy would have to reproduce V11's implemented `Omega_sigma(a)`. This notation adds no coefficient and makes no claim that V11 has already supplied `S_sigma`.

With the covariant-metric sign convention,

`T_m^{mu nu} = (2/sqrt(-g)) delta S_m/delta g_{mu nu}`.

If a proposed medium field `chi^A` determines the physical metric, the chain rule gives

`delta S_m/delta chi^A = (sqrt(-g)/2) T_m^{mu nu} partial_A g_{mu nu} + (delta S_m/delta chi^A)|_g`.

Therefore, defining the medium source on the right-hand side as `J_A=-(1/sqrt(-g))delta S_m/delta chi^A`,

`J_A = -(1/2) T_m^{mu nu} partial_A g_{mu nu} - O_A`,

where `O_A` is the fixed-metric direct matter operator. The strict standard-GR/minimal-coupling reading of V11 sets `O_A=0`. This is the most general parameter-free source formula available. Its normalization is fixed by the definitions and the action, not by a freely chosen multiplier.

Derivative or nonlocal maps `g[chi]` replace `partial_A g` by the corresponding functional kernel. They are mathematically admissible, but V11 does not select one.

## 3. What physical quantity loads the medium?

The unique answer before reduction is **stress-energy**, not rest-mass density alone. A scalar medium equation can receive only a scalar projection of `T_m^{mu nu}`. Which projection occurs is determined by `partial g_{mu nu}/partial chi^A`:

- a conformal response selects the trace `T^mu_mu`;
- a medium-rest-frame response may select energy density and spatial stress separately;
- a tensorial deformation retains corresponding components of `T^{mu nu}`.

For nonrelativistic pressureless matter, several projections collapse numerically to expressions proportional to `rho c^2`. That dust-limit degeneracy cannot establish `rho` as the fundamental operator. Radiation and relativistic pressure distinguish the choices: a trace source can vanish where an energy-density source does not.

## 4. Continuum source and exact gap

Let `P_x^A` denote the normalized reduction from the covariant medium equation to the static MB scalar. Then

`s(x) = P_x^A[J_A] = -P_x^A[(1/2)T_m^{mu nu} partial_A g_{mu nu}]`

for minimal coupling. This is an identity/schema, not a closed prediction. Writing `s(rho)` already assumes more than V11 provides: it discards pressure, momentum flux, anisotropic stress, composition, and possible derivative dependence.

To calculate `s`, PBUF must supply all of the following as one normalized law:

1. the local covariant medium variable `chi` and its dimensions;
2. the physical metric map `g[chi]` (or proof that `chi` is the metric/strain itself);
3. the coarse scalar definition `u=C[chi,g]`, including normalization;
4. the local elastic action `S_sigma`, whose quadratic/static limit derives `K` and `G`;
5. the projection and boundary prescription reducing `J_A` to `s`.

Supplying only `s proportional to rho`, a response coefficient, or a selected V11 alpha would not close this chain and would violate MATTER-001's no-free-parameter rule.

## 5. Multiple admissible actions and symmetry

V11's background equations cannot discriminate among conformal scalar, rest-frame scalar, tensor-strain, derivative, or nonlocal geometric responses. Full spacetime covariance and minimal coupling constrain all of them to use `T^{mu nu}` through metric variation, but do not select the metric-medium map. Locality, isotropy, parity, derivative order, and existence of a preferred medium four-velocity are additional assumptions, not V11 derivations.

If strict Lorentz covariance allows only a scalar conformal perturbation, the trace action is selected conditionally. If the medium has a rest frame, two independent isotropic scalar responses (temporal and spatial) are already possible. If `u` is merely a laboratory scalar projection of a tensor strain, still more responses are hidden. Adopting any of these now would create an unsupported constitutive law even if no symbol were attached to its normalization.

## 6. Conservation-law audit

For a diffeomorphism-invariant, minimally metric-coupled matter action, the matter equations imply `nabla_mu T_m^{mu nu}=0`. If a covariant elastic action exists, the full equations imply conservation of total matter-plus-medium stress-energy. A direct fixed-metric coupling would instead produce a precisely matched exchange term; because V11 supplies no such action, neither its form nor its conservation transfer can be asserted.

MB-001's static equation has the integral compatibility condition

`integral_V s dV = integral_V K u dV - surface_integral G grad(u).n dA`.

The `K u` term is a local recovery term, so MB-001 is not by itself a continuity equation for a conserved charge. Calling `s` a conserved mass source would therefore be unjustified.

## 7. Equation traceability

The complete term-level record is in `equation_traceability.csv`; candidate operators are compared in `candidate_matter_operators.csv`, assumptions in `assumption_audit.csv`, and conservation statements in `conservation_law_audit.csv`. In particular:

- `S_m` and `T^{mu nu}` trace to V11's explicit retention of standard GR/QFT.
- `S_sigma` is required for a local completion but is not given by V11's background `Omega_sigma(a)`.
- `chi`, `g[chi]`, `C`, and `P` are missing definitions, not derived objects.
- `K`, `G`, and `s` remain the conditional MB-001 quantities.
- the alpha hierarchy supplies elastic/cosmological inputs but no matter operator or source projection.

## 8. Recommendation

The next milestone should be **MATTER-002: Covariant medium-to-metric constitutive closure**. It should not run or modify weak lensing. Its acceptance gates should require:

1. one explicit, dimensionally normalized `S_sigma[g,chi]` whose homogeneous limit reproduces the V11 elastic background;
2. an explicit `g[chi]` and `u=C[chi,g]` with no adjustable matter coefficient;
3. derivation of `J_A` by varying the standard matter action;
4. a Noether identity demonstrating total conservation and the minimal-coupling limit;
5. a controlled static/weak-field reduction deriving `K`, `G`, and `s` rather than matching them;
6. discrimination among trace, energy-density, and stress-sensitive loading using theoretical consistency, not observational fitting.

PHOTON closure must remain separate: deriving how matter creates `u` does not derive how light reads it.

## Completion statement

MATTER-001 isolates the missing ingredient without introducing any parameter. Matter carries the stress-energy operator; its action on a medium field is the functional chain-rule projection generated by the metric map. Existing PBUF does not contain that map or the local elastic action, so no unique `s(rho)` can be derived. The milestone is complete under Outcome B/C: the operator-level interaction and conservation law are fixed, while the exact absent constitutive principle is explicitly identified.
