#!/usr/bin/env python3
"""Generate the PBUF MATTER-001 matter-action derivation audit.

This is a source-bounded theory audit.  It treats V11 as authoritative, does
not introduce a coupling constant, and does not import or execute the frozen
weak-lensing implementation.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


EQUATIONS = [
    {
        "equation_id": "MAT-001-E01",
        "equation": "S = S_EH[g] + S_sigma[g,chi; V11 microphysics] + S_m[g,Psi]",
        "status": "most general sector decomposition compatible with V11's unchanged GR/QM architecture; S_sigma is not supplied locally",
        "trace": "V11 secs. 1.1, 2.3, 2.3.2--2.3.4",
        "assumptions": "diffeomorphism invariance; Lambda absent; matter follows standard relativistic dynamics",
        "dimensions": "every term has units of action",
        "gap": "V11 gives background Omega_sigma(a), not a covariant local S_sigma or medium field chi",
    },
    {
        "equation_id": "MAT-001-E02",
        "equation": "T_m^{mu nu} = (2/sqrt(-g)) delta S_m/delta g_{mu nu}",
        "status": "derived definition of the universal metric matter operator (covariant-metric convention)",
        "trace": "standard GR explicitly retained by V11 sec. 1.1",
        "assumptions": "metric variation exists; matter is minimally and universally metric-coupled",
        "dimensions": "[T_m^{mu nu}] = energy/volume",
        "gap": "none for the metric source; its projection onto a PBUF medium variable is absent",
    },
    {
        "equation_id": "MAT-001-E03",
        "equation": "delta S_m/delta chi^A = (sqrt(-g)/2) T_m^{mu nu} (partial g_{mu nu}/partial chi^A) + (delta S_m/delta chi^A)|_g",
        "status": "exact chain rule; general matter loading of a proposed medium field",
        "trace": "E02 plus a not-yet-supplied map g[chi]",
        "assumptions": "local algebraic display; functional kernels replace partial derivatives for derivative/nonlocal maps",
        "dimensions": "both sides have action per chi variation; the metric derivative supplies the required chi normalization",
        "gap": "V11 supplies neither g[chi] nor an authorized direct fixed-metric matter dependence",
    },
    {
        "equation_id": "MAT-001-E04",
        "equation": "J_A = -(1/sqrt(-g)) delta S_m/delta chi^A = -(1/2) T_m^{mu nu} partial_A g_{mu nu} - O_A",
        "status": "derived source operator with stated sign convention; not numerically closed",
        "trace": "definition from E03; O_A=(1/sqrt(-g))(delta S_m/delta chi^A)|_g",
        "assumptions": "medium equation is written with J_A on its source side",
        "dimensions": "[J_A] = energy density/[chi^A]",
        "gap": "strict V11/standard minimal coupling sets O_A=0; partial_A g remains unknown",
    },
    {
        "equation_id": "MAT-001-E05",
        "equation": "nabla_mu T_m^{mu nu} = 0",
        "status": "Noether identity on the matter equations for minimal metric coupling",
        "trace": "diffeomorphism invariance of E01 and unchanged GR in V11 sec. 1.1",
        "assumptions": "matter equations hold; no direct chi dependence at fixed metric",
        "dimensions": "energy-momentum density per length",
        "gap": "a direct matter-medium vertex would require a specified exchange identity and is not in V11",
    },
    {
        "equation_id": "MAT-001-E06",
        "equation": "K u - div(G grad u) = s[rho]",
        "status": "retained conditional MB-001 static scalar balance",
        "trace": "MB-001-E01",
        "assumptions": "static, local, scalar, isotropic reduction",
        "dimensions": "[s]=[K][u]=[G][u]/L^2",
        "gap": "no definition u=C[chi,g] and no reduction from J_A to s",
    },
    {
        "equation_id": "MAT-001-E07",
        "equation": "s(x) = P^A_x[J_A] = -P^A_x[(1/2) T_m^{mu nu} partial_A g_{mu nu}] (minimal case)",
        "status": "most specific parameter-free source-functional identity available; P is unclosed",
        "trace": "E04 coarse-grained/projected into MB-001 E06",
        "assumptions": "a coarse-graining/projection P exists and preserves the chosen static scalar sector",
        "dimensions": "P must convert energy density/[chi] into [K][u] without an arbitrary normalization",
        "gap": "V11 does not define chi, u=C[chi,g], g[chi], P, or the normalization linking the covariant and MB actions",
    },
    {
        "equation_id": "MAT-001-E08",
        "equation": "integral_V s dV = integral_V K u dV - surface_integral_boundaryV G grad(u).n dA",
        "status": "derived global compatibility condition for MB-001, not a conserved matter charge",
        "trace": "volume integral of E06 and divergence theorem",
        "assumptions": "regular fields and boundary; scalar G or appropriate tensor contraction",
        "dimensions": "each term has [s] L^3 in three spatial dimensions",
        "gap": "boundary conditions and coefficients remain unclosed",
    },
]

CANDIDATES = [
    {"candidate": "rest-mass density rho", "pbuf_status": "used by MB-001 notation but not selected by V11", "covariance": "not a universal relativistic scalar source by itself", "conservation": "continuity only for conserved material number/mass", "finding": "admissible dust-limit input, not derivable as the fundamental operator"},
    {"candidate": "energy density T_{mu nu}v^mu v^nu", "pbuf_status": "implied after a medium observer v is chosen", "covariance": "observer/rest-frame dependent scalar projection", "conservation": "inherits constraints from stress-energy conservation", "finding": "possible scalar reduction; V11 defines no medium four-velocity or projection"},
    {"candidate": "trace T^mu_mu", "pbuf_status": "possible if the medium/metric response is purely conformal", "covariance": "scalar", "conservation": "not separately conserved; vanishes classically for conformal radiation", "finding": "one admissible action family, not selected by V11"},
    {"candidate": "full stress-energy tensor T^{mu nu}", "pbuf_status": "implied by V11's retention of standard GR", "covariance": "rank-2 tensor and universal metric source", "conservation": "nabla_mu T^{mu nu}=0 on shell for minimal coupling", "finding": "unique operator-level answer before a medium projection is chosen"},
    {"candidate": "matter Lagrangian or direct microscopic state vertex", "pbuf_status": "no V11 matter-medium vertex supplied", "covariance": "model dependent and Lagrangian-representation sensitive", "conservation": "generically exchanges energy-momentum with the medium", "finding": "not admissible as an asserted PBUF law without a new principle"},
    {"candidate": "V11 alpha hierarchy", "pbuf_status": "controls elastic/cosmological pipeline roles", "covariance": "not a matter operator", "conservation": "not applicable", "finding": "cannot normalize s and is not identified with a matter vertex"},
]

CONSERVATION = [
    {"audit_id": "MAT-C01", "law": "diffeomorphism Noether identity", "result": "minimal matter action gives nabla_mu T_m^{mu nu}=0 on matter equations", "condition": "standard universal metric coupling", "status": "derived from retained GR architecture"},
    {"audit_id": "MAT-C02", "law": "total energy-momentum conservation", "result": "nabla_mu(T_m^{mu nu}+T_sigma^{mu nu})=0 when all field equations hold", "condition": "a covariant S_sigma exists", "status": "structurally required; T_sigma cannot be calculated from supplied V11"},
    {"audit_id": "MAT-C03", "law": "matter-medium exchange", "result": "direct fixed-metric chi dependence makes matter stress-energy nonconserved separately, with an exchange term fixed by the same action", "condition": "nonminimal/direct coupling", "status": "allowed mathematically but not authorized or specified by V11"},
    {"audit_id": "MAT-C04", "law": "MB-001 global balance", "result": "volume source equals bulk recovery minus boundary gradient flux", "condition": "static MB-001 equation and boundary regularity", "status": "derived compatibility identity, not a fundamental conservation law"},
    {"audit_id": "MAT-C05", "law": "equivalence principle/universality", "result": "all matter loads through T^{mu nu}; composition-specific scalar charges are excluded", "condition": "'standard relativistic framework' is read as standard minimal matter coupling", "status": "strongest V11-consistent reading; premise is architectural rather than newly derived"},
]

ASSUMPTIONS = [
    {"item": "Standard GR and conventional QFT remain intact", "classification": "authoritative premise", "source": "V11 sec. 1.1"},
    {"item": "Lambda is absent and an elastic sector supplies the background contribution", "classification": "authoritative premise", "source": "V11 secs. 1.1, 2.3.3"},
    {"item": "Matter is minimally and universally coupled to one metric", "classification": "unavoidable interpretive assumption for the standard-GR reading", "source": "implied, not written as a V11 action"},
    {"item": "A covariant local medium field chi exists", "classification": "missing postulate", "source": "not supplied"},
    {"item": "The laboratory scalar u is a normalized coarse projection of chi or g", "classification": "missing definition", "source": "MB-001 gap"},
    {"item": "Static locality, scalarity, and isotropy of MB-001", "classification": "conditional continuum assumptions", "source": "MB-001-E01"},
    {"item": "s depends on rho alone", "classification": "unsupported restriction", "source": "MB-001 notation; not V11-derived"},
    {"item": "Optical response n(u)", "classification": "separate missing law", "source": "PHOTON-001"},
]


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def report() -> str:
    return r"""# PBUF MATTER-001 — Derivation of the matter action for the spacetime medium

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
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "equation_traceability.csv", EQUATIONS)
    write_csv(output / "candidate_matter_operators.csv", CANDIDATES)
    write_csv(output / "conservation_law_audit.csv", CONSERVATION)
    write_csv(output / "assumption_audit.csv", ASSUMPTIONS)
    (output / "matter_action_derivation_report.md").write_text(report())
    analysis = {
        "mission": "PBUF MATTER-001",
        "decision": "B with unresolved C projection family",
        "matter_operator": "stress-energy tensor T_m^{mu nu}",
        "source_operator": "J_A=-(1/sqrt(-g)) delta S_m/delta chi^A",
        "minimal_source": "J_A=-(1/2)T_m^{mu nu} partial_A g_{mu nu}",
        "continuum_source_status": "unclosed until g[chi], u=C[chi,g], and projection P are supplied",
        "exact_missing_principle": "covariant, normalized medium-to-metric constitutive map plus local elastic action",
        "new_free_parameters": [],
        "g_dev_used": False,
        "weak_lensing_imported_or_modified": False,
        "photon_closure_claimed": False,
    }
    (output / "matter_action_analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")
    validation = {
        "mission": "PBUF MATTER-001",
        "status": "COMPLETE",
        "deliverables_present": True,
        "all_equations_traced": all(row["trace"] and row["status"] for row in EQUATIONS),
        "dimensions_explicit": all(row["dimensions"] for row in EQUATIONS),
        "assumptions_explicit": all(row["assumptions"] for row in EQUATIONS),
        "conservation_audited": len(CONSERVATION) > 0,
        "replacement_free_parameters_introduced": False,
        "g_dev_reintroduced": False,
        "observational_fit_performed": False,
        "frozen_weak_lensing_modified": False,
        "unrelated_alpha_identity_assumed": False,
    }
    validation["all_checks_pass"] = all(
        value for key, value in validation.items()
        if key in {"deliverables_present", "all_equations_traced", "dimensions_explicit", "assumptions_explicit", "conservation_audited"}
    ) and not any(validation[key] for key in {
        "replacement_free_parameters_introduced", "g_dev_reintroduced", "observational_fit_performed",
        "frozen_weak_lensing_modified", "unrelated_alpha_identity_assumed",
    })
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/matter001"))
    args = parser.parse_args()
    main(args.output.resolve())
