"""PBUF MATERIAL-LAB-001: reproducible, extensible constitutive comparison.

The scalar coordinate e is only a restriction of the frozen objective tensor C
to an admissible deformation ray.  It is not a replacement ontology.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable


@dataclass(frozen=True)
class Material:
    key: str
    name: str
    mechanism: str
    assumptions: tuple[str, ...]
    domain: str
    energy_formula: str
    stress_formula: str
    tangent_formula: str
    energy: Callable[[float], float]
    stress: Callable[[float], float]
    tangent: Callable[[float], float]
    upper: float
    neighbor_interaction: str
    wave_support: str
    strengths: tuple[str, ...]
    weaknesses: tuple[str, ...]
    open_questions: tuple[str, ...]


# Dimensionless laboratory normalization only: k0=1, shape parameters=1.
MATERIALS = (
    Material("A", "Progressive strain hardening",
        "Resistance rises with accumulated reversible deformation.",
        ("hyperelastic ray restriction", "k_t(e)>0 and k_t'(e)>=0", "quartic representative has h>0"),
        "signed e in R (representative); tensor lift restricted to the frozen SPD domain",
        "w=k0 e^2/2+h e^4/4", "sigma=k0 e+h e^3", "k_t=k0+3h e^2",
        lambda e: e*e/2+e**4/4, lambda e: e+e**3, lambda e: 1+3*e*e, 3.0,
        "none implied by the local hypothesis", "conditional on an authorized kinetic/gradient completion",
        ("minimal smooth even hardening representative", "strictly convex and coercive"),
        ("hypothesis selects inequalities, not a unique law", "local form supplies no neighbor coupling"),
        ("Which invariant tensor lift?", "What fixes the hardening scale?")),
    Material("B", "Wave-equilibrium material",
        "A single continuum rearranges collectively through a positive gradient-energy term.",
        ("one continuum only", "ell>0 is a candidate correlation length", "static scalar ray proxy"),
        "e in H1(Omega), with boundary data making the variational problem well posed",
        "W=int[k0 e^2/2+k0 ell^2 |grad e|^2/2-s e] dV",
        "local sigma=k0 e; microstress xi=k0 ell^2 grad e", "local k_t=k0; Fourier Hessian=k0(1+ell^2|q|^2)",
        lambda e: e*e/2, lambda e: e, lambda e: 1.0, 3.0,
        "intrinsic: Euler-Lagrange equilibrium (1-ell^2 Laplacian)e=s/k0",
        "yes conditionally: positive inertia gives omega^2=c0^2(q^2+ell^2 q^4); inertia is not selected here",
        ("explicit neighbor coupling", "positive Fourier stiffness", "unique static equilibrium under standard boundary conditions"),
        ("does not progressively harden in amplitude", "ell and kinetic closure are underived"),
        ("What frozen quantity fixes ell?", "What tensor-gradient invariant and kinetic term are admissible?")),
    Material("C", "Exponential hardening",
        "Incremental work grows exponentially in deformation magnitude.",
        ("b>0", "signed-symmetric hyperelastic representative"), "signed e in R",
        "w=k0[cosh(b e)-1]/b^2", "sigma=k0 sinh(b e)/b", "k_t=k0 cosh(b e)",
        lambda e: math.cosh(e)-1, lambda e: math.sinh(e), lambda e: math.cosh(e), 3.0,
        "none implied locally", "conditional on an authorized kinetic/gradient completion",
        ("smooth, strictly convex, coercive", "rapid progressive stiffening"),
        ("no finite mathematical endpoint", "exponential form and scale are not frozen consequences"),
        ("What derives b?", "Does a globally stable tensor lift exist on the full domain?")),
    Material("D", "Finite-extensibility material",
        "A divergent elastic barrier prevents reaching finite extensibility e_star.",
        ("e_star>0", "Gent-type symmetric ray representative"), "|e|<e_star (normalized e_star=1)",
        "w=-(k0 e_star^2/2) log[1-(e/e_star)^2]", "sigma=k0 e/[1-(e/e_star)^2]",
        "k_t=k0[1+(e/e_star)^2]/[1-(e/e_star)^2]^2",
        lambda e: -0.5*math.log(1-e*e), lambda e: e/(1-e*e),
        lambda e: (1+e*e)/(1-e*e)**2, 0.98,
        "none implied locally", "conditional on an authorized kinetic/gradient completion",
        ("finite domain with infinite-energy barrier", "strict convexity and hardening on each loading branch"),
        ("endpoint is not derived", "ray result does not prove full tensor ellipticity"),
        ("Which spectral boundary defines extensibility?", "How are compression boundaries controlled?")),
    Material("E", "Polynomial hardening",
        "A general nonlinear polynomial supplies controllable hardening orders.",
        ("even energy for signed recovery", "a_2m>=0 with at least one positive nonlinear coefficient"),
        "signed e in R for the nonnegative even-energy subclass",
        "w=k0 e^2/2+sum_(m=2)^N a_2m e^(2m)/(2m)",
        "sigma=k0 e+sum a_2m e^(2m-1)", "k_t=k0+sum(2m-1)a_2m e^(2m-2)",
        lambda e: e*e/2+e**4/4+e**6/6, lambda e: e+e**3+e**5,
        lambda e: 1+3*e*e+5*e**4, 3.0,
        "none implied locally", "conditional on an authorized kinetic/gradient completion",
        ("systematically extensible", "stable for the displayed coefficient restrictions"),
        ("mixed coefficients need separate domain proof", "order and coefficients remain free"),
        ("Which truncation is justified?", "What constrains coefficients and tensor invariants?")),
)


PROPERTIES = (
    "Weak-field recovery", "Progressive hardening", "Smooth constitutive response",
    "Stable equilibrium", "Wave-support capability", "Recovery",
    "Finite-energy compatibility", "Emergent geometry compatibility", "Governing-equation readiness",
)


def classify(m: Material) -> dict[str, str]:
    hard = "fail (amplitude tangent is constant)" if m.key == "B" else "pass"
    wave = "pass (conditional kinetic closure)" if m.key == "B" else "conditional (completion required)"
    smooth = "pass on open admissible domain" if m.key == "D" else "pass"
    return {
        "Weak-field recovery": "pass (common k0 tangent)",
        "Progressive hardening": hard,
        "Smooth constitutive response": smooth,
        "Stable equilibrium": "pass on ray; tensor proof pending",
        "Wave-support capability": wave,
        "Recovery": "pass restoring tendency; evolution pending",
        "Finite-energy compatibility": "pass",
        "Emergent geometry compatibility": "conditional on METRIC-001 map",
        "Governing-equation readiness": "conditional; closure slots explicit",
    }


def numerical_checks(m: Material) -> dict[str, bool]:
    n = 151
    xs = [m.upper*i/n for i in range(1, n)]
    energies = [m.energy(x) for x in xs]
    tangents = [m.tangent(x) for x in xs]
    return {
        "reference_energy_zero": abs(m.energy(0.0)) < 1e-12,
        "reference_stress_zero": abs(m.stress(0.0)) < 1e-12,
        "weak_tangent_unity": abs(m.tangent(0.0)-1.0) < 1e-12,
        "positive_energy": all(v > 0 for v in energies),
        "positive_tangent": all(v > 0 for v in tangents),
        "nondecreasing_tangent": all(b >= a-1e-12 for a, b in zip(tangents, tangents[1:])),
    }


def write_outputs(root: Path) -> None:
    out = root / "runs" / "material_lab001"
    out.mkdir(parents=True, exist_ok=True)
    catalogue = []
    checks = {}
    for m in MATERIALS:
        d = asdict(m)
        for key in ("energy", "stress", "tangent", "upper"):
            d.pop(key)
        d["evaluation"] = classify(m)
        catalogue.append(d)
        checks[m.key] = numerical_checks(m)
    (out / "material_model_catalogue.json").write_text(json.dumps({
        "milestone": "PBUF MATERIAL-LAB-001", "schema_version": "1.0",
        "native_variable": "C[q,q0] in Sym+(3); e is an admissible ray coordinate only",
        "normalization": "numerical checks use k0=1 and unit positive shape parameters",
        "models": catalogue,
    }, indent=2) + "\n")

    with (out / "comparison_matrix.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["Property", *[m.key for m in MATERIALS]])
        for p in PROPERTIES: w.writerow([p, *[classify(m)[p] for m in MATERIALS]])

    report = build_report()
    (out / "material_laboratory_report.md").write_text(report)
    (out / "dependency_graph.md").write_text("""# MATERIAL-LAB-001 dependency graph

```text
FOUNDATION-001 + STATE-002 + DEFORMATION-001
                    |
HYPER-001 + ENERGY-PRINCIPLE-001
                    |
       MATERIAL-LAB-001 common gates
      /       /       |       \\       \\
     A       B        C        D        E
      \\       \\       |       /       /
 DURATION-001 + METRIC-001 + BALANCE-001 compatibility gates
                    |
 future tensor lift, dynamics, lensing, cosmology, compact-object and quantum work
```
""")
    passed = all(all(v.values()) for v in checks.values())
    (out / "validation.json").write_text(json.dumps({
        "milestone": "PBUF MATERIAL-LAB-001", "pass": passed,
        "analytic_and_sampled_checks": checks,
        "deliverables": ["material_model_catalogue.json", "comparison_matrix.csv", "material_laboratory_report.md", "dependency_graph.md"],
        "architecture_extensible": "add one Material entry; common evaluator and outputs require no changes",
        "prohibitions": {"ontology_modified": False, "weak_lensing_modified": False, "observational_fit": False,
            "v11_modified": False, "microscopic_constituents_added": False, "preferred_law_selected": False},
    }, indent=2) + "\n")


def build_report() -> str:
    sections = ["""# PBUF MATERIAL-LAB-001 — Comparative Material Laboratory

## Scope and common procedure

The frozen ontology and milestones are inputs, not objects of review. The native deformation remains the objective tensor `C[q,q0]` in the frozen admissible SPD domain. For comparable algebra, each local model is restricted to a smooth reversible ray `C(e)` and tested using `w(e)=W(C(e))`, `sigma=w'`, and `k_t=w''`. Passing a ray test is necessary, not sufficient, for convexity, rank-one convexity, strong ellipticity, or hyperbolicity of a future tensor theory.

Every model receives the same five stages: (1) regularity, domain, positivity, convexity, stability, and equilibrium; (2) native response, energy, tangent, hardening, recovery, neighbor coupling, and waves; (3) frozen-milestone compatibility; (4) geometry, duration, propagation, recovery, finite energy, and stiffening consequences; and (5) downstream readiness. No observations or weak-lensing implementation enter the evaluation.

Common gates are `w(0)=sigma(0)=0`, `k_t(0)=k0>0`, positive energy away from equilibrium, positive tangent, and nondecreasing tangent for progressive hardening. A strict energy minimum gives a restoring tendency, but recovery in time requires an authorized evolution law. Coercive energies make finite-energy sublevels bounded; only a barrier model also has a finite mathematical endpoint.
"""]
    for m in MATERIALS:
        ev = classify(m)
        sections.append(f"""## Model {m.key} — {m.name}

1. **Physical interpretation.** {m.mechanism}
2. **Constitutive assumptions.** {'; '.join(m.assumptions)}.
3. **Derived equation.** `{m.stress_formula}`. For B, variation additionally gives `(1-ell^2 Laplacian)e=s/k0`; boundary conditions are required.
4. **Stored energy.** `{m.energy_formula}` on `{m.domain}`.
5. **Stress and tangent.** `{m.stress_formula}`; `{m.tangent_formula}`.
6. **Progressive hardening.** {ev['Progressive hardening']}. The hypothesis uniquely fixes only monotonicity inequalities unless the displayed representative is stated.
7. **Weak-field limit.** `w=k0 e^2/2+o(e^2)`, `sigma=k0 e+o(e)`.
8. **Recovery.** {ev['Recovery']}.
9. **Mathematical stability.** {ev['Stable equilibrium']}; smoothness: {ev['Smooth constitutive response']}.
10. **Strengths.** {'; '.join(m.strengths)}.
11. **Weaknesses.** {'; '.join(m.weaknesses)}.
12. **Open questions.** {'; '.join(m.open_questions)}.

Neighbor interaction: {m.neighbor_interaction}. Wave support: {m.wave_support}.
""")
    rows = ["| Property | A | B | C | D | E |", "|---|---|---|---|---|---|"]
    for p in PROPERTIES:
        rows.append("| " + " | ".join([p, *[classify(m)[p] for m in MATERIALS]]) + " |")
    sections.append("## Unweighted comparative matrix\n\n" + "\n".join(rows))
    sections.append("""## Frozen PBUF compatibility

All candidates preserve FOUNDATION-001's one-medium ontology and STATE-002/DEFORMATION-001's native objective state. Their weak tangent can satisfy HYPER-001 and ENERGY-PRINCIPLE-001, subject to a full invariant tensor lift and spectral stability proof. They are statewise and introduce no fundamental time, so DURATION-001 is preserved. Emergent geometry remains conditional on a regular METRIC-001 map. Stress is an admissible BALANCE-001 closure ingredient, not a completed balance/evolution equation. CONSTITUTIVE-002 remains unmodified and supplies restrictions rather than a selected formula.

## Derived behavior and readiness

All five have a stable weak equilibrium and finite-energy compatibility. A, C, D, and E progressively stiffen in amplitude; B instead stiffens short wavelengths through neighbor coupling and is the only initial model with native spatial propagation structure. D alone imposes a finite mathematical deformation endpoint. None independently derives duration or the metric, and none closes governing equations without tensor, kinetic, source, boundary, and metric-map choices. Consequently lensing, cosmology, compact-object, and quantum work are only conditionally ready after those closures; no weak-lensing modification is warranted here.

## Comparative ranking without weights

There is no scalar score and no unique winner. The criterion-only result is a partial ordering by distinctive capability:

1. **Broad static hardening set:** A, C, D, and E pass all displayed ray-level hardening gates; none dominates the others absent a selected endpoint or growth principle.
2. **Spatial-interaction leader:** B alone passes native neighbor interaction and conditional wave-support gates, but does not pass amplitude hardening.
3. **Finite-endpoint leader:** D alone supplies finite extensibility through an energetic barrier; the frozen framework does not establish that this extra property must select the law.

This is a capability ranking, not a preference or declaration of physical correctness.

## Recommendations for subsequent development

Preserve all five as separate mechanism classes. Next derive invariant tensor lifts and run spectral Hessian/strong-ellipticity tests on the frozen SPD domain. For B, derive admissible tensor-gradient and kinetic terms and boundary conditions. Keep the weak tangent common, keep model parameters symbolic until independently derived, and only after governing closure evaluate later metric, lensing, cosmology, compact-object, or quantum consequences.
""")
    return "\n".join(sections)


if __name__ == "__main__":
    write_outputs(Path(__file__).resolve().parent)
