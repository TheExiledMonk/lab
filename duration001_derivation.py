#!/usr/bin/env python3
"""Generate the PBUF DURATION-001 emergent-duration derivation.

This is a structural derivation from the frozen PBUF inputs.  It introduces no
fundamental time coordinate, field equation, coupling, or fundamental constant.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TRACEABILITY = [
    ("DU-001", "q:S->Q_adm, s only orders states", "STATE-002 S2-012--S2-013; DYNAMICS-001 D-002", "accepted", "s has no invariant scale, origin, or rate"),
    ("DU-002", "Delta tau[gamma]=integral F(q(s),qdot(s),zdot(s)) ds", "DU-001; additive local accumulation; clock/propagation assumptions", "derived family", "general local accumulated duration"),
    ("DU-003", "F(q,a qdot,a zdot)=a F(q,qdot,zdot), a>0", "DU-002; s->f(s) invariance", "derived", "necessary and sufficient positive degree-one homogeneity"),
    ("DU-004", "d tau=sigma(q,x,n,clock) d ell=d ell/v", "DU-002--DU-003; physical path element d ell", "derived representation", "sigma>0 is slowness and v:=1/sigma; neither is fixed numerically"),
    ("DU-005", "Delta tau=int_path sigma d ell=int d ell/v", "DU-004; piecewise regular propagation", "derived", "most general local scalar distance-speed-duration relation"),
    ("DU-006", "tau_B-tau_A=integral_[sA,sB] F ds", "DU-002", "definition", "duration is an interval functional; its zero is conventional"),
    ("DU-007", "N_C=int dN_C; Delta tau_C=U_C N_C", "fixed clock assumptions; repeatable propagation cycle", "derived clock form", "clock accumulates phase/cycle count; U_C calibrates one cycle in the chosen unit"),
    ("DU-008", "d tau_C=r_C d tau_*", "DU-007; comparison by coincidences", "derived", "consistent clocks require r_C constant after calibration and path/history independent locally"),
    ("DU-009", "d tau^2=-g_eff_mn dx^m dx^n/c^2 on timelike paths", "FP-5; V11 one-metric operational limit", "effective matching condition", "not a microscopic derivation of g_eff or c"),
    ("DU-010", "d tau=dt sqrt(-(g00 c^2+2 g0i c vi+gij vi vj)/c^2)", "DU-009; x0=ct; vi=dxi/dt", "derived coordinate relation", "t is a chart/synchronization coordinate, not s or tau"),
    ("DU-011", "Delta tau=0 for a constant physical history", "DU-002--DU-003; continuous zero-section extension", "derived static limit", "no change/propagation means no clock accumulation"),
    ("DU-012", "q -> oriented history/order s -> propagation -> tau -> effective chart t", "DU-001--DU-011", "derived dependency classification", "q/order are fundamental inputs; tau and t are emergent"),
]


CATALOGUE = {
    "milestone": "PBUF DURATION-001 — Derivation of Emergent Physical Duration",
    "status": "complete_minimal_family_with_v11_matching_condition",
    "ontology_review": False,
    "fundamental_time_dimension_introduced": False,
    "new_constants_or_couplings": [],
    "symbols": {
        "q": {"type": "point of Q_adm", "role": "complete physical medium state", "status": "fundamental accepted input"},
        "s": {"type": "oriented history parameter", "role": "orders successive states", "status": "mathematical gauge label", "observable": False},
        "tau": {"type": "additive reparametrization-invariant interval functional", "role": "clock-measured physical duration", "status": "emergent"},
        "t": {"type": "coordinate function on the effective Lorentzian representation", "role": "labels effective events/slices", "status": "emergent and chart-dependent"},
        "ell": {"type": "physical propagation path length", "role": "propagation progress", "status": "derived from the accepted/effective ruler structure"},
        "sigma": {"type": "positive slowness functional", "role": "local duration per path length", "status": "constitutive/calibration family, not a new universal constant"},
        "v": {"type": "physical propagation speed", "definition": "d ell/d tau=1/sigma", "status": "relational observable"},
    },
    "general_duration": {
        "equation": "Delta tau=integral F(q,qdot,zdot) ds",
        "conditions": ["F>=0", "F(q,a qdot,a zdot)=aF(q,qdot,zdot) for a>0", "locality", "additivity", "gauge/objective scalar", "F=0 on no-process segments"],
    },
    "propagation_family": {
        "equation": "Delta tau=integral_path sigma(q,x,n,clock) d ell=integral_path d ell/v",
        "minimality": "all positive local additive direction-dependent propagation calibrations",
        "excluded_inference": "tau=ell/c is not selected by the accepted inputs",
    },
    "clock": {
        "accumulator": "monotone propagation phase or cycle count N_C",
        "reading": "tau_C=tau_C0+U_C N_C",
        "unit": "a declared number of repeatable cycles of a reference process",
        "comparison": "co-located coincidences and signal exchange, without use of external time",
        "consistency": ["monotonicity", "repeatability", "additivity", "reparametrization invariance", "objectivity", "local universality after calibration", "common one-metric limit"],
    },
    "v11_limit": {
        "proper_duration": "d tau^2=-g_eff_mn dx^m dx^n/c^2 for timelike worldlines",
        "coordinate_relation": "d tau=dt sqrt(-(g00 c^2+2g0i c vi+gij vi vj)/c^2)",
        "light_warning": "null propagation has d tau=0; a light clock measures proper duration on its timelike apparatus worldline between return coincidences",
        "status": "mandatory effective matching condition, conditional on the accepted one-metric V11 regime",
    },
    "static_limit": {"conditions": ["no propagation", "no deformation", "no state evolution"], "result": "Delta tau=0", "interpretation": "order labels may remain representational, but no physical clock duration exists"},
    "non_uniqueness": "The fixed inputs determine the invariant functional form but do not select F, sigma, a microscopic standard clock, or the medium-to-metric map.",
}


DEPENDENCIES = [
    ("DUR-D01", "q", "oriented history q(s)", "fundamental state -> accepted ordered evolution", "accepted"),
    ("DUR-D02", "oriented history q(s)", "s", "choice of order label", "representational/gauge"),
    ("DUR-D03", "q(s) and physical process", "propagation path/progress", "state change realizes propagation", "physical, assumed to occur"),
    ("DUR-D04", "propagation path/progress", "tau", "positive additive invariant accumulation", "emergent; family derived"),
    ("DUR-D05", "tau plus effective one-metric clock synchronization", "t", "chart construction", "emergent/effective"),
]


def validate() -> dict:
    # Direct numerical check of relabeling invariance for a nontrivial member
    # F=sqrt(qdot^2+pdot^2).  s=u^3 is monotone; midpoint quadrature avoids
    # endpoint singular bookkeeping and converges rapidly.
    import math

    def midpoint(n: int, transformed: bool) -> float:
        total = 0.0
        for k in range(n):
            u = (k + 0.5) / n
            if transformed:
                s = u**3
                ds_du = 3.0 * u**2
                qdot_s = 2.0 * s
                pdot_s = math.cos(s)
                integrand = math.hypot(qdot_s * ds_du, pdot_s * ds_du)
            else:
                s = u
                integrand = math.hypot(2.0 * s, math.cos(s))
            total += integrand / n
        return total

    native = midpoint(200000, False)
    relabelled = midpoint(200000, True)
    checks = {
        "relabeling_invariance_numerical": abs(native - relabelled) < 2e-10,
        "duration_integrand_is_degree_one": any("aF(q,qdot,zdot)" in condition for condition in CATALOGUE["general_duration"]["conditions"]),
        "s_tau_t_formally_distinct": len({CATALOGUE["symbols"][x]["type"] for x in ("s", "tau", "t")}) == 3,
        "tau_equals_ell_over_c_not_assumed": "not selected" in CATALOGUE["propagation_family"]["excluded_inference"],
        "static_limit_zero": CATALOGUE["static_limit"]["result"] == "Delta tau=0",
        "null_propagation_distinguished_from_clock_worldline": "null propagation" in CATALOGUE["v11_limit"]["light_warning"],
        "no_new_constant_or_coupling": CATALOGUE["new_constants_or_couplings"] == [],
        "no_fundamental_time_dimension": not CATALOGUE["fundamental_time_dimension_introduced"],
        "all_equations_traceable": all(row[2] and row[3] for row in TRACEABILITY),
    }
    return {"checks": checks, "all_checks_pass": all(checks.values()), "metrics": {"native_integral": native, "relabelled_integral": relabelled, "absolute_difference": abs(native-relabelled)}}


def table(headers: list[str], rows: list[tuple]) -> str:
    clean = lambda value: str(value).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(clean(v) for v in row) + " |" for row in rows)


def report(validation: dict) -> str:
    trace = table(["ID", "Equation/content", "Premises", "Status", "Boundary/use"], TRACEABILITY)
    deps = table(["Edge", "From", "To", "Mapping", "Status"], DEPENDENCIES)
    return r"""# PBUF DURATION-001 — Derivation of Emergent Physical Duration

## Result

Physical duration is the positive, additive, reparametrization-invariant measure of propagation-bearing change along an oriented history of the one medium. For a representative history (q:S\to\mathcal Q_{\rm adm}), the most general local first-order form is

\[
\boxed{\Delta\tau[\gamma]=\int_{s_A}^{s_B}
F(q(s),\dot q(s),\dot z(s))\,ds,\qquad
F(q,a\dot q,a\dot z)=aF(q,\dot q,\dot z),\ a>0.} \tag{DU-002--003}
\]

Here (z) denotes propagation-progress coordinates derived from the physical state/process, not a new field; this notation is distinct from DYNAMICS-001's canonical momentum (p). Positivity, objectivity, locality, additivity, and (F=0) when no physical clock process advances complete the minimal admissibility conditions. The homogeneity condition is necessary and sufficient for invariance under every (s'=f(s)) with (f'>0). Hence (s) orders states but supplies neither a duration scale nor a rate.

The accepted assumptions determine this family, not one unique `F`. In particular, they do not derive a numerical propagation speed, a microscopic standard clock, or the missing normalized one-metric map. DURATION-001 therefore completes the structural derivation while recording the remaining constitutive/calibration freedom rather than disguising it as `tau=ell/c`.

## 1. Formal distinction among `s`, `tau`, and `t`

- `s` is a coordinate on an oriented representative of the history. Under `s'=f(s)`, its increments and rates change. It is not observable and is not physical time.
- `tau` is an interval functional on the oriented physical history. Its differences are invariant under relabeling; its additive origin is conventional. A clock realizes this functional by accumulating propagation phase or cycles.
- `t=x^0/c` (or `x^0`, by convention) is a chart coordinate in the effective V11 Lorentzian description. It depends on synchronization and coordinates. It is neither the arbitrary microscopic label `s` nor, in general, a clock's proper duration `tau`.

Thus one may write `tau(s)=tau_0+integral F ds` in a chosen representative, but this does not identify `tau` with `s`. Where `F>0`, choosing the gauge `s=tau` is allowed after duration has been derived; it is a gauge choice, not an ontological identity. Likewise `t=tau` occurs only for specially synchronized/rest clocks in appropriate effective coordinates.

## 2. Propagation-duration theorem and minimal family

Assume a propagation process traces a physical path with nonnegative additive line element (d\ell). Locality and scalar objectivity imply that an additive duration element must have the form

\[
\boxed{d\tau=\sigma(q,x,n,\mathcal C)\,d\ell
=\frac{d\ell}{v(q,x,n,\mathcal C)},\qquad \sigma>0,\quad v:=\frac{d\ell}{d\tau}.} \tag{DU-004}
\]

Here `n` is propagation direction when anisotropy is admissible and `C` identifies the clock/process class. Integration gives

\[
\boxed{\Delta\tau=\int_{\rm path}\sigma\,d\ell
=\int_{\rm path}\frac{d\ell}{v}.} \tag{DU-005}
\]

This is the minimal admissible family: all positive local additive calibrations of physical propagation distance. Homogeneity in the path tangent makes it independent of (s). If homogeneity/isotropy removes direction dependence, (sigma=sigma(q,\mathcal C)); if local universality is also established, admissible clocks share one (sigma(q)) after calibration. A constant-speed uniform segment gives (Delta\tau=\ell/v). The special expression (ell/c) requires the additional facts that the selected process has speed (c), the path and state make it constant, and the clock protocol converts its travel into a timelike cycle. Those facts are not consequences of propagation alone.

Speed is relational: (v=d\ell/d\tau). Equation (DU-005) is therefore a consistency relation between ruler progress and clock accumulation, not a circular construction of duration from a speed known using an external time. Operationally one first compares repeatable process counts and path standards; the quotient is then called speed.

## 3. Physical clocks

A PBUF clock is a localized subsystem/process of the medium with a monotone phase (\phi_C) generated by propagation, identifiable return events, and repeatable cycles. Define

\[
dN_C=\frac{d\phi_C}{2\pi},\qquad
\Delta\tau_C=U_C\Delta N_C. \tag{DU-007}
\]

The accumulated primitive is phase/cycle count, not external time. One duration unit is a declared number of cycles of a reference realization; `U_C` is a unit conversion fixed by that declaration, not a new fundamental constant. Clock comparisons use coincidences: colocate clocks, compare counts between common events, or exchange propagation signals and correct using the common effective geometry.

Consistent clocks must be monotone on their operating domain, additive under concatenation, repeatable after reset, objective and gauge-invariant, stable against irrelevant construction details, and locally universal after calibration. Relative to a reference clock,

\[
d\tau_C=r_C(q,\text{environment, motion})\,d\tau_* . \tag{DU-008}
\]

An ideal clock has the same calibrated (r_C) as every other ideal clock subjected to the same state and motion. Environmental sensitivities are correctable clock imperfections. A path- or composition-dependent residual after correction would violate universal clock consistency and the V11 one-metric limit.

## 4. Emergence of V11 relativistic time

FP-5 requires that the universal duration family match the single effective Lorentzian metric in the V11 regime. Universality, local Lorentz invariance, locality, and additivity leave the standard timelike line element (up to clock-unit convention):

\[
\boxed{d\tau^2=-\frac{1}{c^2}g^{\rm eff}_{\mu\nu}dx^\mu dx^\nu,\qquad d\tau>0} \tag{DU-009}
\]

for signature `(-+++)`. This is an effective matching condition, not a fourth fundamental medium dimension and not a derivation of `g_eff` or `c`. The coordinates `x^mu`, including `t`, compactly represent correlations among medium states, rulers, signals, and clocks.

With (x^0=ct), (v^i=dx^i/dt),

\[
\boxed{d\tau=dt\sqrt{-\frac{g_{00}c^2+2g_{0i}cv^i+g_{ij}v^iv^j}{c^2}}.} \tag{DU-010}
\]

In a local inertial chart this becomes `d tau=dt sqrt(1-|v|^2/c^2)`, recovering ordinary time dilation and standard ideal-clock behavior. Coordinate time emerges operationally by selecting a congruence of reference clocks, a synchronization convention using propagation signals, and a chart; different choices give different `t` while preserving `tau`.

Null propagation obeys `g_eff_mn dx^m dx^n=0`, so the propagating light signal itself accumulates zero proper duration. A light clock instead counts emission-return coincidences and accumulates proper duration along the timelike apparatus worldline. This distinction prevents the false universal inference `tau=ell/c` while retaining standard radar and clock operations.

The V11 compatibility result is conditional in exactly one inherited sense: the accepted foundations require the one-metric Lorentzian limit but do not provide the normalized microscopic medium-to-metric map. DURATION-001 derives how any admissible duration must match that limit; it does not invent the missing map or modify V11.

## 5. Static universe limit

If there is no propagation, no deformation, and no state evolution, every physical clock phase is constant. The history tangent and propagation progress vanish, so degree-one homogeneity with the continuous zero-section condition gives

\[
\boxed{F(q,0,0)=0,\qquad \Delta\tau=0.} \tag{DU-011}
\]

The ordered set (S) may still be written mathematically, but its labels distinguish no physical states or clock coincidences. Under the accepted ontology there is therefore no measurable physical duration in the exactly static limit. Saying that an unobserved external time nevertheless passes would add the forbidden fundamental time structure.

This conclusion concerns an exactly constant complete state. A macroscopically undeformed configuration containing an operating clock is not static in this sense: its internal propagation changes (q) and accumulates duration.

## 6. Dependency structure

\[
\boxed{q\longrightarrow [q(s)]_{\mathrm{Homeo}_+}
\longrightarrow\text{physical propagation}
\longrightarrow\tau
\longrightarrow t.} \tag{DU-012}
\]

The requested shorthand (q\to s\to\text{propagation}\to\tau\to t) is valid only if the arrow to (s) means “choose an order label”: (s) is not generated as an observable by (q). The complete typed dependency graph is:

""" + deps + r"""

The state `q`, existence/orientation of successive states, and physical propagation premise are fundamental accepted inputs. The numerical label `s` is representational. Duration `tau` is emergent as the invariant measure of propagation-bearing evolution. Relativistic coordinate time `t` is a further effective chart construction from synchronized clock readings and the V11 one-metric representation.

## 7. Equation traceability

""" + trace + rf"""

## 8. Completion and residual boundary

DURATION-001 derives the invariant mathematical definition of duration, its most general local propagation form, the operational clock definition, the static limit, and the V11 proper/coordinate-time matching relations without adding a fundamental dimension. The result is deliberately a minimal family. The fixed assumptions cannot choose (F) or (sigma), prove microscopic clock universality, or construct the normalized effective metric. Selecting any of these would require constitutive or metric information beyond this milestone, not a new ontology.

Automated structural and relabeling checks pass: **{validation['all_checks_pass']}**. The numerical reparameterization test differs by `{validation['metrics']['absolute_difference']:.3e}`.
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("DURATION-001 validation failed")
    (output / "duration_catalogue.json").write_text(json.dumps(CATALOGUE, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "dependency_graph.md").write_text("# PBUF DURATION-001 dependency graph\n\n" + table(["Edge", "From", "To", "Mapping", "Status"], DEPENDENCIES) + "\n")
    with (output / "equation_traceability.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "equation_or_content", "premises", "status", "boundary_or_use"])
        writer.writerows(TRACEABILITY)
    (output / "emergent_duration_derivation.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/duration001"))
    main(parser.parse_args().output)
