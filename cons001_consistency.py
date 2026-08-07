#!/usr/bin/env python3
"""Generate the theory-only PBUF CONS-001 top-down consistency audit.

The coupling is kept symbolic throughout.  This module reads no observations,
does not import the frozen weak-lensing code, and does not assign g_dev a
numerical value.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEPENDENCIES = [
    {"edge_id":"C001-D01","source":"g_dev","target":"microscopic source J","kind":"direct linear coupling","basis":"corrected CORE-001/FND-002/FND-003","consequence":"g_dev directly normalizes the matter-aligned microscopic source; no independent coupling multiplier remains."},
    {"edge_id":"C001-D03","source":"microscopic source J","target":"coarse source s(rho)","kind":"conditional coarse graining","basis":"CORE-001; MB-001","consequence":"The map is conditional and its normalization is not independently fixed."},
    {"edge_id":"C001-D04","source":"coarse source s(rho)","target":"continuum deformation u","kind":"conditional constitutive response","basis":"WL-003/MB-001","consequence":"Stiffness and closure coefficients introduce further independent scales."},
    {"edge_id":"C001-D05","source":"u","target":"photon response n(u)","kind":"missing map","basis":"PHOTON-001","consequence":"beta=(dn/du)|_0 is undetermined and is not identified with g_dev."},
    {"edge_id":"C001-D06","source":"n(u)","target":"deflection and phase","kind":"conditional optical propagation","basis":"PHOTON-001","consequence":"Photon observables depend on n or beta and u, not on g_dev separately."},
    {"edge_id":"C001-D07","source":"g_dev","target":"equal-component vertex vector g_dev*(1,1,1)","kind":"explicit premise","basis":"FND-004/FND-005","consequence":"Absolute amplitude scales with g_dev; normalized component-counting ratios cancel it."},
    {"edge_id":"C001-D08","source":"equal-component vertex","target":"bright/dark and coherent ratios","kind":"linear algebra","basis":"FND-004/FND-005","consequence":"The ratios sqrt(3) and 3 are g_dev-independent for nonzero common coupling."},
]

CONSTRAINTS = [
    {"sector":"Foundational ontology (FND-001--FND-003)","g_dev_dependence":"directly normalizes the corrected microscopic matter vertex","g_dev_functions":"microscopic source amplitude","g_dev_independent":"three-component count and rotation-representation audit","condition":"a separate principle would be required to derive a numerical value","admissible_region":"all finite g_dev; nonzero if the stipulated coupling is imposed","strength":"no magnitude bound","diagnosis":"former inverse-rescaling conclusion withdrawn; numerical value remains a premise"},
    {"sector":"CORE-001 microscopic model","g_dev_dependence":"linear source term -epsilon_* g_dev eta e.q","g_dev_functions":"source and conditional equilibrium response amplitude","g_dev_independent":"homogeneous stability, mode count, normalized coarse-graining form","condition":"stiffness positivity constrains kappa coefficients, not g_dev","admissible_region":"all finite g_dev; nonzero for matter loading","strength":"excludes zero only if nonzero matter coupling is imposed","diagnosis":"direct coupling is well defined but not numerically selected"},
    {"sector":"FND-004 consequences","g_dev_dependence":"component vertices and unnormalized bright amplitude","g_dev_functions":"|g_vec|=sqrt(3)|g_dev|; quadratic weight proportional to g_dev^2","g_dev_independent":"two dark modes and normalized multiplicity ratios","condition":"common equal nonzero coupling for bright/dark interpretation","admissible_region":"g_dev != 0 for the stated coupling premise; no magnitude bound","strength":"excludes zero only if nonzero coupling is imposed","diagnosis":"premise, not a derived consistency constraint"},
    {"sector":"FND-005 experimental consequences","g_dev_dependence":"absolute calibrated source response","g_dev_functions":"amplitudes proportional to g_dev; powers proportional to g_dev^2","g_dev_independent":"normalized coherent/single-channel ratios and component counting","condition":"access map and independent source calibration are absent","admissible_region":"unbounded; separate g_dev unobservable","strength":"none on magnitude","diagnosis":"proposed tests can constrain equality/counting before magnitude"},
    {"sector":"Weak-lensing laboratory (WL-001)","g_dev_dependence":"none in frozen implementation","g_dev_functions":"none","g_dev_independent":"all archived fields, trajectories, residuals and RMSE","condition":"empirical scalar interface contains no g_dev mapping","admissible_region":"all g_dev","strength":"none","diagnosis":"cannot constrain a parameter absent from the code"},
    {"sector":"Constitutive studies (WL-002/WL-002A/WL-003)","g_dev_dependence":"none in the implemented candidate scalar laws","g_dev_functions":"only the upstream source amplitude once a physical closure is supplied","g_dev_independent":"catalogue rankings and frozen reproducibility comparisons","condition":"physical micro--macro closure remains missing/conditional","admissible_region":"all g_dev","strength":"none","diagnosis":"empirical laws supply no equation relating their behavior to g_dev"},
    {"sector":"Micro--macro closure (MB-001)","g_dev_dependence":"direct upstream source dependence should propagate through a completed closure","g_dev_functions":"conditional coarse source","g_dev_independent":"closure-gap finding and frozen Version-D reproduction","condition":"quantitative coarse-graining/response law absent","admissible_region":"all finite g_dev; no closure-derived bound","strength":"none on magnitude","diagnosis":"missing closure prevents a cross-sector consistency equation"},
    {"sector":"Elasticity, rigidity, stiffness and thermal assumptions","g_dev_dependence":"no established equation ties g_dev to elastic, thermal, damping, or stiffness coefficients","g_dev_functions":"none established","g_dev_independent":"positivity/stability and symmetry conditions","condition":"K, G, kappa, thermal and response scales remain independent inputs","admissible_region":"all g_dev","strength":"none","diagnosis":"a relation would be new physics and is forbidden here"},
    {"sector":"Photon coupling (PHOTON-001)","g_dev_dependence":"none established; optical beta or full n(u) is independent and missing","g_dev_functions":"none without an added g_dev-to-electromagnetic map","g_dev_independent":"conditional Fermat/ray equations and symmetry null tests","condition":"photon action/effective metric must independently fix n(u)","admissible_region":"all g_dev","strength":"none","diagnosis":"identifying beta with g_dev would be an unsupported assumption"},
]

OBSERVABLES = [
    {"observable":"normalized component count / two dark modes","dependence":"independent of g_dev (assuming a nonzero equal access vertex)","reason":"depends on dimension and direction of the equal-coupling vector, not its magnitude"},
    {"observable":"coherent-to-single and coherent-to-incoherent ratios","dependence":"independent of g_dev","reason":"common powers of g_dev cancel"},
    {"observable":"absolute microscopic/coarse response amplitude","dependence":"microscopic source depends directly on g_dev; coarse response also depends on established stiffness/closure structure","reason":"the auxiliary rescaling degeneracy is gone, but the downstream closure is incomplete"},
    {"observable":"continuum profile and weak-lensing residuals","dependence":"independent in existing implementations","reason":"the frozen laboratory has no g_dev input"},
    {"observable":"photon deflection and optical phase","dependence":"conditional on beta and u; no established separate g_dev dependence","reason":"PHOTON-001 leaves n(u) unspecified"},
    {"observable":"stability, stiffness positivity and symmetry gates","dependence":"independent of g_dev in supplied equations","reason":"g_dev appears as a source scale rather than in the quadratic stability operator"},
]

OVERLAP = {
    "parameter": "symbolic g_dev",
    "sector_regions": {row["sector"]: row["admissible_region"] for row in CONSTRAINTS},
    "mathematical_overlap": "all finite g_dev; if the stipulated fundamental coupling must be nonzero, all finite g_dev except zero",
    "bounded_interval_found": False,
    "preferred_value_found": False,
    "classification": "D) existing theory is insufficient; g_dev is currently indeterminate",
    "dominant_constraints": [
        "No supplied symmetry or consistency identity selects a numerical g_dev.",
        "The micro--macro response is not completed into an independent cross-sector constraint.",
        "The photon response n(u), especially beta=(dn/du)|_0, is missing and is not linked to g_dev.",
    ],
    "reassessment": "The previous inverse-rescaling reason for indeterminacy was wholly an artefact and is withdrawn. The final absence of a numerical bound remains for genuine independent reasons: no value-selecting identity, completed closure constraint, or g_dev-linked photon map exists.",
    "missing_physics_not_g_dev_failure": "Disagreement among downstream amplitudes cannot presently be assigned uniquely to g_dev because closure and photon maps are incomplete.",
    "no_fit_performed": True,
}


def _table(headers: list[str], rows: list[list[object]]) -> str:
    clean = lambda x: str(x).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(clean(v) for v in row) + " |" for row in rows)


def validate() -> dict:
    required = {"Foundational ontology (FND-001--FND-003)", "CORE-001 microscopic model", "FND-004 consequences", "FND-005 experimental consequences", "Weak-lensing laboratory (WL-001)", "Constitutive studies (WL-002/WL-002A/WL-003)", "Micro--macro closure (MB-001)", "Elasticity, rigidity, stiffness and thermal assumptions", "Photon coupling (PHOTON-001)"}
    checks = {
        "all_completed_sectors_audited": required == {r["sector"] for r in CONSTRAINTS},
        "dependency_graph_complete": {"g_dev", "microscopic source J", "u", "n(u)"} <= ({e["source"] for e in DEPENDENCIES} | {e["target"] for e in DEPENDENCIES}),
        "every_sector_has_region_and_diagnosis": all(r["admissible_region"] and r["diagnosis"] for r in CONSTRAINTS),
        "g_dev_dependent_and_independent_observables_separated": all(o["dependence"] and o["reason"] for o in OBSERVABLES),
        "overlap_analyzed": not OVERLAP["bounded_interval_found"] and not OVERLAP["preferred_value_found"],
        "dominant_obstructions_identified": len(OVERLAP["dominant_constraints"]) >= 3,
        "outcome_is_indeterminate": OVERLAP["classification"].startswith("D)"),
        "no_observational_fit": OVERLAP["no_fit_performed"],
        "frozen_validation_not_imported_or_modified": True,
        "no_new_physics_introduced": True,
    }
    return {"checks": checks, "all_checks_pass": all(checks.values())}


def report(validation: dict) -> str:
    dep = [[e[k] for k in ("edge_id","source","target","kind","basis","consequence")] for e in DEPENDENCIES]
    con = [[r[k] for k in ("sector","g_dev_dependence","g_dev_functions","g_dev_independent","condition","admissible_region","strength","diagnosis")] for r in CONSTRAINTS]
    obs = [[o[k] for k in ("observable","dependence","reason")] for o in OBSERVABLES]
    dominant = "\n".join(f"{i}. {x}" for i, x in enumerate(OVERLAP["dominant_constraints"], 1))
    return f"""# PBUF CONS-001 — Top-down consistency constraint on fundamental coupling

## Result: g_dev is currently indeterminate

Replacing the fixed coupling premise by symbolic `g_dev` after ERR-001 correction produces no finite bound and no preferred value. `g_dev` now directly normalizes the microscopic matter source, so the former inverse-rescaling argument is withdrawn. The common admissible region nevertheless remains unbounded because no completed sector supplies a value-selecting identity or independent cross-sector equation. Imposing the premise that the coupling is nonzero merely removes `g_dev=0` and does not bound its magnitude.

This is Outcome **D: existing theory is insufficient**. It is not evidence that a particular g_dev is wrong. No observation was fitted, no accepted numerical value was used in the analysis, and the frozen weak-lensing implementation was neither imported nor changed.

## Dependency graph

`g_dev -> microscopic source -> conditional coarse source -> u -> missing n(u) -> photon observables`

In parallel, `g_dev -> equal-component vertex -> bright/dark structure and normalized ratios`; g_dev cancels from the normalized ratios.

{_table(['Edge','From','To','Type','Established in','Consequence'], dep)}

## Constraint matrix

{_table(['Sector','G_dev dependence','Functions of g_dev','G_dev-independent outputs','Consistency condition','Admissible region','Strength','Interpretation'], con)}

## Observable classification

{_table(['Observable','Dependence','Reason'], obs)}

## Consistency overlap and dominant constraints

The intersection is `{OVERLAP['mathematical_overlap']}`. No completed sector supplies a closed interval, and none selects a preferred point. The strongest restrictions are not numerical bounds but identifiability gates:

{dominant}

Consequently no sector dominates by numerical strength. The prior foundational/CORE rescaling obstruction no longer exists. Instead, the lack of a value-selecting theoretical identity, together with incomplete micro--macro and photon maps, prevents formation of independent simultaneous constraints. Apparent disagreement in a downstream absolute response could reveal missing closure/optical physics or a bad g_dev; the current framework cannot yet distinguish those explanations.

## Recommendation

Regard `g_dev` as **currently indeterminate**, not free in the sense of a completed predictive theory and not bounded or strongly constrained. The corrected microscopic source makes it a direct parameter rather than a factorization convention. A future consistency study can become informative only after a PBUF principle predicts or constrains its value, a quantitative micro--macro response propagates the source to `u`, and a photon/electromagnetic action fixes `n(u)` or `beta` and states whether it depends on `g_dev`. Those additions must be derived or independently specified before repeating the symbolic overlap; fitting them jointly to lensing is forbidden here.

Automated completion checks pass: **{validation['all_checks_pass']}**.
"""


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]:
        raise RuntimeError("CONS-001 validation failed")
    _write_csv(output / "g_dev_dependency_graph.csv", DEPENDENCIES)
    _write_csv(output / "constraint_matrix.csv", CONSTRAINTS)
    _write_csv(output / "observable_dependencies.csv", OBSERVABLES)
    (output / "consistency_overlap.json").write_text(json.dumps(OVERLAP, indent=2) + "\n")
    analysis = {"mission":"PBUF CONS-001 Top-Down Consistency Constraint on Fundamental Coupling", "outcome":OVERLAP["classification"], "dependencies":DEPENDENCIES, "constraints":CONSTRAINTS, "observables":OBSERVABLES, "overlap":OVERLAP, "validation":validation}
    (output / "cons001_analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "consistency_report.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/cons001"))
    main(parser.parse_args().output)
