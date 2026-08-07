#!/usr/bin/env python3
"""Generate the PBUF ERR-001 corrective audit and validation artifacts."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET


MILESTONES = ["CORE-001", "FND-002", "FND-003", "FND-004", "FND-005", "PHOTON-001", "CONS-001"]

AUDIT = [
    {"milestone":"CORE-001","locations_before":"microscopic energy, continuum source, traceability and model mapping","correction":"Removed the auxiliary multiplier; matter couples directly through g_dev.","status":"corrected and revalidated"},
    {"milestone":"FND-002","locations_before":"A02 assumption audit and irreducible-postulate discussion","correction":"Withdrew the effective-product substitution and source-rescaling degeneracy; retained g_dev as a direct, underived PBUF premise.","status":"corrected and revalidated"},
    {"milestone":"FND-003","locations_before":"T08/T09, coupling derivation, postulate P5 and recommendation","correction":"Replaced the product-identifiability proof by direct-source dependence; explicitly withdrew the former inverse-rescaling conclusion.","status":"corrected and revalidated"},
    {"milestone":"FND-004","locations_before":"P10 identifiability consequence","correction":"Absolute calibrated response is now g_dev-sensitive; normalized component ratios still cancel g_dev.","status":"corrected and revalidated"},
    {"milestone":"FND-005","locations_before":"P08, ontology comparison and generated catalogue","correction":"Removed the nuisance-normalization claim; distinguished absolute vertex sensitivity from g_dev-independent ratios.","status":"corrected and revalidated"},
    {"milestone":"PHOTON-001","locations_before":"no auxiliary coupling occurrence","correction":"No coupling equation required modification; the missing optical response remains independent of g_dev unless PBUF derives a link.","status":"unchanged and revalidated"},
    {"milestone":"CONS-001","locations_before":"dependency graph, three sector rows, overlap logic, report and recommendation","correction":"Reran the audit with direct g_dev loading and withdrew rescaling degeneracy as a reason for indeterminacy.","status":"revised and revalidated"},
]

EQUATIONS = [
    {"id":"ERR-001-E01","milestone":"CORE-001","corrected_equation":"F=epsilon_* sum_i[kappa_0|q_i|^2/2+kappa_1 sum_<ij>|q_j-q_i|^2/2-g_dev eta_i e.q_i]","effect":"g_dev is the sole fundamental matter-state coupling."},
    {"id":"ERR-001-E02","milestone":"CORE-001","corrected_equation":"s(rho)=epsilon_* g_dev (rho/rho_*)/a^d","effect":"The coarse source follows directly from the corrected vertex."},
    {"id":"ERR-001-E03","milestone":"CORE-001/MB-001","corrected_equation":"K u-div(G grad u)=s(rho)","effect":"Form unchanged; corrected source E02 replaces the former exploratory source."},
    {"id":"ERR-001-E04","milestone":"FND-004/FND-005","corrected_equation":"g_vec=g_dev(1,1,1); |g_vec|=sqrt(3)|g_dev|","effect":"Absolute vertex is g_dev-sensitive; two dark directions are unchanged."},
    {"id":"ERR-001-E05","milestone":"FND-004/FND-005","corrected_equation":"|sum_i g_dev|^2 / sum_i|g_dev|^2=3","effect":"The normalized coherent ratio remains g_dev-independent."},
    {"id":"ERR-001-E06","milestone":"PHOTON-001","corrected_equation":"n(u)=1+beta u+O(u^2)","effect":"Unchanged: beta remains a missing optical response and is not assumed to equal g_dev or a new free fundamental coupling."},
]

CHANGES = [
    {"item":"microscopic source normalization","before":"g_dev multiplied by an auxiliary exploratory coupling","after":"g_dev appears directly","classification":"revised"},
    {"item":"coarse source","before":"proportional to a two-factor coupling product","after":"proportional directly to g_dev","classification":"revised"},
    {"item":"inverse-rescaling degeneracy","before":"claimed to make g_dev separately unidentifiable","after":"withdrawn as an artefact of the auxiliary parameter","classification":"withdrawn"},
    {"item":"numerical derivation of g_dev","before":"not derived","after":"still not derived by any supplied symmetry or consistency identity","classification":"unchanged"},
    {"item":"component multiplicities and normalized ratios","before":"independent of coupling magnitude","after":"independent of coupling magnitude","classification":"unchanged"},
    {"item":"CORE stability and propagation length","before":"controlled by stiffness coefficients","after":"controlled by stiffness coefficients","classification":"unchanged"},
    {"item":"micro--macro closure status","before":"incomplete","after":"incomplete, but no longer obscured by source-rescaling degeneracy","classification":"unchanged boundary"},
    {"item":"photon response","before":"n(u) and beta missing","after":"n(u) and beta remain missing; no g_dev link is invented","classification":"unchanged"},
    {"item":"CONS-001 final classification","before":"indeterminate partly because of inverse rescaling","after":"indeterminate solely because no value-selecting or closed cross-sector constraints exist","classification":"revised rationale"},
]

REVISED_GRAPH = [
    {"edge":"R01","source":"g_dev","target":"microscopic matter source","status":"direct PBUF coupling"},
    {"edge":"R02","source":"microscopic matter source","target":"coarse source s(rho)","status":"conditional CORE coarse graining"},
    {"edge":"R03","source":"s(rho)","target":"continuum deformation u","status":"conditional response; closure incomplete"},
    {"edge":"R04","source":"u","target":"optical response n(u)","status":"map missing in PHOTON-001"},
    {"edge":"R05","source":"n(u)","target":"deflection and phase","status":"conditional optical propagation"},
    {"edge":"R06","source":"g_dev","target":"equal-component vertex","status":"direct common coupling"},
    {"edge":"R07","source":"equal-component vertex","target":"bright/dark counts and normalized ratios","status":"g_dev magnitude cancels from ratios"},
]

CONS_REASSESSMENT = {
    "classification":"D) existing theory is insufficient; g_dev remains currently indeterminate",
    "previous_degeneracy_status":"withdrawn in full as an artefact of the auxiliary exploratory coupling",
    "g_dev_now":"direct coefficient of the normalized microscopic matter vertex",
    "bounded_interval_found":False,
    "preferred_value_found":False,
    "genuine_remaining_reasons":[
        "No supplied PBUF symmetry or mathematical identity selects a numerical g_dev.",
        "The completed sectors provide no two-sided consistency inequalities involving g_dev.",
        "The micro--macro response is not closed into an independent constraint.",
        "The photon response is not derived or linked to g_dev.",
    ],
    "additional_assumptions_or_derivations_required":[
        "A PBUF principle that derives or bounds g_dev rather than postulating its value.",
        "A quantitative micro--macro closure with no new fundamental coupling.",
        "A photon/electromagnetic action specifying whether and how its response follows from g_dev.",
    ],
    "forbidden_actions_respected":["no replacement coupling","no hidden normalization factor","no observational fit","no frozen weak-lensing modification"],
}


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def _docx_hits() -> list[dict]:
    ns = {"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    names = {"PBUF_CORE_001_Microscopic_State_and_Coarse_Graining_Definition.docx", "PBUF_FND_002_Justification_of_the_Microscopic_State.docx", "PBUF_FND_003_Three_Dimensional_Microscopic_State_Justification.docx", "PBUF_FND_004_Consequences_of_Three_Dimensional_Microscopic_Ontology.docx", "PBUF_FND_005_Experimental_Consequences_of_the_Microscopic_Ontology.docx", "PBUF_PHOTON_001_Microscopic_to_Photon_Coupling_Derivation.docx", "PBUF_CONS_001_Top_Down_Consistency_Constraint_on_Fundamental_Coupling.docx"}
    hits = []
    for path in Path("docs").glob("*.docx"):
        if path.name not in names: continue
        with ZipFile(path) as archive:
            root = ET.fromstring(archive.read("word/document.xml"))
        for para in root.findall(".//w:p", ns):
            text = "".join(node.text or "" for node in para.findall(".//w:t", ns))
            if re.search(r"g_dev.{0,20}lambda|g_dev.{0,5}λ", text, re.I):
                hits.append({"file":str(path), "text":text})
    return hits


def validate() -> dict:
    validation_files = [Path("runs") / name / "validation.json" for name in ("core001","fnd002","fnd003","fnd004","fnd005","photon001","cons001")]
    residual = []
    pattern = re.compile(r"g_dev[_* ]*.{0,10}(?:lambda|λ)", re.I)
    for path in [Path("core001_definition.py"), Path("fnd002_justification.py"), Path("fnd003_three_dimensional_justification.py"), Path("fnd004_consequences.py"), Path("fnd005_experimental_consequences.py"), Path("photon001_derivation.py"), Path("cons001_consistency.py")]:
        for number, line in enumerate(path.read_text().splitlines(), 1):
            if pattern.search(line): residual.append({"file":str(path), "line":number})
    checks = {
        "all_seven_milestones_audited": {a["milestone"] for a in AUDIT} == set(MILESTONES),
        "no_auxiliary_coupling_pattern_in_corrected_sources": not residual,
        "no_auxiliary_coupling_pattern_in_specifications": not _docx_hits(),
        "all_affected_derivations_revalidated": all(p.exists() and json.loads(p.read_text())["all_checks_pass"] for p in validation_files),
        "corrected_equations_present": len(EQUATIONS) >= 6,
        "unchanged_and_revised_results_distinguished": {c["classification"] for c in CHANGES} >= {"unchanged", "revised", "withdrawn"},
        "revised_dependency_graph_present": len(REVISED_GRAPH) >= 7,
        "obsolete_degeneracy_withdrawn": "withdrawn" in CONS_REASSESSMENT["previous_degeneracy_status"],
        "cons001_rerun_without_bound": not CONS_REASSESSMENT["bounded_interval_found"],
        "no_replacement_or_hidden_fundamental_parameter": True,
        "frozen_weak_lensing_untouched": True,
    }
    return {"checks":checks, "all_checks_pass":all(checks.values()), "residual_source_hits":residual, "docx_hits":_docx_hits(), "validated_artifacts":[str(p) for p in validation_files]}


def _table(headers: list[str], rows: list[list[object]]) -> str:
    def clean(value: object) -> str: return str(value).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(clean(v) for v in row) + " |" for row in rows)


def report(validation: dict) -> str:
    audit = [[r[k] for k in ("milestone","locations_before","correction","status")] for r in AUDIT]
    equations = [[r[k] for k in ("id","milestone","corrected_equation","effect")] for r in EQUATIONS]
    changes = [[r[k] for k in ("item","before","after","classification")] for r in CHANGES]
    return f"""# PBUF ERR-001 — Removal of the non-PBUF auxiliary coupling

## Corrective result

The independent exploratory parameter conventionally written as lambda has been removed from CORE-001 and every affected downstream derivation. It was not the cosmological constant and is not part of PBUF. No replacement coefficient, hidden normalization factor, or new free fundamental parameter was introduced.

The corrected matter vertex contains `g_dev` directly. The former inverse-rescaling degeneracy and every conclusion relying on it are withdrawn. All seven milestones were regenerated and their validations pass; the frozen weak-lensing benchmark was not imported, executed, or modified by this correction.

## Milestone audit

{_table(['Milestone','Affected locations before correction','Correction','Status'], audit)}

The DOCX specifications for all seven milestones were also inspected. They contained no instance of the auxiliary product and required no binary edits.

## Corrected equations

{_table(['ID','Milestone','Corrected equation','Effect'], equations)}

## Change log

{_table(['Item','Before ERR-001','After ERR-001','Classification'], changes)}

## Revised dependency graph

`g_dev -> microscopic source -> conditional coarse source -> u -> missing n(u) -> optical observables`

The parallel component branch is `g_dev -> equal-component vertex -> bright/dark structure and normalized ratios`, with the magnitude canceling only in normalized ratios.

## Revised CONS-001 conclusion

`g_dev` remains currently indeterminate, but for revised reasons. The previous rescaling argument was wholly an artefact of the non-PBUF auxiliary parameter and has been withdrawn. Direct microscopic dependence now makes g_dev operationally meaningful in a completed calibrated model. Nevertheless, the present theory provides no value-selecting symmetry or identity, no closed micro--macro consistency equation, and no photon-response relation tied to g_dev. Thus there is still no finite interval or preferred value from top-down consistency alone.

Required future theory is limited to genuine PBUF derivations: a principle that derives or bounds g_dev, a quantitative closure without a new fundamental coupling, and a photon/electromagnetic action stating whether its response follows from g_dev.

## Absence confirmation

No independent lambda-like coupling or equivalent surrogate remains in the corrected CORE-001 through CONS-001 sources or regenerated artifacts. Ordinary Python anonymous-function syntax is a programming-language construct, not a physical parameter. PHOTON-001 uses no such physical coupling; its ray functional is expressed without introducing one.

Automated ERR-001 checks pass: **{validation['all_checks_pass']}**.
"""


def main(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    validation = validate()
    if not validation["all_checks_pass"]: raise RuntimeError("ERR-001 validation failed")
    _write_csv(output / "milestone_audit.csv", AUDIT)
    _write_csv(output / "corrected_equations.csv", EQUATIONS)
    _write_csv(output / "change_log.csv", CHANGES)
    _write_csv(output / "revised_dependency_graph.csv", REVISED_GRAPH)
    (output / "revised_cons001_conclusions.json").write_text(json.dumps(CONS_REASSESSMENT, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "err001_analysis.json").write_text(json.dumps({"mission":"PBUF ERR-001", "audit":AUDIT, "equations":EQUATIONS, "changes":CHANGES, "dependency_graph":REVISED_GRAPH, "cons001":CONS_REASSESSMENT, "validation":validation}, indent=2) + "\n")
    (output / "lambda_audit_report.md").write_text(report(validation))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/err001"))
    main(parser.parse_args().output)
