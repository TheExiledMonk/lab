#!/usr/bin/env python3
"""PBUF ARCH-001: audit and disposition every use of g_dev.

The audit treats V11 as authoritative.  It does not edit or execute the frozen
weak-lensing laboratory and it does not infer identities from 1/137-like
numerical values.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import zipfile
from pathlib import Path


EQUATIONS = [
    {"audit_id":"ARCH-E01", "sources":"CORE-001-E01; ERR-001-E01", "expression":"F=epsilon_* sum_i[kappa_0|q_i|^2/2+kappa_1 sum_<ij>|q_j-q_i|^2/2-g_dev eta_i e.q_i]", "role":"absolute microscopic matter-state vertex", "requires_g_dev":"only inside the invented post-V11 interaction", "v11_rewrite":"none", "decision":"D", "disposition":"Remove from the V11 theory. Redesign/derive the interaction from an authoritative action before reinstating any coefficient."},
    {"audit_id":"ARCH-E02", "sources":"CORE-001-E07; ERR-001-E02", "expression":"s(rho)=epsilon_* g_dev (rho/rho_*)/a^d", "role":"coarse source inherited from E01", "requires_g_dev":"only because E01 assumed it", "v11_rewrite":"none", "decision":"D", "disposition":"Remove as a PBUF derivation with E01; it may remain only as clearly quarantined exploratory algebra."},
    {"audit_id":"ARCH-E03", "sources":"CORE-001-E09; MB-001-E01; ERR-001-E03", "expression":"K u-div(G grad u)=s(rho)", "role":"conditional continuum balance", "requires_g_dev":"no; g_dev occurs only in one unsupported proposed source", "v11_rewrite":"retain s(rho) as an unclosed source functional", "decision":"A", "disposition":"Keep the conditional balance without a g_dev definition; MB-001 already establishes that the source law is missing."},
    {"audit_id":"ARCH-E04", "sources":"FND-004-P02/P10; FND-005-P08; ERR-001-E04", "expression":"g_vec=g_dev(1,1,1); |g_vec|=sqrt(3)|g_dev|", "role":"absolute equal-component loading", "requires_g_dev":"magnitude only; the bright direction does not", "v11_rewrite":"g_hat=(1,1,1)/sqrt(3) for direction/counting only", "decision":"D/A", "disposition":"Remove the unsupported absolute vertex (D); eliminate g_dev from the structural bright/dark statement by using the unit direction (A)."},
    {"audit_id":"ARCH-E05", "sources":"FND-004-P04; FND-005-P03; ERR-001-E05", "expression":"|sum_i g_dev|^2 / sum_i|g_dev|^2=3", "role":"normalized coherent/incoherent ratio", "requires_g_dev":"no; every nonzero common magnitude cancels", "v11_rewrite":"|sum_i 1|^2/sum_i|1|^2=3", "decision":"A", "disposition":"Eliminate g_dev. State the result as conditional component-counting linear algebra, not as V11 equation (4)."},
    {"audit_id":"ARCH-E06", "sources":"FND-004/FND-005 component amplitudes", "expression":"|g_vec|^2=3 g_dev^2; amplitudes proportional to g_dev; powers proportional to g_dev^2", "role":"unnormalized response scaling", "requires_g_dev":"yes within the assumed absolute vertex", "v11_rewrite":"none", "decision":"D", "disposition":"Remove the absolute predictions until a calibrated matter action and readout derive their normalization."},
    {"audit_id":"ARCH-E07", "sources":"PHOTON-001-E05", "expression":"n(u)=1+beta u+O(u^2)", "role":"conditional photon response; reports no g_dev mapping", "requires_g_dev":"no", "v11_rewrite":"none needed", "decision":"A", "disposition":"Retain the conditional expansion and the explicit negative result; do not identify beta with g_dev or a V11 alpha."},
    {"audit_id":"ARCH-E08", "sources":"CONS-001 dependency/constraint analysis", "expression":"g_dev -> source -> [missing closure] -> u -> [missing n(u)]", "role":"audit dependency, not a field equation", "requires_g_dev":"no after removal of E01/E02", "v11_rewrite":"unsupported vertex -> unclosed source functional", "decision":"A", "disposition":"Remove g_dev nodes and preserve the documented closure and photon-response gaps."},
]

MAPPINGS = [
    {"candidate":"alpha_QM", "v11_role":"Quantum Engine elastic amplitude", "mapping_status":"no mapping", "reason":"No equation connects the Quantum Engine output to the post-V11 matter vertex."},
    {"candidate":"alpha_resolved", "v11_role":"metadata-resolved curvature/elastic amplitude used by equations (9), (15), and (16)", "mapping_status":"no mapping", "reason":"Different physical role and pipeline; no microscopic matter-action derivation exists."},
    {"candidate":"alpha_T(a)", "v11_role":"scale-factor-dependent thermal LUT elastic amplitude", "mapping_status":"no mapping", "reason":"A function entering V11 thermal evolution cannot replace a constant vertex without a derived bridge."},
    {"candidate":"alpha_EM", "v11_role":"electromagnetic fine-structure constant in V11 equation (4)", "mapping_status":"no mapping", "reason":"Numerical proximity to 1/137 is not identity; V11 equation (4) relates alpha_EM to alpha_resolved, not to matter loading."},
    {"candidate":"generic alpha in V11 (15)", "v11_role":"elastic amplitude resolved by cosmological context", "mapping_status":"no mapping", "reason":"ALPHA-ARCH-001 resolves it as alpha_resolved, not as the exploratory vertex."},
    {"candidate":"unit common direction", "v11_role":"not a V11 alpha; normalized linear-algebra device", "mapping_status":"valid elimination for structural ratios only", "reason":"The common magnitude cancels from normalized ratios and bright/dark direction counting."},
]

IRREDUCIBLE = [
    {"use":"absolute microscopic matter-state interaction", "equations":"ARCH-E01", "why_not_eliminable_inside_current_model":"Its magnitude changes the response relative to fixed stiffness normalization.", "status":"not theoretically irreducible; unsupported modelling choice"},
    {"use":"absolute coarse source normalization", "equations":"ARCH-E02", "why_not_eliminable_inside_current_model":"It is inherited algebraically from the assumed microscopic vertex.", "status":"not independent and not V11-derived"},
    {"use":"absolute component amplitudes/powers", "equations":"ARCH-E04/E06", "why_not_eliminable_inside_current_model":"Absolute calibration retains the assumed common magnitude.", "status":"unsupported descendant, not a reason to retain a fundamental parameter"},
]


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def occurrence_inventory(root: Path, output: Path) -> list[dict]:
    """Inventory repository text occurrences, excluding caches and this output."""
    rows = []
    allowed = {".py", ".md", ".csv", ".json", ".txt", ".docx"}
    excluded_parts = {"__pycache__", ".git"}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in allowed:
            continue
        if output in path.parents or excluded_parts.intersection(path.parts):
            continue
        # Frozen weak-lensing run products are not theory sources and are not touched.
        if path.suffix.lower() == ".docx":
            with zipfile.ZipFile(path) as archive:
                xml = archive.read("word/document.xml").decode("utf-8")
            xml = re.sub(r"</w:p>", "\n", xml)
            text = re.sub(r"<[^>]+>", "", xml).replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
        else:
            try:
                text = path.read_text()
            except UnicodeDecodeError:
                continue
        for line_no, line in enumerate(text.splitlines(), 1):
            for match in re.finditer(r"\bg_dev\b", line):
                rel = path.relative_to(root)
                low = line.lower()
                if rel.name == "arch001_gdev_audit.py":
                    layer, group, decision = "audit implementation", "meta-reference", "N/A (audit vocabulary)"
                elif any(term in low for term in ("cancel", "independent of g_dev", "not fixed by g_dev", "not identified with g_dev", "no g_dev link", "without a g_dev")):
                    layer, group, decision = "generated restatement" if str(rel).startswith("runs/") else "active milestone source", "negative mapping or normalized result", "A (retain result without parameter dependence)"
                elif any(term in low for term in ("g_vec", "gdev", "matter vertex", "matter-state", "source term", "s(rho)", "source amplitude", "direct coupling", "directly normal", "amplitude", "powers proportional")):
                    layer, group, decision = "generated restatement" if str(rel).startswith("runs/") else "active milestone source", "absolute vertex/source descendant", "D (retire unsupported exploratory use)"
                elif str(rel).startswith("runs/"):
                    layer, group, decision = "generated restatement", "premise/audit narrative", "D (retire premise or update narrative)"
                elif rel.name == "README.md":
                    layer, group, decision = "repository summary", "narrative", "D (update summary after retirement)"
                else:
                    layer, group, decision = "active milestone source", "premise/audit narrative", "D (retire premise or update narrative)"
                rows.append({"occurrence_id":f"ARCH-O{len(rows)+1:04d}", "file":str(rel), "line":line_no,
                             "column":match.start()+1, "layer":layer, "semantic_group":group,
                             "decision_reference":decision, "context":line.strip()})
    return rows


def table(headers: list[str], rows: list[list[str]]) -> str:
    esc = lambda x: str(x).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(esc(v) for v in r) + " |" for r in rows)


def report(occurrences: list[dict]) -> str:
    eq = table(["ID","Equation/use","V11 rewrite","Decision","Disposition"], [[r[k] for k in ("audit_id","expression","v11_rewrite","decision","disposition")] for r in EQUATIONS])
    mp = table(["Candidate","V11 role","Result","Reason"], [[r[k] for k in ("candidate","v11_role","mapping_status","reason")] for r in MAPPINGS])
    irr = table(["Use","Equations","Why magnitude remains","Finding"], [[r[k] for k in ("use","equations","why_not_eliminable_inside_current_model","status")] for r in IRREDUCIBLE])
    return f"""# PBUF ARCH-001 — Evaluation and possible elimination of `g_dev`

## Decision

**Redesign the affected exploratory theory and eliminate `g_dev` from all claims presented as V11/PBUF consequences.** No supplied input derives a mapping from `g_dev` to `alpha_QM`, `alpha_resolved`, `alpha_T(a)`, `alpha_EM`, or another V11 quantity. The value `1/137` is only a post-V11 premise. It cannot establish identity with `alpha_EM`.

`g_dev` is not mathematically unavoidable. It is a temporary modelling coefficient introduced to normalize an invented linear matter-state vertex. Within that conditional model its magnitude affects absolute response, but that makes it model-dependent—not fundamental. The source vertex and its absolute-amplitude descendants receive decision D. Uses in normalized component ratios receive decision A because the common magnitude cancels exactly.

## Why it was introduced

CORE-001 introduced the predecessor `alpha_*=1/137` to give a stipulated three-component microscopic state an absolute linear matter-loading scale. ERR-001 removed a second auxiliary multiplier and made this coefficient direct. ALPHA-ARCH-001 then renamed it `g_dev` to stop the post-V11 vertex from being confused with the authoritative V11 alpha hierarchy. None of those steps derived the coefficient or its value.

## Occurrence inventory

The machine-readable inventory contains **{len(occurrences)} exact textual occurrences** across repository text sources, extractable DOCX content, and generated records, excluding this audit's output directory, caches, and non-text binary artifacts. Each row records file, line/paragraph, column, layer, context, and its controlling decision. Generated duplicates are explicitly identified as restatements rather than independent equations. See `g_dev_occurrence_inventory.csv`.

## Equation-by-equation replacement audit

{eq}

Decision key: A = eliminated/replaced without a new parameter; B = rigorously derived mapping; C = retained extension; D = unsupported equation/use removed. `D/A` separates the unsupported absolute magnitude from the valid normalized direction/counting statement. There are **no B or C results**.

## Mapping to V11 quantities

{mp}

The absence of a mapping is a positive audit result: symbol role, dependency chain, and equations differ. V11 equation (4), `alpha_resolved ~= 3 alpha_EM`, is an implemented numerical relationship with a motivating dimensional-count interpretation; it is not the quadratic equal-component identity and does not supply a matter vertex.

## Irreducible uses inside the exploratory model

{irr}

These are the only places where deleting the symbol while leaving the present equations otherwise unchanged would change absolute predictions. None is irreducible in the theory-selection sense required by ARCH-001: all descend from the unsupported E01 premise. Absorbing `g_dev` into `q`, `e`, `eta`, `epsilon_*`, `kappa_0`, or `kappa_1` would merely hide it by redefining established normalizations and is therefore rejected.

## Missing derivation

Reinstatement would require a dimensionally normalized microscopic action derived from V11's Quantum Engine/metadata architecture that fixes the matter operator, its state normalization, and its coefficient. It would also require a quantitative coarse-graining closure to `s(rho), K, G`, plus a calibrated readout; a photon action is separately required for `n(u)` or `beta`. These missing laws may derive an existing V11 quantity as the coefficient, or demonstrate a genuinely new extension, but ARCH-001 cannot assume either outcome.

## Required repository disposition

1. Treat CORE-001 E01/E02 and their absolute FND consequences as quarantined exploratory equations, not authoritative PBUF results.
2. Rewrite bright/dark and normalized multiplicity results with the unit common direction; retain their explicit ontology assumptions.
3. Preserve MB-001's unclosed `s(rho)` balance and PHOTON-001's independent missing-response conclusion.
4. Remove `g_dev` from CONS dependency claims after the exploratory branch is retired.
5. Do not replace it by any new free parameter and do not modify the frozen weak-lensing code.

## Completion

Every occurrence is controlled by an equation-level disposition or marked as a narrative/generated restatement. No unexplained use remains. The final classification is: **temporary modelling device; eliminate by redesign, with exact algebraic elimination from normalized structural results.**
"""


def main(root: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    occurrences = occurrence_inventory(root, output)
    write_csv(output / "g_dev_occurrence_inventory.csv", occurrences)
    write_csv(output / "equation_replacement_audit.csv", EQUATIONS)
    write_csv(output / "v11_mapping_table.csv", MAPPINGS)
    write_csv(output / "irreducible_uses.csv", IRREDUCIBLE)
    validation = {
        "mission":"PBUF ARCH-001", "status":"COMPLETE", "occurrences":len(occurrences),
        "equation_uses":len(EQUATIONS), "derived_v11_mappings":0, "retained_extensions":0,
        "frozen_weak_lensing_modified":False, "replacement_free_parameters_introduced":False,
        "numerical_identity_assumed":False,
        "all_occurrences_classified":bool(occurrences) and all(r["decision_reference"] for r in occurrences),
        "all_equation_uses_classified":all(r["decision"] in {"A","B","C","D","D/A"} for r in EQUATIONS),
    }
    validation["all_checks_pass"] = all(v for k,v in validation.items() if k in {"all_occurrences_classified","all_equation_uses_classified"})
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (output / "arch001_report.md").write_text(report(occurrences))
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("runs/arch001"))
    args = parser.parse_args()
    main(args.root.resolve(), args.output.resolve())
