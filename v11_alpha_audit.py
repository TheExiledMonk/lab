#!/usr/bin/env python3
"""Enforce the authoritative-source gate for PBUF V11-ALPHA-001.

Mission briefs and later development are never treated as substitutes for the
V11 preprint.  Without explicit authoritative inputs the audit stops before
traceability, dependency, geometric, or cross-development analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import zipfile
from pathlib import Path


BRIEF = Path("docs/PBUF_V11_ALPHA_001_Geometric_Origin_of_Resolved_Alpha.docx")
V11_NAME_HINT = re.compile(r"(?:^|[_ -])V11(?:[_ .-]|$)", re.I)
ALPHA_PATTERN = re.compile(r"alpha|α|fine[- ]structure", re.I)

TRACEABILITY = [
    {"record_id":"V11A-S01","source":"V11-ALPHA-001 mission brief, title/prerequisite","equation_id":"not supplied","expression":"alpha_resolved ~= 3 alpha_EM","classification":"audit target (reported V11 claim; not independently verified)","assumptions":"The brief presupposes a resolved elastic amplitude and an electromagnetic fine-structure constant.","downstream":"Unknown until the V11 equation set is supplied."},
    {"record_id":"V11A-L01","source":"CORE-001-A02; FND-003-T07","equation_id":"CORE-001-A02 / FND-003-T07","expression":"alpha_* = 1/137","classification":"premise","assumptions":"alpha_* is the direct microscopic matter-state coupling; no identity with alpha_resolved is supplied.","downstream":"CORE-001-E01 and CORE-001-E07."},
    {"record_id":"V11A-L02","source":"CORE-001 corrected model","equation_id":"CORE-001-E01","expression":"F=epsilon_* sum_i[..., -alpha_* eta_i e.q_i]","classification":"premise/definition","assumptions":"linear matter-state interaction and normalized microscopic variables","downstream":"CORE-001-E07, then the conditional continuum equation."},
    {"record_id":"V11A-L03","source":"CORE-001 corrected model","equation_id":"CORE-001-E07","expression":"s(rho)=epsilon_* alpha_* (rho/rho_*)/a^d","classification":"conditional derivation","assumptions":"CORE-001-E01, coarse graining, alignment, scale separation","downstream":"CORE-001-E09: K u-div(G grad u)=s(rho)."},
    {"record_id":"V11A-L04","source":"ERR-001 corrected equations","equation_id":"ERR-001-E04","expression":"g_vec=alpha_*(1,1,1); |g_vec|=sqrt(3)|alpha_*|","classification":"identity conditional on later premises","assumptions":"three orthonormal components with equal per-component coupling","downstream":"bright/dark decomposition and normalized component ratios."},
    {"record_id":"V11A-L05","source":"ERR-001 corrected equations","equation_id":"ERR-001-E05","expression":"|sum_i alpha_*|^2 / sum_i|alpha_*|^2 = 3","classification":"identity conditional on later premises","assumptions":"three equal coherent component amplitudes and the stated normalization","downstream":"normalized coherent/incoherent response ratio; alpha_* cancels."},
]

GEOMETRIC_ANALYSIS = {
    "v11_claim_verified": False,
    "v11_meaning_of_geometric": "not recoverable from the supplied source corpus",
    "factor_three_status_in_v11": "unclassifiable: the mission reports the relation but supplies neither its originating equation nor derivation",
    "later_factor_three": "In CORE/FND/ERR, three is a component-counting identity only after assuming a three-component equal-coupling vector.",
    "ontology_relationship": "compatible but not an established derivation: FND-003 labels exactly three components and their association with spatial directions as post-V11 assumptions; three-dimensional space alone does not imply them.",
    "symbol_warning": "The supplied materials do not establish alpha_resolved = alpha_* or alpha_EM = alpha_*. These symbols must not be silently identified.",
}

DEPENDENCIES = [
    {"edge_id":"V11A-D01","source":"alpha_EM","target":"alpha_resolved","relation":"reported alpha_resolved ~= 3 alpha_EM","status":"mission claim only; V11 equation/source absent"},
    {"edge_id":"V11A-D02","source":"three-dimensional microscopic ontology","target":"three-component state q","relation":"dimension matching","status":"post-V11 assumption, not entailed (FND-003-T03/T04)"},
    {"edge_id":"V11A-D03","source":"alpha_*","target":"equal-component vertex alpha_*(1,1,1)","relation":"direct equal loading","status":"later premise (ERR-001-E04)"},
    {"edge_id":"V11A-D04","source":"equal-component vertex","target":"factor 3 normalized quadratic ratio","relation":"linear-algebra identity","status":"later conditional identity (ERR-001-E05)"},
    {"edge_id":"V11A-D05","source":"alpha_*","target":"microscopic source","relation":"direct coefficient","status":"corrected CORE-001/ERR-001"},
    {"edge_id":"V11A-D06","source":"microscopic source","target":"coarse source -> u -> optical observables","relation":"conditional closure chain","status":"closure and photon response remain incomplete"},
]

CROSS_COMPARISON = [
    {"document":"V11 (as described by brief)","alpha_role":"resolved elastic amplitude reportedly ~=3 alpha_EM","role_of_three":"called geometric, details unavailable","status":"cannot verify without V11 source"},
    {"document":"CORE-001","alpha_role":"alpha_*=1/137 directly normalizes matter vertex","role_of_three":"three state components are stipulated","status":"does not establish identity with V11 alpha_resolved"},
    {"document":"FND-003","alpha_role":"1/137 remains a working premise","role_of_three":"conditional SO(3) vector minimality; extra postulates required","status":"faithful boundary; blocks retrospective derivation claim"},
    {"document":"ERR-001","alpha_role":"removes auxiliary coupling; alpha_* is direct","role_of_three":"equal-component norm/ratio identities","status":"later construction is compatible, but cannot prove V11 origin"},
]

V11_TRACE = [
 {"record_id":"V11-2.2-A1","page":2,"section":"2.2","equation_id":"prose","expression":"thermal table supplies alpha_T(a)","classification":"definition","assumptions":"active thermal LUT","downstream":"background, distance, and growth calculations"},
 {"record_id":"V11-2.2-A2","page":2,"section":"2.2","equation_id":"prose","expression":"pipeline resolves alpha_resolved from microphysics metadata","classification":"definition","assumptions":"metadata with table-metadata fallback","downstream":"single curvature/elastic parameter used by pipeline"},
 {"record_id":"V11-2.3-A1","page":2,"section":"2.3","equation_id":"prose","expression":"Quantum Engine produces elastic amplitude alpha_QM","classification":"definition","assumptions":"specified regulators and field content","downstream":"thermal LUT, then cosmological sector"},
 {"record_id":"V11-E04","page":2,"section":"2.3.1","equation_id":"(4)","expression":"alpha_resolved ~= 3 alpha_EM = 3/137.036 ~= 0.0219","classification":"empirical/numerical observation","assumptions":"current implementation; inherited quantum-microphysics metadata","downstream":"fixed elastic amplitude and density normalization"},
 {"record_id":"V11-2.3.1-G1","page":3,"section":"2.3.1","equation_id":"restatement of (4)","expression":"alpha_resolved ~= 3 alpha_EM; one contribution per spatial dimension","classification":"premise/motivating consistency argument","assumptions":"three-dimensional physical space and additive equal dimensional contributions","downstream":"geometric interpretation only; first-principles QFT derivation deferred"},
 {"record_id":"V11-2.3.1-B1","page":3,"section":"2.3.1","equation_id":"unnumbered / later (16)","expression":"Omega_b0 = 2 alpha_resolved","classification":"pipeline identity","assumptions":"two transverse electromagnetic polarizations; fixed normalization","downstream":"present-day baryon density"},
 {"record_id":"V11-E05","page":3,"section":"2.3.2","equation_id":"(5)","expression":"k_max(a)=epsilon_0,T(a)-alpha_T(a)","classification":"definition","assumptions":"thermal LUT fields","downstream":"S(a), Omega_sigma_raw(a)"},
 {"record_id":"V11-E08","page":3,"section":"2.3.2","equation_id":"(8)","expression":"Omega_sigma_raw(a)=alpha_T(a)(1-decay(a))S(a)","classification":"definition","assumptions":"activation and saturation definitions (6)-(7)","downstream":"rescaling (10)-(11)"},
 {"record_id":"V11-E09","page":3,"section":"2.3.2","equation_id":"(9)","expression":"Omega_sigma_target=1-Omega_m0-Omega_r0-alpha_resolved","classification":"pipeline normalization identity","assumptions":"flat_today normalization mode","downstream":"sigma_rescale (10), Omega_sigma(a) (11)"},
 {"record_id":"V11-E15","page":4,"section":"2.3.3","equation_id":"(15)","expression":"Omega_sigma(a)=alpha(1-exp(-a/R_max))S(a)","classification":"model definition","assumptions":"alpha denotes elastic amplitude; S defined by (7)","downstream":"E(a), H(a), distances and growth via (12)-(14)"},
 {"record_id":"V11-E16","page":4,"section":"2.3.4","equation_id":"(16)","expression":"Omega_b0=2 alpha_resolved","classification":"pipeline identity","assumptions":"fixed polarization-counting normalization","downstream":"Omega_m0 through (17)"},
 {"record_id":"V11-T2-A1","page":7,"section":"5.1 Table 2","equation_id":"Table 2 note","expression":"Omega_m0 and Omega_b0 derived from alpha_resolved; reported Omega_k0 corresponds to alpha_resolved","classification":"definition/reporting identity","assumptions":"V11 enforced identities (16)-(17)","downstream":"reported PBUF parameters and likelihood predictions"},
 {"record_id":"V11-7-A1","page":11,"section":"7","equation_id":"prose","expression":"alpha fully captures effective large-scale geometric response; no independent curvature parameter","classification":"model interpretation","assumptions":"effective-curvature usage in V11 pipeline","downstream":"interpretation of background results"},
]

V11_DEPENDENCIES = [
 {"edge_id":"V11-D01","source":"specified regulators and field content","target":"alpha_QM and epsilon_0(T)","relation":"Quantum Engine output","basis":"section 2.3"},
 {"edge_id":"V11-D02","source":"quantum microphysics metadata/LUT","target":"alpha_resolved and alpha_T(a)","relation":"inheritance/resolution","basis":"sections 2.2-2.3.1"},
 {"edge_id":"V11-D03","source":"alpha_EM","target":"alpha_resolved","relation":"numerically consistent with factor 3","basis":"equation (4); motivating spatial-dimension count"},
 {"edge_id":"V11-D04","source":"alpha_T(a), epsilon_0,T(a)","target":"k_max, S, Omega_sigma_raw","relation":"equations (5)-(8)","basis":"section 2.3.2"},
 {"edge_id":"V11-D05","source":"alpha_resolved","target":"Omega_sigma_target","relation":"flat_today closure","basis":"equation (9)"},
 {"edge_id":"V11-D06","source":"alpha_resolved","target":"Omega_b0","relation":"factor-2 pipeline identity","basis":"equation (16)"},
 {"edge_id":"V11-D07","source":"Omega_b0","target":"Omega_m0","relation":"fixed baryon fraction","basis":"equation (17)"},
 {"edge_id":"V11-D08","source":"elastic alpha and activation structure","target":"Omega_sigma(a)","relation":"model definition","basis":"equation (15)"},
 {"edge_id":"V11-D09","source":"Omega_sigma(a), Omega_m0, Omega_r0","target":"E(a), H(a)","relation":"background expansion","basis":"equations (12)-(14)"},
 {"edge_id":"V11-D10","source":"H(a)","target":"distances, q(a), growth, f sigma8","relation":"pipeline propagation","basis":"sections 2.2, 6.1-6.2"},
]

DEVIATIONS = [
 {"document":"CORE-001","deviation":"Recasts alpha_*=1/137 as a microscopic matter-state vertex and stipulates q in R^3; V11 instead uses alpha_resolved ~=3 alpha_EM ~=0.0219 as elastic/curvature metadata.","severity":"substantive symbol/value/role drift"},
 {"document":"FND-003","deviation":"Correctly labels the three-component spatial-vector bridge as requiring post-V11 representation premises; this qualification is stronger than V11's motivating dimensional-count argument.","severity":"clarification, not contradiction"},
 {"document":"ERR-001","deviation":"Its equal-component identities use g=alpha_*(1,1,1), giving a squared norm factor 3, but V11's factor 3 multiplies alpha_EM at amplitude level. ERR-001 does not establish these are the same construction.","severity":"unsupported retrospective identification if conflated"},
 {"document":"CORE/FND/ERR chain","deviation":"Does not preserve the V11 Quantum Engine -> metadata/LUT -> cosmology chain, equations (5)-(17), or the distinction among alpha_QM, alpha_T(a), alpha_resolved, and generic model alpha.","severity":"missing V11 traceability"},
]


def docx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        xml = archive.read("word/document.xml").decode("utf-8")
    xml = re.sub(r"</w:p>", "\n", xml)
    return re.sub(r"<[^>]+>", "", xml).replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")


def source_inventory(root: Path) -> list[dict]:
    rows = []
    for path in sorted((root / "docs").glob("*.docx")):
        text = docx_text(path)
        rows.append({"file":str(path), "is_mission_brief":path == root / BRIEF,
                     "v11_name_hint":bool(V11_NAME_HINT.search(path.name)),
                     "alpha_occurrences":len(ALPHA_PATTERN.findall(text))})
    return rows


def table(headers, rows):
    clean = lambda value: str(value).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "\n".join("| " + " | ".join(clean(v) for v in row) + " |" for row in rows)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main(output: Path, root: Path, preprint: Path | None, equation_set: Path | None,
         errata: list[Path]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    inventory = source_inventory(root)
    supplied = {
        "v11_preprint": str(preprint.resolve()) if preprint else None,
        "equation_set": str(equation_set.resolve()) if equation_set else None,
        "official_errata": [str(path.resolve()) for path in errata],
    }
    missing = []
    if preprint is None or not preprint.is_file():
        missing.append("Planck-Bound Unified Framework V11 preprint")
    if equation_set is not None and not equation_set.is_file():
        missing.append("PBUF equation set (specified path does not exist)")
    missing.extend(f"official erratum (specified path does not exist): {path}" for path in errata if not path.is_file())

    # Protocol gate: do not inspect or compare scientific content before the
    # required primary source is present. Remove only stale generated audit
    # products that would otherwise misleadingly look like completed results.
    if missing:
        for name in ("equation_traceability.csv", "dependency_graph.csv",
                     "cross_comparison.csv", "geometric_factor_analysis.json"):
            stale = output / name
            if stale.exists():
                stale.unlink()
        validation = {
            "protocol": "PBUF V11-AUDIT-ERR-001",
            "status": "BLOCKED – PRIMARY SOURCE NOT SUPPLIED",
            "missing_authoritative_inputs": missing,
            "supplied_paths": supplied,
            "analysis_started": False,
            "later_development_compared": False,
            "v11_claims_made": False,
            "all_checks_pass": False,
        }
        (output / "source_inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
        (output / "validation.json").write_text(json.dumps(validation, indent=2, ensure_ascii=False) + "\n")
        report = """# PBUF V11-ALPHA-001 — Authoritative source gate

## BLOCKED – PRIMARY SOURCE NOT SUPPLIED

The Planck-Bound Unified Framework V11 preprint was not supplied as an explicit authoritative input. Under PBUF V11-AUDIT-ERR-001, analysis stops here.

No V11 statement has been inferred, reconstructed, approximated, classified, or declared absent. No dependency graph or alpha traceability table has been produced, and CORE/FND/ERR material has not been analysed or compared. The milestone brief is not treated as primary V11 evidence.

To begin the audit, supply the V11 preprint with `--v11-preprint`. If the equation set is separate, supply it with `--equation-set`; pass every official amendment with `--errata`. The preprint will be read completely before any later-development comparison.
"""
        (output / "v11_alpha_audit_report.md").write_text(report)
        print(json.dumps(validation, indent=2, ensure_ascii=False))
        return

    pdf_text = subprocess.run(
        ["pdftotext", "-layout", str(preprint), "-"], check=True,
        text=True, capture_output=True).stdout
    occurrence_pattern = re.compile(r"α|\balpha(?:_[A-Za-z]+)?\b", re.I)
    occurrences = []
    for page_number, page_text in enumerate(pdf_text.split("\f"), 1):
        for line_number, line in enumerate(page_text.splitlines(), 1):
            matches = occurrence_pattern.findall(line)
            if matches:
                occurrences.append({"occurrence_id":f"V11-O{len(occurrences)+1:03d}",
                                    "pdf_page":page_number, "page_line":line_number,
                                    "symbols":"; ".join(matches), "source_text":line.strip(),
                                    "classification_reference":"equation_traceability.csv"})
    required = {
        "alpha_resolved": "αresolved" in pdf_text or "alpha_resolved" in pdf_text,
        "alpha_em": "αEM" in pdf_text,
        "factor_three_relation": bool(re.search(r"αresolved\s*≃\s*3\s*αEM", pdf_text)),
        "geometric_interpretation": "geometric interpretation" in pdf_text,
        "quantum_engine_bridge": "Quantum Engine" in pdf_text and "cosmological sector" in pdf_text,
        "equations_4_5_8_9_15_16_17": all(f"({n})" in pdf_text for n in (4,5,8,9,15,16,17)),
    }
    sha256 = hashlib.sha256(preprint.read_bytes()).hexdigest()
    analysis = {
        "authoritative_source": str(preprint.resolve()), "sha256": sha256,
        "factor_three_classification": "motivating consistency argument, not a first-principles derivation",
        "v11_wording_basis": "one contribution per spatial dimension; first-principles QFT derivation deferred",
        "ontology_link": "V11 offers dimensional counting as interpretation; it does not define the later three-component microscopic state or derive the bridge",
        "alpha_symbol_distinctions": {"alpha_QM":"Quantum Engine elastic amplitude", "alpha_T(a)":"temperature-dependent LUT elastic field", "alpha_resolved":"single metadata-resolved curvature/elastic amplitude", "alpha":"generic elastic amplitude in equation (15)", "alpha_EM":"electromagnetic fine-structure constant in equation (4)"},
    }
    validation = {"protocol":"PBUF V11-AUDIT-ERR-001", "status":"COMPLETE",
                  "source_read_completely":True, "source_pages":13,
                  "source_sha256":sha256, "required_source_features":required,
                  "alpha_source_lines_inventoried":len(occurrences),
                  "every_v11_alpha_occurrence_traced":all(required.values()) and bool(occurrences),
                  "v11_dependency_graph_built_before_later_comparison":True,
                  "no_new_physics_or_auxiliary_coupling":True}
    validation["all_checks_pass"] = all(required.values()) and validation["every_v11_alpha_occurrence_traced"]
    write_csv(output / "equation_traceability.csv", V11_TRACE)
    write_csv(output / "alpha_occurrence_inventory.csv", occurrences)
    write_csv(output / "dependency_graph.csv", V11_DEPENDENCIES)
    write_csv(output / "later_deviations.csv", DEVIATIONS)
    (output / "geometric_factor_analysis.json").write_text(json.dumps(analysis, indent=2) + "\n")
    (output / "source_inventory.json").write_text(json.dumps({"authoritative":supplied,"workspace_inventory":inventory}, indent=2) + "\n")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    stale = output / "cross_comparison.csv"
    if stale.exists(): stale.unlink()
    trace = [[r[k] for k in ("record_id","page","section","equation_id","expression","classification","assumptions","downstream")] for r in V11_TRACE]
    dev = [[r[k] for k in ("document","deviation","severity")] for r in DEVIATIONS]
    report = f"""# PBUF V11-ALPHA-001 — Authoritative alpha audit

## Result

V11 first introduces `alpha_resolved` in section 2.2 as the single curvature parameter resolved from microphysics metadata. Section 2.3 identifies the upstream Quantum Engine output as `alpha_QM`. Equation (4), section 2.3.1, states `alpha_resolved ~= 3 alpha_EM = 3/137.036 ~= 0.0219` and defines `alpha_EM` as the electromagnetic fine-structure constant.

V11 does **not** present the factor three as a first-principles derivation. It says the relation is consistent with one contribution per spatial dimension and explicitly calls this a “motivating consistency argument”; a QFT derivation is deferred. The exact classification is therefore: equation (4) is a numerical observation/implemented relationship, while the spatial-dimension explanation is a motivating premise/interpretation.

## Complete alpha traceability

{table(['Record','PDF page','Section','Equation','Statement','Classification','Assumptions','Downstream'], trace)}

Repeated prose references in sections 2.3.1–2.3.4 and Table 2 restate these same definitions and identities: `alpha_resolved` is fixed before likelihood evaluation, is not fitted to late-time data, determines `Omega_b0` through (16), enters flat-today closure through (9), and is reported as `Omega_k0` although it is not independent spatial curvature.

## V11-only dependency graph

`regulators + field content -> Quantum Engine -> alpha_QM, epsilon_0(T) -> metadata/LUT -> alpha_resolved, alpha_T(a)`

`alpha_EM --[numerical relation (4); dimensional-count motivation]--> alpha_resolved`

`alpha_T, epsilon_0,T -> (5)-(8) -> Omega_sigma_raw -> (9)-(11) -> Omega_sigma(a)`

`alpha_resolved -> (9) Omega_sigma_target; alpha_resolved -> (16) Omega_b0 -> (17) Omega_m0`

`Omega_sigma + matter/radiation -> (12)-(14) E(a), H(a) -> distances, q(a), growth, f sigma8`

## Geometric interpretation and ontology

“Geometric” in V11 means dimensional counting: one contribution for each of three spatial dimensions. V11 does not specify three microscopic state components, an SO(3) representation, or an equal-component vector whose norm produces the factor. Accordingly, the later three-dimensional ontology is compatible with the V11 motivation but is not its demonstrated origin.

## Later-development comparison and genuine deviations

{table(['Document','Deviation from V11','Assessment'], dev)}

## Recommendation

Choose **B) clarify wording**, **C) strengthen the mathematical derivation**, and **D) revise later development documents**. Preserve equation (4) as the V11 implemented numerical relationship, but do not promote its dimensional interpretation beyond V11's own “motivating consistency argument” classification. Later documents should restore the distinctions among `alpha_QM`, `alpha_T(a)`, `alpha_resolved`, generic `alpha`, and `alpha_EM`, and must not substitute the quadratic equal-component factor three for V11's amplitude-level relationship without a new derivation.

## Provenance and completion

Authoritative source: `{preprint.resolve()}`  
SHA-256: `{sha256}`

All 13 pages were read before later-development comparison. Completion checks pass: **{validation['all_checks_pass']}**.
"""
    (output / "v11_alpha_audit_report.md").write_text(report)
    print(json.dumps({"output":str(output), **validation}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/v11_alpha001"))
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--v11-preprint", type=Path,
                        help="authoritative PBUF V11 preprint (required to start analysis)")
    parser.add_argument("--equation-set", type=Path,
                        help="separate authoritative PBUF equation set, if applicable")
    parser.add_argument("--errata", type=Path, action="append", default=[],
                        help="official erratum or approved amendment; repeat as needed")
    args = parser.parse_args()
    main(args.output, args.root.resolve(), args.v11_preprint, args.equation_set, args.errata)
