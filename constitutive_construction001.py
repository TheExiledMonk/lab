"""Build the audit artifacts for PBUF CONSTITUTIVE-CONSTRUCTION-001.

This is a symbolic construction audit.  It deliberately leaves the frozen
weak-field tangent and the previously unselected admissible domain symbolic.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs/constitutive_construction001"

SOURCES = {
    "FOUNDATION-001": "runs/foundation001/foundational_ontology.md",
    "STATE-002": "runs/state002/primitive_medium_state.md",
    "DEFORMATION-001": "runs/deformation001/deformation_measure_report.md",
    "HYPER-001": "runs/hyper001/stored_energy_derivation.md",
    "ENERGY-PRINCIPLE-001": "runs/energy_principle001/energy_selection_derivation.md",
    "DURATION-001": "runs/duration001/emergent_duration_derivation.md",
    "METRIC-001": "runs/metric001/effective_metric_derivation.md",
    "BALANCE-001": "runs/balance001/native_balance_laws.md",
    "LOCALITY-001": "runs/locality001/locality_report.md",
    "ENERGY-SEARCH-001": "runs/energy_search001/energy_search_report.md",
}

CHECKLIST = [
    ("R01", "one scalar stored-energy law for the one medium; no independent sector or freely adjustable constant"),
    ("R02", "state dependence only on the frozen objective rank-three SPD tensor C[q,q0]"),
    ("R03", "dimensionless-deformation dependence invariant under orthogonal similarity; equivalently Phi(I1,I2,I3)"),
    ("R04", "local, parity-even, isotropic, single-valued, statewise hyperelastic and rate/history independent"),
    ("R05", "D_C is path-connected and permutation invariant, contains 1, and has compact closure inside SPD"),
    ("R06", "W is C1 in int(D_C), C2 near 1, and lower semicontinuous when extended-valued"),
    ("R07", "W(C)>=0 and W(1)=0"),
    ("R08", "DW(1)=0"),
    ("R09", "D2W(1)=A0, the frozen strictly positive isotropic weak-field tangent (K0>0, mu0>0)"),
    ("R10", "stress is the variational derivative and is acoustically positive on the declared propagation domain"),
    ("R11", "one authorized endpoint: hard extended value, complete interior blow-up, or finite regular endpoint plus the frozen state constraint"),
    ("R12", "the operational branch remains compatible with a regular V11 metric/cone completion; W alone does not prove that completion"),
    ("R13", "no intrinsic gradient, kernel, hidden state, dissipation, fitted coefficient, or new ontology"),
]

CANDIDATES = [
    {
        "candidate": "A",
        "definition": "W_A=Q on D_C and +infinity outside cl(D_C)",
        "endpoint": "hard extended-value",
        "terms": 2,
        "nonlinear_operations": 0,
        "extra_assumptions": 0,
        "regularity": "analytic interior; lower-semicontinuous extension when D_C is taken closed at the constraint",
        "status": "unconditionally constructible from frozen symbolic data",
    },
    {
        "candidate": "B",
        "definition": "W_B=Q on cl(D_C), with continuation forbidden by the independently frozen state constraint",
        "endpoint": "finite regular constrained endpoint",
        "terms": 2,
        "nonlinear_operations": 0,
        "extra_assumptions": 0,
        "regularity": "analytic wherever Q is evaluated; finite one-sided endpoint",
        "status": "unconditionally constructible from frozen symbolic data",
    },
    {
        "candidate": "C",
        "definition": "W_C=Q+(Q^2/K0)b, where b>=0 is an invariant smooth complete boundary barrier",
        "endpoint": "complete interior barrier",
        "terms": 3,
        "nonlinear_operations": 1,
        "extra_assumptions": 1,
        "regularity": "C1 interior and C2 near reference, conditional on the chosen barrier remainder",
        "status": "class construction; no canonical member exists because the frozen domain and barrier profile are unselected",
    },
]


def main() -> None:
    missing = [v for v in SOURCES.values() if not (ROOT / v).is_file()]
    if missing:
        raise FileNotFoundError(missing)
    OUT.mkdir(parents=True, exist_ok=True)

    with (OUT / "frozen_requirement_checklist.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "frozen_requirement"])
        writer.writerows(CHECKLIST)
    with (OUT / "comparative_complexity.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CANDIDATES[0].keys())
        writer.writeheader()
        writer.writerows(CANDIDATES)

    audit = []
    for c in CANDIDATES:
        for rid, _ in CHECKLIST:
            result = "pass"
            note = "by construction"
            if c["candidate"] == "C" and rid in {"R06", "R10"}:
                result = "conditional"
                note = "requires an admissible zero-2-jet complete barrier remainder and an elliptic declared propagation domain"
            if rid == "R12":
                result = "gate open"
                note = "full metric compatibility cannot be established by a stored energy alone"
            audit.append({"candidate": c["candidate"], "requirement": rid, "result": result, "note": note})
    with (OUT / "frozen_compatibility_audit.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=audit[0].keys())
        writer.writeheader()
        writer.writerows(audit)

    result = {
        "milestone": "PBUF CONSTITUTIVE-CONSTRUCTION-001",
        "pass": True,
        "minimal_interior_energy": "Q(C)=1/2 A0[E,E], E=(C-1)/2 = K0/2 tr(E)^2 + mu0 E_TF:E_TF",
        "minimal_candidates": [c["candidate"] for c in CANDIDATES],
        "lexicographic_winner": "A and B tie in interior formula complexity; A is preferred for self-contained enforcement of capacity",
        "unique_minimum": False,
        "reason": "A and B have the identical minimal interior polynomial and differ only in an already-authorized endpoint semantics; no frozen ordering ranks those semantics.",
        "barrier_member_selected": False,
        "sources": SOURCES,
    }
    (OUT / "validation.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
