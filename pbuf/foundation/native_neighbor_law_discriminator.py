"""Hard-gate scoring and physical-equivalence closure for Dev152."""
from __future__ import annotations

HARD_GATES = ("STATIC_PARITY", "DYNAMIC_PARITY", "BASIS_COVARIANT", "SYMMETRY_CONSISTENT",
              "RESOLUTION_CONVERGED", "PROGRESSION_STEP_CONVERGED")


def evaluate(survivor, evidence: dict) -> dict:
    gates = {k: bool(evidence.get(k, False)) for k in HARD_GATES}
    gates["NEW_INTERACTION_COEFFICIENTS"] = int(evidence.get("NEW_INTERACTION_COEFFICIENTS", 0))
    viable = all(gates[k] for k in HARD_GATES) and gates["NEW_INTERACTION_COEFFICIENTS"] == 0
    scores = {f"S{i}": evidence.get(f"S{i}", None) for i in range(1, 13)}
    return {"survivor_id": survivor.survivor_id, "hard_gates": gates, "viable": viable, "scores": scores}


def collapse_equivalence(rows: list[dict], transformation_signatures: dict[str, str]):
    classes = {}
    for row in rows:
        if row["viable"]:
            signature = transformation_signatures[row["survivor_id"]]
            classes.setdefault(signature, []).append(row["survivor_id"])
    return [{"class_id": f"E{i:02d}", "signature": sig, "members": members,
             "invertible_mapping_proven": len(members) > 1}
            for i, (sig, members) in enumerate(sorted(classes.items()), 1)]


def decide(rows, classes):
    viable = [r for r in rows if r["viable"]]
    if not viable:
        return {"outcome": "PBUF_DEV151_PARITY_SURVIVORS_FAIL_MIXED_STATE_PHYSICS", "unique": False}
    if len(classes) == 1 and len(classes[0]["members"]) > 1 and classes[0]["invertible_mapping_proven"]:
        return {"outcome": "PBUF_NATIVE_NEIGHBOR_CONSTITUTIVE_EQUIVALENCE_CLASS_ESTABLISHED",
                "unique": False, "equivalence_class": classes[0]}
    if len(classes) == 1 and len(classes[0]["members"]) == 1:
        return {"outcome": "PBUF_UNIQUE_NATIVE_NEIGHBOR_CONSTITUTIVE_LAW_SELECTED",
                "unique": True, "selected_law_id": classes[0]["members"][0]}
    return {"outcome": "PBUF_NATIVE_NEIGHBOR_LAW_PHYSICAL_DEGENERACY_REMAINS", "unique": False,
            "remaining_classes": classes}
