"""Validation shared by the registry CLI and build tool."""
from __future__ import annotations

import datetime as _dt

RESULTS = {"FAILED", "PARTIAL", "FULL", "BLOCKED", "NOT_EVALUATED", "DIAGNOSTIC_ONLY"}
STATUSES = {"CANONICAL", "ACTIVE", "OPEN", "PARTIAL", "SUPERSEDED", "REJECTED", "HISTORICAL_ONLY", "BLOCKED", "INFRASTRUCTURE", "DIAGNOSTIC", "FORBIDDEN"}


def validate(registry: dict) -> list[str]:
    errors = []
    if registry.get("registry_schema_version") != 1:
        errors.append("registry_schema_version must be 1")
    targets = registry.get("targets", [])
    attempts = registry.get("attempts", [])
    tids = [x.get("target_id") for x in targets]
    aids = [x.get("attempt_id") for x in attempts]
    if len(tids) != len(set(tids)): errors.append("duplicate target_id")
    if len(aids) != len(set(aids)): errors.append("duplicate attempt_id")
    tset, aset = set(tids), set(aids)
    for target in targets:
        for key in ("canonical_name", "plain_language_question", "aliases", "keywords", "domain", "current_status", "do_not_rederive", "reopen_condition"):
            if key not in target: errors.append(f"target {target.get('target_id')} missing {key}")
        if target.get("current_status") not in STATUSES: errors.append(f"target {target.get('target_id')} invalid current_status")
        for aid in target.get("attempt_ids", []):
            if aid not in aset: errors.append(f"target {target.get('target_id')} dangling attempt {aid}")
        for aid in target.get("canonical_solution_ids", []):
            if aid not in aset: errors.append(f"target {target.get('target_id')} dangling canonical {aid}")
    for attempt in attempts:
        aid = attempt.get("attempt_id")
        for key in ("target_id", "name", "why_attempted", "result", "result_reason", "current_status", "equations", "evidence", "confidence", "reopen_condition", "do_not_repeat_reason"):
            if key not in attempt: errors.append(f"attempt {aid} missing {key}")
        if attempt.get("target_id") not in tset: errors.append(f"attempt {aid} dangling target")
        if attempt.get("result") not in RESULTS: errors.append(f"attempt {aid} invalid result")
        if attempt.get("current_status") not in STATUSES: errors.append(f"attempt {aid} invalid current_status")
        for key in ("superseded_by", "supersedes", "equivalent_to", "derived_from", "ancestor_of", "descendant_of", "related_attempts"):
            for ref in attempt.get(key, []):
                if ref not in aset: errors.append(f"attempt {aid} dangling {key}: {ref}")
        for key in ("date_started", "date_completed"):
            value = attempt.get(key)
            if value:
                try: _dt.date.fromisoformat(value)
                except ValueError: errors.append(f"attempt {aid} invalid {key}")
    for relation in registry.get("equivalences", []):
        if relation.get("source") not in aset or relation.get("target") not in aset:
            errors.append(f"equivalence {relation.get('relation_id')} has dangling endpoint")
    # A selector is a development authorization, not merely narrative metadata.
    # Keep this generic so future gates cannot silently select their own blocked
    # operation.
    for selector in registry.get("development_gate_selectors", []):
        for key in ("gate", "gate_value", "blocked_values", "selected_test", "blocked_operations"):
            if key not in selector:
                errors.append(f"development gate selector missing {key}")
        if (selector.get("gate_value") in selector.get("blocked_values", [])
                and selector.get("selected_test") in selector.get("blocked_operations", [])):
            errors.append(
                f"blocked operation selected: {selector.get('gate')}="
                f"{selector.get('gate_value')} selects {selector.get('selected_test')}"
            )
    return errors
