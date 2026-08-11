"""Case-insensitive semantic registry lookup."""
from __future__ import annotations
import json
from pathlib import Path


def load(path=Path("docs/PBUF_MECHANISM_REGISTRY.json")):
    return json.loads(Path(path).read_text())


def haystack(record: dict) -> str:
    values = []
    def visit(x):
        if isinstance(x, dict):
            for v in x.values(): visit(v)
        elif isinstance(x, list):
            for v in x: visit(v)
        elif x is not None: values.append(str(x))
    visit(record)
    return " ".join(values).lower()


def search(registry: dict, query: str):
    phrase = query.lower().strip()
    words = phrase.split()
    # A known historical phrase/alias is stronger evidence than coincidental
    # occurrence of individual words in a long record.  Fall back to token
    # matching only when the phrase has no hit at all.
    targets = [x for x in registry['targets'] if phrase in haystack(x)]
    if not targets:
        targets = [x for x in registry['targets'] if all(w in haystack(x) for w in words)]
    target_ids = {x['target_id'] for x in targets}
    attempts = [x for x in registry['attempts'] if phrase in haystack(x) or x['target_id'] in target_ids]
    if not attempts:
        attempts = [x for x in registry['attempts'] if all(w in haystack(x) for w in words)]
    return targets, attempts
