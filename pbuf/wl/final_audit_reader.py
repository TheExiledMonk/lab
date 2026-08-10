"""Read and validate the frozen Dev Doc 110/112 root audit log."""

from __future__ import annotations

import json
from pathlib import Path


class FinalAuditError(ValueError):
    """Raised when the frozen audit evidence is absent or incomplete."""


def load_final_audit(path: str | Path = "final.log") -> dict:
    """Return the final structured ``RESULT_JSON`` object from *path*.

    The file is opened read-only as UTF-8.  Logs may contain more than one
    marker; the last complete object is the authoritative final result.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"DEV113_REQUIRED_FINAL_LOG_MISSING: {source}")
    text = source.read_text(encoding="utf-8")
    lines = text.splitlines()
    parsed: list[dict] = []
    for index, line in enumerate(lines):
        if line.strip() != "RESULT_JSON":
            continue
        for candidate in lines[index + 1 :]:
            if not candidate.strip():
                continue
            try:
                value = json.loads(candidate)
            except json.JSONDecodeError as exc:
                raise FinalAuditError(
                    f"malformed RESULT_JSON after line {index + 1}: {exc}"
                ) from exc
            if not isinstance(value, dict):
                raise FinalAuditError("RESULT_JSON must contain a JSON object")
            parsed.append(value)
            break
    if not parsed:
        raise FinalAuditError("final.log has no complete RESULT_JSON block")
    return parsed[-1]
