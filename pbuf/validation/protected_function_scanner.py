"""Protected-function scanner.

A source-tree scanner that fails if a lab defines a protected function
locally. Labs may import them only.
"""
from __future__ import annotations
import ast
import json
from pathlib import Path

__all__ = ["scan_protected_functions", "ProtectedFunctionViolation"]


class ProtectedFunctionViolation(Exception):
    pass


def _load_protected_functions(path: Path) -> list:
    with open(path) as f:
        data = json.load(f)
    return list(data["protected_functions"])


def scan_protected_functions(labs_dir: Path, registry_path: Path):
    """Scan every Python file under ``labs_dir`` for local definitions
    of protected functions.

    Returns a list of violations; raises ProtectedFunctionViolation if
    any violation is found.
    """
    protected = _load_protected_functions(registry_path)
    violations = []
    for path in sorted(labs_dir.rglob("*.py")):
        try:
            source = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if node.name in protected:
                    violations.append({
                        "file": str(path),
                        "function": node.name,
                        "line": node.lineno,
                    })
    return violations


def scan_for_file(labs_file: Path, registry_path: Path):
    """Scan a single lab file."""
    protected = _load_protected_functions(registry_path)
    violations = []
    try:
        source = labs_file.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return violations
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return violations
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name in protected:
                violations.append({
                    "file": str(labs_file),
                    "function": node.name,
                    "line": node.lineno,
                })
    return violations


if __name__ == "__main__":
    import sys
    repo_root = Path(__file__).resolve().parents[2]
    labs_dir = repo_root / "pbuf" / "labs"
    registry = Path(__file__).resolve().parent / "protected_functions.json"
    v = scan_protected_functions(labs_dir, registry)
    if v:
        print(f"VIOLATIONS: {len(v)}")
        for vio in v:
            print(f"  {vio['file']}:{vio['line']}: {vio['function']}")
        sys.exit(1)
    else:
        print("No protected-function violations found")