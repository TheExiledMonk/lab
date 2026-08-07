"""Dependency-free structural checks for the HYPER-001 deliverables."""

import csv
import json
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main() -> None:
    with (ROOT / "equation_catalogue.json").open() as stream:
        catalogue = json.load(stream)
    with (ROOT / "validation.json").open() as stream:
        validation = json.load(stream)

    ids = [item["id"] for item in catalogue["equations"]]
    assert len(ids) == len(set(ids)) == 23
    assert ids == [f"H-{index:03d}" for index in range(1, 24)]
    assert validation["checks"]["primitive_rank"] == 3

    with (ROOT / "equation_traceability.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert [row["equation_id"] for row in rows] == ids

    # H-011: for diagonal epsilon, the omitted determinant term is cubic.
    for _ in range(100):
        x, y, z = (random.uniform(-0.05, 0.05) for _ in range(3))
        c = (1 + 2 * x, 1 + 2 * y, 1 + 2 * z)
        t = x + y + z
        s2 = x * x + y * y + z * z
        i1 = sum(c)
        i2 = sum(c[a] * c[b] for a in range(3) for b in range(a + 1, 3))
        i3 = c[0] * c[1] * c[2]
        assert abs(i1 - (3 + 2 * t)) < 1e-12
        assert abs(i2 - (3 + 4 * t + 2 * (t * t - s2))) < 1e-12
        assert abs(i3 - (1 + 2 * t + 2 * (t * t - s2)) - 8 * x * y * z) < 1e-12

    print("HYPER-001 catalogue, traceability, and H-011 checks passed")


if __name__ == "__main__":
    main()
