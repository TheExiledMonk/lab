from pathlib import Path
def test_generated_parameter_inventory_has_active_and_retired():
    import json
    p=Path(__file__).parents[1]/"runs/native_microphysics_reconstruction001/parameter_inventory.json"
    if p.exists():
      rows=json.loads(p.read_text()); assert any(x["status"]=="active" for x in rows) and any("retired" in x["status"] for x in rows)
