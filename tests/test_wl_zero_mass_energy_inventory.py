from pbuf.wl.native_zero_mass_energy import energy_state_inventory,momentum_state_inventory,candidate_registry,energy_driver_registry

def test_inventory_keeps_medium_and_mode_energy_distinct():
    rows=energy_state_inventory()
    assert any(r["name"].startswith("bounded-strain") and r["status"]=="RELATION_ONLY" for r in rows)
    assert any(r["name"]=="zero-mass mode energy" and r["status"]=="MISSING_NATIVE_STATE" for r in rows)

def test_all_candidate_families_registered():
    assert [r["candidate_id"] for r in candidate_registry()]==[f"Q{i:02d}" for i in range(1,36)]
    assert [r["candidate_id"] for r in energy_driver_registry()]==[f"E{i:02d}" for i in range(1,17)]
    assert momentum_state_inventory()[0]["available"] and not momentum_state_inventory()[1]["available"]
