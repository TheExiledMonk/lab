from pbuf.audit.native_equation_registry import equations
def test_frozen_mechanisms_present():
    names=" ".join(x["name"] for x in equations()); assert "excitation norm" in names and "frame" in names and "equilibrium" in names
