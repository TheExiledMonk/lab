from pbuf.audit.native_equation_registry import equations
def test_equations_are_unique_and_provenanced():
    e=equations(); assert len({x["equation_id"] for x in e})==len(e)
    assert all(x["code_location"] and x["latex"] and x["status"] for x in e)
