from pbuf.quantum.native_emission_absorption import fractional_incident_controls
def test_required_fractional_controls():
    rows=fractional_incident_controls(4); assert [r['incident_norm'] for r in rows]==[1,2,3,4,5,6,8]
