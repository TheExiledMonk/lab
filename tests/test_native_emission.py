from pbuf.quantum.native_emission_absorption import emission_audit
def test_no_packet_is_injected_or_claimed():
    r=emission_audit(); assert not r['emitted_packet_injected'] and not r['emission_established']
