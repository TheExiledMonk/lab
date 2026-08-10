from pbuf.quantum.native_emission_absorption import selection_audit
def test_basis_and_spin_guards():
    r=selection_audit(); assert not r['basis_dependence_introduced'] and not r['handedness_labelled_spin']
