from pbuf.excitation.native_excitation_n6 import execute_operator,gaussian_packet
def test_dev148_candidates_are_not_registry_only():
    x=gaussian_packet((8,8,8))
    for i in range(1,11):
      y=execute_operator(f"O{i:02d}",x); assert isinstance(y,dict) or hasattr(y,"shape")
