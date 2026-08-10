from pbuf.wl.native_source_depth import scale_free_ratios
def test_coordinate_rescaling_invariance():
 a=scale_free_ratios(0,2,5,.5,2)
 for s in (.5,1,2,4):
  b=scale_free_ratios(0,2*s,5*s,.5*s,2*s)
  for k in ("D_LS_over_D_OL","D_OS_over_D_OL","D_LS_over_D_OS","R_source_over_R_lens"):assert a[k]==b[k]
