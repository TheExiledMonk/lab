from pbuf.wl.native_spatial_wave_evolution import scale_cancellation

def test_scale_free_integral_exactly_cancels_L0():
    for a in (.5,1,2,4):
        x=scale_cancellation(.25,.75,a); assert x["native"]==x["physical"]
