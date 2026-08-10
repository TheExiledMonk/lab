from pbuf.foundation.native_loaded_link_response import execute
def test_zero_load_all_candidates():
    for i in range(1,21):
        r=execute(f"T{i:02d}","LOAD00","EX03"); assert abs(r["norm_out"]-r["norm_in"])<1e-12 and r["transfer_fraction"]==1
