from pbuf.foundation.native_loaded_link_response import execute
def test_full_family_norm_conservation():
    for t in range(1,21):
      for l in range(11):
       r=execute(f"T{t:02d}",f"LOAD{l:02d}","EX04"); assert abs(r["norm_out"]-r["norm_in"])<1e-11
