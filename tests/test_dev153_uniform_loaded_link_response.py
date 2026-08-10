from pbuf.foundation.native_loaded_link_response import execute
def test_uniform_load_does_not_create_unlicensed_response():
    for i in range(1,21): assert execute(f"T{i:02d}","LOAD03","EX01")["progression_ratio"]==execute(f"T{i:02d}","LOAD00","EX01")["progression_ratio"]
