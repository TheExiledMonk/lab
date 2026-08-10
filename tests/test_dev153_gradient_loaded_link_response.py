from pbuf.foundation.native_loaded_link_response import execute
def test_gradient_does_not_create_unlicensed_response():
    for i in range(1,21): assert execute(f"T{i:02d}","LOAD07","EX06")["transfer_fraction"]==1
