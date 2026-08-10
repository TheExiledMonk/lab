import numpy as np
from pbuf.wl.native_incremental_elastic_energy import incremental_elastic_energy,integrate_packet

def test_spreading_quadratic_packet_integral():
    x=np.linspace(-20,20,20001)
    vals=[]
    for width in (1,2,4):
        # L2-normalized shapes model redistribution, not a new propagation law.
        de=.03/np.sqrt(width)*np.exp(-.5*(x/width)**2)
        vals.append(integrate_packet(incremental_elastic_energy(0,de),cell_volume=x[1]-x[0])["signed"])
    assert np.std(vals)/np.mean(vals)<5e-4

