import numpy as np
def allocation_conserved(before, after, atol=1e-12): return bool(np.isclose(np.sum(before),np.sum(after),atol=atol))
