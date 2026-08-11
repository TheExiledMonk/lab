import numpy as np
def allocate(relations):
    w=1+np.abs(np.asarray(relations,float)); return w/w.sum()
