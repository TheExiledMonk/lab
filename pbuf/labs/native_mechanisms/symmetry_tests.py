import numpy as np
def permutation_covariant(rule, relations):
    p=np.array([2,3,4,5,0,1]); return bool(np.allclose(rule(np.asarray(relations)[p]),rule(relations)[p]))
