"""N6 fixtures for Dev165."""
import numpy as np

def directed_differences(field, index=None):
    q=np.asarray(field,float); index=index or tuple(s//2 for s in q.shape)
    z,y,x=index; c=q[z,y,x]
    return np.array([q[z,y,min(x+1,q.shape[2]-1)]-c,q[z,y,max(x-1,0)]-c,
      q[z,min(y+1,q.shape[1]-1),x]-c,q[z,max(y-1,0),x]-c,
      q[min(z+1,q.shape[0]-1),y,x]-c,q[max(z-1,0),y,x]-c])

def asymmetric_sample(field):
    q=np.asarray(field,float); best=None
    for z in range(1,q.shape[0]-1):
      for y in range(1,q.shape[1]-1):
       for x in range(1,q.shape[2]-1):
        d=directed_differences(q,(z,y,x)); score=float(np.std(np.abs(d)))
        if best is None or score>best[0]: best=(score,(z,y,x),d)
    return best[1],best[2]
