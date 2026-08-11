"""Represent a projected detector field as a relative native 2D constraint."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Native2DSourceConstraint:
    amplitude: np.ndarray
    uncertainty: np.ndarray
    support: np.ndarray
    channels: tuple[str, ...]
    normalization: str = "L1_RELATIVE_PER_CHANNEL"
    depth_assigned: bool = False
    physical_mass_scale: bool = False

    def __post_init__(self):
        if self.amplitude.ndim != 3 or self.amplitude.shape != self.uncertainty.shape:
            raise ValueError("amplitude and uncertainty must be channel,y,x arrays")
        if self.support.shape != self.amplitude.shape or len(self.channels) != self.amplitude.shape[0]:
            raise ValueError("support/channels do not match amplitude")
        if self.depth_assigned or self.physical_mass_scale:
            raise ValueError("this bridge permits neither depth nor a physical mass scale")


def native_2d_constraint(images, uncertainties, channels):
    amplitudes=[]; sigmas=[]; supports=[]
    for image, uncertainty in zip(images, uncertainties):
        support=np.isfinite(image)&np.isfinite(uncertainty)&(uncertainty>0)&(image>0)
        positive=np.where(support,image,0.0); scale=float(np.sum(positive))
        if not scale > 0: raise ValueError("channel has no positive background-subtracted support")
        amplitudes.append(positive/scale); sigmas.append(np.where(support,uncertainty/scale,np.inf)); supports.append(support)
    return Native2DSourceConstraint(np.stack(amplitudes),np.stack(sigmas),np.stack(supports),tuple(channels))


def support_diagnostics(field):
    rows=[]
    for channel, amplitude, support in zip(field.channels,field.amplitude,field.support):
        yy,xx=np.indices(amplitude.shape); total=float(amplitude.sum())
        cx=float((amplitude*xx).sum()/total); cy=float((amplitude*yy).sum()/total)
        dx=xx-cx;dy=yy-cy; cov=np.array([[(amplitude*dx*dx).sum(),(amplitude*dx*dy).sum()],[(amplitude*dx*dy).sum(),(amplitude*dy*dy).sum()]])/total
        eig=np.linalg.eigvalsh(cov)[::-1]; peaks=np.argpartition(amplitude.ravel(),-min(5,amplitude.size))[-5:]
        rows.append({"channel":channel,"support_pixels":int(support.sum()),"support_fraction":float(support.mean()),
          "centroid_xy":[cx,cy],"rms_radius":float(np.sqrt(np.trace(cov))),"principal_axis_variances":eig.tolist(),
          "anisotropy":float((eig[0]-eig[1])/(eig.sum()+np.finfo(float).eps)),
          "peak_locations_xy":[[int(p%amplitude.shape[1]),int(p//amplitude.shape[1])] for p in peaks[np.argsort(amplitude.ravel()[peaks])[::-1]]],
          "connected_components":"NOT_COMPUTED_WITHOUT_MORPHOLOGY_THRESHOLD"})
    return rows
