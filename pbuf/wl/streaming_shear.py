"""Bounded-memory sufficient statistics for the frozen Jacobian/TSC shear readout."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


MOMENT_NAMES = ("count", "mean_u0", "mean_v0", "mean_uf", "mean_vf",
                "cxx00", "cxx01", "cxx11", "cxy00", "cxy01", "cxy10", "cxy11")


@dataclass
class JacobianTSCAccumulator:
    """Incremental equivalent of construct_local_primitives + TSC deposition.

    ``transfer[target, source]`` holds the accumulated TSC weight from rays in
    one hard reception cell to each weighted observer cell.  The Jacobian is
    constant within a source cell, so this is an exact sufficient statistic.
    """

    bins: int
    extent: float
    moments: np.ndarray
    transfer: np.ndarray
    occupancy: np.ndarray

    @classmethod
    def empty(cls, bins: int, extent: float):
        cells = bins * bins
        return cls(bins, extent, np.zeros((len(MOMENT_NAMES), cells), np.float64),
                   np.zeros((cells, cells), np.float64), np.zeros(cells, np.int64))

    def add(self, u0, v0, uf, vf) -> None:
        u0, v0, uf, vf = (np.asarray(x, np.float64) for x in (u0, v0, uf, vf))
        width = 2.0 * self.extent / self.bins
        col = np.floor((uf + self.extent) / width).astype(np.int64)
        row = np.floor((vf + self.extent) / width).astype(np.int64)
        valid = np.isfinite(u0 + v0 + uf + vf) & (row >= 0) & (row < self.bins) & (col >= 0) & (col < self.bins)
        source = row[valid] * self.bins + col[valid]
        # Stable parallel covariance merge.  Centered cross-products avoid the
        # catastrophic cancellation of raw sum-of-squares and reproduce the
        # SVD least-squares path to near machine precision.
        all_values = np.column_stack((u0[valid], v0[valid], uf[valid], vf[valid]))
        for key in np.unique(source):
            values = all_values[source == key]; nb = float(len(values)); mb = values.mean(axis=0)
            z = values-mb; cxx = z[:,:2].T@z[:,:2]; cxy = z[:,:2].T@z[:,2:]
            na = self.moments[0,key]; ma = self.moments[1:5,key].copy(); total=na+nb
            delta=mb-ma
            if na:
                factor=na*nb/total
                self.moments[5,key] += cxx[0,0]+factor*delta[0]*delta[0]
                self.moments[6,key] += cxx[0,1]+factor*delta[0]*delta[1]
                self.moments[7,key] += cxx[1,1]+factor*delta[1]*delta[1]
                self.moments[8,key] += cxy[0,0]+factor*delta[0]*delta[2]
                self.moments[9,key] += cxy[0,1]+factor*delta[0]*delta[3]
                self.moments[10,key] += cxy[1,0]+factor*delta[1]*delta[2]
                self.moments[11,key] += cxy[1,1]+factor*delta[1]*delta[3]
                self.moments[1:5,key] = ma+delta*(nb/total)
            else:
                self.moments[1:5,key]=mb
                self.moments[5:8,key]=(cxx[0,0],cxx[0,1],cxx[1,1])
                self.moments[8:12,key]=(cxy[0,0],cxy[0,1],cxy[1,0],cxy[1,1])
            self.moments[0,key]=total
        self.occupancy += np.bincount(source, minlength=self.bins*self.bins).astype(np.int64)

        # Existing TSC3x3 definition and edge renormalization, accumulated by
        # (target cell, hard source cell) rather than retaining individual rays.
        x, y = uf[valid], vf[valid]
        qx, qy = (x + self.extent) / width - .5, (y + self.extent) / width - .5
        cx, cy = np.floor(qx + .5).astype(np.int64), np.floor(qy + .5).astype(np.int64)
        def weight(d):
            a = np.abs(d)
            return np.where(a <= .5, .75-a*a, np.where(a <= 1.5, .5*(1.5-a)**2, 0.))
        xis, yis = (cx-1, cx, cx+1), (cy-1, cy, cy+1)
        xws, yws = tuple(weight(qx-i) for i in xis), tuple(weight(qy-i) for i in yis)
        norm = np.zeros(x.size)
        for xi, xw in zip(xis, xws):
            for yi, yw in zip(yis, yws):
                ok = (xi >= 0) & (xi < self.bins) & (yi >= 0) & (yi < self.bins)
                norm[ok] += xw[ok]*yw[ok]
        for xi, xw in zip(xis, xws):
            for yi, yw in zip(yis, yws):
                ok = (xi >= 0) & (xi < self.bins) & (yi >= 0) & (yi < self.bins)
                target = yi[ok]*self.bins + xi[ok]
                np.add.at(self.transfer, (target, source[ok]), xw[ok]*yw[ok]/norm[ok])

    def finalize(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        m = {name: self.moments[i] for i, name in enumerate(MOMENT_NAMES)}
        n = m["count"]
        q1 = np.full(n.size, np.nan); q2 = np.full(n.size, np.nan)
        for k in np.flatnonzero(n >= 6):
            xtx = np.array([[m["cxx00"][k],m["cxx01"][k]],[m["cxx01"][k],m["cxx11"][k]]])
            xty = np.array([[m["cxy00"][k],m["cxy01"][k]],[m["cxy10"][k],m["cxy11"][k]]])
            try: a = np.linalg.solve(xtx, xty)
            except np.linalg.LinAlgError: continue
            q1[k], q2[k] = a[0, 0]-a[1, 1], a[0, 1]+a[1, 0]
        good = np.isfinite(q1) & np.isfinite(q2)
        count = self.transfer[:, good].sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            g1 = (self.transfer[:, good] @ q1[good]) / count
            g2 = (self.transfer[:, good] @ q2[good]) / count
        g1[count == 0] = np.nan; g2[count == 0] = np.nan
        shape = (self.bins, self.bins)
        return g1.reshape(shape), g2.reshape(shape), count.reshape(shape)

    def state(self) -> dict[str, np.ndarray]:
        return {"moments": self.moments, "transfer": self.transfer, "occupancy": self.occupancy}

    @classmethod
    def from_state(cls, bins, extent, state):
        return cls(bins, extent, state["moments"], state["transfer"], state["occupancy"])
