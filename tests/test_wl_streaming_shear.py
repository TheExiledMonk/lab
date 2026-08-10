import unittest
import numpy as np

from pbuf.wl.shear_readout import build_shear_candidates, construct_local_primitives, evaluate_candidate
from pbuf.wl.streaming_shear import JacobianTSCAccumulator


class StreamingShearTests(unittest.TestCase):
    def test_streamed_matches_full_path(self):
        rng = np.random.default_rng(115)
        u0 = rng.uniform(-.99, .99, 20000); v0 = rng.uniform(-.99, .99, 20000)
        uf = u0 + .03*u0 + .02*v0 + rng.normal(0, .002, u0.size)
        vf = v0 - .01*u0 - .02*v0 + rng.normal(0, .002, u0.size)
        rays = {"u0": u0, "v0": v0, "uf": uf, "vf": vf,
                "dx": np.zeros_like(u0), "dy": np.zeros_like(u0)}
        rays.update(construct_local_primitives(rays, bins=8, extent=1.0))
        spec = next(x for x in build_shear_candidates(("tsc_3x3",)) if x.family == "D_jacobian")
        reference = evaluate_candidate(spec, rays, bins=8, extent=1.0)
        acc = JacobianTSCAccumulator.empty(8, 1.0)
        for sl in (slice(0, 7000), slice(7000, 13000), slice(13000, None)):
            acc.add(u0[sl], v0[sl], uf[sl], vf[sl])
        actual = acc.finalize()[:2]
        np.testing.assert_allclose(actual, reference, rtol=1e-12, atol=1e-12, equal_nan=True)


if __name__ == "__main__": unittest.main()
