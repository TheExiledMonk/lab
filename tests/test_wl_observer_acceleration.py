import unittest
import numpy as np

import observable_lab001 as OLD
from pbuf.wl.backends.vulkan_kde import CpuExactKDE, make_kde_backend
from pbuf.wl.backends.vulkan_runtime import vulkan_available
from pbuf.wl.observer_cache import ObserverPrimitiveCache, ObserverStateId
from pbuf.wl.observer_dependencies import describe_observer_dependencies, KDE_UNIFORM_TRANSLATION_INVARIANT
from pbuf.wl.observer_profile import ObserverProfile


class ObserverAccelerationTests(unittest.TestCase):
    def test_dependency_inventory(self):
        inventory=describe_observer_dependencies();self.assertEqual(len(inventory),45);self.assertEqual(len(set(inventory)),45)
        self.assertTrue(KDE_UNIFORM_TRANSLATION_INVARIANT)
    def test_cpu_translation_invariance_and_self(self):
        u=np.array([-.3,.1,.8,1.2,-.7,0.]);v=np.array([.4,-.2,.9,-.5,.3,.1]);cpu=CpuExactKDE();base=cpu.evaluate(u,v)
        for du,dv in ((1e-12,0),(-1e-12,0),(0,1e-12),(0,-1e-12),(1e-12,1e-12)):
            np.testing.assert_allclose(base,cpu.evaluate(u+du,v+dv),rtol=1e-13,atol=1e-15)
        h=OLD._gaussian_kde_bandwidth(np.vstack((u,v)));self.assertTrue(np.all(base >= 1/(2*np.pi*h[0]*h[1]*len(u))))
    def test_cache_state_backend_and_deposition_independence(self):
        p=ObserverProfile();c=ObserverPrimitiveCache(p);calls=[]
        def calc():calls.append(1);return np.array([3.])
        a=ObserverStateId("cpu_base",backend="cpu");b=ObserverStateId("cpu_other",backend="cpu");g=ObserverStateId("cpu_base",backend="vulkan")
        key=lambda s,m="hard":c.key("pairwise_kde",s,coordinates=(np.arange(3.),),parameters=(1,),translation_invariant=True)
        c.get_or_compute(key(a),calc,"pairwise_kde");c.get_or_compute(key(a,"six methods"),calc,"pairwise_kde")
        c.get_or_compute(key(b),calc,"pairwise_kde");c.get_or_compute(key(g),calc,"pairwise_kde")
        self.assertEqual(len(calls),3);self.assertEqual(p.describe()["pairwise_kde"]["cache_hit_count"],1)
    def test_unknown_backend(self):
        with self.assertRaisesRegex(ValueError,"unsupported KDE backend: wat"):make_kde_backend("wat")

@unittest.skipUnless(vulkan_available(),"float64 Vulkan unavailable")
class VulkanExactKDETests(unittest.TestCase):
    def test_sizes_tail_repeatability(self):
        rng=np.random.default_rng(111)
        from pbuf.wl.backends.vulkan_kde import VulkanExactKDE
        with VulkanExactKDE() as gpu:
            for n in (1,2,7,31,257,1000):
                u=rng.normal(size=n);v=rng.normal(size=n);ref=CpuExactKDE().evaluate(u,v);a=gpu.evaluate(u,v);b=gpu.evaluate(u,v)
                np.testing.assert_allclose(a,ref,rtol=1e-11,atol=1e-13);self.assertTrue(np.array_equal(a,b));self.assertTrue(np.all(np.isfinite(a)))
                self.assertTrue(gpu.last_timing["no_n_squared_buffer_allocation"])

if __name__=="__main__":unittest.main()
