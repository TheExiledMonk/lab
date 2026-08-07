# BALANCE-001 dependency graph

```text
FOUNDATION-001 + STATE-002
          |
          v
 q -> C[q,q0], oriented unparameterized histories
          |
          +------------------------------+
          |                              |
          v                              v
 HYPER/ENERGY: W, P_C=DW       DYNAMICS: degree-one action,
          |                    p=D_v L, <p,v>-L=0
          +---------------+--------------+
                          v
              BALANCE-001 templates
        d_tau rho + Div J = sigma
        nabla_eff . J = Sigma
                          |
                          v
                constitutive closure
   (stored energy selection, inertia, flux/source maps,
       optional dissipation, boundary work/loading)
                          |
             DURATION-001 + METRIC-001
             (tau calibration and G[q,C;D])
                          |
                          v
            native/effective source projection
                          |
                          v
       governing evolution equations and constraints
```

The required shorthand is therefore

`q -> Balance -> Constitutive Law -> Evolution Equations`,

with duration, metric, source projection, and data as explicit closure gates.
