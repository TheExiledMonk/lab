# DYNAMICS-001 dependency graph

```text
FOUNDATION-001 FP-1..FP-6          V11 operational limits
              |                            |
              v                            v
STATE-002: q in Q_adm          covariance / effective-limit gates
       |              \
       |               -> DEFORMATION-001: C[q,q0]
       |                              |
       |                              -> HYPER-001 + ENERGY-PRINCIPLE-001
       |                                     W=Phi(I1,I2,I3)
       v
oriented C1 histories / Diff_+(S)
       |
       -> v in T_q Q_adm, modulo positive rescaling
       |
       -> L(q,a v)=a L(q,v)
                    |
                    -> S[q]=integral L ds
                              |
                              -> delta S=0 on admissible variations

Missing before FIELD-001:
realization + function spaces
kinetic/action selection or derived relational clock
inertial normalization
boundary/endpoint data
one-metric and clock/ruler map
source projection and optional matter action
well-posedness + V11 effective-limit proof
```

The chain required by the brief is

```text
q -> S[q] -> delta S
```

No arrow in this graph denotes a field equation.
