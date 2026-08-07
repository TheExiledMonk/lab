# CONTINUUM-STATE-001 dependency graph

```mermaid
flowchart TD
    V["V11:<br/>effective GR and Lorentz invariance"]
    O["PBUF premise:<br/>real physical continuum"]
    K["KINEMATICS-001:<br/>current/reference configuration,<br/>objective relative deformation"]
    C["CONSTITUTIVE/NATURE:<br/>recoverable stored energy,<br/>state-dependent response"]
    S["Complete constitutive state<br/>on an admissible state slice"]
    U["Mechanical single-occupancy:<br/>one configuration/deformation<br/>per element at that state position"]
    E["Ordered succession of<br/>complete continuum states"]
    T["TIME-001:<br/>candidate emergent ordering"]
    MT["Missing temporal identification:<br/>relational clock -> metric proper time"]
    R["V11-compatible effective<br/>relativistic time"]
    A["Missing cross-history<br/>actuality-and-slicing principle"]
    P1["Single-actual-state interpretation"]
    P2["Growing-history interpretation"]
    P3["Whole-history interpretation"]
    D["Missing dynamics:<br/>successors, reversibility,<br/>causality, arrow, memory"]

    O --> K
    K --> S
    C --> S
    S --> U
    U --> E
    E --> T
    MT --> R
    T --> R
    A --> P1
    A --> P2
    A --> P3
    E --> P1
    E --> P2
    E --> P3
    D --> E
    R --> V
    P1 -. "compatibility gate" .-> V
    P2 -. "compatibility gate" .-> V
    P3 -. "compatibility gate" .-> V
```

Mechanical single-occupancy constrains every actuality branch: no branch permits two incompatible complete deformations of one element on the same state slice. The branch point occurs only at cross-history actuality, where continuum mechanics supplies no selection rule. Dashed arrows are required V11 compatibility tests, not derived identities.

