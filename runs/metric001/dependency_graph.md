# METRIC-001 dependency graph

```mermaid
flowchart LR
  F["FOUNDATION-001"] --> Q["q in Q_adm"]
  S["STATE-002"] --> Q
  Q --> C["C[q,q0] in D_C"]
  D["DEFORMATION-001"] --> C
  H["HYPER-001 + ENERGY-PRINCIPLE-001<br/>admissible finite-bound domain"] --> C
  Q --> G["admissible operational constitutive map G"]
  C --> G
  DU["DURATION-001<br/>clock functional and V11 matching"] --> G
  V["V11 fixed Lorentzian target"] --> G
  G --> M["g_eff = G[q,C]"]
  M --> O["clock durations, ruler lengths,<br/>signal cones, synchronization"]
  M -. "future input" .-> FI["FIELD-001"]
```

The requested core chain is

\[
\boxed{q\longrightarrow C[q,q_0]\longrightarrow
G\in\mathfrak G\longrightarrow g^{\rm eff}_{\mu\nu}.}
\]

The direct edge \(q\to G\) is essential: accepted rank-three \(C\) does not
contain clock normalization or temporal–spatial information.  The dashed edge
does not derive or authorize field equations.
