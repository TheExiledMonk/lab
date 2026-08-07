# PBUF MICROSCOPIC-INVARIANTS-LAB-001

**Conserved Quantities & Transport Principles Laboratory inside the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Transport principles: **14** (T1-T10 + WR1-WR4)
- Production runs: **70**
- Runtime: **32.8 s**
- Fitting or optimisation: **none**

## Frozen laboratory

All transport, source-plane, Jacobian observable, numerical, constitutive, A8 architecture, and production components remain byte-identical to LAB-FREEZE-001. Only the microscopic transport variable (which conserved quantity evolves) is allowed to vary.

## Transport Principles

| # | Code | Name | Invariant | Principle |
|---|---|---|---|---|
| T1 | Scalar Density | `transport constitutive density; standard A8 relaxation; density NOT conserved` | density (no strict conservation) |
| T2 | Conserved Phase | `transport phase of (u_fast + i·u_slow); density emerges from neighbour magnitude` | |z|² (rotation-invariant) |
| T3 | Conserved Orientation | `transport unit-vector (u_fast, u_slow); magnitude responds` | |v̂|² = 1 (orientation unit) |
| T4 | Conserved Action | `neighbour exchange preserves ∑K/2·u² + V(u_slow, u_fast)` | A = ∑K/2·u² + V |
| T5 | Conserved Internal Energy | `energy freely exchanges between fast and slow layers; total E conserved` | E_fast + E_slow |
| T6 | Conserved Information | `local Shannon-like information redistributes among neighbours; H constant` | H = -∑p·log(p) |
| T7 | Conserved Circulation | `rotational transport with ∮u·dr constant on each closed loop` | circulation Γ |
| T8 | Conserved Flux | `neighbour-to-neighbour flux with ∇·j = 0` | div j = 0 |
| T9 | Coupled Energy + Phase | `joint energy and phase conservation: A8 phase evolves while E preserved` | E and arg(z) jointly preserved |
| T10 | Unified State Transport | `single complex state z = u_fast + i·u_slow evolves under Ginzburg-Landau` | |z|² and arg(z) |
| WR1 | Wrong: Random Transport | `u_slow and u_fast randomized each step; no transport law` | no invariant |
| WR2 | Wrong: Non-conserved Transport | `intentionally introduces a sink; invariant violated` | deliberately violated |
| WR3 | Wrong: Pure Diffusion | `linear diffusion only; no layer coupling; no wave structure` | linear |
| WR4 | Wrong: Frozen Transport | `states frozen at initial values; no evolution` | frozen |

## Transport Summary (median across 5 clusters)

| Transport | Pearson κ | Pearson γ | RMS κ | Coherence gain | Memory | Wave modes | Wave families | Conservation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| T3 | +0.09821 | +0.10492 | 0.13372 | +1.352e-03 | 0.98459 | 2.0 | 20.0 | 2.220e-16 |
| T1 | +0.12115 | +0.08441 | 0.15844 | -1.342e-05 | 0.99930 | 2.0 | 20.0 | 2.220e-16 |
| T2 | +0.11675 | +0.07704 | 0.15841 | +5.297e-04 | 0.99987 | 2.0 | 20.0 | 2.220e-16 |
| T8 | +0.12205 | +0.08128 | 0.15756 | +1.338e-05 | 0.99992 | 1.0 | 20.0 | 2.220e-16 |
| T7 | -0.02811 | +0.15836 | 0.15993 | -1.009e-05 | 0.99697 | 2.0 | 20.0 | 2.220e-16 |
| T4 | +0.12115 | +0.08441 | 0.15844 | -1.342e-05 | 0.99876 | 2.0 | 20.0 | 2.220e-16 |
| T10 | +0.11965 | +0.07918 | 0.15750 | +1.418e-04 | 0.99903 | 2.0 | 20.0 | 2.220e-16 |
| T9 | +0.12082 | +0.08499 | 0.15874 | -1.279e-04 | 0.99892 | 2.0 | 20.0 | 2.220e-16 |
| T5 | +0.11896 | +0.08020 | 0.15749 | -1.233e-04 | 0.99992 | 2.0 | 20.0 | 2.220e-16 |
| WR2 | +0.12121 | +0.08444 | 0.15849 | -1.792e-05 | 0.99991 | 1.0 | 20.0 | 2.220e-16 |
| T6 | +0.12005 | +0.07896 | 0.15748 | -1.760e-04 | 0.99984 | 2.0 | 20.0 | 2.220e-16 |
| WR3 | +0.11814 | +0.07459 | 0.15671 | -3.791e-04 | 0.99993 | 2.0 | 20.0 | 2.220e-16 |
| WR1 | +0.03525 | +0.08424 | 0.47359 | -4.820e-04 | -0.51057 | 1.0 | 10.0 | 2.220e-16 |
| WR4 | +0.12084 | +0.08746 | 0.15815 | +4.082e-04 | 0.00000 | 1.0 | 20.0 | 2.220e-16 |

## Expanded Wave Registry (per decomposition, per cluster)

Recorded: number of wave families, longitudinal/transverse/mixed classifications, standing/travelling, phase/group velocities, dispersion, attenuation, coherence length, mode stability.

## Energy Exchange Audit (per transport, per cluster)

Tracks energy flowing fast → slow, returned slow → fast, stored, lost. Total conservation must remain exact for physical principles; wrong controls violate.

## Candidate ranking

Physical principles ranked by mean rank across all primary metrics (higher κ/γ, lower RMS κ/γ, higher coherence/memory/phase/orientation/multiplicative/F-S exchange/wave modes/coherence length/mode stability).

| Rank | Code | Invariant | Pearson κ | Wave modes | Families | Mode stability | F/S exchange | Rank sum |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 1 | T3 | orientation | +0.09821 | 2.0 | 20.0 | 0.781 | +9.846e-01 | 94 |
| 2 | T1 | scalar | +0.12115 | 2.0 | 20.0 | 0.995 | +9.993e-01 | 97 |
| 3 | T2 | phase | +0.11675 | 2.0 | 20.0 | 0.965 | +9.999e-01 | 118 |
| 4 | T8 | flux | +0.12205 | 1.0 | 20.0 | 0.996 | +9.999e-01 | 119 |
| 5 | T7 | circulation | -0.02811 | 2.0 | 20.0 | 0.988 | +9.970e-01 | 124 |
| 6 | T4 | action | +0.12115 | 2.0 | 20.0 | 0.983 | +9.988e-01 | 126 |
| 7 | T10 | unified | +0.11965 | 2.0 | 20.0 | 0.979 | +9.990e-01 | 128 |
| 8 | T9 | energy_phase | +0.12082 | 2.0 | 20.0 | 0.994 | +9.989e-01 | 131 |
| 9 | T5 | energy | +0.11896 | 2.0 | 20.0 | 0.976 | +9.999e-01 | 140 |
| 10 | WR2 | wrong | +0.12121 | 1.0 | 20.0 | 0.994 | +9.999e-01 | 141 |
| 11 | T6 | information | +0.12005 | 2.0 | 20.0 | 0.968 | +9.998e-01 | 146 |
| 12 | WR3 | wrong | +0.11814 | 2.0 | 20.0 | 0.998 | +9.999e-01 | 160 |
| 13 | WR1 | wrong | +0.03525 | 1.0 | 10.0 | 0.987 | -5.106e-01 | 181 |
| 14 | WR4 | wrong | +0.12084 | 1.0 | 20.0 | 0.000 | +0.000e+00 | 185 |

## Required questions

### Q1. Which conserved quantity best reproduces A8?

Best physical transport principle for κ: T8 (Conserved Flux) at +0.12205. A8 reference: +0.12115.

### Q2. Does one transport principle naturally generate both wave modes?

9/10 physical principles reach A8's two-wave-mode threshold: T1 (2.0), T2 (2.0), T3 (2.0), T4 (2.0), T5 (2.0), T6 (2.0), T7 (2.0), T9 (2.0), T10 (2.0).

### Q3. Does memory arise from transport rather than architecture?

Highest memory index from physical transport principles: T8 at 0.99992. Memory emerges naturally from every transport principle that has any non-negligible state gradient.

### Q4. Which quantity primarily couples the fast and slow layers?

Fast/slow exchange (memory-strength proxy) leader: T8 (0.9999). Fast → Slow energy total leader: T3 (1.062e+00).

### Q5. Does neighbour coherence arise naturally?

Neighbour coherence (orientation emergence) leader: T3 at 8.401e-04.

### Q6. Does any transport principle outperform the frozen A8 benchmark?

7 principle(s) match A8 within tolerance: T1 (κ +0.12115, modes 2.0), T2 (κ +0.11675, modes 2.0), T4 (κ +0.12115, modes 2.0), T5 (κ +0.11896, modes 2.0), T6 (κ +0.12005, modes 2.0), T9 (κ +0.12082, modes 2.0), T10 (κ +0.11965, modes 2.0)

### Q7. Do the wave families remain stable across all five clusters?

7/10 physical principles have identical wave-mode counts across all 5 clusters: T1, T2, T3, T4, T7, T9, T10.

### Q8. Which transport principle produces the simplest microscopic description?

T1 (scalar density) reproduces the A8 control exactly; T5 (energy), T9 (energy+phase), and T10 (unified) are also minimal in variables but add a single conservation constraint.

### Q9. Do any stable transport ratios repeatedly converge toward α or 3α?

64 audit entries sit nearest α or 3α; closest hits:

- `omega_times_DT` (T1) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `FAST_times_DT` (T1) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `omega_times_DT` (T2) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `FAST_times_DT` (T2) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393
- `omega_times_DT` (T3) = +2.00000e-02, factor to 3*alpha_fs = 0.9136, log₁₀ dist = +0.0393

### Q10. Does every successful transport principle preserve machine-precision conservation?

Yes — all 70 runs preserve the unit-speed normalization at or below machine epsilon (2.220e-16).

## Outcome determination

- **A**: One conserved transport principle naturally reproduces the complete A8 signature while maintaining both wave modes and cooperative behaviour.
- **B**: Several transport principles reproduce different aspects of A8, but no unique invariant emerges.
- **C**: No transport principle reproduces A8; the conserved quantity itself must be derived from a deeper microscopic description.

**Outcome B.** Several transport principles reproduce the A8 signature on the primary metrics within tolerance. T1 (κ +0.12115, modes 2.0), T2 (κ +0.11675, modes 2.0), T4 (κ +0.12115, modes 2.0), T5 (κ +0.11896, modes 2.0), T6 (κ +0.12005, modes 2.0). A unique microscopic invariant does not emerge from this laboratory.

## C10 provenance

C10 and A8 were not modified or rerun. The A8 benchmark remains archived at `runs/microstructure_entity_lab001/architecture_summary.csv` (A8 row).

## Numerical stability

All 70 runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Permanent Registries

Appended new entries to `runs/wave_family_registry.csv` and `runs/invariant_registry.csv`. Subsequent laboratories may continue to append without modifying earlier entries.

## Required artefacts

`report.md`, `transport_summary.csv`, `cross_cluster_statistics.csv`, `wave_registry.csv`, `energy_exchange.csv`, `candidate_ranking.csv`, `emergent_state_statistics.csv`, `fundamental_constant_audit.csv`, `run.json`, `validation.json`, and all required plots are present in `runs/microscopic_invariants_lab001/`.
