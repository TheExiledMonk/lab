# PBUF MICROSCOPIC-TRANSPORT-EQUIVALENCE-LAB-001

**Transport Representation and Equivalence-Class Laboratory inside the frozen Version 1 weak-lensing laboratory (LAB-FREEZE-001).**

## Status

- Frozen hash verification: **PASS**
- Successful candidates evaluated: **6** (T1, T4, T5, T6, T9, T10)
- Production runs: **30** (6 candidates × 5 clusters)
- Wrong-control runs: **8** (4 controls × 2 clusters)
- Runtime: **100.7 s**
- Fitting, lag optimisation, transport modification: **none**

## Frozen baseline

Reproduces `runs/microscopic_invariants_lab001/` with byte-identical transport equations for the six successful candidates. Hashes verified before execution.

## Sampling protocol

1024 cells per cluster, deterministic (seed = 42 + cluster hash): 256 high-density, 256 medium-density, 256 low-density, 256 uniform-grid. Same coordinates across all six candidates (see `sampled_cells.csv`).

## Required questions

### Q1. Are T1, T4, T5, T6, T9, and T10 only macroscopically equivalent, or dynamically equivalent?

15 pairs reach Level 3-4 (dynamical/representation), 0 pairs reach Level 2 (mode), 0 pairs reach Level 0-1 (macroscopic/observable).

### Q2. Which candidate pairs reach each equivalence level?

T1-T4: Level 4 (representation equivalence)
T1-T5: Level 3 (dynamical equivalence)
T1-T6: Level 3 (dynamical equivalence)
T1-T9: Level 4 (representation equivalence)
T1-T10: Level 3 (dynamical equivalence)
T4-T5: Level 3 (dynamical equivalence)
T4-T6: Level 3 (dynamical equivalence)
T4-T9: Level 4 (representation equivalence)
T4-T10: Level 3 (dynamical equivalence)
T5-T6: Level 3 (dynamical equivalence)
T5-T9: Level 3 (dynamical equivalence)
T5-T10: Level 3 (dynamical equivalence)
T6-T9: Level 3 (dynamical equivalence)
T6-T10: Level 3 (dynamical equivalence)
T9-T10: Level 3 (dynamical equivalence)

### Q3. Does one fixed representation transformation map any candidate pair across all five clusters?

3 pair(s) reach Level 4 (representation equivalence) under a fixed transformation across all five clusters.

### Q4. Do successful principles share the same fast/slow exchange cycle?

3 of 15 pairs share an equivalent fast/slow exchange cycle.

### Q5. Are the two wave modes equivalent across transport principles?

Wave-mode equivalence: see `wave_mode_equivalence.csv`. 15 pair(s) have equivalent longitudinal and transverse wave modes.

### Q6. Do similar final κ values arise from similar or different microscopic trajectories?

15 pairs match κ but also share trajectories; 0 pairs match κ but have divergent trajectories.

### Q7. Does conserved action provide any internal behaviour not already present in scalar-density transport?

T1 vs T4: Level 4 (representation equivalence), trajectory correlation 1.000, derivative correlation 1.000.

### Q8. Are energy, information, and unified-state transport distinguishable from scalar-density transport?

T1 vs T4: trajectory corr 1.000, level Level 4 (representation equivalence)
T1 vs T5: trajectory corr 0.992, level Level 3 (dynamical equivalence)
T1 vs T6: trajectory corr 0.992, level Level 3 (dynamical equivalence)
T1 vs T9: trajectory corr 1.000, level Level 4 (representation equivalence)
T1 vs T10: trajectory corr 0.978, level Level 3 (dynamical equivalence)

### Q9. Does coupled energy-plus-phase transport introduce a genuinely new state-space geometry?

T1 vs T9 geometry difference: {'path_length': 0.00015145168337671986}

### Q10. Do wrong controls correctly destroy dynamical equivalence while preserving selected marginal statistics?

WR1: trajectory destroyed at 0.354, derivative at -0.002
WR2: trajectory destroyed at 0.355, derivative at -0.002
WR3: trajectory destroyed at 0.496, derivative at 0.483
WR4: trajectory destroyed at 0.762, derivative at 0.744
WR1: trajectory destroyed at 0.301, derivative at 0.003
WR2: trajectory destroyed at 0.302, derivative at -0.003
WR3: trajectory destroyed at 0.494, derivative at 0.482
WR4: trajectory destroyed at 0.746, derivative at 0.742

### Q11. Do any independently generated dimensionless quantities recur near α, 3α, or α⁻¹?

108 audit entries nearest α or 3α; 0 nearest inverse α.

### Q12. Is the laboratory observing one transport equivalence class or several physically distinct classes?

Of 15 pairs: 15 reach dynamical/representation equivalence, 0 reach mode equivalence, 0 reach only macroscopic equivalence.

## Outcome determination

- **A**: Most successful principles reach Level 3 or Level 4 equivalence.
- **B**: Multiple equivalence subclasses; some principles dynamically equivalent, others only macroscopically equivalent.
- **C**: Macroscopic degeneracy only; weak-lensing cannot uniquely determine the microscopic invariant.
- **D**: No stable classification; equivalence assignments vary strongly by cluster or normalization.

**Outcome A.** 15 of 15 pairs reach Level 3 or 4 equivalence: transport principles are largely representation-equivalent descriptions of one underlying flow.

## Numerical stability

All 30 production runs preserve the frozen unit-speed normalization at or below machine epsilon (2.220e-16).

## Required artefacts

`report.md`, `transport_pair_summary.csv`, `trajectory_equivalence.csv`, `derivative_equivalence.csv`, `state_space_geometry.csv`, `wave_mode_equivalence.csv`, `fast_slow_exchange_equivalence.csv`, `representation_transform_results.csv`, `wrong_control_results.csv`, `fundamental_constant_audit.csv`, `sampled_cells.csv`, `candidate_classification.csv`, `run.json`, `validation.json`, and all required plots are present in `runs/microscopic_transport_equivalence_lab001/`.
