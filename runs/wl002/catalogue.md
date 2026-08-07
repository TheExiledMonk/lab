# WL-002 constitutive-equation catalogue

All candidates used identical WL-001 inputs, geometry, propagation, reconstruction, benchmarks, residual calculations, and validation gates.

| Version | Local/propagating | Linear/nonlinear | Stiffness | PBUF RMSE | Gates |
|---|---|---|---|---:|---|
| A | local | local linear | constant | 0.0003822288849 | PASS |
| B | local | local nonlinear | constant | 0.0001921056932 | FAIL |
| C | local | local nonlinear | increases with loading | 0.0005291126983 | PASS |
| D | propagating | nonlocal linear propagation of nonlinear loading | constant propagation rigidity | 0.0001921056932 | PASS |

## Findings

- A is the linear local WL-001 baseline.
- B suppresses deformation in dilute outskirts. It improves the broad, horizontally shifted residual lobes because this synthetic reconstruction responds only to the spatial mean deformation. It fails Gate 4: the frozen path sampler does not encounter a large enough gradient to exceed the required trajectory-change threshold.
- C softens the normalized response at low loading while preserving the peak. Its broader deformation worsens those same residual lobes.
- D propagates B's loading through an elastic Helmholtz response. It changes the deformation and gradient maps, but preserves the source mean, so this reconstruction gives it the same RMSE as B (up to numerical precision).
- No new fitted constant was introduced: D's propagation length is the observed baryonic width already present in the fixed geometry.

## Recommendation

Carry Version D into WL-003 provisionally (RMSE 0.0001921056932, versus 0.0003822288849 for A). D is the strongest accepted law: it passes every gate, implements spatial propagation, and retains B's RMSE improvement; B is rejected despite the numerical tie because it fails Gate 4.

## Limitation exposed

The frozen WL-001 image reconstruction uses only `deformation.mean()` rather than the gradient field or photon landing positions. Consequently it cannot use RMSE to distinguish constitutive laws with equal mean deformation, and spatial residual-pattern claims beyond the global image shift are not identifiable. This is a pipeline limitation, not a reason to tune the constitutive laws or modify the frozen reconstruction during WL-002.

## Candidate dossiers

### 1. Version D — PREFERRED

**Equation:** `(1 - sigma_rho^2 Laplacian) u = u0 (rho/rho_max)^2`

**Motivation:** Distributed elastic recovery balances local loading; the observed baryonic width supplies the only length scale.

**Assumptions:** scalar isotropic medium; periodic numerical boundary; recovery is distributed.

**Strengths:** propagates deformation; stable positive spectral response; no fitted length.

**Weaknesses:** scalar proxy cannot represent shear stress; periodic boundary is a laboratory approximation.

**Evidence:** RMSE=0.0001921056932; deformation range=3.72575e-08..0.0606551; gradient RMS=0.00615742; gradient max=0.0440664; photon max deviation=0.000336594; finite=True; cost=0.5895s; gates=PASS.

### 2. Version A — RETAIN FOR COMPARISON

**Equation:** `u = u0 rho/rho_max`

**Motivation:** Minimal scalar, isotropic local response and the WL-001 limiting law.

**Assumptions:** deformation is scalar; response is instantaneous and local.

**Strengths:** fewest assumptions; exact WL-001 baseline.

**Weaknesses:** cannot propagate deformation.

**Evidence:** RMSE=0.0003822288849; deformation range=4.64025e-55..0.18; gradient RMS=0.0197473; gradient max=0.144808; photon max deviation=9.26395e-06; finite=True; cost=0.6859s; gates=PASS.

### 3. Version C — RETAIN FOR COMPARISON

**Equation:** `u = 2 u0 q/(1+q), q=rho/rho_max`

**Motivation:** A local compliance law tests a medium that stiffens with loading.

**Assumptions:** local equilibrium; K/K0=(1+q)/2.

**Strengths:** bounded; recovers the baseline at zero and peak loading.

**Weaknesses:** rigidity interpolation is postulated; no propagation.

**Evidence:** RMSE=0.0005291126983; deformation range=9.2805e-55..0.18; gradient RMS=0.0214573; gradient max=0.132898; photon max deviation=1.85154e-05; finite=True; cost=0.5877s; gates=PASS.

### 4. Version B — REJECT

**Equation:** `u = u0 (rho/rho_max)^2`

**Motivation:** Tests nonlinear loading without changing the medium model.

**Assumptions:** deformation is scalar; quadratic loading is physically selected.

**Strengths:** stable; suppresses dilute loading.

**Weaknesses:** quadratic exponent is not derived from PBUF; no propagation.

**Evidence:** RMSE=0.0001921056932; deformation range=1.19622e-108..0.18; gradient RMS=0.019711; gradient max=0.203925; photon max deviation=6.36741e-09; finite=True; cost=0.5894s; gates=FAIL.

## Families rejected on physical/interface grounds

### diffusion-type

**Equation:** `partial_t u = D Laplacian(u) + S(rho)`

**Motivation:** Local deformation spreads down spatial gradients.

**Assumptions:** a physical evolution time and diffusivity exist.

**Strengths:** causal evolution can be represented; smooths small-scale structure.

**Weaknesses:** the frozen static laboratory supplies neither time nor diffusivity.

**Outcome:** REJECTED BEFORE NUMERICAL RANKING: would require an arbitrary constant and an unfrozen time model.

**Recommendation:** Revisit only when PBUF defines a time scale.

### elastic-medium tensor

**Equation:** `div C:epsilon(u_vec) = f(rho)`

**Motivation:** A genuine elastic medium carries directional displacement and shear stress.

**Assumptions:** deformation is vector/tensor valued; elastic moduli and boundary tractions are known.

**Strengths:** represents shear and anisotropy; has conservation-law structure.

**Weaknesses:** WL-001 accepts one scalar deformation field; PBUF supplies no elastic moduli or boundary data.

**Outcome:** REJECTED BEFORE NUMERICAL RANKING: incompatible with the frozen scalar interface and underdetermined.

**Recommendation:** Use as the next generalization after observables constrain tensor components.

