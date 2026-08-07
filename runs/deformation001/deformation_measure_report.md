# PBUF DEFORMATION-001 — Native deformation measure of the spacetime medium

## Foundation status and decision

FOUNDATION-001 and FP-1 through FP-6 are fixed inputs. No ontology is revised here, and no constitutive equation, stored-energy functional, material constant, field equation, or observational fit is introduced.

**Decision:** the native deformation variable is fixed uniquely at the abstract kinematic level as a **dimensionless objective relative-deformation endomorphism**

\[
C[q,q_0]:V_0\rightarrow V_0,\qquad C[q_0,q_0]=\mathbf 1.
\]

Here `q` is the current complete medium configuration, `q_0` is an unloaded reference representative, and `V_0` is the relevant reference tangent space. Physical deformation is the gauge-invariant spectral content of `C`, not its coordinate components.

A scalar or vector is insufficient because it cannot encode generic volumetric and shear deformation. A raw deformation gradient `F` contains an unphysical rigid-rotation part. When a material or coframe map exists, the simplest objective representative is therefore the right Cauchy–Green-type tensor

\[
C=F^{\sharp}F,
\]

where `sharp` is the adjoint defined by the unloaded metric. On the admissible material branch, `C` is symmetric positive-definite and dimensionless. Equivalent one-to-one strain coordinates, such as

\[
E=\tfrac12(C-\mathbf1),\qquad H=\tfrac12\log C,
\]

do not constitute different physical deformation variables; they are reparametrizations of `C` on their domains.

The accepted inputs do not close one remaining realization choice. For instantaneous three-direction material kinematics, `C` has rank three. A four-direction relative coframe/tensor representation is also mathematically admissible if it is the structure that carries emergent clocks and rulers. Thus the concrete field is reduced to the smallest unavoidable family:

1. material coordinates: `C^I_J=B^IK kappa_KJ`, with `B^IJ=g^mu nu partial_mu phi^I partial_nu phi^J`;
2. relative coframe: `F^a_b=E^a_mu (Ebar^-1)^mu_b`, followed by `C=F^sharp F`; or
3. relative symmetric tensor: `C^mu_nu=qbar^mu alpha q_alpha nu` on an admissible real spectral branch.

These are realizations of the same abstract relative-endomorphism requirement, but are not interchangeable physical theories until the clock/ruler and one-metric identification is supplied.

## Reference configuration: undeformed spacetime

“Undeformed spacetime” is not a zero tensor and not a preferred coordinate chart. It is the equivalence class

\[
\mathcal R_0=\{q_0\;/\;(\mathrm{Diff}\times G_{\rm int})\mid C[q_0,q_0]=\mathbf1\}.
\]

Every representative must be nondegenerate, orientation- and signature-admissible, locally compatible with the V11 Minkowski limit in a freely falling frame, and homogeneous and isotropic when used as the V11 cosmological background. Rigid coordinate, frame, or material-label changes do not create deformation.

The corpus does not determine whether `R_0` is a fixed material configuration, an instantaneous FLRW natural configuration, or a temperature-dependent family `R_0(T,a)`. This is physically consequential because it determines whether cosmological expansion is counted as deformation. V11's `alpha`, `epsilon_0`, `Omega_sigma`, saturation history, and `Rmax` retain their original meanings and do not define the reference state or a limiting stretch.

## Complete isotropic invariant catalogue

Let `lambda_A>0` be the three eigenvalues of a rank-three `C`. For a local parity-even isotropic medium, a complete algebraically independent set is

\[
I_1=\operatorname{tr}C,
\]

\[
I_2=\tfrac12\left[(\operatorname{tr}C)^2-\operatorname{tr}(C^2)\right],
\]

\[
I_3=\det C.
\]

Their unloaded values are `(3,3,1)`. Equivalently, one may use the unordered principal values `{lambda_1,lambda_2,lambda_3}`, or a volumetric coordinate `J=sqrt(I_3)` plus two independent invariants of `Cbar=I_3^(-1/3)C`. These are coordinate changes on invariant space, not extra physical inputs.

For a genuine rank-four clock/ruler comparison, the complete generic set is the four elementary symmetric polynomials `e_n(C)`, `n=1,...,4`, with unloaded values `(4,6,4,1)`. Three invariants suffice only after the temporal eigenvalue is fixed or gauged by an additional identification principle.

Parity-odd pseudoscalars, derivatives `nabla C`, curvature, two-point invariants, and invariants involving an independent four-velocity are not part of the minimal local isotropic set. They become admissible only if corresponding orientation, gradient, nonlocal, or congruence structure is separately established.

## Small-deformation limit

Write

\[
C=\mathbf1+2\varepsilon+O(\varepsilon^2),\qquad
\varepsilon=\tfrac12(C-\mathbf1)+O((C-\mathbf1)^2).
\]

Then

\[
\varepsilon=\tfrac13(\operatorname{tr}\varepsilon)\mathbf1+\varepsilon_{\rm TF}.
\]

The trace is the infinitesimal volumetric channel and the trace-free symmetric part is the infinitesimal shear channel. For a displacement realization `F=1+grad u`,

\[
\varepsilon=\operatorname{sym}(\nabla u)+O(|\nabla u|^2),
\]

so rigid infinitesimal rotations drop out. For a relative metric/coframe realization with `q=q_0+delta q`,

\[
\delta C=q_0^{-1}\delta q+O(\delta q^2),
\]

and hence `epsilon` is linearly related to the usual weak metric perturbation after the gauge quotient. This is exactly the tensorial linear kinematics required by V11's retained weak-field GR regime. It does not derive Einstein dynamics or a normalization map from medium strain to the physical metric; those remain downstream gates.

The invariants expand as

\[
I_1=3+2\operatorname{tr}\varepsilon+O(\varepsilon^2),
\]

\[
I_2=3+4\operatorname{tr}\varepsilon+O(\varepsilon^2),
\]

\[
I_3=1+2\operatorname{tr}\varepsilon+O(\varepsilon^2).
\]

Shear first appears through trace-free information at quadratic invariant order, while remaining present linearly in the tensor `epsilon`. This is why a scalar alone cannot be the complete native variable.

## Behaviour near a finite elastic bound

The deformation measure itself supplies no saturation law. Let an admissible path approach a finite elastic endpoint with principal values

\[
\lambda_A(s)\rightarrow\lambda_A^*,\qquad 0<\lambda_A^*<\infty.
\]

Then `C`, `E`, the elementary invariants, and `H=1/2 log C` all approach finite values continuously. The boundary is a boundary of the admissible spectral domain, not a singularity forced by the deformation measure.

If an eigenvalue tends to zero, `C` can remain finite but loses invertibility; `log C`, reciprocal/Eulerian measures, and determinant-normalized shear coordinates become singular. If an eigenvalue diverges, the elementary invariants diverge. Accordingly, the regularity-safe native description is the principal spectrum or elementary invariants on an explicitly nondegenerate, orientation-preserving branch. A bounded “saturation coordinate” cannot be selected without adding the prohibited constitutive normalization and bound map.

V11's finite activation/saturation history does not by itself identify any `lambda_A^*`; no inversion from that homogeneous history to local strain is authorized.

## Compatibility audit

- **Rotational symmetry/objectivity:** common rigid rotations alter `F` but cancel from `C=F^sharp F`; only principal values and invariant combinations are physical.
- **Covariance:** `C` is constructed tensorially from `q` and `q_0`; diffeomorphisms change representatives, not scalar invariants. Internal frame or material relabellings act by similarity and leave the spectrum unchanged.
- **Isotropy:** isotropy permits dependence only on symmetric functions of the principal values. It does not reduce generic deformation to one scalar and does not remove shear.
- **Emergent time:** `C` is defined on each complete instantaneous medium state and compared along an ordered, reparametrization-equivalent history. No fundamental time component or rate is built into the deformation measure. A four-direction clock/ruler realization is an effective covariant representation and requires the separate temporal-identification principle before being called native.
- **V11 relativistic compatibility:** the unloaded state admits the local Minkowski limit, and the tangent perturbation is a symmetric tensor compatible with weak-field metric kinematics. The missing one-metric map and dynamics must still be matched downstream; DEFORMATION-001 does not modify V11.

## Dependency graph toward HYPER-001

```text
FOUNDATION-001 (FP-1--FP-6) + V11 retained limits
                         |
                         v
       current configuration q + unloaded class R_0
                         |
                         v
       objective relative endomorphism C[q,R_0]
                  /              \
                 v                v
       rank-3 material branch   rank-4 clock/ruler branch
                 \                /
                  v              v
     admissible signature/orientation/spectral domain
                         |
                         v
       characteristic invariants {I_A} or {e_n}
                         |
                         v
   closure gate: primitive realization + one-metric identification
                         |
                         v
                     HYPER-001
       (future invariant functional; not specified here)
```

HYPER-001 may begin only after choosing the primitive realization, the reference-state family, rank three versus four, the internal gauge quotient, and the relation to the single effective metric. It may then use only the corresponding complete invariant set as local isotropic arguments.

## Completion statement

The milestone is complete under the stated success criterion. A bare scalar, vector, raw displacement, or raw deformation gradient is excluded as the complete native deformation. The unique abstract measure is the objective dimensionless relative endomorphism `C`; its physical content is its spectrum. The remaining admissible concrete family is irreducibly split between rank-three material and rank-four clock/ruler realizations until PBUF supplies the already-identified reference-and-identification principle. No constitutive content has been inferred.
