# PBUF OBSERVATION-BRIDGE-001

Audit-only milestone. No PBUF modification. No parameter fitting. No cosmological scaling introduced.

## Scope

Determine whether the observational products used in
WEAK-LENSING-OBSERVATION-001 correspond to the same physical
stage of the lensing pipeline as the frozen Version A outputs.

## Version A chain (frozen, identical parameters)

Every stage is documented in `version_a_chain.md`. Summary:

- **Matter density ρ(X)**: Input scalar field declared to be a non-negative dimensionless matter density on the pipeline grid....
- **Constitutive field C(X) = 0.18 · ρ(X) / ρ_max**: Local linear scalar proxy of the medium's deformation by matter. C is the output of Version A's constitutive equation. No physical dimension is attached; C exis...
- **Gradient field ∇C**: Gradient of the constitutive scalar field. Has the magnitude of a response per dimensionless length....
- **Response field r = (rx, ry)**: Vector field representing the local transverse response of the medium. The 90-degree rotation is the frozen 'transport' choice; the response amplitude equals th...
- **Photon propagation**: Iterative ray-tracing through the response field. Photons start at x = -8 with v = (1, 0); the pipeline runs for steps = 80 iterations with step = 0.06, so the ...
- **Predicted convergence κ**: Local photon-count ratio in a 64 x 64 histogram of (x_f, y_f) after propagation versus the initial (x_0, y_0) histogram. The pipeline only produces finite κ on ...
- **Predicted shear γ_1, γ_2**: Components of the 2 x 2 Jacobian of the photon-displacement field, evaluated on the same 64 x 64 grid. No convergence-to-shear correction is applied; no reduced...
- **Predicted magnification μ**: Magnification derived from the standard lensing identity, using the predicted κ and γ. NaN wherever the denominator is non-positive or one of the inputs is unde...

## Published Frontier Fields chain

Every stage is documented in `published_chain.md`. Summary:

- **Weak-lensing shear catalogue**: Direct image-ellipticity measurements from Subaru / VLT / HST / ESO-WFI imaging, passed through photo-z cuts to select background galaxies. The ellipticity is t...
- **Strong-lensing multiple-image systems**: Direct observational input: spectroscopically and photometrically confirmed multiply-imaged sources with redshift estimates. Each system provides positional con...
- **SaWLens parametric model (Merten et al. 2009, 2011)**: Bayesian inversion that fits both the convergence κ and the shear components g_1, g_2 (or γ_1, γ_2 depending on parameterisation) to the combined weak-shear + s...
- **kappa.fits**: Surface-mass-density map of the lens, in units of the critical density Σ_crit(z_l, z_s). Includes dark matter, gas and stellar contributions. Depends explicitly...
- **gamma1.fits, gamma2.fits**: Components of the complex shear field at the same source redshift. The Frontier Fields / SaWLens outputs are most commonly interpreted as the *reduced shear* co...
- **gamma.fits**: Scalar magnitude of the (reduced) shear, supplied as an internal-consistency check. Always non-negative....
- **jacdet.fits, magnification.fits**: Determinant of the lens mapping Jacobian and the corresponding magnification, derived from the reconstructed κ and γ....

## Pipeline Comparison Diagram

![Pipeline comparison]
(pipeline_comparison_diagram.png)

The two pipelines share field symbols (κ, γ, μ) but originate
from physically different stages: the published pipeline
starts from real astronomical observations (galaxy
ellipticities and confirmed multiple-image systems) and
ends at reconstructed posterior-mean fields at a chosen
source redshift Z_S = 9 with explicit cosmological scaling;
Version A starts from an internal dimensionless scalar
input and ends at dimensionless lensing-like observables
with no cosmology, no redshift, and no physical length.

## Product Classification

| File | Class |
|---|---|
| `kappa.fits` | Reconstruction |
| `gamma.fits` | Derived observable |
| `gamma1.fits` | Reconstruction |
| `gamma2.fits` | Reconstruction |

Full justification in `product_classification.md`.

## Physical Mapping Table

| Published | Version A | Comparable | Reason |
|---|---|---|---|
| kappa.fits | Predicted κ (from photon-count ratio) | PARTIALLY | Same mathematical symbol κ but Version A κ is a local photon-density distortion (dimensionless Cartesian units, no cosmology, no Σ_crit); published κ = Σ / Σ_crit is the surface-mass-density ratio and depends explicitly on cosmological distance ratios and source redshift Z_S. |
| gamma1.fits | Predicted γ_1 (from deflection-gradient) | PARTIALLY | Same symbol γ_1 but different physical process. Version A γ_1 is derived from the displacement of photons through a PBUF response field; published γ_1 (or g_1) is the posterior-mean shear component from a SaWLens inversion of observed galaxy ellipticities. Even if the published quantity is the reduced shear g_1 = γ_1 / (1 - κ), Version A does not perform the (1 - κ) division and the two γ_1 share only the field name. |
| gamma2.fits | Predicted γ_2 (from deflection-gradient) | PARTIALLY | Same caveats as gamma1.fits; same conclusion. |
| gamma.fits | Predicted |γ| | PARTIALLY | |γ| = sqrt(γ_1^2 + γ_2^2) by construction in both pipelines, so the magnitude is internally consistent, but the underlying γ_1, γ_2 carry different physical content as documented above. |
| jacdet.fits / magnification.fits | Predicted μ | NO | Published μ = 1 / ((1 - κ)^2 - |γ|^2) is computed from the SaWLens-reconstructed κ and γ at Z_S = 9; Version A μ uses the Version A κ and γ and is not rescaled by any cosmological factor. A direct comparison would require matching Z_S and adding the missing cosmological bridge. |

## Unit Table

| Quantity | Published Units | Version A Units | Compatible |
|---|---|---|---|
| Convergence κ | dimensionless (Σ / Σ_crit); depends on D_ls(z_l, z_s) / D_s(z_s); scaled to Z_S = 9 | dimensionless; internal pipeline units; no cosmology; no distance ratio; no Σ_crit | NO |
| Shear γ_1 | dimensionless; reduced-shear component g_1 (γ_1 / (1 - κ)) from observed galaxy ellipticity; scaled to Z_S = 9 | dimensionless; raw shear component γ_1 from deflection gradient; no (1 - κ) division; no cosmology | NO |
| Shear γ_2 | dimensionless; reduced-shear component g_2 (γ_2 / (1 - κ)) from observed galaxy ellipticity; scaled to Z_S = 9 | dimensionless; raw shear component γ_2 from deflection gradient; no (1 - κ) division; no cosmology | NO |
| Shear magnitude |γ| | dimensionless; magnitude of the reduced shear; scaled to Z_S = 9 | dimensionless; magnitude of the raw shear; no cosmology | NO |
| Magnification μ | dimensionless; 1 / ((1 - κ)^2 - |γ|^2) from reconstructed fields | dimensionless; 1 / ((1 - κ)^2 - |γ|^2) from Version A fields | NO |
| Deflection α | dimensionless in arcsec-units? actually NOT supplied in the benchmark (no deflection maps published) | dimensionless in pipeline grid units [-8, 8] | NO (no published deflection to compare against) |
| Spatial coordinate x, y | RA / Dec on WCS grid; CDELT in deg / pixel; origin at CRVAL; pixel scale 6.25-11.36 arcsec / pixel per cluster | Dimensionless Cartesian on [-8, 8]; no WCS; no angular scale; no RA/Dec; origin at pipeline centre | NO (irreconcilable without external angular scaling) |

## Matter Input Audit

**Claim:** WEAK-LENSING-OBSERVATION-001 used ρ = max(κ, 0) / max(κ) as the matter-density input to the frozen Version A constitutive law.

**Verdict:** `approximation`

Full justification in `matter_input_audit.md`. In short,
treating max(κ, 0) as a bare matter density ρ is an
*approximation*. The published κ is a Σ/Σ_crit map, not a
mass density. It depends on cosmology through Σ_crit and
includes dark matter; the normalisation and any baryonic/
dark partition are lost in the substitution. The non-negativity
clamp and peak normalisation further suppress the published
field's negative tails.

## Coordinate Audit

Published products use equatorial RA/Dec on a WCS grid
(TAN projection) with explicit angular pixel scales
(6.25-11.36 arcsec/pixel). Version A uses a dimensionless
Cartesian grid on [-8, 8]. The two cannot be aligned without
imposing an external angular scale. The bilinear resampling
performed in WEAK-LENSING-OBSERVATION-001 mapped
pixel index (0, N-1) -> (-extent, +extent) and discarded all
angular information, including the WCS handedness sign
(CD1_1 < 0).

## Comparison Layer Audit

The current comparison is a direct numerical comparison of
Version A κ, γ_1, γ_2, |γ| against the published κ, γ_1,
γ_2, |γ| after bilinear resampling onto the pipeline grid.
This layer is physically incomplete:

1. The two field families carry different physical content
   (raw shear vs. reduced shear, dimensionless pipeline
   units vs. Σ/Σ_crit, no Z_S vs. Z_S = 9).
2. The comparison ignores the cosmological bridge (Σ_crit,
   D_ls/D_s) that the published products already encode.
3. The comparison ignores the unit/angular bridge that the
   WCS encodes.

Full discussion in `comparison_layer_audit.md`. The four
alternative layers considered are tabulated there.

## Required Conclusion

**Are the currently compared quantities physically equivalent?**

**PARTIALLY**.

The field symbols (κ, γ_1, γ_2, |γ|, μ) are mathematically
related to the lensing observables in both pipelines, but the
physical content is different:

- **Convergence κ**: same symbol, different physical
  quantity. Version A κ is a photon-density distortion;
  published κ is a Σ/Σ_crit mass map at Z_S = 9. They are
  not the same physical object.
- **Shear γ_1, γ_2**: same symbol, different physical
  quantity. Version A γ is the 2x2 Jacobian of the photon
  displacement; published γ (most likely the reduced shear
  g) is a posterior-mean reconstruction from image
  ellipticities. Even if both stored γ rather than g, the
  underlying deflection potential is generated by a
  different physical mechanism (PBUF transport vs.
  gravitational potential).
- **Magnification μ**: same symbol, different inputs.

**First point of physical divergence**

Stage 0: the two pipelines do not even share an input. The
published pipeline takes real astronomical observations
(galaxy ellipticities + strong-lens positions) and the
Version A pipeline takes an internal dimensionless scalar
field ρ(X). There is no point in either pipeline at which a
direct physical quantity is interchangeable.

**Can the discrepancy be removed by unit conversion,
coordinate conversion, or normalisation?**

**No.** The discrepancies are *not* a matter of units,
coordinates, or normalisation. They are a matter of physical
content:

- A *cosmological bridge* (Σ_crit at the given z_l, z_s) is
  required to convert Version A's dimensionless lensing-like
  outputs into a Σ/Σ_crit map. The frozen laboratory does not
  select a cosmology or compute Σ_crit.
- A *reduced-shear bridge* (γ -> g = γ / (1 - κ)) is
  required if the published maps store g rather than γ.
  Version A does not compute this division.
- A *source-redshift scaling bridge* is required if Version A
  were to produce outputs at a specific z_s. The frozen
  Version A has no source redshift at all.
- A *physical-unit bridge* is required to convert Version A's
  dimensionless Cartesian grid into an angular grid. The
  frozen grid [-8, 8] carries no physical length.

These bridges require an explicit observational product and
explicit cosmological inputs that are not frozen. Under the
frozen laboratory the comparison performed in
WEAK-LENSING-OBSERVATION-001 is therefore **not a comparison
of like with like**. The large disagreement recorded in that
milestone is, in part, an artefact of the missing bridges
rather than a property of the frozen Version A pipeline.

## Identical-pipeline verification (SHA-256)

| File | SHA-256 |
|---|---|
| `observation_bridge001.py` | `73ee7256bd0c4c6170a42ec4edf3ce5c22be2499c25807bd52ef11e8b9448b71` |
| `constitutive_equations.py` | `e2c789d19fd559753519704c6668c7a2879c53eb61a315604ec81af6795aca9f` |
| `weak_lensing_observation001.py` | `a5c3632fec9adfc2659d5c283d07c599db6db7edb2e34e4aebd84e35434642bc` |

## Notes

- The frozen Version A pipeline (constitutive + transport +
  response + observables) is unchanged in this milestone.
- No cosmological scaling, no fitting, no parameter change,
  no reinterpretation of Version A has been performed.
- All five benchmark datasets are read-only; the FITS
  products and the WEAK-LENSING-OBSERVATION-001 outputs are
  consumed for comparison and metadata recording only.
- Total execution time: 0.44 s.
