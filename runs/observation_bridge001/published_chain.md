# Published Frontier Fields Pipeline - Stage-by-Stage Audit

Reconstruction method: SaWLens (Merten et al. 2009, 2011), as documented in the benchmark README files. All five clusters (Abell 2744, MACS J0416, MACS J1149, Abell S1063, Abell 370) are Frontier Fields lensing reconstructions at `Z_S = 9.0`. The weak-lensing inputs are Subaru / VLT / HST / ESO-WFI shear catalogues; the strong-lensing inputs are confirmed multiple-image systems.

## Weak-lensing shear catalogue

- Mathematical quantity: `Tabulated ellipticities e_1, e_2 of background galaxies at effective source redshift z_s, on-sky positions (RA, Dec)`
- Physical meaning: Direct image-ellipticity measurements from Subaru / VLT / HST / ESO-WFI imaging, passed through photo-z cuts to select background galaxies. The ellipticity is the *reduced shear* estimator: <e> ≈ g = γ / (1 - κ).
- Units: dimensionless (ellipticity components)
- Assumptions: Background selection by colour / photo-z; calibration of point-spread function and shear-estimation bias are external; the observed ellipticity is an estimator of the reduced shear only after correction.

## Strong-lensing multiple-image systems

- Mathematical quantity: `Sets of (RA_i, Dec_i, z_i) for each confirmed multiple image`
- Physical meaning: Direct observational input: spectroscopically and photometrically confirmed multiply-imaged sources with redshift estimates. Each system provides positional constraints on the deflection potential.
- Units: RA, Dec in degrees; redshift dimensionless
- Assumptions: Image identifications, redshifts and pairing assignments are supplied by the ST Frontier Fields map makers; these are external inputs to SaWLens.

## SaWLens parametric model (Merten et al. 2009, 2011)

- Mathematical quantity: `Multi-scale adaptive-mesh lensing inversion with a joint weak + strong likelihood; the field is parameterised on a three-level grid: low (full field), medium, high (cluster core).`
- Physical meaning: Bayesian inversion that fits both the convergence κ and the shear components g_1, g_2 (or γ_1, γ_2 depending on parameterisation) to the combined weak-shear + strong-position likelihood. Output is the posterior mean of the convergence field and of the shear components at a chosen source redshift.
- Units: dimensionless convergence and shear, scaled to the chosen z_S
- Assumptions: Parametric form (multi-resolution grid); uniform prior on convergence at each level; the published maps are the posterior-mean reconstructions at the chosen source redshift Z_S (here Z_S = 9.0, effectively an infinite-source approximation). The reconstruction is model-dependent.

## kappa.fits

- Mathematical quantity: `Convergence κ = Σ / Σ_crit at the source plane`
- Physical meaning: Surface-mass-density map of the lens, in units of the critical density Σ_crit(z_l, z_s). Includes dark matter, gas and stellar contributions. Depends explicitly on cosmology through the angular-diameter-distance ratio D_ls / D_s.
- Units: dimensionless (Σ / Σ_crit), scaled to z_S = 9
- Assumptions: Reconstructed posterior mean from SaWLens; Z_S = 9; lens redshift Z_L per cluster; standard ΛCDM distance ratios.

## gamma1.fits, gamma2.fits

- Mathematical quantity: `Components of the shear: g_1, g_2 (or γ_1, γ_2)`
- Physical meaning: Components of the complex shear field at the same source redshift. The Frontier Fields / SaWLens outputs are most commonly interpreted as the *reduced shear* components g_1 = γ_1 / (1 - κ), g_2 = γ_2 / (1 - κ), since the observable is the galaxy ellipticity, which is a direct estimator of g.
- Units: dimensionless (g or γ), scaled to z_S = 9
- Assumptions: Same SaWLens reconstruction; same Z_S; same cosmology. Whether the published map stores γ or g is determined by the SaWLens internal parameterisation and is not explicitly disambiguated in the supplied README.

## gamma.fits

- Mathematical quantity: `|γ| or |g| = sqrt(γ_1^2 + γ_2^2) (or sqrt(g_1^2 + g_2^2))`
- Physical meaning: Scalar magnitude of the (reduced) shear, supplied as an internal-consistency check. Always non-negative.
- Units: dimensionless
- Assumptions: Same as gamma1/gamma2 above.

## jacdet.fits, magnification.fits

- Mathematical quantity: `J = (1 - κ)^2 - |γ|^2;    magnification μ = 1/J`
- Physical meaning: Determinant of the lens mapping Jacobian and the corresponding magnification, derived from the reconstructed κ and γ.
- Units: dimensionless
- Assumptions: Derived observables; not directly observed but computed from the reconstructed κ and γ.

