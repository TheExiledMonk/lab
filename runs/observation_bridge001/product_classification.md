# Product Classification

Classification classes (per milestone spec):
- Direct observation
- Derived observable
- Reconstruction
- Inversion product
- Model-dependent quantity

| File | Class | Justification |
|---|---|---|
| kappa.fits | Reconstruction | Posterior-mean convergence map from a parametric SaWLens inversion of weak-shear + strong-lens catalogues. Not directly observed; inferred from a model fit. |
| gamma.fits | Derived observable | Magnitude of the (reduced) shear, derived from the gamma1 / gamma2 components of the SaWLens reconstruction. Computed from the reconstructed field, not measured directly. |
| gamma1.fits | Reconstruction | First component of the (reduced) shear from the SaWLens inversion. The observable is galaxy ellipticity; the published map is the posterior mean of the field component. |
| gamma2.fits | Reconstruction | Second component of the (reduced) shear, same provenance as gamma1.fits. |
