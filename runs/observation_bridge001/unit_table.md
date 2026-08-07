# Unit Table

| Quantity | Published Units | Version A Units | Compatible |
|---|---|---|---|
| Convergence κ | dimensionless (Σ / Σ_crit); depends on D_ls(z_l, z_s) / D_s(z_s); scaled to Z_S = 9 | dimensionless; internal pipeline units; no cosmology; no distance ratio; no Σ_crit | NO |
| Shear γ_1 | dimensionless; reduced-shear component g_1 (γ_1 / (1 - κ)) from observed galaxy ellipticity; scaled to Z_S = 9 | dimensionless; raw shear component γ_1 from deflection gradient; no (1 - κ) division; no cosmology | NO |
| Shear γ_2 | dimensionless; reduced-shear component g_2 (γ_2 / (1 - κ)) from observed galaxy ellipticity; scaled to Z_S = 9 | dimensionless; raw shear component γ_2 from deflection gradient; no (1 - κ) division; no cosmology | NO |
| Shear magnitude |γ| | dimensionless; magnitude of the reduced shear; scaled to Z_S = 9 | dimensionless; magnitude of the raw shear; no cosmology | NO |
| Magnification μ | dimensionless; 1 / ((1 - κ)^2 - |γ|^2) from reconstructed fields | dimensionless; 1 / ((1 - κ)^2 - |γ|^2) from Version A fields | NO |
| Deflection α | dimensionless in arcsec-units? actually NOT supplied in the benchmark (no deflection maps published) | dimensionless in pipeline grid units [-8, 8] | NO (no published deflection to compare against) |
| Spatial coordinate x, y | RA / Dec on WCS grid; CDELT in deg / pixel; origin at CRVAL; pixel scale 6.25-11.36 arcsec / pixel per cluster | Dimensionless Cartesian on [-8, 8]; no WCS; no angular scale; no RA/Dec; origin at pipeline centre | NO (irreconcilable without external angular scaling) |
