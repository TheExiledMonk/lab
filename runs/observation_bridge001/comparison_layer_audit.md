# Comparison Layer Audit

## Current comparison layer

Direct numerical comparison between Version A κ, γ_1, γ_2, |γ| and the published κ, γ_1, γ_2, |γ| after a bilinear resampling onto the pipeline grid.

## Physical problem

The two fields share field names but represent different physical stages. The Version A outputs are internal dimensionless lensing-like observables on a dimensionless Cartesian grid; the published maps are reconstructed posterior-mean convergence and shear at a chosen source redshift on a real angular WCS.

## Alternative comparison layers

| Candidate | Status | Reason |
|---|---|---|
| Compare Version A's photon deflection direction to the published shear orientation | POSSIBLE BUT UNDERSPECIFIED | Both pipelines have a 'deflection-like' field; the mapping between Version A's 90-degree-transverse photon response and the GR deflection potential is not established by the frozen theory. |
| Compare Version A's response field direction to the published κ gradient | NOT APPLICABLE | Version A's response is the 90-degree rotation of ∇C; the published κ gradient has a different meaning (mass-density gradient on the sky). The two are not the same quantity. |
| Compare Version A's κ_pred as a qualitative convergence pattern to the published κ | PARTIALLY | As a qualitative *shape* comparison (peak location, asymmetry, extent) the two share the same dimensionless convergence label, but no quantitative agreement is expected because the matter input is itself a unit-bearing proxy and the response model is different. |
| Forward-pipe through a cosmological bridge | REQUIRED | The physically correct bridge is: ρ -> κ(Σ / Σ_crit, z_l, z_s, cosmology) -> shear -> observables. This bridge is NOT supplied by the frozen laboratory; it would require an explicit Σ_crit evaluation and a chosen cosmology, neither of which is frozen. |
