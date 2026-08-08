# C_STATE Bridge Propagation Shape 001

This lab advances the supported native accumulation bridge into a synthetic propagation-shape audit without assigning a physical lensing amplitude.

Frozen native chain:

`rho -> existing A8 transport -> raw c_state -> six-neighbor bounded-strain equilibrium -> accumulated u`

For a static weak metric, first-order null-ray bending depends on the transverse gradient of the Weyl combination `Phi+Psi`. The current bridge supplies one accumulated scalar `u`, but the physical coefficient and tensor split are not yet derived. Therefore the lab evaluates only the coefficient-independent unit-channel shape

`alpha_tilde(b) = integral dz d_x u`.

Any constant map `Phi+Psi = C_W u` multiplies all values by the same unresolved `C_W`, leaving the tested symmetries and log-log exponents unchanged.

Predeclared structural checks:

- zero source gives zero propagation response;
- central transverse response vanishes for a centered source;
- reflection gives `alpha(+b) = -alpha(-b)`;
- for a positive unit Weyl-channel coefficient, the gradient points toward the source;
- weak-regime response scales linearly with source mass;
- impact-parameter magnitude is near `b^-1` without inserting a `1/r` or `1/b` law.

This lab does not claim the physical map from `u` to `Phi+Psi`, does not split `Phi` and `Psi`, does not compute a calibrated deflection angle, and does not compare to kappa, shear, HST, or any observed lensing target.

Valid outcomes:

- `C_STATE_BRIDGE_PROPAGATION_SHAPE_SUPPORTED`
- `C_STATE_BRIDGE_PROPAGATION_SHAPE_PARTIAL_SUPPORT`
- `C_STATE_BRIDGE_PROPAGATION_SHAPE_NOT_SUPPORTED`

Run:

`PYTHONPATH=. python pbuf/labs/foundation/c_state_bridge_propagation_shape001.py`
