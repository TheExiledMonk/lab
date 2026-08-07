# PBUF MATERIAL-DISCOVERY-001 — Systematic Material-Mechanism Survey

## Scope, method, and non-selection rule

This catalogue compares established material mechanisms; it does not claim that spacetime literally consists of chains, cells, particles, fields, resonators, fluids, or fibres. FOUNDATION-001 through MATERIAL-LAB-002 are frozen inputs. The native state remains the objective relative deformation `C[q,q0]` on its admissible SPD domain. Scalar formulas, displacement fields, lattices, and microstructural stories are comparison devices only.

Every row is tested against the same laboratory contract: a local elastic candidate must admit an objective invariant lift `W(C)=Phi(I1,I2,I3)`, the unloaded state must be a stable minimum, and the weak tangent must agree with HYPER-001. A local energy supplies stress but not neighbour communication. A gradient, nonlocal, lattice, or adjacency operator supplies neighbour communication but not inertia. Wave claims therefore remain conditional until the frozen balance/duration architecture supplies a positive kinetic closure. Recovery means energetic restoring tendency unless an evolution law is explicitly present. Ray stability never proves tensor rank-one convexity, polyconvexity, strong ellipticity, or hyperbolicity.

Compatibility labels have strict meanings: **compatible** would require no missing gate; **conditional** can be represented using the frozen state after stated closures/derivations; **ontology-change** requires an independent field, constituent, dimensional carrier, or sector absent from the frozen state; **incompatible** fails a required behavior such as shear storage, stability, or recovery; **unknown** is reserved for insufficient mathematical specification. No mechanism is selected or ranked by taste.

## E01 — Hookean linear elasticity

| Required item | Evaluation |
|---|---|
| Physical mechanism | Spring-like resistance proportional to infinitesimal strain. |
| Why it exists in nature | Atomic/network bonds linearized near equilibrium. |
| Typical representation | `sigma=lambda tr(eps)I+2mu eps` |
| Stored energy | W=lambda(tr eps)^2/2+mu eps:eps |
| Interaction behavior | Local stress divergence after balance closure. |
| Neighbour coupling | Only through the continuum balance operator. |
| Wave support | Supports nondispersive longitudinal and shear waves with positive inertia. |
| Recovery | Exact elastic return in its infinitesimal domain. |
| Progressive hardening | None; tangent is constant. |
| Finite deformation | Not objective as a complete finite-strain law. |
| Weak field | Required Lamé-type tangent when moduli are positive. |
| Large deformation | Unbounded linear response becomes physically incomplete. |
| Mathematical strengths | Simple, differentiable, analytically ready. |
| Mathematical weaknesses | Fails finite-deformation objectivity and progressive hardening. |
| PBUF compatibility | **conditional** — Admissible only as the frozen weak-field tangent. |
| Potential PBUF role | weak-field limit |
| Sources | S08 |

## E02 — Saint-Venant-Kirchhoff

| Required item | Evaluation |
|---|---|
| Physical mechanism | Extends Hooke elasticity by making energy quadratic in Green strain. |
| Why it exists in nature | Finite kinematics with a linear material stress-strain relation. |
| Typical representation | `E=(C-I)/2; S=lambda tr(E)I+2mu E` |
| Stored energy | W=lambda(tr E)^2/2+mu E:E |
| Interaction behavior | Local hyperelastic stress plus balance. |
| Neighbour coupling | No intrinsic length or nonlocal coupling. |
| Wave support | Conditional elastic waves where strong ellipticity holds. |
| Recovery | Reversible on a stable branch. |
| Progressive hardening | Not reliably progressive; can soften under finite loading. |
| Finite deformation | Objective finite kinematics but unstable/nonphysical in parts of large compression. |
| Weak field | Matches Hooke elasticity. |
| Large deformation | Loss of ellipticity and poor large-strain behavior possible. |
| Mathematical strengths | Very simple tensor lift. |
| Mathematical weaknesses | Global admissibility is weak. |
| PBUF compatibility | **conditional** — Fits the native C but fails a global stability guarantee. |
| Potential PBUF role | baseline/rejected complete law |
| Sources | S03 |

## E03 — Compressible neo-Hookean

| Required item | Evaluation |
|---|---|
| Physical mechanism | An isotropic rubber-like network resists distortional and volumetric change. |
| Why it exists in nature | Gaussian-chain elasticity plus a volumetric penalty. |
| Typical representation | `W=mu(I1bar-3)/2+U(J)` |
| Stored energy | Same as representation. |
| Interaction behavior | Local first variation of invariant energy. |
| Neighbour coupling | No intrinsic neighbour term. |
| Wave support | Finite elastic waves if the tangent is strongly elliptic. |
| Recovery | Reversible to the unique energy minimum. |
| Progressive hardening | Usually mild; no finite-chain stiffening. |
| Finite deformation | Objective and finite-strain capable on J>0. |
| Weak field | Positive Lamé tangent with suitable U. |
| Large deformation | No finite extension barrier; volumetric choice controls extremes. |
| Mathematical strengths | Minimal invariant finite-strain model. |
| Mathematical weaknesses | May miss strong hardening; polyconvexity depends on U. |
| PBUF compatibility | **conditional** — Direct invariant lift is possible; full spectral gates remain to prove. |
| Potential PBUF role | constitutive law |
| Sources | S01, S03 |

## E04 — Mooney-Rivlin

| Required item | Evaluation |
|---|---|
| Physical mechanism | Two invariant channels describe rubber-like shear response. |
| Why it exists in nature | Phenomenological/network elasticity beyond one-invariant neo-Hooke. |
| Typical representation | `W=C10(I1bar-3)+C01(I2bar-3)+U(J)` |
| Stored energy | Same as representation. |
| Interaction behavior | Local energetic stress. |
| Neighbour coupling | No intrinsic neighbour term. |
| Wave support | Conditional on positive acoustic tensor and inertia. |
| Recovery | Reversible on a stable parameter branch. |
| Progressive hardening | Parameter-dependent; not necessarily monotone. |
| Finite deformation | Good moderate finite strain. |
| Weak field | Can share the common positive tangent. |
| Large deformation | No finite stretch barrier and stability is coefficient/domain dependent. |
| Mathematical strengths | Simple two-channel invariant response. |
| Mathematical weaknesses | Coefficients are underived and global stability is nonautomatic. |
| PBUF compatibility | **conditional** — Compatible functional form on C; parameter restrictions required. |
| Potential PBUF role | constitutive law |
| Sources | S01, S02 |

## E05 — Ogden principal-stretch series

| Required item | Evaluation |
|---|---|
| Physical mechanism | Power-law terms in principal stretches flexibly represent nonlinear elasticity. |
| Why it exists in nature | Phenomenological spectral representation of elastomer response. |
| Typical representation | `W=sum_p mu_p/alpha_p(sum_i lambdabar_i^alpha_p-3)+U(J)` |
| Stored energy | Same as representation. |
| Interaction behavior | Local energetic stress. |
| Neighbour coupling | No intrinsic neighbour term. |
| Wave support | Conditional; acoustic tensor can lose positivity for some parameters. |
| Recovery | Reversible on stable branches. |
| Progressive hardening | Can harden or soften depending on exponents and coefficients. |
| Finite deformation | Excellent finite-deformation expressiveness. |
| Weak field | Parameters can be constrained to the frozen tangent. |
| Large deformation | Extrapolation and global ellipticity are delicate. |
| Mathematical strengths | Spectral and highly flexible. |
| Mathematical weaknesses | Many parameters; admissibility not transparent. |
| PBUF compatibility | **conditional** — Uses the frozen spectrum directly but needs strict coefficient/domain audit. |
| Potential PBUF role | constitutive family |
| Sources | S01, S02 |

## E06 — Yeoh/reduced polynomial

| Required item | Evaluation |
|---|---|
| Physical mechanism | Higher powers of the first distortional invariant increase resistance. |
| Why it exists in nature | Empirical nonlinear elasticity of filled elastomers. |
| Typical representation | `W=sum_n C_n0(I1bar-3)^n+U(J)` |
| Stored energy | Same as representation. |
| Interaction behavior | Local energetic stress. |
| Neighbour coupling | No intrinsic neighbour term. |
| Wave support | Conditional on tangent positivity and kinetic closure. |
| Recovery | Reversible if energy remains single-well. |
| Progressive hardening | Natural for nonnegative higher coefficients. |
| Finite deformation | Objective finite-strain model. |
| Weak field | Quadratic tangent is selectable. |
| Large deformation | Leading positive terms are coercive but no finite endpoint. |
| Mathematical strengths | Systematic polynomial hierarchy. |
| Mathematical weaknesses | One-invariant form can miss deformation modes; mixed signs destabilize. |
| PBUF compatibility | **conditional** — Matches MATERIAL-LAB polynomial class after tensor stability checks. |
| Potential PBUF role | constitutive law |
| Sources | S02 |

## E07 — Gent finite extensibility

| Required item | Evaluation |
|---|---|
| Physical mechanism | Resistance diverges as an invariant approaches a limiting chain extensibility. |
| Why it exists in nature | Finite extensibility of polymer networks. |
| Typical representation | `W=-(mu Jm/2)log(1-(I1bar-3)/Jm)+U(J)` |
| Stored energy | Logarithmic barrier. |
| Interaction behavior | Local energetic stress. |
| Neighbour coupling | No intrinsic neighbour term. |
| Wave support | Conditional before the barrier. |
| Recovery | Reversible within the open elastic domain. |
| Progressive hardening | Strong asymptotic hardening. |
| Finite deformation | Finite invariant endpoint; tensor boundary coverage is incomplete alone. |
| Weak field | Reduces to neo-Hooke at small strain. |
| Large deformation | Infinite energy at the selected limit. |
| Mathematical strengths | Simple barrier and correct weak limit. |
| Mathematical weaknesses | Jm is underived; one invariant does not guard all SPD boundaries. |
| PBUF compatibility | **conditional** — Natural barrier representative, never selected or parameterized here. |
| Potential PBUF role | constitutive law |
| Sources | S04, S05 |

## E08 — Arruda-Boyce eight-chain

| Required item | Evaluation |
|---|---|
| Physical mechanism | An affine network of finite chains stiffens as chains approach full extension. |
| Why it exists in nature | Statistical mechanics of an eight-chain polymer cell. |
| Typical representation | `Inverse-Langevin network energy, often expanded in powers of I1bar` |
| Stored energy | Finite-chain stored energy. |
| Interaction behavior | Local energetic stress after homogenization. |
| Neighbour coupling | Microscopic network is motivation, not a native PBUF constituent. |
| Wave support | Conditional elastic waves. |
| Recovery | Reversible below chain limit. |
| Progressive hardening | Progressive finite-chain hardening. |
| Finite deformation | Finite extensibility through the inverse Langevin response. |
| Weak field | Neo-Hookean leading limit. |
| Large deformation | Divergent/very steep limiting response. |
| Mathematical strengths | Mechanistic link between network geometry and barrier behavior. |
| Mathematical weaknesses | Imports a chain picture and chain count that PBUF cannot assume. |
| PBUF compatibility | **conditional** — Only its macroscopic invariant response may be compared. |
| Potential PBUF role | constitutive exemplar |
| Sources | S03, S05 |

## E09 — Fung-Demiray exponential

| Required item | Evaluation |
|---|---|
| Physical mechanism | Collagenous/soft tissues exhibit an exponentially rising tangent. |
| Why it exists in nature | Recruitment and alignment of biological microstructure under stretch. |
| Typical representation | `W=A(exp[B(I1bar-3)]-1)+U(J)` |
| Stored energy | Exponential invariant energy. |
| Interaction behavior | Local energetic stress. |
| Neighbour coupling | No intrinsic neighbour term. |
| Wave support | Conditional on strong ellipticity. |
| Recovery | Reversible idealization. |
| Progressive hardening | Strong smooth progressive hardening. |
| Finite deformation | Objective finite deformation. |
| Weak field | Expansion supplies a quadratic weak tangent. |
| Large deformation | Coercive without a finite endpoint. |
| Mathematical strengths | Smooth, simple rapid hardening. |
| Mathematical weaknesses | Growth scale and biological rationale are not PBUF derivations. |
| PBUF compatibility | **conditional** — Macroscopic exponential class is admissible after normalization. |
| Potential PBUF role | constitutive law |
| Sources | S02 |

## E10 — Hencky logarithmic elasticity

| Required item | Evaluation |
|---|---|
| Physical mechanism | Energy is quadratic in logarithmic principal strain. |
| Why it exists in nature | Multiplicative stretches become additive in log strain. |
| Typical representation | `H=log(C)/2; W=K(tr H)^2/2+mu devH:devH` |
| Stored energy | Quadratic logarithmic energy. |
| Interaction behavior | Local energetic stress. |
| Neighbour coupling | No intrinsic neighbour term. |
| Wave support | Conditional where ellipticity holds. |
| Recovery | Reversible on SPD domain. |
| Progressive hardening | Constant in H; not progressive along log-strain rays. |
| Finite deformation | Naturally handles large stretch ratios while C stays SPD. |
| Weak field | Matches linear elasticity. |
| Large deformation | Diverges at zero/infinite stretch but lacks finite endpoint. |
| Mathematical strengths | Clean volumetric/deviatoric split. |
| Mathematical weaknesses | Global convexity/ellipticity is domain-sensitive. |
| PBUF compatibility | **conditional** — H is an authorized reparametrization of C. |
| Potential PBUF role | constitutive coordinate/baseline |
| Sources | S03 |

## E11 — Incompressible hyperelastic constraint

| Required item | Evaluation |
|---|---|
| Physical mechanism | Volume preservation is enforced by a pressure multiplier while shear stores energy. |
| Why it exists in nature | Near-incompressibility of rubbers and liquids. |
| Typical representation | `J=1; sigma=-pI+sigma_dev` |
| Stored energy | W=What(I1bar,I2bar)+indicator_{J=1} |
| Interaction behavior | Pressure constraint plus local deviatoric stress. |
| Neighbour coupling | No neighbour term by itself. |
| Wave support | Shear waves possible; compressional mode becomes constrained. |
| Recovery | Reversible if the deviatoric law is elastic. |
| Progressive hardening | Inherited from What. |
| Finite deformation | Finite deformation on the J=1 manifold. |
| Weak field | Eliminates volumetric compliance, which may conflict with the frozen tangent. |
| Large deformation | Singular bulk limit. |
| Mathematical strengths | Mathematically clear constraint. |
| Mathematical weaknesses | Cannot represent a generic volumetric PBUF channel. |
| PBUF compatibility | **conditional** — Usable only if later derivation removes the volumetric mode. |
| Potential PBUF role | restricted constitutive subclass |
| Sources | S02 |

## W01 — Strain-gradient elasticity

| Required item | Evaluation |
|---|---|
| Physical mechanism | Energy penalizes rapid spatial variations as well as local strain. |
| Why it exists in nature | Finite microstructural correlation lengths in solids. |
| Typical representation | `E=int[W(C)+ell^2|Grad C|^2/2]dV` |
| Stored energy | Local plus positive gradient energy. |
| Interaction behavior | Euler-Lagrange operator contains higher spatial derivatives. |
| Neighbour coupling | Intrinsic nearest-neighbour continuum coupling. |
| Wave support | Supports dispersive elastic waves with positive inertia. |
| Recovery | Conservative recovery; boundary conditions matter. |
| Progressive hardening | Inherited from W; gradient stiffens short wavelengths, not amplitudes. |
| Finite deformation | Finite local strain plus smoothness regularization. |
| Weak field | Local tangent remains common; q-dependent stiffness appears. |
| Large deformation | High-wave-number stiffness; extra boundary data required. |
| Mathematical strengths | Explicit conservative communication and regularization. |
| Mathematical weaknesses | Introduces length scale ell and higher-order boundary conditions. |
| PBUF compatibility | **conditional** — Matches LAB-002 interaction slot if ell is derived under FP-6. |
| Potential PBUF role | interaction/hybrid mechanism |
| Sources | S23 |

## W02 — Integral nonlocal elasticity

| Required item | Evaluation |
|---|---|
| Physical mechanism | Stress/energy at one point samples a finite neighbourhood. |
| Why it exists in nature | Long-range bonds or homogenized microstructure. |
| Typical representation | `E=1/4 int int (C(x)-C(y)):K(x,y):(C(x)-C(y)) dxdy` |
| Stored energy | Positive symmetric kernel energy. |
| Interaction behavior | Nonlocal integral force operator. |
| Neighbour coupling | Explicit finite-horizon or decaying neighbour coupling. |
| Wave support | Dispersive waves for a positive Fourier symbol and inertia. |
| Recovery | Conservative for symmetric kernels. |
| Progressive hardening | Kernel nonlinearity may harden; linear kernels do not. |
| Finite deformation | Handles finite deformation only with objective two-point measures. |
| Weak field | Kernel moments can recover local elasticity. |
| Large deformation | Boundary/horizon effects and kernel positivity dominate. |
| Mathematical strengths | Direct neighbour communication without a lattice ontology. |
| Mathematical weaknesses | Kernel, horizon and objectivity are underived. |
| PBUF compatibility | **conditional** — Admissible interaction form after objective tensor construction. |
| Potential PBUF role | interaction mechanism |
| Sources | S12 |

## W03 — Bond-based peridynamics

| Required item | Evaluation |
|---|---|
| Physical mechanism | Material points interact across a finite horizon through pairwise bond forces. |
| Why it exists in nature | Continuum idealization of long-range cohesive interactions and fracture. |
| Typical representation | `rho u_ddot=int_H f(x'-x,u'-u)dV'+b` |
| Stored energy | Microelastic versions possess a bond potential. |
| Interaction behavior | Integral bond-force operator. |
| Neighbour coupling | Intrinsic finite-horizon integral coupling. |
| Wave support | Supports dispersive waves; stability follows positive bond stiffness. |
| Recovery | Elastic bonds recover; bond-breaking versions do not. |
| Progressive hardening | Nonlinear bond force can harden or soften. |
| Finite deformation | Large motion possible with objective bond stretch. |
| Weak field | Local elasticity recovered asymptotically, with bond-based Poisson restrictions. |
| Large deformation | Fracture/bond failure creates irreversible large-deformation behavior. |
| Mathematical strengths | No spatial derivatives and explicit nonlocality. |
| Mathematical weaknesses | Pairwise form restricts elastic constants; horizon is extra data. |
| PBUF compatibility | **conditional** — Elastic state-based limit may serve interaction; breaking is incompatible with recovery. |
| Potential PBUF role | interaction mechanism |
| Sources | S13 |

## W04 — Micropolar/Cosserat elasticity

| Required item | Evaluation |
|---|---|
| Physical mechanism | Independent microrotation and couple stress capture rotational cell mechanics. |
| Why it exists in nature | Lattices, foams and granular media with rotational units. |
| Typical representation | `W(epsilon(u,phi),kappa=Grad phi)` |
| Stored energy | Strain and curvature energy. |
| Interaction behavior | Force-stress and couple-stress divergences. |
| Neighbour coupling | Intrinsic through rotation gradients. |
| Wave support | Additional rotational/optic wave branches. |
| Recovery | Elastic recovery when energy is positive. |
| Progressive hardening | Usually linear; nonlinear versions may harden. |
| Finite deformation | Finite micropolar kinematics exist. |
| Weak field | Adds rotational modes beyond ordinary elasticity. |
| Large deformation | Size effects and extra branches persist. |
| Mathematical strengths | Captures lattice rotations and dispersion. |
| Mathematical weaknesses | Requires an independent rotation field and moduli. |
| PBUF compatibility | **ontology-change** — Independent microrotation is absent from the frozen state unless derived from C/q. |
| Potential PBUF role | rejected pending derived internal field |
| Sources | S11 |

## W05 — Nonlinear elastic wave continuum

| Required item | Evaluation |
|---|---|
| Physical mechanism | A nonlinear hyperelastic stress and inertia transport finite-amplitude disturbances. |
| Why it exists in nature | Elastic solids at finite wave amplitude. |
| Typical representation | `rho u_ddot=Div P(F), P=dW/dF` |
| Stored energy | Kinetic plus hyperelastic energy. |
| Interaction behavior | Spatial divergence of nonlinear stress. |
| Neighbour coupling | Intrinsic local-neighbour continuum coupling through gradients of displacement. |
| Wave support | Supports amplitude-dependent speeds, steepening and shocks. |
| Recovery | Conservative smooth solutions can recover; shocks require entropy/dissipation. |
| Progressive hardening | Inherited from W; can be progressive. |
| Finite deformation | Designed for finite deformation. |
| Weak field | Linearized acoustic tensor recovers weak waves. |
| Large deformation | Shock formation can destroy smooth classical solutions. |
| Mathematical strengths | Directly joins energy, interaction and waves. |
| Mathematical weaknesses | Requires authorized kinetic closure and hyperbolicity proof. |
| PBUF compatibility | **conditional** — A governing-equation template, not a selected PBUF law. |
| Potential PBUF role | wave/hybrid mechanism |
| Sources | S24 |

## F01 — Tension-only cable/rope network

| Required item | Evaluation |
|---|---|
| Physical mechanism | Slender members carry tension, go slack in compression and align with load. |
| Why it exists in nature | Ropes, tendons and cable nets. |
| Typical representation | `W=sum_a w_a(lambda_a) with zero/low compressive branch` |
| Stored energy | Fibre stretch energy. |
| Interaction behavior | Node/member force balance. |
| Neighbour coupling | Explicit network adjacency. |
| Wave support | Tension waves along taut paths; slack regions do not transmit them. |
| Recovery | Elastic if fibres do not slip or yield. |
| Progressive hardening | Geometric recruitment and fibre law can harden. |
| Finite deformation | Large rotations and stretches are natural. |
| Weak field | An isotropic random network may homogenize to an elastic tangent. |
| Large deformation | Nonsmooth slack-taut transitions and anisotropic localization. |
| Mathematical strengths | Clear neighbour topology and load paths. |
| Mathematical weaknesses | Discrete fibres/topology are not frozen ontology; compression support is poor. |
| PBUF compatibility | **conditional** — Only homogenized isotropic network response may be compared. |
| Potential PBUF role | topology-driven exemplar |
| Sources | S07 |

## F02 — Worm-like-chain network

| Required item | Evaluation |
|---|---|
| Physical mechanism | Semiflexible chains straighten entropically and stiffen sharply near contour length. |
| Why it exists in nature | Biopolymers such as actin and DNA. |
| Typical representation | `f(x)~(kBT/Lp)[1/(4(1-x/L)^2)-1/4+x/L]` |
| Stored energy | Chain free energy integrated from force. |
| Interaction behavior | Network nodes transmit chain tension. |
| Neighbour coupling | Explicit network coupling before homogenization. |
| Wave support | Elastic waves after network homogenization and inertia. |
| Recovery | Entropic reversible recovery in ideal conditions. |
| Progressive hardening | Strong finite-extensibility hardening. |
| Finite deformation | Finite contour-length endpoint. |
| Weak field | Small extension has a finite tangent after prestress/reference choice. |
| Large deformation | Force diverges near contour length. |
| Mathematical strengths | Mechanistic finite-extensible response. |
| Mathematical weaknesses | Temperature, chain and network ontology cannot be imported into PBUF. |
| PBUF compatibility | **conditional** — Macroscopic barrier structure only; microscopic story is excluded. |
| Potential PBUF role | constitutive exemplar |
| Sources | S07 |

## F03 — Holzapfel-Gasser-Ogden fibre composite

| Required item | Evaluation |
|---|---|
| Physical mechanism | An isotropic matrix is reinforced by dispersed fibre families that engage in tension. |
| Why it exists in nature | Arterial walls and collagenous tissues. |
| Typical representation | `W=Wiso(I1)+sum_a k1/(2k2)(exp[k2 Ea^2]-1)+U(J)` |
| Stored energy | Matrix plus anisotropic exponential energy. |
| Interaction behavior | Local energetic stress. |
| Neighbour coupling | No spatial coupling unless fibres are modeled as a network. |
| Wave support | Direction-dependent elastic waves. |
| Recovery | Reversible idealization. |
| Progressive hardening | Fibre recruitment gives strong hardening. |
| Finite deformation | Finite deformation and dispersion of orientations. |
| Weak field | Can share an isotropic tangent if orientations average isotropically. |
| Large deformation | Strong anisotropy/recruitment at large strain. |
| Mathematical strengths | Mature invariant fibre formulation. |
| Mathematical weaknesses | Needs structural tensors/fibre directions absent from minimal isotropic state. |
| PBUF compatibility | **ontology-change** — Independent preferred directions violate the current minimal isotropic closure unless derived. |
| Potential PBUF role | rejected pending anisotropy derivation |
| Sources | S06 |

## M01 — Prestretched membrane

| Required item | Evaluation |
|---|---|
| Physical mechanism | In-plane tension of a thin sheet restores transverse displacement. |
| Why it exists in nature | Drums, films and biological membranes. |
| Typical representation | `E_s=int_A W_s(a) dA; linearized transverse operator -T Delta w` |
| Stored energy | Surface stretching energy. |
| Interaction behavior | Surface stress divergence. |
| Neighbour coupling | Intrinsic within the two-dimensional sheet. |
| Wave support | Strong transverse membrane waves under positive tension. |
| Recovery | Elastic return while tension remains positive. |
| Progressive hardening | Geometric stiffening can occur with deflection. |
| Finite deformation | Large surface deformation possible. |
| Weak field | Small waves depend on prestress, not solely material modulus. |
| Large deformation | Wrinkling under compression and no bulk volumetric mode. |
| Mathematical strengths | Simple wave-support mechanism. |
| Mathematical weaknesses | A 2D carrier/embedding and prestress conflict with a generic 3D medium unless derived. |
| PBUF compatibility | **ontology-change** — Dimensional reduction is not authorized as the complete spacetime medium. |
| Potential PBUF role | rejected as complete law |
| Sources | S25 |

## M02 — Elastic plate/shell bending

| Required item | Evaluation |
|---|---|
| Physical mechanism | Stretching and curvature resistance jointly restore a thin surface. |
| Why it exists in nature | Plates, shells, capsules and cell walls. |
| Typical representation | `E=int_A[Ws(a)+B|b-b0|^2/2]dA` |
| Stored energy | Surface plus bending energy. |
| Interaction behavior | Fourth-order bending and second-order membrane operators. |
| Neighbour coupling | Intrinsic surface neighbour coupling. |
| Wave support | Dispersive flexural waves plus membrane waves. |
| Recovery | Elastic on a stable branch. |
| Progressive hardening | Geometric nonlinearities may harden or buckle/soften. |
| Finite deformation | Finite shell kinematics available. |
| Weak field | Flexural omega scales as q^2 without tension. |
| Large deformation | Buckling, multiple equilibria and dimension-specific response. |
| Mathematical strengths | Explicit curvature interaction and rich waves. |
| Mathematical weaknesses | Requires surface geometry, thickness and bending scale. |
| PBUF compatibility | **ontology-change** — Not a complete 3D isotropic medium without a derived shell reduction. |
| Potential PBUF role | rejected as complete law |
| Sources | S25 |

## V01 — Ideal/barotropic fluid

| Required item | Evaluation |
|---|---|
| Physical mechanism | Isotropic pressure depends on density; no static shear resistance. |
| Why it exists in nature | Liquids and gases near local equilibrium. |
| Typical representation | `sigma=-p(rho)I; p=rho^2 d(e/rho)/drho` |
| Stored energy | Internal energy depends on density. |
| Interaction behavior | Pressure gradient in momentum balance. |
| Neighbour coupling | Local continuum coupling through density gradients. |
| Wave support | Supports longitudinal sound; no elastic shear waves. |
| Recovery | Returns density via pressure but not shear shape. |
| Progressive hardening | Equation of state may stiffen in compression. |
| Finite deformation | Large flow kinematics natural; not elastic shape storage. |
| Weak field | Sound speed squared dp/drho must be positive. |
| Large deformation | Shocks possible; shear modulus is zero. |
| Mathematical strengths | Simple conservative propagation. |
| Mathematical weaknesses | Cannot reproduce generic frozen shear tangent or deformation memory. |
| PBUF compatibility | **incompatible** — Fails the required shear channel of objective C as a complete material. |
| Potential PBUF role | rejected complete law |
| Sources | S08 |

## V02 — Kelvin-Voigt viscoelastic solid

| Required item | Evaluation |
|---|---|
| Physical mechanism | An elastic spring and viscous dashpot act in parallel. |
| Why it exists in nature | Polymers, tissues and damped solids. |
| Typical representation | `sigma=E epsilon+eta epsilon_dot` |
| Stored energy | Elastic W=E epsilon^2/2 plus Rayleigh dissipation. |
| Interaction behavior | Local stress divergence. |
| Neighbour coupling | Via balance; no intrinsic length. |
| Wave support | Damped/dispersive elastic waves with inertia. |
| Recovery | Creep is bounded and unloading recovers asymptotically. |
| Progressive hardening | Linear version has no progressive hardening. |
| Finite deformation | Finite variants require objective rates. |
| Weak field | Recovers Hooke elasticity at low rate/static limit. |
| Large deformation | High-frequency response depends on clock rate and viscosity. |
| Mathematical strengths | Stable damping and simple recovery. |
| Mathematical weaknesses | Needs fundamental rate/relaxation data absent before emergent duration calibration. |
| PBUF compatibility | **conditional** — Could be downstream dissipation only after DURATION-001-compatible derivation. |
| Potential PBUF role | dissipative hybrid |
| Sources | S08 |

## V03 — Maxwell viscoelastic fluid

| Required item | Evaluation |
|---|---|
| Physical mechanism | A spring and dashpot in series relax stress and permit permanent flow. |
| Why it exists in nature | Polymer melts and stress-relaxing fluids. |
| Typical representation | `sigma_dot/E+sigma/eta=epsilon_dot` |
| Stored energy | Elastic spring energy plus viscous dissipation. |
| Interaction behavior | Local stress divergence. |
| Neighbour coupling | Via balance only. |
| Wave support | Frequency-dependent waves may be overdamped at low frequency. |
| Recovery | Does not recover total strain after unloading. |
| Progressive hardening | Linear model has no progressive hardening. |
| Finite deformation | Finite variants require objective stress rates. |
| Weak field | Has elastic response at high frequency. |
| Large deformation | Stress relaxes to zero and strain can remain. |
| Mathematical strengths | Canonical relaxation model. |
| Mathematical weaknesses | Fails full recovery and introduces a relaxation time. |
| PBUF compatibility | **incompatible** — Permanent flow conflicts with the reversible stored-deformation role. |
| Potential PBUF role | rejected complete law |
| Sources | S08 |

## V04 — Standard linear solid

| Required item | Evaluation |
|---|---|
| Physical mechanism | Parallel elastic and Maxwell branches combine instantaneous and relaxed stiffness. |
| Why it exists in nature | Broad relaxation spectra approximated by spring-dashpot networks. |
| Typical representation | `sigma+tau_sigma sigma_dot=E0(epsilon+tau_epsilon epsilon_dot)` |
| Stored energy | Recoverable elastic energy plus dissipation. |
| Interaction behavior | Local stress divergence. |
| Neighbour coupling | Via balance only. |
| Wave support | Damped waves with bounded low/high frequency moduli. |
| Recovery | Recovers to zero strain for zero stress. |
| Progressive hardening | Linear unless nonlinear springs are used. |
| Finite deformation | Objective finite variants exist. |
| Weak field | Positive relaxed modulus supplies weak recovery. |
| Large deformation | Additional relaxation scales and history. |
| Mathematical strengths | Better recovery than Maxwell and causal frequency response. |
| Mathematical weaknesses | Clock-dependent parameters are not frozen consequences. |
| PBUF compatibility | **conditional** — Possible downstream dissipative completion, not native energy law. |
| Potential PBUF role | dissipative hybrid |
| Sources | S08 |

## V05 — Oldroyd-B viscoelastic fluid

| Required item | Evaluation |
|---|---|
| Physical mechanism | A Newtonian solvent couples to a convected elastic conformation tensor. |
| Why it exists in nature | Dilute polymer solutions modeled as Hookean dumbbells. |
| Typical representation | `tau+lambda upper_convected(tau)=2eta_p D; sigma=-pI+2eta_sD+tau` |
| Stored energy | Conformation free energy plus solvent/polymer dissipation. |
| Interaction behavior | Stress divergence and advected internal tensor. |
| Neighbour coupling | Continuum coupling plus advection. |
| Wave support | Complex elastic/viscous modes; often diffusive or damped. |
| Recovery | Conformation relaxes, but material elements flow. |
| Progressive hardening | Hookean chains do not finitely harden. |
| Finite deformation | Objective rate supports large flow. |
| Weak field | Positive zero-shear viscosity, elastic high-rate response. |
| Large deformation | Unbounded chain stretch causes extensional singular behavior. |
| Mathematical strengths | Established objective viscoelastic tensor model. |
| Mathematical weaknesses | Adds an independent conformation tensor, solvent split and relaxation clock. |
| PBUF compatibility | **ontology-change** — Multiple sectors/internal state are unauthorized and total deformation is not recovered. |
| Potential PBUF role | rejected complete law |
| Sources | S09 |

## V06 — Two-fluid superfluid hydrodynamics

| Required item | Evaluation |
|---|---|
| Physical mechanism | Interpenetrating inviscid superfluid and viscous normal components carry distinct velocities. |
| Why it exists in nature | Quantum fluids below a superfluid transition. |
| Typical representation | `Mass, entropy and two momentum balances; quantized circulation constraints` |
| Stored energy | Thermodynamic internal energy, not a solid strain energy. |
| Interaction behavior | Pressure, entropy and mutual-friction couplings. |
| Neighbour coupling | Continuum neighbour coupling. |
| Wave support | First and second sound plus vortical excitations. |
| Recovery | No generic recovery of shear deformation. |
| Progressive hardening | Equation of state, not progressive elastic hardening. |
| Finite deformation | Large flows and vortices possible. |
| Weak field | Sound modes exist with positive thermodynamic derivatives. |
| Large deformation | Requires two velocities and quantum circulation structure. |
| Mathematical strengths | Rich lossless propagation. |
| Mathematical weaknesses | Contradicts one minimal state sector and lacks static shear storage. |
| PBUF compatibility | **ontology-change** — Interpenetrating components and quantum postulates are not authorized. |
| Potential PBUF role | rejected complete law |
| Sources | S10 |

## V07 — Gross-Pitaevskii/BEC continuum

| Required item | Evaluation |
|---|---|
| Physical mechanism | A complex order parameter has phase stiffness, interaction pressure and quantum-gradient energy. |
| Why it exists in nature | Dilute weakly interacting Bose condensates. |
| Typical representation | `i hbar psi_t=[-hbar^2 Delta/(2m)+g|psi|^2]psi` |
| Stored energy | E=int[hbar^2|Grad psi|^2/(2m)+g|psi|^4/2+V|psi|^2]dV |
| Interaction behavior | Hamiltonian field interaction. |
| Neighbour coupling | Intrinsic Laplacian/gradient coupling. |
| Wave support | Bogoliubov sound and dispersive short waves. |
| Recovery | Hamiltonian but not recovery of an elastic C state. |
| Progressive hardening | Nonlinear density pressure; not finite-strain hardening. |
| Finite deformation | Supports vortices and large phase gradients. |
| Weak field | Linearization gives acoustic dispersion. |
| Large deformation | Dispersion dominates at short scales; no shear modulus. |
| Mathematical strengths | Clean conservative wave medium. |
| Mathematical weaknesses | Imports complex field, hbar, mass and coupling; wrong native variable. |
| PBUF compatibility | **incompatible** — Would replace rather than instantiate the frozen C-based ontology. |
| Potential PBUF role | rejected |
| Sources | S10 |

## L01 — Harmonic crystal/phonon lattice

| Required item | Evaluation |
|---|---|
| Physical mechanism | Discrete sites coupled by quadratic neighbour springs generate collective phonons. |
| Why it exists in nature | Crystalline solids near equilibrium. |
| Typical representation | `E=sum_i p_i^2/(2m)+1/2 sum_ij u_i K_ij u_j` |
| Stored energy | Quadratic spring energy. |
| Interaction behavior | Discrete dynamical matrix. |
| Neighbour coupling | Explicit adjacency/range coupling. |
| Wave support | Acoustic and optical dispersive branches. |
| Recovery | Exact in the harmonic idealization. |
| Progressive hardening | None; harmonic tangent is constant. |
| Finite deformation | Invalid at large displacement without nonlinear potentials. |
| Weak field | Continuum long-wave limit is Hookean. |
| Large deformation | Brillouin-zone dispersion and lattice anisotropy. |
| Mathematical strengths | Transparent origin of waves and neighbour coupling. |
| Mathematical weaknesses | Inventing sites/atoms is forbidden; introduces spacing and preferred lattice. |
| PBUF compatibility | **ontology-change** — Only its homogenized continuum interaction structure is comparable. |
| Potential PBUF role | rejected microscopic mechanism |
| Sources | S11 |

## L02 — Fermi-Pasta-Ulam nonlinear lattice

| Required item | Evaluation |
|---|---|
| Physical mechanism | Neighbour springs include cubic/quartic nonlinearities, producing amplitude-dependent waves. |
| Why it exists in nature | Anharmonic crystals and canonical nonlinear-lattice dynamics. |
| Typical representation | `E=sum_i[p_i^2/2m+V(u_{i+1}-u_i)], V=kx^2/2+alpha x^3/3+beta x^4/4` |
| Stored energy | Anharmonic bond energy. |
| Interaction behavior | Nearest-neighbour forces. |
| Neighbour coupling | Explicit discrete adjacency. |
| Wave support | Nonlinear dispersion, mode coupling, solitons and recurrences. |
| Recovery | Conservative; recovery depends on phase/mode distribution. |
| Progressive hardening | Positive quartic terms harden at large amplitude. |
| Finite deformation | Large bond strain limited by potential validity. |
| Weak field | Harmonic long-wave limit. |
| Large deformation | Strong nonlinearity; cubic terms can destabilize one branch. |
| Mathematical strengths | Joins neighbour coupling and progressive hardening minimally. |
| Mathematical weaknesses | Discrete sites and coefficients violate ontology/FP-6 if taken literally. |
| PBUF compatibility | **ontology-change** — Its continuum gradient-hyperelastic limit is an admissible analogy only. |
| Potential PBUF role | hybrid exemplar |
| Sources | S14 |

## L03 — Central pair-potential continuum

| Required item | Evaluation |
|---|---|
| Physical mechanism | Distance-dependent conservative pair forces create equilibrium spacing and elastic moduli. |
| Why it exists in nature | Molecular crystals, colloids and particle networks. |
| Typical representation | `E=1/2 sum_{i!=j} V(|x_i-x_j|) or continuum double integral` |
| Stored energy | Pair potential energy. |
| Interaction behavior | Force is the distance derivative of V. |
| Neighbour coupling | Explicit distance-dependent coupling. |
| Wave support | Phonon-like waves about a stable minimum. |
| Recovery | Conservative return locally. |
| Progressive hardening | Depends on anharmonic curvature of V. |
| Finite deformation | Large separation/compression response potential-dependent. |
| Weak field | Hessian at equilibrium gives elastic constants. |
| Large deformation | May soften, fracture or collapse outside stable well. |
| Mathematical strengths | Clear conservative origin of neighbour forces. |
| Mathematical weaknesses | Requires constituents, pair distance and interaction scale absent from frozen ontology. |
| PBUF compatibility | **ontology-change** — Only the symmetric nonlocal-kernel limit can be retained. |
| Potential PBUF role | interaction analogy |
| Sources | S12 |

## C01 — Open-cell/bending-dominated foam

| Required item | Evaluation |
|---|---|
| Physical mechanism | Cell ribs bend, buckle and densify under compression. |
| Why it exists in nature | Polymer, metal and biological cellular solids. |
| Typical representation | `Homogenized W_eff(C;relative density, topology)` |
| Stored energy | Effective energy from beam/cell deformation. |
| Interaction behavior | Network force and moment balance. |
| Neighbour coupling | Explicit cell connectivity before homogenization. |
| Wave support | Elastic waves with band/dispersion effects. |
| Recovery | Recoverable only before buckling damage/plasticity. |
| Progressive hardening | Densification produces strong compressive hardening after a plateau. |
| Finite deformation | Very large compressive strain through cell collapse. |
| Weak field | Low effective moduli controlled by relative density. |
| Large deformation | Nonconvex plateau, buckling and hysteresis common. |
| Mathematical strengths | Topology controls response and permits densification hardening. |
| Mathematical weaknesses | Cells, voids and damage are extra structure; tension/compression asymmetry. |
| PBUF compatibility | **conditional** — Homogenized reversible branch only; literal foam ontology excluded. |
| Potential PBUF role | topology-driven exemplar |
| Sources | S15 |

## C02 — Affine entropic network

| Required item | Evaluation |
|---|---|
| Physical mechanism | A connected random network deforms approximately affinely and stores energy in strands. |
| Why it exists in nature | Rubber networks, gels and cytoskeletal networks. |
| Typical representation | `W=<w_chain(|F a|)> over orientation distribution` |
| Stored energy | Orientation-averaged strand energy. |
| Interaction behavior | Network node forces homogenize to stress. |
| Neighbour coupling | Explicit connectivity before homogenization. |
| Wave support | Elastic waves after homogenization. |
| Recovery | Reversible ideal network. |
| Progressive hardening | Finite chains or nonaffinity can harden. |
| Finite deformation | Finite deformation natural. |
| Weak field | Isotropic orientation average yields standard elastic tangent. |
| Large deformation | May become nonaffine, lose stability or approach chain limits. |
| Mathematical strengths | Explains convergence of many network models to invariant hyperelasticity. |
| Mathematical weaknesses | Microscopic strands and distribution are not PBUF primitives. |
| PBUF compatibility | **conditional** — Use only the homogenized invariant energy class. |
| Potential PBUF role | constitutive exemplar |
| Sources | S05, S07 |

## A01 — Rate-independent plasticity with work hardening

| Required item | Evaluation |
|---|---|
| Physical mechanism | Yielding produces irreversible strain; the yield surface expands with accumulated plastic work. |
| Why it exists in nature | Dislocation motion and irreversible rearrangement in solids. |
| Typical representation | `F=Fe Fp; yield f(tau,kappa)<=0; kappa_dot>=0` |
| Stored energy | Elastic energy W(Fe) plus dissipation, not total-strain potential. |
| Interaction behavior | Local stress plus plastic flow rule. |
| Neighbour coupling | No intrinsic neighbour coupling in classical plasticity. |
| Wave support | Elastic precursor waves; plastic waves/shocks possible. |
| Recovery | Does not recover plastic deformation. |
| Progressive hardening | Yes, through evolving yield strength, not reversible tangent hardening. |
| Finite deformation | Designed for large permanent deformation. |
| Weak field | Elastic tangent before yield. |
| Large deformation | History, path dependence and dissipation dominate. |
| Mathematical strengths | Mature hardening framework. |
| Mathematical weaknesses | Fails reversible recovery and requires internal history variables. |
| PBUF compatibility | **incompatible** — Work hardening is not the frozen progressive elastic hardening requirement. |
| Potential PBUF role | rejected |
| Sources | S08 |

## A02 — Continuum damage/phase-field fracture

| Required item | Evaluation |
|---|---|
| Physical mechanism | A damage field progressively degrades stiffness and regularizes cracks. |
| Why it exists in nature | Microcracking, decohesion and fracture. |
| Typical representation | `E=int[g(d)W(C)+Gc(d^2/(2ell)+ell|Grad d|^2/2)]dV` |
| Stored energy | Degraded elastic plus crack-surface energy. |
| Interaction behavior | Gradient damage equation couples neighbours. |
| Neighbour coupling | Intrinsic through Grad d. |
| Wave support | Waves attenuate/scatter as damage grows. |
| Recovery | Irreversibility constraint d_dot>=0 prevents recovery. |
| Progressive hardening | Softening, not progressive elastic hardening. |
| Finite deformation | Handles cracks and large deformation in variants. |
| Weak field | Undamaged tangent can match elasticity. |
| Large deformation | Loss of stiffness and localization; finite fracture energy. |
| Mathematical strengths | Mathematically controlled fracture regularization. |
| Mathematical weaknesses | Adds damage field, length and irreversible evolution. |
| PBUF compatibility | **incompatible** — Violates recovery and minimal frozen state as a complete mechanism. |
| Potential PBUF role | rejected |
| Sources | S19 |

## A03 — Shape-memory phase-transforming solid

| Required item | Evaluation |
|---|---|
| Physical mechanism | Stress/temperature switches among crystal variants, permitting recoverable transformation strain. |
| Why it exists in nature | Shape-memory alloys and martensitic transformations. |
| Typical representation | `Free energy with multiple wells in strain, temperature and phase fractions` |
| Stored energy | Multiwell thermoelastic energy plus dissipation. |
| Interaction behavior | Local stress; gradient interfacial terms in regularized models. |
| Neighbour coupling | Conditional through phase gradients/interfaces. |
| Wave support | Complex dispersive/attenuated waves across phases. |
| Recovery | Recovery requires thermal or stress cycling; hysteretic. |
| Progressive hardening | Plateaus and subsequent hardening, not monotone single-well response. |
| Finite deformation | Large recoverable strains possible. |
| Weak field | One well has an elastic tangent. |
| Large deformation | Multiple equilibria and hysteresis. |
| Mathematical strengths | Established large recoverable transformation mechanism. |
| Mathematical weaknesses | Requires temperature, phases and internal variables not in frozen minimal state. |
| PBUF compatibility | **ontology-change** — Cannot enter without authorized phase/thermal state enlargement. |
| Potential PBUF role | rejected |
| Sources | S20 |

## A04 — Jammed granular/contact network

| Required item | Evaluation |
|---|---|
| Physical mechanism | Repulsive contacts form a load-bearing network above a jamming threshold. |
| Why it exists in nature | Dense grains, emulsions and foams. |
| Typical representation | `E=sum_contacts V(overlap) for frictionless soft particles` |
| Stored energy | Contact energy; friction adds nonconservative forces. |
| Interaction behavior | Contact force balance. |
| Neighbour coupling | Explicit evolving adjacency. |
| Wave support | Acoustic modes emerge above jamming; anomalous soft modes near threshold. |
| Recovery | Frictionless elastic contacts may recover; rearrangements/hysteresis generally do not. |
| Progressive hardening | Hertz contacts harden with overlap; network changes can soften. |
| Finite deformation | Large rearrangements and contact changes. |
| Weak field | Moduli vanish or scale near the jamming point. |
| Large deformation | Non-smooth topology changes and marginal stability. |
| Mathematical strengths | Emergent rigidity from neighbour topology. |
| Mathematical weaknesses | Requires particles, contacts and evolving topology forbidden by the mission. |
| PBUF compatibility | **ontology-change** — Continuum scaling may inspire tests but cannot be imported literally. |
| Potential PBUF role | rejected microscopic mechanism |
| Sources | S21 |

## X01 — Magnetoelastic continuum

| Required item | Evaluation |
|---|---|
| Physical mechanism | Deformation and a magnetic field exchange energy, producing field-tunable stress. |
| Why it exists in nature | Magnetorheological elastomers and ferromagnetic solids. |
| Typical representation | `W=W(C,B) with invariant magnetoelastic couplings` |
| Stored energy | Coupled mechanical-field energy. |
| Interaction behavior | Mechanical stress and Maxwell-type field equations. |
| Neighbour coupling | Field mediates spatial interaction. |
| Wave support | Coupled magnetoacoustic waves. |
| Recovery | Potentially reversible without hysteresis. |
| Progressive hardening | Field can stiffen or soften depending on coupling. |
| Finite deformation | Finite magnetoelastic formulations exist. |
| Weak field | At zero field can reduce to ordinary elasticity. |
| Large deformation | Stability requires joint Hessian positivity; hysteresis common. |
| Mathematical strengths | Established conservative field-mediated coupling. |
| Mathematical weaknesses | Adds an independent field and electromagnetic ontology explicitly disallowed. |
| PBUF compatibility | **ontology-change** — Mathematical comparison only; no EM ontology may be assumed. |
| Potential PBUF role | rejected as native mechanism |
| Sources | S18 |

## X02 — Liquid-crystal elastomer

| Required item | Evaluation |
|---|---|
| Physical mechanism | A polymer network couples strain to an orientational director, enabling soft modes. |
| Why it exists in nature | Cross-linked liquid-crystal polymers. |
| Typical representation | `W=W(C,n)+Frank gradient energy in n` |
| Stored energy | Elastic plus orientational gradient energy. |
| Interaction behavior | Stress and director Euler-Lagrange equations. |
| Neighbour coupling | Director gradients couple neighbours. |
| Wave support | Coupled acoustic and orientational modes. |
| Recovery | May recover but can have domains/hysteresis. |
| Progressive hardening | Often soft or semisoft rather than progressively hardening. |
| Finite deformation | Large spontaneous anisotropic deformation. |
| Weak field | Director relaxation changes the apparent tangent. |
| Large deformation | Nonconvex multiwell behavior and domain formation. |
| Mathematical strengths | Rich coupled internal-field mechanics. |
| Mathematical weaknesses | Independent director violates minimal isotropy/state. |
| PBUF compatibility | **ontology-change** — Excluded unless orientation is derived from the existing q/C field. |
| Potential PBUF role | rejected |
| Sources | S22 |

## T01 — Negative-stiffness stabilized composite

| Required item | Evaluation |
|---|---|
| Physical mechanism | A locally unstable inclusion is stabilized by a positive host, yielding extreme effective response. |
| Why it exists in nature | Buckled structures and phase-transforming inclusions in composites. |
| Typical representation | `W_eff from host plus nonconvex inclusion energy` |
| Stored energy | Composite energy can be stable globally despite negative local tangent. |
| Interaction behavior | Coupling through host matrix. |
| Neighbour coupling | Explicit microstructural coupling. |
| Wave support | Unusual dispersion; stability requires positive total dynamic energy. |
| Recovery | Often metastable and hysteretic. |
| Progressive hardening | Can switch, soften or show negative incremental stiffness. |
| Finite deformation | Large snap-through events. |
| Weak field | Weak effective modulus can be enhanced or approach instability. |
| Large deformation | Nonconvexity and sensitivity to boundary conditions. |
| Mathematical strengths | Demonstrates architected non-standard response. |
| Mathematical weaknesses | Local negative tangent fails PBUF stability/progressive-hardening gates. |
| PBUF compatibility | **incompatible** — Useful elimination control, not an admissible complete law. |
| Potential PBUF role | rejected |
| Sources | S16 |

## T02 — Locally resonant elastic metamaterial

| Required item | Evaluation |
|---|---|
| Physical mechanism | Internal resonators exchange energy with a host to create band gaps and dynamic effective mass. |
| Why it exists in nature | Sonic/phononic architected materials. |
| Typical representation | `Host elasticity coupled to oscillator fields m_r q_ddot+k_r(q-u)=0` |
| Stored energy | Positive host and resonator energy. |
| Interaction behavior | Local resonator-host coupling plus host neighbour interaction. |
| Neighbour coupling | Host continuum plus internal coupling. |
| Wave support | Strongly dispersive waves and stop bands. |
| Recovery | Elastic idealization recovers; damping can be added. |
| Progressive hardening | Not amplitude hardening unless resonators are nonlinear. |
| Finite deformation | Finite motion variants exist but usually linearized. |
| Weak field | Long-wave effective parameters are frequency dependent. |
| Large deformation | Resonances can yield negative effective dynamic parameters. |
| Mathematical strengths | Powerful programmable wave control. |
| Mathematical weaknesses | Adds resonator degrees, masses, frequencies and microstructure. |
| PBUF compatibility | **ontology-change** — Cannot be native without deriving internal modes from existing state. |
| Potential PBUF role | rejected complete mechanism |
| Sources | S17 |

## H01 — Gradient hyperelastic hybrid

| Required item | Evaluation |
|---|---|
| Physical mechanism | A progressively hardening local invariant energy is combined with positive gradient interaction. |
| Why it exists in nature | Established strain-gradient nonlinear elasticity. |
| Typical representation | `E=int[Phi(I1,I2,I3)+ell^2|Grad C|^2/2]dV` |
| Stored energy | Local nonlinear plus gradient energy. |
| Interaction behavior | Fourth/second-order variational operator depending on chosen gradient measure. |
| Neighbour coupling | Intrinsic conservative neighbour coupling. |
| Wave support | Dispersive nonlinear elastic waves with positive inertia. |
| Recovery | Reversible on a stable branch. |
| Progressive hardening | Inherited from Phi; can be polynomial, exponential or barrier. |
| Finite deformation | Finite deformation if gradient term is objective. |
| Weak field | Common PBUF tangent plus q-dependent correction. |
| Large deformation | Combines large-amplitude hardening with short-scale stiffening. |
| Mathematical strengths | Fills both LAB-002 constitutive slots in one variational functional. |
| Mathematical weaknesses | Length scale, objective gradient and boundary data remain underived. |
| PBUF compatibility | **conditional** — A family for future testing, not a preferred model. |
| Potential PBUF role | hybrid mechanism |
| Sources | S23, S24 |

## H02 — Fibre-reinforced matrix hybrid

| Required item | Evaluation |
|---|---|
| Physical mechanism | A soft isotropic matrix shares load with hardening oriented fibres. |
| Why it exists in nature | Biological tissues and composites. |
| Typical representation | `W=Wmatrix(C)+Wfibres(C;A_a)` |
| Stored energy | Additive coupled hyperelastic energy. |
| Interaction behavior | Local stress; network extensions add coupling. |
| Neighbour coupling | None intrinsically in homogenized local form. |
| Wave support | Anisotropic waves. |
| Recovery | Reversible idealization. |
| Progressive hardening | Recruitment/finite fibres progressively harden. |
| Finite deformation | Finite deformation mature. |
| Weak field | Isotropic tangent possible only for isotropic orientation ensemble. |
| Large deformation | Strong anisotropy and locking at large stretch. |
| Mathematical strengths | Explains independent convergence to exponential/barrier responses. |
| Mathematical weaknesses | Structural tensors are additional state data. |
| PBUF compatibility | **ontology-change** — Conditional only if fibre directions emerge from C/q without ontology change. |
| Potential PBUF role | hybrid exemplar |
| Sources | S06, S07 |

## H03 — Visco-hyperelastic solid

| Required item | Evaluation |
|---|---|
| Physical mechanism | A finite-strain hyperelastic equilibrium network combines with relaxing overstress branches. |
| Why it exists in nature | Rubbers and soft tissues with rate dependence. |
| Typical representation | `sigma=dW_eq/dC+sum overstress_a; Qdot_a=G_a(C,Q_a)` |
| Stored energy | Equilibrium energy plus branch free energies and dissipation. |
| Interaction behavior | Local stress divergence. |
| Neighbour coupling | No intrinsic length unless gradient terms are added. |
| Wave support | Damped nonlinear elastic waves. |
| Recovery | Equilibrium branch recovers; transient branches relax. |
| Progressive hardening | Can inherit progressive hardening from W_eq. |
| Finite deformation | Finite-deformation objective formulations exist. |
| Weak field | Equilibrium tangent can match PBUF. |
| Large deformation | Rate-dependent hysteresis and multiple relaxation times. |
| Mathematical strengths | Separates equilibrium hardening from damping. |
| Mathematical weaknesses | Adds internal variables and clock scales. |
| PBUF compatibility | **conditional** — Only as a later duration-compatible dissipative completion. |
| Potential PBUF role | hybrid mechanism |
| Sources | S08, S09 |

## Family clustering

- **classical elastic:** E01
- **classical hyperelastic:** E02
- **isotropic hyperelastic:** E03, E04, E05
- **progressive hyperelastic:** E06, E09
- **finite-extensible hyperelastic:** E07, E08
- **strain-measure hyperelastic:** E10
- **constraint class:** E11
- **gradient/nonlocal:** W01, W02, W03
- **generalized continuum:** W04
- **wave/nonlinear:** W05
- **fibre/network:** F01, F02
- **anisotropic hyperelastic:** F03
- **membrane/shell:** M01, M02
- **fluid:** V01
- **viscoelastic:** V02, V03, V04
- **viscoelastic fluid:** V05
- **quantum/fluid:** V06, V07
- **lattice/crystal:** L01, L02
- **field-mediated/network:** L03
- **cellular/topological:** C01, C02
- **inelastic/history:** A01, A02
- **phase-transforming:** A03
- **contact/topological:** A04
- **field-coupled:** X01
- **field-coupled/anisotropic:** X02
- **metamaterial:** T01, T02
- **hybrid:** H01, H02, H03
- **naturally C-based local energies:** E02, E03, E04, E05, E06, E07, E08, E09, E10
- **intrinsic neighbour-coupled:** E02, E03, E04, E05, E06, E07, E09, E10, W01, W02, W03, W04, W05, F01, F02, M01, M02, V02, V07, L01, L02, L03, C01, C02, A01, A02, A04, T01, H01, H02, H03
- **finite-extensible/barrier:** E02, E03, E04, E08, F02, V05, V07, C02, H02
- **irreversible/history-dependent:** V01, V03, V07, A01, A02, T01

## Elimination analysis

- **Immediately incompatible (V01, V03, V07, A01, A02, T01):** these fail recovery, stable positive tangent, generic shear storage, or the frozen native-variable requirement. They remain controls in the matrix rather than being silently removed.
- **Require ontology/state changes (W04, F03, M01, M02, V05, V06, L01, L02, L03, A03, A04, X01, X02, T02, H02):** these require independent directors, microrotations, particles, phases, fields, resonators, fibres, surfaces, or multiple fluid sectors. Their homogenized mathematics may still illuminate a conditional continuum class, but their literal mechanisms are not admissible.
- **Conditional/naturally representable (E01, E02, E03, E04, E05, E06, E07, E08, E09, E10, E11, W01, W02, W03, W05, F01, F02, V02, V04, C01, C02, H01, H03):** these can be written as local invariant energies, objective gradient/nonlocal operators, or later dissipative completions. Every one still has named gates—tensor stability, parameter derivation, objective gradient structure, boundary data, or duration-compatible kinetics.
- **Unconditionally compatible (none):** none. The frozen milestones deliberately do not select a formula, length, kinetic operator, endpoint, or parameter set.
- **Mathematically incomplete/unknown:** no catalogue row is left formula-free, but all wave claims are incomplete until kinetic closure; all scalar/ray claims are incomplete until an invariant tensor lift passes spectral stability. This is reported as conditional rather than hidden.

## Convergence analysis

Independent physical stories repeatedly collapse to a small number of mathematical structures:

1. **Invariant hyperelasticity:** rubber phenomenology, affine networks, biological matrices, and homogenized cellular systems converge to `W=Phi(I1,I2,I3)` and `P_C=DW`. Their microscopic stories do not survive as PBUF ontology.
2. **Common quadratic tangent:** Hooke, neo-Hooke, Mooney-Rivlin, Ogden, polynomial, exponential, finite-chain, and multiwell single-branch models all reduce near the unloaded state to the same positive volumetric/shear quadratic form. This is a universality class, not model selection.
3. **Superquadratic or barrier hardening:** polynomial anharmonicity, exponential fibre recruitment, Gent/Arruda-Boyce/worm-like-chain finite extensibility, and contact densification all produce an increasing tangent. Two distinct subclasses remain: coercive growth on an unbounded domain and divergence at a finite admissible boundary.
4. **Positive spatial operator:** gradient elasticity, nonlocal kernels, peridynamics, lattices, membranes/shells, and field energies all generate neighbour communication through a positive self-adjoint operator (or a stable discrete dynamical matrix). Its generic continuum signature is a positive `q`-dependent stiffness. This does not derive inertia.
5. **Internal-variable relaxation:** Maxwell, standard-linear-solid, Oldroyd-B, visco-hyperelastic, phase-transforming, and damage models use extra internal variables and clock scales. This shared structure explains attenuation and memory but is not contained in a state-local `W(C)`.
6. **Topology-controlled response:** fibres, lattices, foams, jammed contacts, and metamaterials obtain response from adjacency and architecture. Homogenization can yield ordinary invariant energy, gradient terms, anisotropy, or resonant internal modes; which limit occurs is additional constitutive data.

The deepest common admissible structure exposed by the survey is therefore not a named material: it is a stable invariant local energy plus a separately testable positive neighbour operator, embedded in an authorized kinetic/balance closure. Hardening is an inequality on the local tangent; spatial dispersion is an inequality on the interaction symbol. They are independent axes.

## Recommendations for future constitutive development

Retain a portfolio rather than a winner. Tensor-lift the local polynomial, exponential, principal-stretch, and barrier families and test stress-free normalization, spectral Hessian positivity, rank-one convexity/strong ellipticity, coercivity, and every SPD boundary. In parallel, test objective gradient and symmetric nonlocal interaction families for positivity, well-posed boundary conditions, and their long-wave limits. Only then combine representatives factorially (local energy × neighbour operator × authorized kinetic closure) so that hardening and dispersion are not confounded. Keep every scale and coefficient symbolic until independently derived under FP-6. Do not import microscopic stories, modify V11 or weak lensing, or fit observations.

## Sources and reproducibility

The source register records established provenance for each family. It is not observational evidence for PBUF. `material_discovery001.py` is the canonical generator for the JSON catalogue, CSV matrix, clustering, eliminations, and this report.

- **S01:** [Ogden-model historical review](https://pmc.ncbi.nlm.nih.gov/articles/PMC9421375/)
- **S02:** [Review of isotropic incompressible hyperelastic models](https://pubmed.ncbi.nlm.nih.gov/26087063/)
- **S03:** [Finite elasticity of elastomeric materials](https://doi.org/10.1093/oso/9780198864721.003.0031)
- **S04:** [Gent model review](https://doi.org/10.1016/j.ijmecsci.2014.05.010)
- **S05:** [Limited-stretch rubber models](https://arxiv.org/abs/2005.09648)
- **S06:** [Holzapfel-Gasser-Ogden arterial model](https://doi.org/10.1016/S0021-9290(00)00053-4)
- **S07:** [Worm-like chain model](https://doi.org/10.1115/1.2798296)
- **S08:** [Continuum mechanics and thermodynamics excerpt](https://assets.cambridge.org/97811070/89952/excerpt/9781107089952_excerpt.pdf)
- **S09:** [Oldroyd-B constitutive description](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/fast-flow-of-an-oldroydb-model-fluid-through-a-narrow-slowly-varying-contraction/FA6DC7D4141EE96A72E1D7C83677A400)
- **S10:** [Gross-Pitaevskii superfluid dynamics](https://doi.org/10.1103/PhysRevA.97.013627)
- **S11:** [Lattice material homogenization review](https://pmc.ncbi.nlm.nih.gov/articles/PMC8778170/)
- **S12:** [Nonlocal continuum theories review](https://pmc.ncbi.nlm.nih.gov/articles/PMC9863499/)
- **S13:** [Peridynamic continuum theory](https://doi.org/10.1016/S0022-5096(99)00029-0)
- **S14:** [Fermi-Pasta-Ulam nonlinear lattice](https://doi.org/10.2172/4376203)
- **S15:** [Cellular solids mechanics](https://doi.org/10.1017/CBO9781139878326)
- **S16:** [Negative-stiffness inclusions](https://doi.org/10.1103/PhysRevLett.86.2890)
- **S17:** [Locally resonant sonic materials](https://doi.org/10.1126/science.289.5485.1734)
- **S18:** [Magnetoelasticity framework](https://doi.org/10.1098/rspa.2013.0058)
- **S19:** [Phase-field fracture](https://doi.org/10.1016/S0022-5096(00)00028-0)
- **S20:** [Shape-memory alloy constitutive review](https://doi.org/10.1016/S0079-6425(03)00012-7)
- **S21:** [Jamming review](https://doi.org/10.1103/RevModPhys.82.2633)
- **S22:** [Liquid-crystal elastomer theory](https://doi.org/10.1103/PhysRevE.65.061710)
- **S23:** [Gradient elasticity overview](https://doi.org/10.1007/978-3-319-06385-3)
- **S24:** [Nonlinear elastic waves](https://doi.org/10.1016/j.ijnonlinmec.2014.06.007)
- **S25:** [Shell theory review](https://doi.org/10.1016/j.apm.2008.10.015)
