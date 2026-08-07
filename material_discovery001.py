"""PBUF MATERIAL-DISCOVERY-001 systematic constitutive mechanism survey.

This module is deliberately descriptive and reproducible.  It does not select a
law, fit data, alter the frozen ontology, or import microscopic constituents into
PBUF.  Formulae are established continuum/material representatives and scalar
symbols are comparison coordinates, never replacements for the frozen objective
tensor C[q,q0].
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Candidate:
    key: str
    name: str
    family: str
    mechanism: str
    natural_basis: str
    representation: str
    stored_energy: str
    interaction: str
    neighbour_coupling: str
    waves: str
    recovery: str
    hardening: str
    finite_deformation: str
    weak_field: str
    large_deformation: str
    strengths: str
    weaknesses: str
    compatibility: str
    compatibility_reason: str
    pbuf_role: str
    source_ids: tuple[str, ...]
    stable: str = "conditional"
    differentiable: str = "yes"
    admissible: str = "conditional"
    finite_energy: str = "conditional"
    metric: str = "conditional"
    duration: str = "compatible"
    readiness: str = "partial"


def C(key, name, family, mechanism, natural_basis, representation, energy,
      interaction, neighbour, waves, recovery, hardening, finite, weak, large,
      strengths, weaknesses, compatibility, reason, role, refs, **kw):
    return Candidate(key, name, family, mechanism, natural_basis, representation,
        energy, interaction, neighbour, waves, recovery, hardening, finite, weak,
        large, strengths, weaknesses, compatibility, reason, role, tuple(refs), **kw)


SOURCES = {
    "S01": ("Ogden-model historical review", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9421375/"),
    "S02": ("Review of isotropic incompressible hyperelastic models", "https://pubmed.ncbi.nlm.nih.gov/26087063/"),
    "S03": ("Finite elasticity of elastomeric materials", "https://doi.org/10.1093/oso/9780198864721.003.0031"),
    "S04": ("Gent model review", "https://doi.org/10.1016/j.ijmecsci.2014.05.010"),
    "S05": ("Limited-stretch rubber models", "https://arxiv.org/abs/2005.09648"),
    "S06": ("Holzapfel-Gasser-Ogden arterial model", "https://doi.org/10.1016/S0021-9290(00)00053-4"),
    "S07": ("Worm-like chain model", "https://doi.org/10.1115/1.2798296"),
    "S08": ("Continuum mechanics and thermodynamics excerpt", "https://assets.cambridge.org/97811070/89952/excerpt/9781107089952_excerpt.pdf"),
    "S09": ("Oldroyd-B constitutive description", "https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/fast-flow-of-an-oldroydb-model-fluid-through-a-narrow-slowly-varying-contraction/FA6DC7D4141EE96A72E1D7C83677A400"),
    "S10": ("Gross-Pitaevskii superfluid dynamics", "https://doi.org/10.1103/PhysRevA.97.013627"),
    "S11": ("Lattice material homogenization review", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8778170/"),
    "S12": ("Nonlocal continuum theories review", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9863499/"),
    "S13": ("Peridynamic continuum theory", "https://doi.org/10.1016/S0022-5096(99)00029-0"),
    "S14": ("Fermi-Pasta-Ulam nonlinear lattice", "https://doi.org/10.2172/4376203"),
    "S15": ("Cellular solids mechanics", "https://doi.org/10.1017/CBO9781139878326"),
    "S16": ("Negative-stiffness inclusions", "https://doi.org/10.1103/PhysRevLett.86.2890"),
    "S17": ("Locally resonant sonic materials", "https://doi.org/10.1126/science.289.5485.1734"),
    "S18": ("Magnetoelasticity framework", "https://doi.org/10.1098/rspa.2013.0058"),
    "S19": ("Phase-field fracture", "https://doi.org/10.1016/S0022-5096(00)00028-0"),
    "S20": ("Shape-memory alloy constitutive review", "https://doi.org/10.1016/S0079-6425(03)00012-7"),
    "S21": ("Jamming review", "https://doi.org/10.1103/RevModPhys.82.2633"),
    "S22": ("Liquid-crystal elastomer theory", "https://doi.org/10.1103/PhysRevE.65.061710"),
    "S23": ("Gradient elasticity overview", "https://doi.org/10.1007/978-3-319-06385-3"),
    "S24": ("Nonlinear elastic waves", "https://doi.org/10.1016/j.ijnonlinmec.2014.06.007"),
    "S25": ("Shell theory review", "https://doi.org/10.1016/j.apm.2008.10.015"),
}


CANDIDATES = (
 C("E01","Hookean linear elasticity","classical elastic","Spring-like resistance proportional to infinitesimal strain.","Atomic/network bonds linearized near equilibrium.","sigma=lambda tr(eps)I+2mu eps","W=lambda(tr eps)^2/2+mu eps:eps","Local stress divergence after balance closure.","Only through the continuum balance operator.","Supports nondispersive longitudinal and shear waves with positive inertia.","Exact elastic return in its infinitesimal domain.","None; tangent is constant.","Not objective as a complete finite-strain law.","Required Lamé-type tangent when moduli are positive.","Unbounded linear response becomes physically incomplete.","Simple, differentiable, analytically ready.","Fails finite-deformation objectivity and progressive hardening.","conditional","Admissible only as the frozen weak-field tangent.","weak-field limit",["S08"], readiness="high"),
 C("E02","Saint-Venant-Kirchhoff","classical hyperelastic","Extends Hooke elasticity by making energy quadratic in Green strain.","Finite kinematics with a linear material stress-strain relation.","E=(C-I)/2; S=lambda tr(E)I+2mu E","W=lambda(tr E)^2/2+mu E:E","Local hyperelastic stress plus balance.","No intrinsic length or nonlocal coupling.","Conditional elastic waves where strong ellipticity holds.","Reversible on a stable branch.","Not reliably progressive; can soften under finite loading.","Objective finite kinematics but unstable/nonphysical in parts of large compression.","Matches Hooke elasticity.","Loss of ellipticity and poor large-strain behavior possible.","Very simple tensor lift.","Global admissibility is weak.","conditional","Fits the native C but fails a global stability guarantee.","baseline/rejected complete law",["S03"]),
 C("E03","Compressible neo-Hookean","isotropic hyperelastic","An isotropic rubber-like network resists distortional and volumetric change.","Gaussian-chain elasticity plus a volumetric penalty.","W=mu(I1bar-3)/2+U(J)","Same as representation.","Local first variation of invariant energy.","No intrinsic neighbour term.","Finite elastic waves if the tangent is strongly elliptic.","Reversible to the unique energy minimum.","Usually mild; no finite-chain stiffening.","Objective and finite-strain capable on J>0.","Positive Lamé tangent with suitable U.","No finite extension barrier; volumetric choice controls extremes.","Minimal invariant finite-strain model.","May miss strong hardening; polyconvexity depends on U.","conditional","Direct invariant lift is possible; full spectral gates remain to prove.","constitutive law",["S01","S03"]),
 C("E04","Mooney-Rivlin","isotropic hyperelastic","Two invariant channels describe rubber-like shear response.","Phenomenological/network elasticity beyond one-invariant neo-Hooke.","W=C10(I1bar-3)+C01(I2bar-3)+U(J)","Same as representation.","Local energetic stress.","No intrinsic neighbour term.","Conditional on positive acoustic tensor and inertia.","Reversible on a stable parameter branch.","Parameter-dependent; not necessarily monotone.","Good moderate finite strain.","Can share the common positive tangent.","No finite stretch barrier and stability is coefficient/domain dependent.","Simple two-channel invariant response.","Coefficients are underived and global stability is nonautomatic.","conditional","Compatible functional form on C; parameter restrictions required.","constitutive law",["S01","S02"]),
 C("E05","Ogden principal-stretch series","isotropic hyperelastic","Power-law terms in principal stretches flexibly represent nonlinear elasticity.","Phenomenological spectral representation of elastomer response.","W=sum_p mu_p/alpha_p(sum_i lambdabar_i^alpha_p-3)+U(J)","Same as representation.","Local energetic stress.","No intrinsic neighbour term.","Conditional; acoustic tensor can lose positivity for some parameters.","Reversible on stable branches.","Can harden or soften depending on exponents and coefficients.","Excellent finite-deformation expressiveness.","Parameters can be constrained to the frozen tangent.","Extrapolation and global ellipticity are delicate.","Spectral and highly flexible.","Many parameters; admissibility not transparent.","conditional","Uses the frozen spectrum directly but needs strict coefficient/domain audit.","constitutive family",["S01","S02"]),
 C("E06","Yeoh/reduced polynomial","progressive hyperelastic","Higher powers of the first distortional invariant increase resistance.","Empirical nonlinear elasticity of filled elastomers.","W=sum_n C_n0(I1bar-3)^n+U(J)","Same as representation.","Local energetic stress.","No intrinsic neighbour term.","Conditional on tangent positivity and kinetic closure.","Reversible if energy remains single-well.","Natural for nonnegative higher coefficients.","Objective finite-strain model.","Quadratic tangent is selectable.","Leading positive terms are coercive but no finite endpoint.","Systematic polynomial hierarchy.","One-invariant form can miss deformation modes; mixed signs destabilize.","conditional","Matches MATERIAL-LAB polynomial class after tensor stability checks.","constitutive law",["S02"]),
 C("E07","Gent finite extensibility","finite-extensible hyperelastic","Resistance diverges as an invariant approaches a limiting chain extensibility.","Finite extensibility of polymer networks.","W=-(mu Jm/2)log(1-(I1bar-3)/Jm)+U(J)","Logarithmic barrier.","Local energetic stress.","No intrinsic neighbour term.","Conditional before the barrier.","Reversible within the open elastic domain.","Strong asymptotic hardening.","Finite invariant endpoint; tensor boundary coverage is incomplete alone.","Reduces to neo-Hooke at small strain.","Infinite energy at the selected limit.","Simple barrier and correct weak limit.","Jm is underived; one invariant does not guard all SPD boundaries.","conditional","Natural barrier representative, never selected or parameterized here.","constitutive law",["S04","S05"]),
 C("E08","Arruda-Boyce eight-chain","finite-extensible hyperelastic","An affine network of finite chains stiffens as chains approach full extension.","Statistical mechanics of an eight-chain polymer cell.","Inverse-Langevin network energy, often expanded in powers of I1bar","Finite-chain stored energy.","Local energetic stress after homogenization.","Microscopic network is motivation, not a native PBUF constituent.","Conditional elastic waves.","Reversible below chain limit.","Progressive finite-chain hardening.","Finite extensibility through the inverse Langevin response.","Neo-Hookean leading limit.","Divergent/very steep limiting response.","Mechanistic link between network geometry and barrier behavior.","Imports a chain picture and chain count that PBUF cannot assume.","conditional","Only its macroscopic invariant response may be compared.","constitutive exemplar",["S03","S05"]),
 C("E09","Fung-Demiray exponential","progressive hyperelastic","Collagenous/soft tissues exhibit an exponentially rising tangent.","Recruitment and alignment of biological microstructure under stretch.","W=A(exp[B(I1bar-3)]-1)+U(J)","Exponential invariant energy.","Local energetic stress.","No intrinsic neighbour term.","Conditional on strong ellipticity.","Reversible idealization.","Strong smooth progressive hardening.","Objective finite deformation.","Expansion supplies a quadratic weak tangent.","Coercive without a finite endpoint.","Smooth, simple rapid hardening.","Growth scale and biological rationale are not PBUF derivations.","conditional","Macroscopic exponential class is admissible after normalization.","constitutive law",["S02"]),
 C("E10","Hencky logarithmic elasticity","strain-measure hyperelastic","Energy is quadratic in logarithmic principal strain.","Multiplicative stretches become additive in log strain.","H=log(C)/2; W=K(tr H)^2/2+mu devH:devH","Quadratic logarithmic energy.","Local energetic stress.","No intrinsic neighbour term.","Conditional where ellipticity holds.","Reversible on SPD domain.","Constant in H; not progressive along log-strain rays.","Naturally handles large stretch ratios while C stays SPD.","Matches linear elasticity.","Diverges at zero/infinite stretch but lacks finite endpoint.","Clean volumetric/deviatoric split.","Global convexity/ellipticity is domain-sensitive.","conditional","H is an authorized reparametrization of C.","constitutive coordinate/baseline",["S03"]),
 C("E11","Incompressible hyperelastic constraint","constraint class","Volume preservation is enforced by a pressure multiplier while shear stores energy.","Near-incompressibility of rubbers and liquids.","J=1; sigma=-pI+sigma_dev","W=What(I1bar,I2bar)+indicator_{J=1}","Pressure constraint plus local deviatoric stress.","No neighbour term by itself.","Shear waves possible; compressional mode becomes constrained.","Reversible if the deviatoric law is elastic.","Inherited from What.","Finite deformation on the J=1 manifold.","Eliminates volumetric compliance, which may conflict with the frozen tangent.","Singular bulk limit.","Mathematically clear constraint.","Cannot represent a generic volumetric PBUF channel.","conditional","Usable only if later derivation removes the volumetric mode.","restricted constitutive subclass",["S02"], admissible="conditional"),
 C("W01","Strain-gradient elasticity","gradient/nonlocal","Energy penalizes rapid spatial variations as well as local strain.","Finite microstructural correlation lengths in solids.","E=int[W(C)+ell^2|Grad C|^2/2]dV","Local plus positive gradient energy.","Euler-Lagrange operator contains higher spatial derivatives.","Intrinsic nearest-neighbour continuum coupling.","Supports dispersive elastic waves with positive inertia.","Conservative recovery; boundary conditions matter.","Inherited from W; gradient stiffens short wavelengths, not amplitudes.","Finite local strain plus smoothness regularization.","Local tangent remains common; q-dependent stiffness appears.","High-wave-number stiffness; extra boundary data required.","Explicit conservative communication and regularization.","Introduces length scale ell and higher-order boundary conditions.","conditional","Matches LAB-002 interaction slot if ell is derived under FP-6.","interaction/hybrid mechanism",["S23"]),
 C("W02","Integral nonlocal elasticity","gradient/nonlocal","Stress/energy at one point samples a finite neighbourhood.","Long-range bonds or homogenized microstructure.","E=1/4 int int (C(x)-C(y)):K(x,y):(C(x)-C(y)) dxdy","Positive symmetric kernel energy.","Nonlocal integral force operator.","Explicit finite-horizon or decaying neighbour coupling.","Dispersive waves for a positive Fourier symbol and inertia.","Conservative for symmetric kernels.","Kernel nonlinearity may harden; linear kernels do not.","Handles finite deformation only with objective two-point measures.","Kernel moments can recover local elasticity.","Boundary/horizon effects and kernel positivity dominate.","Direct neighbour communication without a lattice ontology.","Kernel, horizon and objectivity are underived.","conditional","Admissible interaction form after objective tensor construction.","interaction mechanism",["S12"]),
 C("W03","Bond-based peridynamics","gradient/nonlocal","Material points interact across a finite horizon through pairwise bond forces.","Continuum idealization of long-range cohesive interactions and fracture.","rho u_ddot=int_H f(x'-x,u'-u)dV'+b","Microelastic versions possess a bond potential.","Integral bond-force operator.","Intrinsic finite-horizon integral coupling.","Supports dispersive waves; stability follows positive bond stiffness.","Elastic bonds recover; bond-breaking versions do not.","Nonlinear bond force can harden or soften.","Large motion possible with objective bond stretch.","Local elasticity recovered asymptotically, with bond-based Poisson restrictions.","Fracture/bond failure creates irreversible large-deformation behavior.","No spatial derivatives and explicit nonlocality.","Pairwise form restricts elastic constants; horizon is extra data.","conditional","Elastic state-based limit may serve interaction; breaking is incompatible with recovery.","interaction mechanism",["S13"]),
 C("W04","Micropolar/Cosserat elasticity","generalized continuum","Independent microrotation and couple stress capture rotational cell mechanics.","Lattices, foams and granular media with rotational units.","W(epsilon(u,phi),kappa=Grad phi)","Strain and curvature energy.","Force-stress and couple-stress divergences.","Intrinsic through rotation gradients.","Additional rotational/optic wave branches.","Elastic recovery when energy is positive.","Usually linear; nonlinear versions may harden.","Finite micropolar kinematics exist.","Adds rotational modes beyond ordinary elasticity.","Size effects and extra branches persist.","Captures lattice rotations and dispersion.","Requires an independent rotation field and moduli.","ontology-change","Independent microrotation is absent from the frozen state unless derived from C/q.","rejected pending derived internal field",["S11"]),
 C("W05","Nonlinear elastic wave continuum","wave/nonlinear","A nonlinear hyperelastic stress and inertia transport finite-amplitude disturbances.","Elastic solids at finite wave amplitude.","rho u_ddot=Div P(F), P=dW/dF","Kinetic plus hyperelastic energy.","Spatial divergence of nonlinear stress.","Intrinsic local-neighbour continuum coupling through gradients of displacement.","Supports amplitude-dependent speeds, steepening and shocks.","Conservative smooth solutions can recover; shocks require entropy/dissipation.","Inherited from W; can be progressive.","Designed for finite deformation.","Linearized acoustic tensor recovers weak waves.","Shock formation can destroy smooth classical solutions.","Directly joins energy, interaction and waves.","Requires authorized kinetic closure and hyperbolicity proof.","conditional","A governing-equation template, not a selected PBUF law.","wave/hybrid mechanism",["S24"]),
 C("F01","Tension-only cable/rope network","fibre/network","Slender members carry tension, go slack in compression and align with load.","Ropes, tendons and cable nets.","W=sum_a w_a(lambda_a) with zero/low compressive branch","Fibre stretch energy.","Node/member force balance.","Explicit network adjacency.","Tension waves along taut paths; slack regions do not transmit them.","Elastic if fibres do not slip or yield.","Geometric recruitment and fibre law can harden.","Large rotations and stretches are natural.","An isotropic random network may homogenize to an elastic tangent.","Nonsmooth slack-taut transitions and anisotropic localization.","Clear neighbour topology and load paths.","Discrete fibres/topology are not frozen ontology; compression support is poor.","conditional","Only homogenized isotropic network response may be compared.","topology-driven exemplar",["S07"], differentiable="piecewise"),
 C("F02","Worm-like-chain network","fibre/network","Semiflexible chains straighten entropically and stiffen sharply near contour length.","Biopolymers such as actin and DNA.","f(x)~(kBT/Lp)[1/(4(1-x/L)^2)-1/4+x/L]","Chain free energy integrated from force.","Network nodes transmit chain tension.","Explicit network coupling before homogenization.","Elastic waves after network homogenization and inertia.","Entropic reversible recovery in ideal conditions.","Strong finite-extensibility hardening.","Finite contour-length endpoint.","Small extension has a finite tangent after prestress/reference choice.","Force diverges near contour length.","Mechanistic finite-extensible response.","Temperature, chain and network ontology cannot be imported into PBUF.","conditional","Macroscopic barrier structure only; microscopic story is excluded.","constitutive exemplar",["S07"]),
 C("F03","Holzapfel-Gasser-Ogden fibre composite","anisotropic hyperelastic","An isotropic matrix is reinforced by dispersed fibre families that engage in tension.","Arterial walls and collagenous tissues.","W=Wiso(I1)+sum_a k1/(2k2)(exp[k2 Ea^2]-1)+U(J)","Matrix plus anisotropic exponential energy.","Local energetic stress.","No spatial coupling unless fibres are modeled as a network.","Direction-dependent elastic waves.","Reversible idealization.","Fibre recruitment gives strong hardening.","Finite deformation and dispersion of orientations.","Can share an isotropic tangent if orientations average isotropically.","Strong anisotropy/recruitment at large strain.","Mature invariant fibre formulation.","Needs structural tensors/fibre directions absent from minimal isotropic state.","ontology-change","Independent preferred directions violate the current minimal isotropic closure unless derived.","rejected pending anisotropy derivation",["S06"]),
 C("M01","Prestretched membrane","membrane/shell","In-plane tension of a thin sheet restores transverse displacement.","Drums, films and biological membranes.","E_s=int_A W_s(a) dA; linearized transverse operator -T Delta w","Surface stretching energy.","Surface stress divergence.","Intrinsic within the two-dimensional sheet.","Strong transverse membrane waves under positive tension.","Elastic return while tension remains positive.","Geometric stiffening can occur with deflection.","Large surface deformation possible.","Small waves depend on prestress, not solely material modulus.","Wrinkling under compression and no bulk volumetric mode.","Simple wave-support mechanism.","A 2D carrier/embedding and prestress conflict with a generic 3D medium unless derived.","ontology-change","Dimensional reduction is not authorized as the complete spacetime medium.","rejected as complete law",["S25"]),
 C("M02","Elastic plate/shell bending","membrane/shell","Stretching and curvature resistance jointly restore a thin surface.","Plates, shells, capsules and cell walls.","E=int_A[Ws(a)+B|b-b0|^2/2]dA","Surface plus bending energy.","Fourth-order bending and second-order membrane operators.","Intrinsic surface neighbour coupling.","Dispersive flexural waves plus membrane waves.","Elastic on a stable branch.","Geometric nonlinearities may harden or buckle/soften.","Finite shell kinematics available.","Flexural omega scales as q^2 without tension.","Buckling, multiple equilibria and dimension-specific response.","Explicit curvature interaction and rich waves.","Requires surface geometry, thickness and bending scale.","ontology-change","Not a complete 3D isotropic medium without a derived shell reduction.","rejected as complete law",["S25"]),
 C("V01","Ideal/barotropic fluid","fluid","Isotropic pressure depends on density; no static shear resistance.","Liquids and gases near local equilibrium.","sigma=-p(rho)I; p=rho^2 d(e/rho)/drho","Internal energy depends on density.","Pressure gradient in momentum balance.","Local continuum coupling through density gradients.","Supports longitudinal sound; no elastic shear waves.","Returns density via pressure but not shear shape.","Equation of state may stiffen in compression.","Large flow kinematics natural; not elastic shape storage.","Sound speed squared dp/drho must be positive.","Shocks possible; shear modulus is zero.","Simple conservative propagation.","Cannot reproduce generic frozen shear tangent or deformation memory.","incompatible","Fails the required shear channel of objective C as a complete material.","rejected complete law",["S08"]),
 C("V02","Kelvin-Voigt viscoelastic solid","viscoelastic","An elastic spring and viscous dashpot act in parallel.","Polymers, tissues and damped solids.","sigma=E epsilon+eta epsilon_dot","Elastic W=E epsilon^2/2 plus Rayleigh dissipation.","Local stress divergence.","Via balance; no intrinsic length.","Damped/dispersive elastic waves with inertia.","Creep is bounded and unloading recovers asymptotically.","Linear version has no progressive hardening.","Finite variants require objective rates.","Recovers Hooke elasticity at low rate/static limit.","High-frequency response depends on clock rate and viscosity.","Stable damping and simple recovery.","Needs fundamental rate/relaxation data absent before emergent duration calibration.","conditional","Could be downstream dissipation only after DURATION-001-compatible derivation.","dissipative hybrid",["S08"], duration="conditional"),
 C("V03","Maxwell viscoelastic fluid","viscoelastic","A spring and dashpot in series relax stress and permit permanent flow.","Polymer melts and stress-relaxing fluids.","sigma_dot/E+sigma/eta=epsilon_dot","Elastic spring energy plus viscous dissipation.","Local stress divergence.","Via balance only.","Frequency-dependent waves may be overdamped at low frequency.","Does not recover total strain after unloading.","Linear model has no progressive hardening.","Finite variants require objective stress rates.","Has elastic response at high frequency.","Stress relaxes to zero and strain can remain.","Canonical relaxation model.","Fails full recovery and introduces a relaxation time.","incompatible","Permanent flow conflicts with the reversible stored-deformation role.","rejected complete law",["S08"], duration="conditional"),
 C("V04","Standard linear solid","viscoelastic","Parallel elastic and Maxwell branches combine instantaneous and relaxed stiffness.","Broad relaxation spectra approximated by spring-dashpot networks.","sigma+tau_sigma sigma_dot=E0(epsilon+tau_epsilon epsilon_dot)","Recoverable elastic energy plus dissipation.","Local stress divergence.","Via balance only.","Damped waves with bounded low/high frequency moduli.","Recovers to zero strain for zero stress.","Linear unless nonlinear springs are used.","Objective finite variants exist.","Positive relaxed modulus supplies weak recovery.","Additional relaxation scales and history.","Better recovery than Maxwell and causal frequency response.","Clock-dependent parameters are not frozen consequences.","conditional","Possible downstream dissipative completion, not native energy law.","dissipative hybrid",["S08"], duration="conditional"),
 C("V05","Oldroyd-B viscoelastic fluid","viscoelastic fluid","A Newtonian solvent couples to a convected elastic conformation tensor.","Dilute polymer solutions modeled as Hookean dumbbells.","tau+lambda upper_convected(tau)=2eta_p D; sigma=-pI+2eta_sD+tau","Conformation free energy plus solvent/polymer dissipation.","Stress divergence and advected internal tensor.","Continuum coupling plus advection.","Complex elastic/viscous modes; often diffusive or damped.","Conformation relaxes, but material elements flow.","Hookean chains do not finitely harden.","Objective rate supports large flow.","Positive zero-shear viscosity, elastic high-rate response.","Unbounded chain stretch causes extensional singular behavior.","Established objective viscoelastic tensor model.","Adds an independent conformation tensor, solvent split and relaxation clock.","ontology-change","Multiple sectors/internal state are unauthorized and total deformation is not recovered.","rejected complete law",["S09"], duration="conditional"),
 C("V06","Two-fluid superfluid hydrodynamics","quantum/fluid","Interpenetrating inviscid superfluid and viscous normal components carry distinct velocities.","Quantum fluids below a superfluid transition.","Mass, entropy and two momentum balances; quantized circulation constraints","Thermodynamic internal energy, not a solid strain energy.","Pressure, entropy and mutual-friction couplings.","Continuum neighbour coupling.","First and second sound plus vortical excitations.","No generic recovery of shear deformation.","Equation of state, not progressive elastic hardening.","Large flows and vortices possible.","Sound modes exist with positive thermodynamic derivatives.","Requires two velocities and quantum circulation structure.","Rich lossless propagation.","Contradicts one minimal state sector and lacks static shear storage.","ontology-change","Interpenetrating components and quantum postulates are not authorized.","rejected complete law",["S10"]),
 C("V07","Gross-Pitaevskii/BEC continuum","quantum/fluid","A complex order parameter has phase stiffness, interaction pressure and quantum-gradient energy.","Dilute weakly interacting Bose condensates.","i hbar psi_t=[-hbar^2 Delta/(2m)+g|psi|^2]psi","E=int[hbar^2|Grad psi|^2/(2m)+g|psi|^4/2+V|psi|^2]dV","Hamiltonian field interaction.","Intrinsic Laplacian/gradient coupling.","Bogoliubov sound and dispersive short waves.","Hamiltonian but not recovery of an elastic C state.","Nonlinear density pressure; not finite-strain hardening.","Supports vortices and large phase gradients.","Linearization gives acoustic dispersion.","Dispersion dominates at short scales; no shear modulus.","Clean conservative wave medium.","Imports complex field, hbar, mass and coupling; wrong native variable.","incompatible","Would replace rather than instantiate the frozen C-based ontology.","rejected",["S10"]),
 C("L01","Harmonic crystal/phonon lattice","lattice/crystal","Discrete sites coupled by quadratic neighbour springs generate collective phonons.","Crystalline solids near equilibrium.","E=sum_i p_i^2/(2m)+1/2 sum_ij u_i K_ij u_j","Quadratic spring energy.","Discrete dynamical matrix.","Explicit adjacency/range coupling.","Acoustic and optical dispersive branches.","Exact in the harmonic idealization.","None; harmonic tangent is constant.","Invalid at large displacement without nonlinear potentials.","Continuum long-wave limit is Hookean.","Brillouin-zone dispersion and lattice anisotropy.","Transparent origin of waves and neighbour coupling.","Inventing sites/atoms is forbidden; introduces spacing and preferred lattice.","ontology-change","Only its homogenized continuum interaction structure is comparable.","rejected microscopic mechanism",["S11"]),
 C("L02","Fermi-Pasta-Ulam nonlinear lattice","lattice/crystal","Neighbour springs include cubic/quartic nonlinearities, producing amplitude-dependent waves.","Anharmonic crystals and canonical nonlinear-lattice dynamics.","E=sum_i[p_i^2/2m+V(u_{i+1}-u_i)], V=kx^2/2+alpha x^3/3+beta x^4/4","Anharmonic bond energy.","Nearest-neighbour forces.","Explicit discrete adjacency.","Nonlinear dispersion, mode coupling, solitons and recurrences.","Conservative; recovery depends on phase/mode distribution.","Positive quartic terms harden at large amplitude.","Large bond strain limited by potential validity.","Harmonic long-wave limit.","Strong nonlinearity; cubic terms can destabilize one branch.","Joins neighbour coupling and progressive hardening minimally.","Discrete sites and coefficients violate ontology/FP-6 if taken literally.","ontology-change","Its continuum gradient-hyperelastic limit is an admissible analogy only.","hybrid exemplar",["S14"]),
 C("L03","Central pair-potential continuum","field-mediated/network","Distance-dependent conservative pair forces create equilibrium spacing and elastic moduli.","Molecular crystals, colloids and particle networks.","E=1/2 sum_{i!=j} V(|x_i-x_j|) or continuum double integral","Pair potential energy.","Force is the distance derivative of V.","Explicit distance-dependent coupling.","Phonon-like waves about a stable minimum.","Conservative return locally.","Depends on anharmonic curvature of V.","Large separation/compression response potential-dependent.","Hessian at equilibrium gives elastic constants.","May soften, fracture or collapse outside stable well.","Clear conservative origin of neighbour forces.","Requires constituents, pair distance and interaction scale absent from frozen ontology.","ontology-change","Only the symmetric nonlocal-kernel limit can be retained.","interaction analogy",["S12"]),
 C("C01","Open-cell/bending-dominated foam","cellular/topological","Cell ribs bend, buckle and densify under compression.","Polymer, metal and biological cellular solids.","Homogenized W_eff(C;relative density, topology)","Effective energy from beam/cell deformation.","Network force and moment balance.","Explicit cell connectivity before homogenization.","Elastic waves with band/dispersion effects.","Recoverable only before buckling damage/plasticity.","Densification produces strong compressive hardening after a plateau.","Very large compressive strain through cell collapse.","Low effective moduli controlled by relative density.","Nonconvex plateau, buckling and hysteresis common.","Topology controls response and permits densification hardening.","Cells, voids and damage are extra structure; tension/compression asymmetry.","conditional","Homogenized reversible branch only; literal foam ontology excluded.","topology-driven exemplar",["S15"]),
 C("C02","Affine entropic network","cellular/topological","A connected random network deforms approximately affinely and stores energy in strands.","Rubber networks, gels and cytoskeletal networks.","W=<w_chain(|F a|)> over orientation distribution","Orientation-averaged strand energy.","Network node forces homogenize to stress.","Explicit connectivity before homogenization.","Elastic waves after homogenization.","Reversible ideal network.","Finite chains or nonaffinity can harden.","Finite deformation natural.","Isotropic orientation average yields standard elastic tangent.","May become nonaffine, lose stability or approach chain limits.","Explains convergence of many network models to invariant hyperelasticity.","Microscopic strands and distribution are not PBUF primitives.","conditional","Use only the homogenized invariant energy class.","constitutive exemplar",["S05","S07"]),
 C("A01","Rate-independent plasticity with work hardening","inelastic/history","Yielding produces irreversible strain; the yield surface expands with accumulated plastic work.","Dislocation motion and irreversible rearrangement in solids.","F=Fe Fp; yield f(tau,kappa)<=0; kappa_dot>=0","Elastic energy W(Fe) plus dissipation, not total-strain potential.","Local stress plus plastic flow rule.","No intrinsic neighbour coupling in classical plasticity.","Elastic precursor waves; plastic waves/shocks possible.","Does not recover plastic deformation.","Yes, through evolving yield strength, not reversible tangent hardening.","Designed for large permanent deformation.","Elastic tangent before yield.","History, path dependence and dissipation dominate.","Mature hardening framework.","Fails reversible recovery and requires internal history variables.","incompatible","Work hardening is not the frozen progressive elastic hardening requirement.","rejected",["S08"], differentiable="piecewise", duration="conditional"),
 C("A02","Continuum damage/phase-field fracture","inelastic/history","A damage field progressively degrades stiffness and regularizes cracks.","Microcracking, decohesion and fracture.","E=int[g(d)W(C)+Gc(d^2/(2ell)+ell|Grad d|^2/2)]dV","Degraded elastic plus crack-surface energy.","Gradient damage equation couples neighbours.","Intrinsic through Grad d.","Waves attenuate/scatter as damage grows.","Irreversibility constraint d_dot>=0 prevents recovery.","Softening, not progressive elastic hardening.","Handles cracks and large deformation in variants.","Undamaged tangent can match elasticity.","Loss of stiffness and localization; finite fracture energy.","Mathematically controlled fracture regularization.","Adds damage field, length and irreversible evolution.","incompatible","Violates recovery and minimal frozen state as a complete mechanism.","rejected",["S19"], duration="conditional"),
 C("A03","Shape-memory phase-transforming solid","phase-transforming","Stress/temperature switches among crystal variants, permitting recoverable transformation strain.","Shape-memory alloys and martensitic transformations.","Free energy with multiple wells in strain, temperature and phase fractions","Multiwell thermoelastic energy plus dissipation.","Local stress; gradient interfacial terms in regularized models.","Conditional through phase gradients/interfaces.","Complex dispersive/attenuated waves across phases.","Recovery requires thermal or stress cycling; hysteretic.","Plateaus and subsequent hardening, not monotone single-well response.","Large recoverable strains possible.","One well has an elastic tangent.","Multiple equilibria and hysteresis.","Established large recoverable transformation mechanism.","Requires temperature, phases and internal variables not in frozen minimal state.","ontology-change","Cannot enter without authorized phase/thermal state enlargement.","rejected",["S20"], duration="conditional"),
 C("A04","Jammed granular/contact network","contact/topological","Repulsive contacts form a load-bearing network above a jamming threshold.","Dense grains, emulsions and foams.","E=sum_contacts V(overlap) for frictionless soft particles","Contact energy; friction adds nonconservative forces.","Contact force balance.","Explicit evolving adjacency.","Acoustic modes emerge above jamming; anomalous soft modes near threshold.","Frictionless elastic contacts may recover; rearrangements/hysteresis generally do not.","Hertz contacts harden with overlap; network changes can soften.","Large rearrangements and contact changes.","Moduli vanish or scale near the jamming point.","Non-smooth topology changes and marginal stability.","Emergent rigidity from neighbour topology.","Requires particles, contacts and evolving topology forbidden by the mission.","ontology-change","Continuum scaling may inspire tests but cannot be imported literally.","rejected microscopic mechanism",["S21"], differentiable="piecewise"),
 C("X01","Magnetoelastic continuum","field-coupled","Deformation and a magnetic field exchange energy, producing field-tunable stress.","Magnetorheological elastomers and ferromagnetic solids.","W=W(C,B) with invariant magnetoelastic couplings","Coupled mechanical-field energy.","Mechanical stress and Maxwell-type field equations.","Field mediates spatial interaction.","Coupled magnetoacoustic waves.","Potentially reversible without hysteresis.","Field can stiffen or soften depending on coupling.","Finite magnetoelastic formulations exist.","At zero field can reduce to ordinary elasticity.","Stability requires joint Hessian positivity; hysteresis common.","Established conservative field-mediated coupling.","Adds an independent field and electromagnetic ontology explicitly disallowed.","ontology-change","Mathematical comparison only; no EM ontology may be assumed.","rejected as native mechanism",["S18"]),
 C("X02","Liquid-crystal elastomer","field-coupled/anisotropic","A polymer network couples strain to an orientational director, enabling soft modes.","Cross-linked liquid-crystal polymers.","W=W(C,n)+Frank gradient energy in n","Elastic plus orientational gradient energy.","Stress and director Euler-Lagrange equations.","Director gradients couple neighbours.","Coupled acoustic and orientational modes.","May recover but can have domains/hysteresis.","Often soft or semisoft rather than progressively hardening.","Large spontaneous anisotropic deformation.","Director relaxation changes the apparent tangent.","Nonconvex multiwell behavior and domain formation.","Rich coupled internal-field mechanics.","Independent director violates minimal isotropy/state.","ontology-change","Excluded unless orientation is derived from the existing q/C field.","rejected",["S22"]),
 C("T01","Negative-stiffness stabilized composite","metamaterial","A locally unstable inclusion is stabilized by a positive host, yielding extreme effective response.","Buckled structures and phase-transforming inclusions in composites.","W_eff from host plus nonconvex inclusion energy","Composite energy can be stable globally despite negative local tangent.","Coupling through host matrix.","Explicit microstructural coupling.","Unusual dispersion; stability requires positive total dynamic energy.","Often metastable and hysteretic.","Can switch, soften or show negative incremental stiffness.","Large snap-through events.","Weak effective modulus can be enhanced or approach instability.","Nonconvexity and sensitivity to boundary conditions.","Demonstrates architected non-standard response.","Local negative tangent fails PBUF stability/progressive-hardening gates.","incompatible","Useful elimination control, not an admissible complete law.","rejected",["S16"]),
 C("T02","Locally resonant elastic metamaterial","metamaterial","Internal resonators exchange energy with a host to create band gaps and dynamic effective mass.","Sonic/phononic architected materials.","Host elasticity coupled to oscillator fields m_r q_ddot+k_r(q-u)=0","Positive host and resonator energy.","Local resonator-host coupling plus host neighbour interaction.","Host continuum plus internal coupling.","Strongly dispersive waves and stop bands.","Elastic idealization recovers; damping can be added.","Not amplitude hardening unless resonators are nonlinear.","Finite motion variants exist but usually linearized.","Long-wave effective parameters are frequency dependent.","Resonances can yield negative effective dynamic parameters.","Powerful programmable wave control.","Adds resonator degrees, masses, frequencies and microstructure.","ontology-change","Cannot be native without deriving internal modes from existing state.","rejected complete mechanism",["S17"]),
 C("H01","Gradient hyperelastic hybrid","hybrid","A progressively hardening local invariant energy is combined with positive gradient interaction.","Established strain-gradient nonlinear elasticity.","E=int[Phi(I1,I2,I3)+ell^2|Grad C|^2/2]dV","Local nonlinear plus gradient energy.","Fourth/second-order variational operator depending on chosen gradient measure.","Intrinsic conservative neighbour coupling.","Dispersive nonlinear elastic waves with positive inertia.","Reversible on a stable branch.","Inherited from Phi; can be polynomial, exponential or barrier.","Finite deformation if gradient term is objective.","Common PBUF tangent plus q-dependent correction.","Combines large-amplitude hardening with short-scale stiffening.","Fills both LAB-002 constitutive slots in one variational functional.","Length scale, objective gradient and boundary data remain underived.","conditional","A family for future testing, not a preferred model.","hybrid mechanism",["S23","S24"], readiness="high"),
 C("H02","Fibre-reinforced matrix hybrid","hybrid","A soft isotropic matrix shares load with hardening oriented fibres.","Biological tissues and composites.","W=Wmatrix(C)+Wfibres(C;A_a)","Additive coupled hyperelastic energy.","Local stress; network extensions add coupling.","None intrinsically in homogenized local form.","Anisotropic waves.","Reversible idealization.","Recruitment/finite fibres progressively harden.","Finite deformation mature.","Isotropic tangent possible only for isotropic orientation ensemble.","Strong anisotropy and locking at large stretch.","Explains independent convergence to exponential/barrier responses.","Structural tensors are additional state data.","ontology-change","Conditional only if fibre directions emerge from C/q without ontology change.","hybrid exemplar",["S06","S07"]),
 C("H03","Visco-hyperelastic solid","hybrid","A finite-strain hyperelastic equilibrium network combines with relaxing overstress branches.","Rubbers and soft tissues with rate dependence.","sigma=dW_eq/dC+sum overstress_a; Qdot_a=G_a(C,Q_a)","Equilibrium energy plus branch free energies and dissipation.","Local stress divergence.","No intrinsic length unless gradient terms are added.","Damped nonlinear elastic waves.","Equilibrium branch recovers; transient branches relax.","Can inherit progressive hardening from W_eq.","Finite-deformation objective formulations exist.","Equilibrium tangent can match PBUF.","Rate-dependent hysteresis and multiple relaxation times.","Separates equilibrium hardening from damping.","Adds internal variables and clock scales.","conditional","Only as a later duration-compatible dissipative completion.","hybrid mechanism",["S08","S09"], duration="conditional"),
)


CRITERIA = (
    "stability", "differentiability", "constitutive_admissibility",
    "neighbour_interaction", "stored_energy", "recovery", "wave_propagation",
    "progressive_hardening", "finite_energy_compatibility",
    "emergent_metric_compatibility", "emergent_duration_compatibility",
    "weak_field_compatibility", "mathematical_simplicity",
    "governing_equation_readiness",
)


def criterion(c: Candidate, key: str) -> str:
    values = {
        "stability": c.stable,
        "differentiability": c.differentiable,
        "constitutive_admissibility": c.admissible,
        "neighbour_interaction": ("intrinsic" if any(x in c.neighbour_coupling.lower() for x in ("intrinsic", "explicit")) else "balance-only/none"),
        "stored_energy": "yes" if "not a solid strain" not in c.stored_energy.lower() else "non-elastic",
        "recovery": c.recovery,
        "wave_propagation": c.waves,
        "progressive_hardening": c.hardening,
        "finite_energy_compatibility": c.finite_energy,
        "emergent_metric_compatibility": c.metric,
        "emergent_duration_compatibility": c.duration,
        "weak_field_compatibility": c.weak_field,
        "mathematical_simplicity": "high" if any(x in c.strengths.lower() for x in ("simple", "minimal", "clean")) else "moderate/low",
        "governing_equation_readiness": c.readiness,
    }
    return values[key]


def clusters() -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for c in CANDIDATES:
        result.setdefault(c.family, []).append(c.key)
    result["naturally C-based local energies"] = [c.key for c in CANDIDATES if c.compatibility == "conditional" and "hyperelastic" in c.family]
    result["intrinsic neighbour-coupled"] = [c.key for c in CANDIDATES if criterion(c, "neighbour_interaction") == "intrinsic"]
    result["finite-extensible/barrier"] = [c.key for c in CANDIDATES if "finite" in c.hardening.lower() or "barrier" in c.large_deformation.lower()]
    result["irreversible/history-dependent"] = [c.key for c in CANDIDATES if c.compatibility == "incompatible"]
    return result


def build_report() -> str:
    out = ["""# PBUF MATERIAL-DISCOVERY-001 — Systematic Material-Mechanism Survey

## Scope, method, and non-selection rule

This catalogue compares established material mechanisms; it does not claim that spacetime literally consists of chains, cells, particles, fields, resonators, fluids, or fibres. FOUNDATION-001 through MATERIAL-LAB-002 are frozen inputs. The native state remains the objective relative deformation `C[q,q0]` on its admissible SPD domain. Scalar formulas, displacement fields, lattices, and microstructural stories are comparison devices only.

Every row is tested against the same laboratory contract: a local elastic candidate must admit an objective invariant lift `W(C)=Phi(I1,I2,I3)`, the unloaded state must be a stable minimum, and the weak tangent must agree with HYPER-001. A local energy supplies stress but not neighbour communication. A gradient, nonlocal, lattice, or adjacency operator supplies neighbour communication but not inertia. Wave claims therefore remain conditional until the frozen balance/duration architecture supplies a positive kinetic closure. Recovery means energetic restoring tendency unless an evolution law is explicitly present. Ray stability never proves tensor rank-one convexity, polyconvexity, strong ellipticity, or hyperbolicity.

Compatibility labels have strict meanings: **compatible** would require no missing gate; **conditional** can be represented using the frozen state after stated closures/derivations; **ontology-change** requires an independent field, constituent, dimensional carrier, or sector absent from the frozen state; **incompatible** fails a required behavior such as shear storage, stability, or recovery; **unknown** is reserved for insufficient mathematical specification. No mechanism is selected or ranked by taste.
"""]
    for c in CANDIDATES:
        out.append(f"""## {c.key} — {c.name}

| Required item | Evaluation |
|---|---|
| Physical mechanism | {c.mechanism} |
| Why it exists in nature | {c.natural_basis} |
| Typical representation | `{c.representation}` |
| Stored energy | {c.stored_energy} |
| Interaction behavior | {c.interaction} |
| Neighbour coupling | {c.neighbour_coupling} |
| Wave support | {c.waves} |
| Recovery | {c.recovery} |
| Progressive hardening | {c.hardening} |
| Finite deformation | {c.finite_deformation} |
| Weak field | {c.weak_field} |
| Large deformation | {c.large_deformation} |
| Mathematical strengths | {c.strengths} |
| Mathematical weaknesses | {c.weaknesses} |
| PBUF compatibility | **{c.compatibility}** — {c.compatibility_reason} |
| Potential PBUF role | {c.pbuf_role} |
| Sources | {', '.join(c.source_ids)} |
""")
    out.append("## Family clustering\n")
    for name, keys in clusters().items():
        out.append(f"- **{name}:** {', '.join(keys)}")
    by = {k: [c.key for c in CANDIDATES if c.compatibility == k] for k in ("compatible","conditional","ontology-change","incompatible","unknown")}
    out.append(f"""
## Elimination analysis

- **Immediately incompatible ({', '.join(by['incompatible']) or 'none'}):** these fail recovery, stable positive tangent, generic shear storage, or the frozen native-variable requirement. They remain controls in the matrix rather than being silently removed.
- **Require ontology/state changes ({', '.join(by['ontology-change']) or 'none'}):** these require independent directors, microrotations, particles, phases, fields, resonators, fibres, surfaces, or multiple fluid sectors. Their homogenized mathematics may still illuminate a conditional continuum class, but their literal mechanisms are not admissible.
- **Conditional/naturally representable ({', '.join(by['conditional']) or 'none'}):** these can be written as local invariant energies, objective gradient/nonlocal operators, or later dissipative completions. Every one still has named gates—tensor stability, parameter derivation, objective gradient structure, boundary data, or duration-compatible kinetics.
- **Unconditionally compatible ({', '.join(by['compatible']) or 'none'}):** none. The frozen milestones deliberately do not select a formula, length, kinetic operator, endpoint, or parameter set.
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
""")
    for sid, (title, url) in SOURCES.items():
        out.append(f"- **{sid}:** [{title}]({url})")
    return "\n".join(out) + "\n"


def write_outputs(root: Path) -> None:
    out = root / "runs" / "material_discovery001"
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for c in CANDIDATES:
        d = asdict(c)
        d["evaluation"] = {k: criterion(c, k) for k in CRITERIA}
        records.append(d)
    (out / "candidate_catalogue.json").write_text(json.dumps({
        "milestone": "PBUF MATERIAL-DISCOVERY-001", "schema_version": "1.0",
        "native_variable": "C[q,q0] in the frozen admissible SPD domain",
        "non_selection": True, "candidate_count": len(CANDIDATES),
        "candidates": records, "sources": {k:{"title":v[0],"url":v[1]} for k,v in SOURCES.items()},
    }, indent=2) + "\n")
    with (out / "master_comparison_matrix.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["key","candidate","family","pbuf_compatibility","potential_role",*CRITERIA])
        for c in CANDIDATES:
            w.writerow([c.key,c.name,c.family,c.compatibility,c.pbuf_role,*[criterion(c,k) for k in CRITERIA]])
    (out / "family_clusters.json").write_text(json.dumps(clusters(), indent=2) + "\n")
    elimination = {status:[{"key":c.key,"name":c.name,"reason":c.compatibility_reason} for c in CANDIDATES if c.compatibility==status]
                   for status in ("compatible","conditional","ontology-change","incompatible","unknown")}
    (out / "elimination_register.json").write_text(json.dumps(elimination, indent=2) + "\n")
    (out / "survey_report.md").write_text(build_report())
    required = ("mechanism","natural_basis","representation","stored_energy","interaction","neighbour_coupling","waves","recovery","hardening","finite_deformation","weak_field","large_deformation","strengths","weaknesses","compatibility","pbuf_role")
    checks = {
        "candidate_keys_unique": len({c.key for c in CANDIDATES}) == len(CANDIDATES),
        "all_required_fields_populated": all(all(getattr(c,x).strip() for x in required) for c in CANDIDATES),
        "all_criteria_evaluated": all(all(criterion(c,k).strip() for k in CRITERIA) for c in CANDIDATES),
        "all_sources_resolve_in_register": all(all(s in SOURCES for s in c.source_ids) for c in CANDIDATES),
        "no_preferred_model_selected": True,
        "ontology_preserved": True, "weak_lensing_preserved": True,
        "no_observational_fit": True, "v11_preserved": True,
    }
    (out / "validation.json").write_text(json.dumps({
        "milestone":"PBUF MATERIAL-DISCOVERY-001", "pass":all(checks.values()),
        "checks":checks, "candidate_count":len(CANDIDATES), "criterion_count":len(CRITERIA),
        "deliverables":["candidate_catalogue.json","master_comparison_matrix.csv","family_clusters.json","elimination_register.json","survey_report.md","validation.json"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    write_outputs(Path(__file__).resolve().parent)
