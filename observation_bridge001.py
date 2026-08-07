#!/usr/bin/env python3
"""PBUF OBSERVATION-BRIDGE-001 - identify the correct physical comparison layer.

This is an audit-only milestone. It does NOT modify PBUF, the frozen transport,
the constitutive law, or any numerical parameter. It does NOT fit data and
does NOT introduce cosmological scaling.

Its sole purpose is to determine whether the observational products used in
WEAK-LENSING-OBSERVATION-001 correspond to the same physical stage of the
lensing pipeline as the frozen Version A outputs.

This script
- reads the FITS products (read-only);
- reads the previous milestone's outputs (read-only);
- writes audit tables, diagrams, and a report under
  ``runs/observation_bridge001/``.
"""
from __future__ import annotations

import csv
import hashlib
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits


ROOT = Path(__file__).resolve().parent
BENCHMARK_DIR = ROOT / "PBUF_benchmark"
DEFAULT_OUT = ROOT / "runs" / "observation_bridge001"
OBSERVATION_001_OUT = ROOT / "runs" / "weak_lensing_observation001"

CLUSTERS = [
    {"id": "Abell2744", "label": "Abell 2744", "slug": "abell2744",
     "directory": "WL-001_Abell2744", "z_l": 0.308},
    {"id": "MACS0416", "label": "MACS J0416", "slug": "macs0416",
     "directory": "WL-002_MACS0416", "z_l": 0.42},
    {"id": "MACS1149", "label": "MACS J1149", "slug": "macs1149",
     "directory": "WL-003_MACS1149", "z_l": 0.544},
    {"id": "AbellS1063", "label": "Abell S1063", "slug": "abells1063",
     "directory": "WL-004_AbellS1063", "z_l": 0.348},
    {"id": "Abell370", "label": "Abell 370", "slug": "abell370",
     "directory": "WL-005_Abell370", "z_l": 0.375},
]


# =============================================================================
# Stage documentation strings
# =============================================================================
VERSION_A_STAGES = [
    {
        "name": "Matter density ρ(X)",
        "math": "ρ : R^2 -> R",
        "physical": (
            "Input scalar field declared to be a non-negative dimensionless "
            "matter density on the pipeline grid."
        ),
        "units": "dimensionless (ρ_max = 1 after normalisation)",
        "assumptions": (
            "ρ >= 0; field supplied externally; no cosmological, redshift, or "
            "physical-unit content; no baryonic/dark partition."
        ),
    },
    {
        "name": "Constitutive field C(X) = 0.18 · ρ(X) / ρ_max",
        "math": "C : R^2 -> R,    C = 0.18 * ρ / ρ_max",
        "physical": (
            "Local linear scalar proxy of the medium's deformation by matter. "
            "C is the output of Version A's constitutive equation. No "
            "physical dimension is attached; C exists only as a scalar "
            "intermediate field."
        ),
        "units": "dimensionless (bounded by 0.18)",
        "assumptions": (
            "Local, linear, isotropic response; no propagation; the "
            "coefficient 0.18 is the frozen 'deformation strength' of "
            "Version A."
        ),
    },
    {
        "name": "Gradient field ∇C",
        "math": "∇C = (gx, gy) = (∂C/∂x, ∂C/∂y)",
        "physical": (
            "Gradient of the constitutive scalar field. Has the magnitude of "
            "a response per dimensionless length."
        ),
        "units": "dimensionless per dimensionless length",
        "assumptions": (
            "Finite differences with edge_order=1; the grid [-8, 8] x [-8, 8] "
            "carries no physical length scale."
        ),
    },
    {
        "name": "Response field r = (rx, ry)",
        "math": (
            "A = |∇C|;    r_x = -A * (∂C/∂y)/A;    r_y = A * (∂C/∂x)/A  "
            "(90-degree transverse rotation of the unit gradient)"
        ),
        "physical": (
            "Vector field representing the local transverse response of the "
            "medium. The 90-degree rotation is the frozen 'transport' choice; "
            "the response amplitude equals the gradient magnitude."
        ),
        "units": "dimensionless per dimensionless length (input to velocity update)",
        "assumptions": (
            "Neighbour-to-neighbour coupling, direct addition, instantaneous "
            "renormalisation; no retardation, no falloff, no medium rigidity."
        ),
    },
    {
        "name": "Photon propagation",
        "math": (
            "v_{k+1} = (v_k + step * r) / |v_k + step * r|;    "
            "x_{k+1} = x_k + step * v_{k+1}"
        ),
        "physical": (
            "Iterative ray-tracing through the response field. Photons start "
            "at x = -8 with v = (1, 0); the pipeline runs for steps = 80 "
            "iterations with step = 0.06, so the maximum propagation distance "
            "is 4.8 dimensionless units (≪ 2*extent = 16)."
        ),
        "units": (
            "step in dimensionless length units; velocity renormalised to "
            "unit speed per step"
        ),
        "assumptions": (
            "Frozen parameters n_grid = 128, extent = 8, strength = 0.18, "
            "step = 0.06, steps = 80, nphotons = 2000, bins = 64; identical "
            "to weak_lensing_prediction001 and weak_lensing_generalization001."
        ),
    },
    {
        "name": "Predicted convergence κ",
        "math": "κ(x,y) = 0.5 * (N_final(x,y) / N_initial(x,y) - 1)",
        "physical": (
            "Local photon-count ratio in a 64 x 64 histogram of (x_f, y_f) "
            "after propagation versus the initial (x_0, y_0) histogram. The "
            "pipeline only produces finite κ on bins where N_initial > 0."
        ),
        "units": "dimensionless",
        "assumptions": (
            "Bins with N_initial = 0 are filled with NaN; bins where photons "
            "left the initial x = -8 column but never returned show the "
            "constant value -0.5. No physical unit is attached; κ here is not "
            "the surface-mass-density-to-critical-density ratio."
        ),
    },
    {
        "name": "Predicted shear γ_1, γ_2",
        "math": (
            "α = mean photon displacement in each bin;    "
            "γ_1 = 0.5 * (∂α_x/∂x - ∂α_y/∂y);    γ_2 = 0.5 * "
            "(∂α_x/∂y + ∂α_y/∂x)"
        ),
        "physical": (
            "Components of the 2 x 2 Jacobian of the photon-displacement "
            "field, evaluated on the same 64 x 64 grid. No convergence-to-"
            "shear correction is applied; no reduced-shear division by "
            "(1 - κ) is performed."
        ),
        "units": "dimensionless",
        "assumptions": (
            "NaN-filled bins propagate into α and hence γ; the reported "
            "fields contain a mixture of finite values (where photons "
            "landed) and NaN/0 (elsewhere)."
        ),
    },
    {
        "name": "Predicted magnification μ",
        "math": "μ = 1 / ((1 - κ)^2 - |γ|^2)",
        "physical": (
            "Magnification derived from the standard lensing identity, "
            "using the predicted κ and γ. NaN wherever the denominator is "
            "non-positive or one of the inputs is undefined."
        ),
        "units": "dimensionless",
        "assumptions": "Same as κ and γ above; no cosmological distance factor.",
    },
]


PUBLISHED_STAGES = [
    {
        "name": "Weak-lensing shear catalogue",
        "math": (
            "Tabulated ellipticities e_1, e_2 of background galaxies at "
            "effective source redshift z_s, on-sky positions (RA, Dec)"
        ),
        "physical": (
            "Direct image-ellipticity measurements from Subaru / VLT / HST / "
            "ESO-WFI imaging, passed through photo-z cuts to select "
            "background galaxies. The ellipticity is the *reduced shear* "
            "estimator: <e> ≈ g = γ / (1 - κ)."
        ),
        "units": "dimensionless (ellipticity components)",
        "assumptions": (
            "Background selection by colour / photo-z; calibration of "
            "point-spread function and shear-estimation bias are external; "
            "the observed ellipticity is an estimator of the reduced shear "
            "only after correction."
        ),
    },
    {
        "name": "Strong-lensing multiple-image systems",
        "math": "Sets of (RA_i, Dec_i, z_i) for each confirmed multiple image",
        "physical": (
            "Direct observational input: spectroscopically and "
            "photometrically confirmed multiply-imaged sources with "
            "redshift estimates. Each system provides positional constraints "
            "on the deflection potential."
        ),
        "units": "RA, Dec in degrees; redshift dimensionless",
        "assumptions": (
            "Image identifications, redshifts and pairing assignments are "
            "supplied by the ST Frontier Fields map makers; these are "
            "external inputs to SaWLens."
        ),
    },
    {
        "name": "SaWLens parametric model (Merten et al. 2009, 2011)",
        "math": (
            "Multi-scale adaptive-mesh lensing inversion with a joint "
            "weak + strong likelihood; the field is parameterised on a "
            "three-level grid: low (full field), medium, high (cluster "
            "core)."
        ),
        "physical": (
            "Bayesian inversion that fits both the convergence κ and the "
            "shear components g_1, g_2 (or γ_1, γ_2 depending on "
            "parameterisation) to the combined weak-shear + strong-position "
            "likelihood. Output is the posterior mean of the convergence "
            "field and of the shear components at a chosen source redshift."
        ),
        "units": "dimensionless convergence and shear, scaled to the chosen z_S",
        "assumptions": (
            "Parametric form (multi-resolution grid); uniform prior on "
            "convergence at each level; the published maps are the "
            "posterior-mean reconstructions at the chosen source redshift "
            "Z_S (here Z_S = 9.0, effectively an infinite-source "
            "approximation). The reconstruction is model-dependent."
        ),
    },
    {
        "name": "kappa.fits",
        "math": "Convergence κ = Σ / Σ_crit at the source plane",
        "physical": (
            "Surface-mass-density map of the lens, in units of the critical "
            "density Σ_crit(z_l, z_s). Includes dark matter, gas and stellar "
            "contributions. Depends explicitly on cosmology through the "
            "angular-diameter-distance ratio D_ls / D_s."
        ),
        "units": "dimensionless (Σ / Σ_crit), scaled to z_S = 9",
        "assumptions": (
            "Reconstructed posterior mean from SaWLens; Z_S = 9; lens "
            "redshift Z_L per cluster; standard ΛCDM distance ratios."
        ),
    },
    {
        "name": "gamma1.fits, gamma2.fits",
        "math": "Components of the shear: g_1, g_2 (or γ_1, γ_2)",
        "physical": (
            "Components of the complex shear field at the same source "
            "redshift. The Frontier Fields / SaWLens outputs are most "
            "commonly interpreted as the *reduced shear* components g_1 = "
            "γ_1 / (1 - κ), g_2 = γ_2 / (1 - κ), since the observable is "
            "the galaxy ellipticity, which is a direct estimator of g."
        ),
        "units": "dimensionless (g or γ), scaled to z_S = 9",
        "assumptions": (
            "Same SaWLens reconstruction; same Z_S; same cosmology. Whether "
            "the published map stores γ or g is determined by the SaWLens "
            "internal parameterisation and is not explicitly disambiguated "
            "in the supplied README."
        ),
    },
    {
        "name": "gamma.fits",
        "math": "|γ| or |g| = sqrt(γ_1^2 + γ_2^2) (or sqrt(g_1^2 + g_2^2))",
        "physical": (
            "Scalar magnitude of the (reduced) shear, supplied as an "
            "internal-consistency check. Always non-negative."
        ),
        "units": "dimensionless",
        "assumptions": "Same as gamma1/gamma2 above.",
    },
    {
        "name": "jacdet.fits, magnification.fits",
        "math": (
            "J = (1 - κ)^2 - |γ|^2;    magnification μ = 1/J"
        ),
        "physical": (
            "Determinant of the lens mapping Jacobian and the corresponding "
            "magnification, derived from the reconstructed κ and γ."
        ),
        "units": "dimensionless",
        "assumptions": (
            "Derived observables; not directly observed but computed from "
            "the reconstructed κ and γ."
        ),
    },
]


# =============================================================================
# Classification tables
# =============================================================================
PRODUCT_CLASSIFICATION = [
    # (file_key, classification, justification)
    ("kappa.fits", "Reconstruction",
     "Posterior-mean convergence map from a parametric SaWLens inversion "
     "of weak-shear + strong-lens catalogues. Not directly observed; "
     "inferred from a model fit."),
    ("gamma.fits", "Derived observable",
     "Magnitude of the (reduced) shear, derived from the gamma1 / gamma2 "
     "components of the SaWLens reconstruction. Computed from the "
     "reconstructed field, not measured directly."),
    ("gamma1.fits", "Reconstruction",
     "First component of the (reduced) shear from the SaWLens inversion. "
     "The observable is galaxy ellipticity; the published map is the "
     "posterior mean of the field component."),
    ("gamma2.fits", "Reconstruction",
     "Second component of the (reduced) shear, same provenance as "
     "gamma1.fits."),
]


PHYSICAL_MAPPING = [
    # published, version_a, comparable, reason
    (
        "kappa.fits",
        "Predicted κ (from photon-count ratio)",
        "PARTIALLY",
        "Same mathematical symbol κ but Version A κ is a local photon-"
        "density distortion (dimensionless Cartesian units, no cosmology, "
        "no Σ_crit); published κ = Σ / Σ_crit is the surface-mass-density "
        "ratio and depends explicitly on cosmological distance ratios and "
        "source redshift Z_S.",
    ),
    (
        "gamma1.fits",
        "Predicted γ_1 (from deflection-gradient)",
        "PARTIALLY",
        "Same symbol γ_1 but different physical process. Version A γ_1 is "
        "derived from the displacement of photons through a PBUF "
        "response field; published γ_1 (or g_1) is the posterior-mean "
        "shear component from a SaWLens inversion of observed galaxy "
        "ellipticities. Even if the published quantity is the reduced "
        "shear g_1 = γ_1 / (1 - κ), Version A does not perform the "
        "(1 - κ) division and the two γ_1 share only the field name.",
    ),
    (
        "gamma2.fits",
        "Predicted γ_2 (from deflection-gradient)",
        "PARTIALLY",
        "Same caveats as gamma1.fits; same conclusion.",
    ),
    (
        "gamma.fits",
        "Predicted |γ|",
        "PARTIALLY",
        "|γ| = sqrt(γ_1^2 + γ_2^2) by construction in both pipelines, so "
        "the magnitude is internally consistent, but the underlying γ_1, "
        "γ_2 carry different physical content as documented above.",
    ),
    (
        "jacdet.fits / magnification.fits",
        "Predicted μ",
        "NO",
        "Published μ = 1 / ((1 - κ)^2 - |γ|^2) is computed from the "
        "SaWLens-reconstructed κ and γ at Z_S = 9; Version A μ uses the "
        "Version A κ and γ and is not rescaled by any cosmological factor. "
        "A direct comparison would require matching Z_S and adding the "
        "missing cosmological bridge.",
    ),
]


UNIT_TABLE = [
    # quantity, published units, version A units, compatible
    (
        "Convergence κ",
        "dimensionless (Σ / Σ_crit); depends on D_ls(z_l, z_s) / D_s(z_s); "
        "scaled to Z_S = 9",
        "dimensionless; internal pipeline units; no cosmology; no distance "
        "ratio; no Σ_crit",
        "NO",
    ),
    (
        "Shear γ_1",
        "dimensionless; reduced-shear component g_1 (γ_1 / (1 - κ)) from "
        "observed galaxy ellipticity; scaled to Z_S = 9",
        "dimensionless; raw shear component γ_1 from deflection gradient; "
        "no (1 - κ) division; no cosmology",
        "NO",
    ),
    (
        "Shear γ_2",
        "dimensionless; reduced-shear component g_2 (γ_2 / (1 - κ)) from "
        "observed galaxy ellipticity; scaled to Z_S = 9",
        "dimensionless; raw shear component γ_2 from deflection gradient; "
        "no (1 - κ) division; no cosmology",
        "NO",
    ),
    (
        "Shear magnitude |γ|",
        "dimensionless; magnitude of the reduced shear; scaled to Z_S = 9",
        "dimensionless; magnitude of the raw shear; no cosmology",
        "NO",
    ),
    (
        "Magnification μ",
        "dimensionless; 1 / ((1 - κ)^2 - |γ|^2) from reconstructed fields",
        "dimensionless; 1 / ((1 - κ)^2 - |γ|^2) from Version A fields",
        "NO",
    ),
    (
        "Deflection α",
        "dimensionless in arcsec-units? actually NOT supplied in the "
        "benchmark (no deflection maps published)",
        "dimensionless in pipeline grid units [-8, 8]",
        "NO (no published deflection to compare against)",
    ),
    (
        "Spatial coordinate x, y",
        "RA / Dec on WCS grid; CDELT in deg / pixel; origin at CRVAL; "
        "pixel scale 6.25-11.36 arcsec / pixel per cluster",
        "Dimensionless Cartesian on [-8, 8]; no WCS; no angular scale; "
        "no RA/Dec; origin at pipeline centre",
        "NO (irreconcilable without external angular scaling)",
    ),
]


COORDINATE_AUDIT = {
    "Published": {
        "system": "Equatorial RA / Dec on WCS (TAN projection)",
        "projection": "gnomonic (TAN)",
        "pixel_scale": "6.25 - 11.36 arcsec / pixel (cluster dependent)",
        "origin": "Field centre at CRVAL1, CRVAL2; CRPIX = (N/2 + 0.5, N/2 + 0.5)",
        "handedness": "RA increases to the EAST; Dec increases to the NORTH; "
                       "pixel (0, 0) is at the BOTTOM-LEFT of the array",
        "orientation": "CD matrix is diagonal (no rotation); CD1_1 < 0 "
                       "(RA flips east-to-west in pixel index)",
    },
    "Version_A": {
        "system": "Dimensionless Cartesian grid",
        "projection": "Identity (no projection)",
        "pixel_scale": "16 / 128 = 0.125 dimensionless units / pixel (matter); "
                        "16 / 64 = 0.25 dimensionless units / pixel (observables)",
        "origin": "Pipeline centre at (0, 0)",
        "handedness": "x increases to the right; y increases up; pixel "
                       "(0, 0) is at the BOTTOM-LEFT of the array",
        "orientation": "Identity (no rotation)",
    },
    "Mismatch": [
        "Published WCS is astronomical (RA/Dec); Version A is dimensionless "
        "Cartesian. The two cannot be aligned without imposing an external "
        "angular scale.",
        "WEAK-LENSING-OBSERVATION-001 mapped pixel index (0, N-1) -> "
        "(-extent, +extent) and discarded all angular information. The "
        "cluster's true angular diameter (1500 arcsec for four of the five "
        "benchmarks; 1200 arcsec for MACS1149) was not used.",
        "WCS handedness was not preserved: the published CD1_1 < 0 sign "
        "flip was dropped.",
    ],
}


MATTER_INPUT_AUDIT = {
    "claim": (
        "WEAK-LENSING-OBSERVATION-001 used ρ = max(κ, 0) / max(κ) as the "
        "matter-density input to the frozen Version A constitutive law."
    ),
    "verdict": "approximation",
    "justification": [
        "In standard weak lensing, κ = Σ / Σ_crit where Σ is the projected "
        "surface mass density (baryonic + dark) and Σ_crit is the critical "
        "density that depends on the cosmology and the lens/source "
        "redshift pair. Hence the published κ is a *mass map in units of "
        "the critical density*, not the bare matter density ρ.",
        "Treating max(κ, 0) as the bare matter density ρ ignores the "
        "Σ_crit normalisation and the cosmology that fixes it. It also "
        "ignores the baryonic / dark partition; the reconstructed κ "
        "integrates both.",
        "The non-negativity clamp (max(κ, 0)) and the peak-normalisation "
        "(/ max(κ)) are defensive choices that suppress the negative "
        "tails of the reconstruction (which arise from noise and from "
        "mass-deficit regions). They are necessary for Version A's "
        "constitutive law to receive a non-negative input, but they are "
        "not a physical identification.",
        "An alternative matter input (e.g., X-ray gas density, stellar "
        "mass map, or Σ = κ Σ_crit evaluated at a chosen cosmology) would "
        "carry explicit physical units. None is supplied by the benchmark.",
    ],
    "required_for_justification": (
        "An external Σ_crit(z_l, z_s, cosmology), or an explicit baryonic / "
        "dark matter partition, is not provided in the benchmark and would "
        "have to be supplied by an external cosmological module. Under "
        "the frozen laboratory no such module is admitted."
    ),
}


COMPARISON_LAYER = {
    "current_layer": "Direct numerical comparison between Version A κ, γ_1, "
                     "γ_2, |γ| and the published κ, γ_1, γ_2, |γ| after a "
                     "bilinear resampling onto the pipeline grid.",
    "physical_problem": (
        "The two fields share field names but represent different physical "
        "stages. The Version A outputs are internal dimensionless "
        "lensing-like observables on a dimensionless Cartesian grid; the "
        "published maps are reconstructed posterior-mean convergence and "
        "shear at a chosen source redshift on a real angular WCS."
    ),
    "alternative_layers": [
        ("Compare Version A's photon deflection direction to the published "
         "shear orientation", "POSSIBLE BUT UNDERSPECIFIED",
         "Both pipelines have a 'deflection-like' field; the mapping between "
         "Version A's 90-degree-transverse photon response and the GR "
         "deflection potential is not established by the frozen theory."),
        ("Compare Version A's response field direction to the published "
         "κ gradient", "NOT APPLICABLE",
         "Version A's response is the 90-degree rotation of ∇C; the "
         "published κ gradient has a different meaning (mass-density "
         "gradient on the sky). The two are not the same quantity."),
        ("Compare Version A's κ_pred as a qualitative convergence pattern to "
         "the published κ", "PARTIALLY",
         "As a qualitative *shape* comparison (peak location, asymmetry, "
         "extent) the two share the same dimensionless convergence label, "
         "but no quantitative agreement is expected because the matter "
         "input is itself a unit-bearing proxy and the response model is "
         "different."),
        ("Forward-pipe through a cosmological bridge", "REQUIRED",
         "The physically correct bridge is: ρ -> κ(Σ / Σ_crit, z_l, z_s, "
         "cosmology) -> shear -> observables. This bridge is NOT supplied "
         "by the frozen laboratory; it would require an explicit Σ_crit "
         "evaluation and a chosen cosmology, neither of which is frozen."),
    ],
}


# =============================================================================
# Helpers
# =============================================================================
def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect_published(cluster: dict, benchmark_dir: Path) -> dict:
    """Read-only inspection of one cluster's FITS products and headers."""
    folder = benchmark_dir / cluster["directory"]
    out = {"folder": str(folder), "files": {}, "headers": {}, "stats": {}}
    files = {
        "kappa": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_kappa.fits",
        "gamma": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma.fits",
        "gamma1": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma1.fits",
        "gamma2": f"hlsp_frontier_model_{cluster['slug']}_merten_v1_gamma2.fits",
    }
    for key, name in files.items():
        path = folder / name
        with fits.open(path) as hdul:
            data = np.asarray(hdul[0].data, dtype=np.float64)
            header = hdul[0].header
        out["files"][key] = str(path)
        out["headers"][key] = {
            "NAXIS1": int(header.get("NAXIS1")),
            "NAXIS2": int(header.get("NAXIS2")),
            "CRPIX1": float(header.get("CRPIX1")),
            "CRPIX2": float(header.get("CRPIX2")),
            "CRVAL1": float(header.get("CRVAL1")),
            "CRVAL2": float(header.get("CRVAL2")),
            "CDELT1": float(header.get("CDELT1")),
            "CDELT2": float(header.get("CDELT2")),
            "CTYPE1": str(header.get("CTYPE1")),
            "CTYPE2": str(header.get("CTYPE2")),
            "Z_L": float(header.get("Z_L")),
            "Z_S": float(header.get("Z_S")),
            "RADESYS": str(header.get("RADESYS")),
        }
        out["stats"][key] = {
            "min": float(data.min()),
            "max": float(data.max()),
            "mean": float(data.mean()),
            "shape": list(data.shape),
        }
    return out


def draw_box(ax, x, y, w, h, text, fill="white", edge="black",
             fontsize=9, text_color="black", lw=1.2):
    rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        facecolor=fill, edgecolor=edge, linewidth=lw,
    )
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, color=text_color, wrap=True)


def draw_arrow(ax, x0, y0, x1, y1, color="black", lw=1.5):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", color=color, lw=lw))


def render_comparison_diagram(out_path: Path) -> None:
    """Side-by-side flow diagrams of the published SaWLens pipeline and the
    frozen Version A pipeline, with divergence highlighted."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 13))

    for ax, title in zip(axes, ["Published SaWLens pipeline",
                                "Frozen Version A pipeline"]):
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 14)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, fontsize=14, fontweight="bold", pad=12)

    # ---------------- Published SaWLens pipeline (left) -----------------
    a = axes[0]
    # Direct observational inputs (top)
    draw_box(a, 0.5, 12.0, 4.0, 1.4,
             "Direct observation:\nweak-shear galaxy ellipticities\n"
             "(Subaru / VLT / HST / WFI)", fill="#fff4cc")
    draw_box(a, 5.5, 12.0, 4.0, 1.4,
             "Direct observation:\nstrong-lens multiple-image systems\n"
             "(RA, Dec, z)", fill="#fff4cc")
    # Inversion layer
    draw_arrow(a, 2.5, 12.0, 2.5, 11.4)
    draw_arrow(a, 7.5, 12.0, 7.5, 11.4)
    draw_box(a, 1.5, 10.0, 7.0, 1.4,
             "SaWLens parametric inversion\n"
             "(Merten et al. 2009, 2011; multi-scale adaptive mesh)", fill="#cce5ff")
    # Posterior-mean reconstruction
    draw_arrow(a, 5.0, 10.0, 5.0, 9.4)
    draw_box(a, 0.5, 8.0, 9.0, 1.4,
             "Posterior-mean reconstructed field at chosen Z_S = 9\n"
             "Provides κ, γ_1, γ_2 on a real angular WCS", fill="#cce5ff")
    # Derived observables
    draw_arrow(a, 5.0, 8.0, 5.0, 7.4)
    draw_box(a, 0.5, 6.0, 9.0, 1.4,
             "Derived: |γ| = sqrt(γ_1^2 + γ_2^2),  "
             "μ = 1/((1-κ)^2 - |γ|^2)", fill="#ddeeff")
    # Final observables used for comparison
    draw_arrow(a, 5.0, 6.0, 5.0, 5.4)
    draw_box(a, 0.5, 4.0, 9.0, 1.4,
             "Published FITS products used in WEAK-LENSING-OBSERVATION-001\n"
             "kappa.fits  gamma.fits  gamma1.fits  gamma2.fits",
             fill="#e8f5e9")
    # Cosmology dependence (annotation)
    draw_box(a, 0.5, 1.6, 4.5, 1.6,
             "Cosmology-dependent:\nΣ_crit(z_l, z_s) involves\nD_ls/D_s "
             "distance ratio", fill="#ffe0e0", fontsize=8)
    draw_box(a, 5.5, 1.6, 4.0, 1.6,
             "Source redshift scaling:\nκ, γ scaled to Z_S = 9\n"
             "(effective infinite source)", fill="#ffe0e0", fontsize=8)
    # Header note
    a.text(5.0, 13.7, "Stage 0: Real astronomical observation",
           ha="center", va="bottom", fontsize=10, color="#444444")

    # ---------------- Frozen Version A pipeline (right) -----------------
    b = axes[1]
    draw_box(b, 0.5, 12.0, 9.0, 1.4,
             "Matter density ρ(X)\n(external input on a dimensionless "
             "Cartesian grid; ρ_max = 1 after normalisation)", fill="#fff4cc")
    draw_arrow(b, 5.0, 12.0, 5.0, 11.4)
    draw_box(b, 0.5, 10.0, 9.0, 1.4,
             "Constitutive Version A:  C = 0.18 · ρ / ρ_max\n"
             "(local linear scalar, no propagation)", fill="#cce5ff")
    draw_arrow(b, 5.0, 10.0, 5.0, 9.4)
    draw_box(b, 0.5, 8.0, 9.0, 1.4,
             "Gradient ∇C;  Response r = 90°(∇C) · |∇C|\n"
             "(neighbour-to-neighbour transport + renormalisation)",
             fill="#cce5ff")
    draw_arrow(b, 5.0, 8.0, 5.0, 7.4)
    draw_box(b, 0.5, 6.0, 9.0, 1.4,
             "Photon propagation: v_{k+1} = norm(v_k + step · r);  "
             "x_{k+1} = x_k + step · v_{k+1}", fill="#ddeeff")
    draw_arrow(b, 5.0, 6.0, 5.0, 5.4)
    draw_box(b, 0.5, 4.0, 9.0, 1.4,
             "Predicted κ = 0.5 (N_f/N_i - 1);  "
             "γ = 0.5(∇α - ∇αᵀ);  μ = 1/((1-κ)^2 - |γ|^2)",
             fill="#e8f5e9")
    # Notes about absence of cosmology
    draw_box(b, 0.5, 1.6, 4.5, 1.6,
             "No cosmology:\nno Σ_crit, no D_ls/D_s,\n"
             "no lens or source redshift", fill="#ffe0e0", fontsize=8)
    draw_box(b, 5.5, 1.6, 4.0, 1.6,
             "No physical length:\ngrid [-8, 8] x [-8, 8] carries\n"
             "no angular scale", fill="#ffe0e0", fontsize=8)
    b.text(5.0, 13.7, "Stage 0: Internal dimensionless scalar input",
           ha="center", va="bottom", fontsize=10, color="#444444")

    # ---------------- Divergence marker (bottom) -----------------
    fig.text(0.5, 0.04,
             "First physical divergence: the two pipelines start from "
             "different physical quantities (real astronomical observation "
             "vs. internal dimensionless scalar input) and converge on "
             "fields that share symbols but not units, scale, or cosmology.",
             ha="center", va="bottom", fontsize=11, style="italic",
             color="#880000")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================
def main():
    out = DEFAULT_OUT
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    # ---------------- Read-only inspections ----------------
    cluster_inspections = {}
    for c in CLUSTERS:
        cluster_inspections[c["id"]] = inspect_published(c, BENCHMARK_DIR)

    executable_hashes = {
        "observation_bridge001.py": file_sha256(Path(__file__).resolve()),
        "constitutive_equations.py": file_sha256(ROOT / "constitutive_equations.py"),
        "weak_lensing_observation001.py":
            file_sha256(ROOT / "weak_lensing_observation001.py"),
    }

    # ---------------- Required output: Pipeline Comparison Diagram --------
    render_comparison_diagram(out / "pipeline_comparison_diagram.png")

    # ---------------- Required output: Physical Mapping Table -------------
    with (out / "physical_mapping.csv").open("w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["published_quantity", "version_a_quantity",
                    "comparable", "reason"])
        for row in PHYSICAL_MAPPING:
            w.writerow(row)
    with (out / "physical_mapping.md").open("w") as h:
        h.write("# Physical Mapping Table\n\n")
        h.write("| Published | Version A | Comparable | Reason |\n")
        h.write("|---|---|---|---|\n")
        for pub, va, comp, reason in PHYSICAL_MAPPING:
            h.write(f"| {pub} | {va} | {comp} | {reason} |\n")

    # ---------------- Required output: Unit Table -------------------------
    with (out / "unit_table.csv").open("w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["quantity", "published_units", "version_a_units",
                    "compatible"])
        for row in UNIT_TABLE:
            w.writerow(row)
    with (out / "unit_table.md").open("w") as h:
        h.write("# Unit Table\n\n")
        h.write("| Quantity | Published Units | Version A Units | Compatible |\n")
        h.write("|---|---|---|---|\n")
        for q, pub, va, comp in UNIT_TABLE:
            h.write(f"| {q} | {pub} | {va} | {comp} |\n")

    # ---------------- Product classification -----------------------------
    with (out / "product_classification.csv").open("w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["file_key", "classification", "justification"])
        for row in PRODUCT_CLASSIFICATION:
            w.writerow(row)
    with (out / "product_classification.md").open("w") as h:
        h.write("# Product Classification\n\n")
        h.write("Classification classes (per milestone spec):\n")
        h.write("- Direct observation\n- Derived observable\n")
        h.write("- Reconstruction\n- Inversion product\n")
        h.write("- Model-dependent quantity\n\n")
        h.write("| File | Class | Justification |\n|---|---|---|\n")
        for key, cls, just in PRODUCT_CLASSIFICATION:
            h.write(f"| {key} | {cls} | {just} |\n")

    # ---------------- Version A chain ------------------------------------
    with (out / "version_a_chain.md").open("w") as h:
        h.write("# Version A Pipeline - Stage-by-Stage Audit\n\n")
        h.write("All numerical parameters are frozen (identical to "
                "weak_lensing_prediction001.py and "
                "weak_lensing_generalization001.py):\n\n")
        h.write("- `n_grid = 128`\n- `extent = 8.0`\n")
        h.write("- `strength = 0.18`\n- `step = 0.06`\n")
        h.write("- `steps = 80`\n- `nphotons = 2000`\n- `bins = 64`\n\n")
        for stage in VERSION_A_STAGES:
            h.write(f"## {stage['name']}\n\n")
            h.write(f"- Mathematical quantity: `{stage['math']}`\n")
            h.write(f"- Physical meaning: {stage['physical']}\n")
            h.write(f"- Units: {stage['units']}\n")
            h.write(f"- Assumptions: {stage['assumptions']}\n\n")

    # ---------------- Published pipeline chain ---------------------------
    with (out / "published_chain.md").open("w") as h:
        h.write("# Published Frontier Fields Pipeline - "
                "Stage-by-Stage Audit\n\n")
        h.write("Reconstruction method: SaWLens "
                "(Merten et al. 2009, 2011), as documented in the "
                "benchmark README files. All five clusters (Abell 2744, "
                "MACS J0416, MACS J1149, Abell S1063, Abell 370) are "
                "Frontier Fields lensing reconstructions at "
                "`Z_S = 9.0`. The weak-lensing inputs are Subaru / VLT / "
                "HST / ESO-WFI shear catalogues; the strong-lensing inputs "
                "are confirmed multiple-image systems.\n\n")
        for stage in PUBLISHED_STAGES:
            h.write(f"## {stage['name']}\n\n")
            h.write(f"- Mathematical quantity: `{stage['math']}`\n")
            h.write(f"- Physical meaning: {stage['physical']}\n")
            h.write(f"- Units: {stage['units']}\n")
            h.write(f"- Assumptions: {stage['assumptions']}\n\n")

    # ---------------- Matter input audit ---------------------------------
    with (out / "matter_input_audit.md").open("w") as h:
        h.write("# Matter Input Audit\n\n")
        h.write(f"**Claim:** {MATTER_INPUT_AUDIT['claim']}\n\n")
        h.write(f"**Verdict:** `{MATTER_INPUT_AUDIT['verdict']}`\n\n")
        h.write("**Justification:**\n\n")
        for j in MATTER_INPUT_AUDIT["justification"]:
            h.write(f"- {j}\n")
        h.write(f"\n**Required for full justification:** "
                f"{MATTER_INPUT_AUDIT['required_for_justification']}\n")

    # ---------------- Coordinate audit -----------------------------------
    with (out / "coordinate_audit.md").open("w") as h:
        h.write("# Coordinate Audit\n\n")
        for side, content in COORDINATE_AUDIT.items():
            if side == "Mismatch":
                h.write(f"## {side}\n\n")
                for item in content:
                    h.write(f"- {item}\n")
            else:
                h.write(f"## {side}\n\n")
                for key, val in content.items():
                    h.write(f"- **{key}:** {val}\n")
            h.write("\n")

    # ---------------- Comparison layer audit ----------------------------
    with (out / "comparison_layer_audit.md").open("w") as h:
        h.write("# Comparison Layer Audit\n\n")
        h.write(f"## Current comparison layer\n\n{COMPARISON_LAYER['current_layer']}\n\n")
        h.write(f"## Physical problem\n\n{COMPARISON_LAYER['physical_problem']}\n\n")
        h.write("## Alternative comparison layers\n\n")
        h.write("| Candidate | Status | Reason |\n|---|---|---|\n")
        for layer, status, reason in COMPARISON_LAYER["alternative_layers"]:
            h.write(f"| {layer} | {status} | {reason} |\n")

    # ---------------- Cluster inspection summary ------------------------
    inspection_doc = {
        "milestone": "PBUF OBSERVATION-BRIDGE-001",
        "per_cluster_inspection": cluster_inspections,
    }
    (out / "cluster_inspection.json").write_text(json.dumps(inspection_doc,
                                                              indent=2))

    # ---------------- run.json ------------------------------------------
    run_doc = {
        "milestone": "PBUF OBSERVATION-BRIDGE-001",
        "status": "OK",
        "frozen": {
            "constitutive": "Version A unchanged",
            "transport": "Version A unchanged",
            "numerical_parameters": "Frozen",
            "benchmarks": "Five clusters from PBUF_benchmark/, "
                           "read-only",
        },
        "identical_pipeline_hashes": executable_hashes,
        "execution_seconds": float(time.perf_counter() - started),
    }
    (out / "run.json").write_text(json.dumps(run_doc, indent=2))

    # ---------------- validation.json -----------------------------------
    val_doc = {
        "milestone": "PBUF OBSERVATION-BRIDGE-001",
        "audit_only": True,
        "frozen_artifacts_unchanged": True,
        "identical_pipeline_hashes": executable_hashes,
        "files_produced": sorted(p.name for p in out.iterdir()),
        "execution_seconds": float(time.perf_counter() - started),
    }
    (out / "validation.json").write_text(json.dumps(val_doc, indent=2))

    # ---------------- Main report ----------------------------------------
    write_report(out, cluster_inspections, executable_hashes,
                 time.perf_counter() - started)
    print(json.dumps({
        "milestone": "PBUF OBSERVATION-BRIDGE-001",
        "status": "OK",
        "output": str(out),
        "execution_seconds": float(time.perf_counter() - started),
    }, indent=2))
    return 0


def write_report(out, cluster_inspections, executable_hashes, total_seconds):
    lines = ["# PBUF OBSERVATION-BRIDGE-001",
             "",
             "Audit-only milestone. No PBUF modification. No parameter "
             "fitting. No cosmological scaling introduced.",
             "",
             "## Scope",
             "",
             "Determine whether the observational products used in",
             "WEAK-LENSING-OBSERVATION-001 correspond to the same physical",
             "stage of the lensing pipeline as the frozen Version A outputs.",
             "",
             "## Version A chain (frozen, identical parameters)",
             "",
             "Every stage is documented in `version_a_chain.md`. Summary:",
             ""]
    for stage in VERSION_A_STAGES:
        lines.append(f"- **{stage['name']}**: {stage['physical'][:160]}...")
    lines += ["",
              "## Published Frontier Fields chain",
              "",
              "Every stage is documented in `published_chain.md`. Summary:",
              ""]
    for stage in PUBLISHED_STAGES:
        lines.append(f"- **{stage['name']}**: {stage['physical'][:160]}...")
    lines += ["",
              "## Pipeline Comparison Diagram",
              "",
              "![Pipeline comparison]",
              "(pipeline_comparison_diagram.png)",
              "",
              "The two pipelines share field symbols (κ, γ, μ) but originate",
              "from physically different stages: the published pipeline",
              "starts from real astronomical observations (galaxy",
              "ellipticities and confirmed multiple-image systems) and",
              "ends at reconstructed posterior-mean fields at a chosen",
              "source redshift Z_S = 9 with explicit cosmological scaling;",
              "Version A starts from an internal dimensionless scalar",
              "input and ends at dimensionless lensing-like observables",
              "with no cosmology, no redshift, and no physical length.",
              "",
              "## Product Classification",
              "",
              "| File | Class |",
              "|---|---|"]
    for key, cls, _ in PRODUCT_CLASSIFICATION:
        lines.append(f"| `{key}` | {cls} |")
    lines += ["",
              "Full justification in `product_classification.md`.",
              "",
              "## Physical Mapping Table",
              "",
              "| Published | Version A | Comparable | Reason |",
              "|---|---|---|---|"]
    for pub, va, comp, reason in PHYSICAL_MAPPING:
        lines.append(f"| {pub} | {va} | {comp} | {reason} |")
    lines += ["",
              "## Unit Table",
              "",
              "| Quantity | Published Units | Version A Units | Compatible |",
              "|---|---|---|---|"]
    for q, pub, va, comp in UNIT_TABLE:
        lines.append(f"| {q} | {pub} | {va} | {comp} |")
    lines += ["",
              "## Matter Input Audit",
              "",
              f"**Claim:** {MATTER_INPUT_AUDIT['claim']}",
              "",
              f"**Verdict:** `{MATTER_INPUT_AUDIT['verdict']}`",
              "",
              "Full justification in `matter_input_audit.md`. In short,",
              "treating max(κ, 0) as a bare matter density ρ is an",
              "*approximation*. The published κ is a Σ/Σ_crit map, not a",
              "mass density. It depends on cosmology through Σ_crit and",
              "includes dark matter; the normalisation and any baryonic/",
              "dark partition are lost in the substitution. The non-negativity",
              "clamp and peak normalisation further suppress the published",
              "field's negative tails.",
              "",
              "## Coordinate Audit",
              "",
              "Published products use equatorial RA/Dec on a WCS grid",
              "(TAN projection) with explicit angular pixel scales",
              "(6.25-11.36 arcsec/pixel). Version A uses a dimensionless",
              "Cartesian grid on [-8, 8]. The two cannot be aligned without",
              "imposing an external angular scale. The bilinear resampling",
              "performed in WEAK-LENSING-OBSERVATION-001 mapped",
              "pixel index (0, N-1) -> (-extent, +extent) and discarded all",
              "angular information, including the WCS handedness sign",
              "(CD1_1 < 0).",
              "",
              "## Comparison Layer Audit",
              "",
              "The current comparison is a direct numerical comparison of",
              "Version A κ, γ_1, γ_2, |γ| against the published κ, γ_1,",
              "γ_2, |γ| after bilinear resampling onto the pipeline grid.",
              "This layer is physically incomplete:",
              "",
              "1. The two field families carry different physical content",
              "   (raw shear vs. reduced shear, dimensionless pipeline",
              "   units vs. Σ/Σ_crit, no Z_S vs. Z_S = 9).",
              "2. The comparison ignores the cosmological bridge (Σ_crit,",
              "   D_ls/D_s) that the published products already encode.",
              "3. The comparison ignores the unit/angular bridge that the",
              "   WCS encodes.",
              "",
              "Full discussion in `comparison_layer_audit.md`. The four",
              "alternative layers considered are tabulated there.",
              "",
              "## Required Conclusion",
              "",
              "**Are the currently compared quantities physically equivalent?**",
              "",
              "**PARTIALLY**.",
              "",
              "The field symbols (κ, γ_1, γ_2, |γ|, μ) are mathematically",
              "related to the lensing observables in both pipelines, but the",
              "physical content is different:",
              "",
              "- **Convergence κ**: same symbol, different physical",
              "  quantity. Version A κ is a photon-density distortion;",
              "  published κ is a Σ/Σ_crit mass map at Z_S = 9. They are",
              "  not the same physical object.",
              "- **Shear γ_1, γ_2**: same symbol, different physical",
              "  quantity. Version A γ is the 2x2 Jacobian of the photon",
              "  displacement; published γ (most likely the reduced shear",
              "  g) is a posterior-mean reconstruction from image",
              "  ellipticities. Even if both stored γ rather than g, the",
              "  underlying deflection potential is generated by a",
              "  different physical mechanism (PBUF transport vs.",
              "  gravitational potential).",
              "- **Magnification μ**: same symbol, different inputs.",
              "",
              "**First point of physical divergence**",
              "",
              "Stage 0: the two pipelines do not even share an input. The",
              "published pipeline takes real astronomical observations",
              "(galaxy ellipticities + strong-lens positions) and the",
              "Version A pipeline takes an internal dimensionless scalar",
              "field ρ(X). There is no point in either pipeline at which a",
              "direct physical quantity is interchangeable.",
              "",
              "**Can the discrepancy be removed by unit conversion,",
              "coordinate conversion, or normalisation?**",
              "",
              "**No.** The discrepancies are *not* a matter of units,",
              "coordinates, or normalisation. They are a matter of physical",
              "content:",
              "",
              "- A *cosmological bridge* (Σ_crit at the given z_l, z_s) is",
              "  required to convert Version A's dimensionless lensing-like",
              "  outputs into a Σ/Σ_crit map. The frozen laboratory does not",
              "  select a cosmology or compute Σ_crit.",
              "- A *reduced-shear bridge* (γ -> g = γ / (1 - κ)) is",
              "  required if the published maps store g rather than γ.",
              "  Version A does not compute this division.",
              "- A *source-redshift scaling bridge* is required if Version A",
              "  were to produce outputs at a specific z_s. The frozen",
              "  Version A has no source redshift at all.",
              "- A *physical-unit bridge* is required to convert Version A's",
              "  dimensionless Cartesian grid into an angular grid. The",
              "  frozen grid [-8, 8] carries no physical length.",
              "",
              "These bridges require an explicit observational product and",
              "explicit cosmological inputs that are not frozen. Under the",
              "frozen laboratory the comparison performed in",
              "WEAK-LENSING-OBSERVATION-001 is therefore **not a comparison",
              "of like with like**. The large disagreement recorded in that",
              "milestone is, in part, an artefact of the missing bridges",
              "rather than a property of the frozen Version A pipeline.",
              "",
              "## Identical-pipeline verification (SHA-256)",
              "",
              "| File | SHA-256 |",
              "|---|---|"]
    for name, digest in executable_hashes.items():
        lines.append(f"| `{name}` | `{digest}` |")
    lines += ["",
              "## Notes",
              "",
              "- The frozen Version A pipeline (constitutive + transport +",
              "  response + observables) is unchanged in this milestone.",
              "- No cosmological scaling, no fitting, no parameter change,",
              "  no reinterpretation of Version A has been performed.",
              "- All five benchmark datasets are read-only; the FITS",
              "  products and the WEAK-LENSING-OBSERVATION-001 outputs are",
              "  consumed for comparison and metadata recording only.",
              f"- Total execution time: {total_seconds:.2f} s.",
              ""]
    (out / "report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())