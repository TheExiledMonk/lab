"""Dev150 transition registries and conservative discreteness tests."""
from __future__ import annotations
import numpy as np

QL_NAMES=("localized-state quantization","emission transition quantization","absorption transition quantization",
 "source boundary quantization","interaction threshold quantization","topological localized-state quantization",
 "nonlinear localized-state quantization","loading-bound mode quantization","combined loading/excitation quantization",
 "no quantization in current localized physics")
QX_NAMES=("discrete stable localized norms","discrete internal wavelengths","discrete node counts","discrete winding states",
 "discrete loading/excitation composites","transition-only discreteness","absorption thresholds","emission thresholds",
 "state-difference packet norm","state-difference emitted k","mode-selective coupling","polarization-selective transitions",
 "handedness-selective transitions","discrete source burst","subthreshold accumulation","topological transition",
 "nonlinear attractor transitions","continuous localized states but discrete interaction","fully continuous interaction",
 "current physics insufficient")

def quantization_location_registry():
    return [{"id":f"QL{i:02d}","name":n,"attempted":True,
             "status":"ESTABLISHED" if i==10 else "MISSING_INTERACTION_LAW"} for i,n in enumerate(QL_NAMES,1)]
def quantization_family_registry():
    return [{"id":f"QX{i:02d}","name":n,"attempted":True,
             "status":"ESTABLISHED" if i==20 else "MISSING_LOCALIZED_STATE"} for i,n in enumerate(QX_NAMES,1)]

def cluster_final_states(initial_norms, final_norms, *, tolerance=1e-6):
    a=np.asarray(initial_norms,float); b=np.asarray(final_norms,float)
    if a.shape!=b.shape or a.ndim!=1 or a.size<2 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("matching finite one-dimensional sweeps required")
    order=np.sort(b); clusters=[]
    for value in order:
        if not clusters or abs(value-np.mean(clusters[-1])) > tolerance: clusters.append([float(value)])
        else: clusters[-1].append(float(value))
    centers=[float(np.mean(c)) for c in clusters]
    continuous=bool(np.allclose(b,a,rtol=tolerance,atol=tolerance))
    return {"cluster_centers":centers,"cluster_count":len(centers),"input_count":len(a),
            "classification":"CONTINUOUS_FINAL_STATE_FAMILY" if continuous else "DISCRETE_FINAL_STATE_FAMILY" if len(centers)<len(a) else "MIXED",
            "discrete_attractor":bool(len(centers)<len(a) and not continuous)}

def artifact_controls(grid_results, domain_results):
    grid=np.asarray(grid_results,float); domain=np.asarray(domain_results,float)
    converged=bool(grid.size>1 and domain.size>1 and np.allclose(grid,grid[-1],rtol=.02) and np.allclose(domain,domain[-1],rtol=.02))
    return {"grid_quantization_rejected":converged,"box_quantization_rejected":converged,
            "classification":"PHYSICAL_CANDIDATE" if converged else "NUMERICAL_QUANTIZATION_ONLY"}

def transition_graph(states):
    return {"states":list(states),"edges":[],"status":"MISSING_INTERACTION_LAW"}
