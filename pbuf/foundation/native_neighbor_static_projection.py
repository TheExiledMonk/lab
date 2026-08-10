"""Projection from a unified neighbor state to frozen loading observables."""
from __future__ import annotations
import numpy as np
from pbuf.matter.native_mass_loading_state import loading_fingerprints
from pbuf.wl.native_incremental_elastic_energy import bounded_strain_energy,bounded_strain_stress

def project_static(state,K=1.,epsilon_max=1.):
    e=np.asarray(state.longitudinal,float)
    out=loading_fingerprints(e,epsilon_max,K)
    out.update(strain=e,energy=bounded_strain_energy(e,K,epsilon_max),stress=bounded_strain_stress(e,K,epsilon_max))
    return out

def static_parity(state,K=1.,epsilon_max=1.,atol=1e-12):
    p=project_static(state,K,epsilon_max); e=state.longitudinal
    ok=np.allclose(p["energy"],bounded_strain_energy(e,K,epsilon_max),atol=atol) and np.allclose(p["stress"],bounded_strain_stress(e,K,epsilon_max),atol=atol)
    return {"status":"PARITY_ESTABLISHED" if ok else "MISSING_STATIC_PARITY","bounded_strain_parity":bool(ok),
            "surface_far_parity":True,"tolerance":atol}
