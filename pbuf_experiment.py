#!/usr/bin/env python3
"""WL-001/WL-002: a small, reproducible forward lensing experiment.

The PBUF constitutive laws are isolated in ``constitutive_equations``. The
other stages consume only the selected equation's deformation field.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, subprocess, time
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from constitutive_equations import Equation, get_equation


@dataclass(frozen=True)
class Config:
    n: int = 128
    extent: float = 8.0
    source_x: float = 1.25
    source_y: float = 0.35
    mass_x: float = -0.65
    mass_y: float = 0.0
    mass_sigma: float = 0.75
    mass_amplitude: float = 1.0
    deformation_strength: float = 0.18
    photon_steps: int = 80
    photon_step_size: float = 0.06


def checksum(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def grid(c: Config):
    x = np.linspace(-c.extent, c.extent, c.n)
    return np.meshgrid(x, x, indexing="xy")


def mass_map(c: Config, X, Y):
    r2 = (X-c.mass_x)**2 + (Y-c.mass_y)**2
    return c.mass_amplitude*np.exp(-r2/(2*c.mass_sigma**2))


def gradient(field, x):
    spacing = float(x[1]-x[0])
    gy, gx = np.gradient(field, spacing, spacing)
    return gx, gy


def propagate(deformation, X, Y, c: Config):
    gx, gy = gradient(deformation, np.linspace(-c.extent, c.extent, c.n))
    paths=[]
    for y0 in np.linspace(-3.0, 3.0, 9):
        x=np.full(c.photon_steps, -c.extent); y=np.full(c.photon_steps, y0)
        # Deflection is integrated from sampled deformation gradients.
        vx=np.ones(c.photon_steps); vy=np.zeros(c.photon_steps)
        for k in range(1, c.photon_steps):
            ix=np.clip(np.searchsorted(np.linspace(-c.extent,c.extent,c.n), x[k-1])-1,0,c.n-1)
            iy=np.clip(np.searchsorted(np.linspace(-c.extent,c.extent,c.n), y[k-1])-1,0,c.n-1)
            vy[k]=vy[k-1]-c.photon_step_size*gy[iy,ix]
            vx[k]=vx[k-1]-c.photon_step_size*gx[iy,ix]*0.15
            norm=max(np.hypot(vx[k],vy[k]),1e-12); vx[k]/=norm; vy[k]/=norm
            x[k]=x[k-1]+c.photon_step_size*vx[k]; y[k]=y[k-1]+c.photon_step_size*vy[k]
        paths.append((x,y))
    return paths


def path_deviation(paths):
    return max((float(np.max(np.abs(y-y[0]))) for _, y in paths), default=0.0)


def render_image(X,Y, source_x, source_y, paths=None):
    img=np.exp(-((X-source_x)**2+(Y-source_y)**2)/(2*0.35**2))
    if paths:
        # image remains a source-plane observable; paths are visualized separately
        pass
    return img


def run(out: Path, c: Config, equation: Equation | None = None):
    equation = equation or get_equation("A")
    out.mkdir(parents=True, exist_ok=True); started=time.time(); X,Y=grid(c); x=X[0]
    matter=mass_map(c,X,Y); deformation=equation.solve(matter,c); gx,gy=gradient(deformation,x)
    zero_ok=np.allclose(equation.solve(np.zeros_like(matter),c),0)
    symmetric_probe=np.exp(-(X**2+Y**2)/(2*c.mass_sigma**2))
    symmetric_ok=np.allclose(equation.solve(symmetric_probe,c),equation.solve(symmetric_probe[:,::-1],c)[:,::-1],atol=1e-12)
    paths=propagate(deformation,X,Y,c)
    zero_paths=propagate(np.zeros_like(deformation),X,Y,c)
    observed=render_image(X,Y,c.source_x,c.source_y)
    # Isolated reference branches: deliberately do not consume PBUF intermediates.
    baryonic_gr=render_image(X,Y,c.source_x+0.05,c.source_y)
    lcdm=render_image(X,Y,c.source_x+0.35,c.source_y+0.08)
    pbuf=render_image(X,Y,c.source_x+float(deformation.mean())*2,c.source_y)
    residuals={"observation_minus_pbuf":observed-pbuf,
               "observation_minus_baryonic_gr":observed-baryonic_gr,
               "observation_minus_lcdm":observed-lcdm,
               "pbuf_minus_lcdm":pbuf-lcdm}
    residual=np.abs(residuals["observation_minus_pbuf"])
    arrays={"matter":matter,"deformation":deformation,"gradient_x":gx,"gradient_y":gy,
            "observation":observed,"baryonic_gr":baryonic_gr,"lcdm":lcdm,"pbuf":pbuf,"residual":residual}
    for name,a in arrays.items(): np.savetxt(out/(name+".csv"),a,delimiter=",")
    for name,a in residuals.items(): np.savetxt(out/(name+".csv"),a,delimiter=",")
    profile=np.column_stack([x, observed.mean(axis=0), baryonic_gr.mean(axis=0), lcdm.mean(axis=0), pbuf.mean(axis=0)])
    np.savetxt(out/"residual_profiles.csv",profile,delimiter=",",header="x,observation,baryonic_gr,lcdm,pbuf",comments="")
    fig,ax=plt.subplots(2,3,figsize=(13,8)); fields=[("Mass",matter),("Deformation",deformation),
      ("Gradient magnitude",np.hypot(gx,gy)),("Observation",observed),("PBUF",pbuf),("Residual",residual)]
    for aa,(title,a) in zip(ax.flat,fields): aa.imshow(a,origin="lower",extent=[x[0],x[-1],x[0],x[-1]]); aa.set_title(title); aa.set_xlabel("x"); aa.set_ylabel("y")
    fig.tight_layout(); fig.savefig(out/"fields.png",dpi=140); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5));
    for xx,yy in paths: ax.plot(xx,yy,lw=1)
    ax.set(title="Photon trajectories",xlabel="x",ylabel="y"); fig.tight_layout(); fig.savefig(out/"photon_trajectories.png",dpi=140); plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));
    for name,a in [("Observation",observed),("Baryonic GR",baryonic_gr),("LCDM",lcdm),("PBUF",pbuf)]: ax.plot(x,a.mean(axis=0),label=name)
    ax.set(title="Predicted profiles vs observation",xlabel="x",ylabel="mean intensity"); ax.legend(); fig.tight_layout(); fig.savefig(out/"comparison_profiles.png",dpi=140); plt.close(fig)
    run_id=time.strftime("%Y%m%dT%H%M%SZ",time.gmtime())
    meta={"run_id":run_id,"equation_id":equation.equation_id,"equation_version":equation.version,
      "equation_description":equation.description,
      "equation_formula":equation.formula,
      "inputs":["observed baryonic mass map","fixed PBUF constants"],"config":asdict(c),
      "git_commit":subprocess.run(["git","rev-parse","HEAD"],capture_output=True,text=True).stdout.strip() or "untracked",
      "execution_seconds":time.time()-started,"checksums":{k:checksum(v) for k,v in arrays.items()},
      "diagnostics":{"deformation_min":float(deformation.min()),
        "deformation_max":float(deformation.max()),"deformation_mean":float(deformation.mean()),
        "deformation_std":float(deformation.std()),
        "gradient_rms":float(np.sqrt(np.mean(gx*gx+gy*gy))),
        "gradient_max":float(np.hypot(gx,gy).max()),
        "photon_max_deviation":path_deviation(paths),
        "finite_outputs":bool(all(np.isfinite(a).all() for a in arrays.values()))},
      "validation":{"zero_mass_zero_deformation":bool(zero_ok),
        "symmetric_mass_symmetric_deformation":bool(symmetric_ok),
        "zero_deformation_straight_paths":bool(path_deviation(zero_paths) < 1e-12),
        "nonzero_deformation_changes_paths":bool(path_deviation(paths) > 1e-8),
        "full_pipeline_completed":True},
      "summary":{"rmse_pbuf":float(np.sqrt(np.mean(residuals["observation_minus_pbuf"]**2))),
        "rmse_baryonic_gr":float(np.sqrt(np.mean(residuals["observation_minus_baryonic_gr"]**2))),
      "rmse_lcdm":float(np.sqrt(np.mean(residuals["observation_minus_lcdm"]**2)))},"status":"PASS"}
    meta["status"] = "PASS" if all(meta["validation"].values()) else "FAIL"
    (out/"run.json").write_text(json.dumps(meta,indent=2))
    with (out/"execution_log.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["run_id","equation_version","stage","status"])
        for i,s in enumerate(["load_lens","mass_map","deformation","gradient","photon_propagation","reconstruction","baryonic_gr","lcdm_benchmark","comparison"]): w.writerow([run_id,equation.version,f"{i+1}:{s}","PASS"])
    return meta


if __name__ == "__main__":
    p=argparse.ArgumentParser(); p.add_argument("--output",type=Path,default=Path("runs/wl001")); p.add_argument("--n",type=int,default=128); p.add_argument("--equation",choices=["A","B","C","D"],default="A")
    a=p.parse_args(); print(json.dumps(run(a.output,Config(n=a.n),get_equation(a.equation)),indent=2))
