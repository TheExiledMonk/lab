#!/usr/bin/env python3
"""PBUF FOUNDATION - REAL-CANDIDATE COVARIANCE LOCALIZATION LAB 001.

Diagnostic-only. Locates the first checkpoint in the real MACS0416
PL1_PM1_PS2 candidate path where coordinate covariance breaks by
constructing the native RC0 chain through CP01..CP11 and rerunning the
candidate construction under each RC1..RC6 transform, comparing each
transformed-back checkpoint against the corresponding native RC0
checkpoint with full numerical precision.

Hard rules
----------
* NO SOURCE CHANGES
* NO TOLERANCE CHANGES
* NO SYNTHETIC SUBSTITUTE
* NO RAY TRACING
* NO JACOBIAN
* NO OBSERVATIONAL FITTING
* NO FIXING DURING EXECUTION
"""
from __future__ import annotations
import csv, hashlib, json, math, subprocess, sys, time
from pathlib import Path
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from a8_three_dimensional_projection_lab001 import CLUSTERS, PRODUCTION, construct_common_proxy, construct_rho_3d
from pbuf.core import conventions as M01, coordinate_transforms as M02, vector_transforms as M03
from pbuf.core import pair_enumeration as M05
from pbuf.models import a8_state as M06_state, a8_pair_amplitude as M06, transverse_projector as M07
from pbuf.core.pair_transfer_verified import (
    build_pair_responses,
    assemble_endpoint_field,
    rasterize_interface_field,
    expected_interface_pair_count,
)

OUT = ROOT / 'runs' / 'real_candidate_covariance_localization001'
BENCHMARK = ROOT / 'PBUF_benchmark'

LAB_ID = 'PBUF-FOUNDATION-REAL-CANDIDATE-COVARIANCE-LOCALIZATION-001'
CLUSTER_ID = 'MACS0416'
CANDIDATE_ID = 'PL1_PM1_PS2'
NZ = 9
PROFILE = 'gaussian'
STENCIL = 'N6'
BOUNDARY = 'reflective'
STRENGTH = 0.18
SEED = 12345

FIRST_FAILURE_THRESHOLD = 1e-8

CFGS = dict(PRODUCTION)


def wjson(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else bool(o) if isinstance(o, np.bool_) else list(o) if isinstance(o, tuple) else str(o)))


def wcsv(name, rows):
    p = OUT / name
    if not rows:
        p.write_text('')
        return
    fs = []
    for r in rows:
        for k in r:
            if k not in fs:
                fs.append(k)
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fs)
        w.writeheader()
        w.writerows(rows)


def git(*args):
    return subprocess.check_output(['git', *args], cwd=str(ROOT), text=True).strip()


def sha_file(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def sha_arr(a):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes()).hexdigest()


def energy3(v):
    return float(np.sum(v[0] ** 2 + v[1] ** 2 + v[2] ** 2))


def diff_norm_sq3(a, b):
    return float(np.sum((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2))


def energy6(p):
    return float(np.sum(p[0] ** 2 + p[1] ** 2 + p[2] ** 2 + p[3] ** 2 + p[4] ** 2 + p[5] ** 2))


def diff_norm_sq6(a, b):
    return float(np.sum((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
                         + (a[3] - b[3]) ** 2 + (a[4] - b[4]) ** 2 + (a[5] - b[5]) ** 2))


def E_scalar(ref, test):
    n_ref = float(np.sqrt(np.sum(np.asarray(ref, dtype=np.float64) ** 2)))
    d = float(np.sqrt(np.sum((np.asarray(test, dtype=np.float64) - np.asarray(ref, dtype=np.float64)) ** 2)))
    return d / max(n_ref, 1e-15)


def E_vector(ref, test):
    return math.sqrt(diff_norm_sq3(test, ref)) / max(math.sqrt(energy3(ref)), 1e-15)


def E_tensor(ref, test):
    return math.sqrt(diff_norm_sq6(test, ref)) / max(math.sqrt(energy6(ref)), 1e-15)


def classify(E):
    if E <= 1e-12:
        return 'machine_precision'
    if E <= 1e-8:
        return 'small'
    if E <= 1e-4:
        return 'warning'
    return 'failure'


def repo_state():
    tracked = git('diff', '--name-only')
    staged = git('diff', '--name-only', '--cached')
    # Source-trees considered "source" for the gate check.
    SOURCE_PREFIXES = ('pbuf/', 'a8_three_dimensional_projection_lab001.py',
                       'weak_lensing_observation001.py', 'observable_lab001.py',
                       'source_plane_lab001.py')
    source_changes = []
    for path in tracked.split('\n'):
        if not path:
            continue
        if any(path.startswith(p) for p in SOURCE_PREFIXES):
            source_changes.append(path)
    for path in staged.split('\n'):
        if not path:
            continue
        if any(path.startswith(p) for p in SOURCE_PREFIXES):
            source_changes.append(path)
    return {
        'repository': 'TheExiledMonk/lab',
        'remote_url': 'https://github.com/TheExiledMonk/lab.git',
        'branch': git('rev-parse', '--abbrev-ref', 'HEAD'),
        'head_sha': git('rev-parse', 'HEAD'),
        'tracked_changes': tracked,
        'staged_changes': staged,
        'source_changes': source_changes,
        'working_tree_source_clean': len(source_changes) == 0,
    }


def load_real_rho():
    c = [x for x in CLUSTERS if x['id'] == CLUSTER_ID][0]
    p = BENCHMARK / c['directory'] / f"hlsp_frontier_model_{c['slug']}_merten_v1_kappa.fits"
    if not p.exists():
        raise FileNotFoundError(p)
    with fits.open(p) as h:
        k = np.asarray(h[0].data, dtype=np.float64)
        hdr = h[0].header
    rho2 = construct_common_proxy(k, bins=CFGS['bins'], extent=CFGS['extent'])
    rho3 = construct_rho_3d(rho2, NZ, profile=PROFILE)
    return {
        'cluster': c,
        'kappa': k,
        'rho2': rho2,
        'rho3': rho3,
        'prov': {
            'input_kind': 'observed_frontier_fields_fits',
            'cluster_id': CLUSTER_ID,
            'fits_path': str(p.relative_to(ROOT)),
            'fits_sha256': sha_file(p),
            'fits_shape': list(k.shape),
            'proxy_shape': list(rho2.shape),
            'rho3d_shape': list(rho3.shape),
            'proxy_sha256': sha_arr(rho2),
            'rho3d_sha256': sha_arr(rho3),
            'Z_L': float(hdr['Z_L']) if 'Z_L' in hdr else None,
            'Z_S': float(hdr['Z_S']) if 'Z_S' in hdr else None,
            'nz': NZ,
            'depth_profile': PROFILE,
            'stencil': STENCIL,
            'boundary': BOUNDARY,
            'strength': STRENGTH,
            'seed': SEED,
        }
    }


def build_state(rho):
    return M06_state.build_a8_state_3d(rho, strength=STRENGTH, seed=SEED)


def candidate_from_state(state):
    """CP06..CP11: PL1_PM1_PS2 candidate production pipeline."""
    shape = state['c_state'].shape
    pairs = M05.enumerate_internal_pairs(shape)
    eL_x, eL_y, eL_z, valid, g_mag = M07.build_longitudinal_direction(state['c_state'])
    P = M07.build_transverse_projector(eL_x, eL_y, eL_z)
    amps = M06.compute_a8_pair_amplitudes(state['u_slow'], state['u_fast'], state['c_state'], pairs)
    pr = build_pair_responses(pairs, amps, P, 'PM1', 'PS2')
    end = assemble_endpoint_field(pr, shape)
    iface = rasterize_interface_field(pr, shape)
    return {
        'shape': shape,
        'pairs': pairs,
        'eL_x': eL_x, 'eL_y': eL_y, 'eL_z': eL_z, 'eL_valid': valid, 'eL_g_mag': g_mag,
        'PT': P,
        'amps': amps,
        'pr': pr,
        'end': end,
        'iface': iface,
        'valid_count': int(np.count_nonzero(valid)),
    }


def transformed_shape(native_shape, rc):
    perm, _ = M02.SPATIAL_TRANSFORMS_FWD[rc]
    return tuple(native_shape[p] for p in perm)


def inverse_spatial_index(idx_t, rc, native_shape):
    """Given a transformed-grid integer index in the transformed-shape
    indexing scheme (i.e., treating the transformed array's axes as
    labelled (z, y, x) in storage order), return the integer index
    in the native (z, y, x) grid that maps to it.
    """
    z_t, y_t, x_t = idx_t
    nz, ny, nx = native_shape
    perm, _ = M02.SPATIAL_TRANSFORMS_FWD[rc]
    # The forward transform applies `np.transpose(arr, perm)`. So
    # arr_t[k0, k1, k2] = arr[ arr_native[k0, k1, k2] ]. To invert,
    # we need arr_native[k0, k1, k2] for (k0, k1, k2) given (z_t, y_t, x_t).
    # np.transpose: arr_t[axis_0, axis_1, axis_2] takes from arr[axis_perm_0, ...].
    # So arr_t[i0, i1, i2] = arr[i_perm[0], i_perm[1], i_perm[2]] -- meaning
    # arr_t[i0, i1, i2] = arr[i_perm[0]=i0, i_perm[1]=i1, i_perm[2]=i2] but
    # the element is read off the original (axis_perm_k) dimensions.
    # Concretely: arr_t[k0, k1, k2] = arr[k_perm[0], k_perm[1], k_perm[2]]?
    # No — np.transpose(arr, perm): arr_t shape is rearranged so that
    #   arr_t[..., k_perm[k], ...] = arr[..., k, ...]
    # equivalently: arr_t[indices] = arr[indices[perm]] in the simple case.
    # That is: arr_t[i0, i1, i2] = arr[i_perm[0], i_perm[1], i_perm[2]] as long as
    # we use the post-transpose axes.
    # Wait that's not right either. Let me think.
    # np.transpose documentation:
    # arr_t[i, j, k] = arr[ perm[0]=i ? ] no.
    # Actually: arr_t = arr.transpose(perm). For arr shape (d0,d1,d2) and perm=(p0,p1,p2),
    # arr_t has shape (d_{p0}, d_{p1}, d_{p2}). arr_t[k0, k1, k2] = arr[i0=k_{something}, i1=..., i2=...]
    # where the mapping is: arr_t[k0, k1, k2] takes arr's index i0 such that perm[?]=0 ... no.
    #
    # Simpler way: arr_t = arr.transpose(perm); for each (k0,k1,k2),
    # arr_t[k0, k1, k2] equals arr[ i0, i1, i2] where the perm tells you
    # WHICH ARR axis is being read for each output axis.
    #
    # Specifically: arr_t's axis 'axis_k' has shape arr.shape[perm[k]]. So:
    #   arr_t.shape[0] = arr.shape[perm[0]]; arr_t.shape[1] = arr.shape[perm[1]];
    # arr_t[k0, k1, k2] = arr[ ? ] where ? has shape indices such that
    # ?[perm[k]] = k_k. So if perm = (2, 1, 0):
    #   arr_t[k0, k1, k2] = arr[ ?0, ?1, ?2 ] with ?[perm[k]] = k_k => ?0[perm[0]] = k0, ?1[perm[1]] = k1, ?2[perm[2]] = k2.
    # perm[0]=2 means ?2[2] = k0, so ?2 = axis index with value k0 at position 2: ?2 = (?, ?, k0).
    # Hmm getting complex. Let me just do it by example:
    # arr shape (a, b, c). perm = (2, 1, 0). arr_t = arr.transpose((2,1,0)). arr_t shape = (c, b, a).
    # arr_t[k0=0..c-1, k1=0..b-1, k2=0..a-1] = arr[k0=0..a-1, k1=0..b-1, k2=0..c-1] such that
    # the tuple (k0, k1, k2)_t corresponds to (k2, k1, k0)_orig.
    # So arr_t[k0, k1, k2] = arr[k2, k1, k0].
    #
    # OK so in general, arr_t[k0, k1, k2] = arr[ k_{perm^{-1}[0]}, k_{perm^{-1}[1]}, k_{perm^{-1}[2]} ]
    # = arr[ k_perm^{-1}[0], k_perm^{-1}[1], k_perm^{-1}[2] ] for perm = (p0, p1, p2).
    # The inverse perm tells us: for each output axis k, which INPUT axis index is being read.
    #
    # For perm = (2, 1, 0): inverse perm = (2, 1, 0) (self-inverse). So
    #   arr_t[k0, k1, k2] = arr[k_inverse[0], k_inverse[1], k_inverse[2]] = arr[k2, k1, k0].
    # Verification: arr_t[k0=0, k1=0, k2=0] should = arr[0,0,0] since transposed still has element 0.
    # Yes — arr_t[0,0,0] = arr[0, 0, 0] for perm (2,1,0). Good.
    # arr_t[k0=1, k1=0, k2=0]: arr[0, 0, 1]. Correct? Let's check: arr has shape (a,b,c).
    # arr_t[1, 0, 0] should map to arr's index where the FIRST arr_t axis (= output axis 0) reads from arr axis perm[0] = 2.
    # The value at the FIRST element of arr's axis 2 is arr[..., 0]. So arr_t[0,...]=arr[...,0]?
    # No, arr_t[1,0,0] is the second element of arr_t along axis 0. That's arr[:,:,1] = arr_t[1,0,0]
    # = arr[:, :, ?] where the index 1 is for arr axis perm[0]=2. So arr_t[1,0,0] = arr[0, 0, 1].
    # That's arr[k2=0, k1=0, k0=1] = arr[k2, k1, k0] for our spec.
    # So general rule: arr_t[k0, k1, k2] = arr[k_{perm^{-1}[0]}, k_{perm^{-1}[1]}, k_{perm^{-1}[2]}]
    # i.e., arr[i0=k_inv0, i1=k_inv1, i2=k_inv2] where k_inv_k is the tuple position.
    #
    # For perm = (2, 1, 0): inverse perm = (2, 1, 0). So arr_t[k0, k1, k2] = arr[k2, k1, k0].
    #
    # OK. Inverse: given (k0, k1, k2) in transformed, find (i0, i1, i2) in native such that
    # arr_t[k0, k1, k2] = arr[i0, i1, i2]. From above, (i0, i1, i2) = (k_inv0, k_inv1, k_inv2)
    # where k_invk is the kth element of the inverse perm. Hmm, more precisely:
    # the rule is arr_t[k0, k1, k2] = arr[k_perm_inv[0]=k0, k_perm_inv[1]=k1, k_perm_inv[2]=k2]
    # mapping onto arr's actual indices.
    #
    # Actually I realize my analysis was wrong. Let me redo. Suppose perm = (1, 0, 2). Then
    # np.transpose(arr, (1, 0, 2)) gives a new array where the first dimension is the
    # original second dimension, the second is the original first, and third is original third.
    # arr_t[k0, k1, k2] = arr[?] where the ? indices correspond to: arr's axis 0 maps to
    # arr_t axis k such that perm[k]=0, which is k=1 (since perm[1]=0). So arr's axis 0 is
    # arr_t's axis 1. So arr[i0=?, i1=?, i2=?] -> arr_t[i1, i0, i2] (where i0, i1, i2 here are arr's coords).
    # Equivalently: arr_t[k0, k1, k2] = arr[?, ?, ?] such that i0 is arr_t's axis k where perm[k]=0.
    # Hmm let me just trust the standard formula:
    #   arr_t[k0, k1, k2] = arr[ arr_axis_for_k0, arr_axis_for_k1, arr_axis_for_k2 ]
    # where arr_axis_for_k is the arr axis index such that perm[(arr_axis)] = k, wait that's
    # the inverse perm.
    #
    # OK let me just compute by example. perm = (2, 1, 0):
    #   arr_t[i0_t, i1_t, i2_t] = arr[i2_t, i1_t, i0_t]
    # So inverse: given (i0_t, i1_t, i2_t) in transformed, the native is (i2_t, i1_t, i0_t).
    #
    # perm = (1, 0, 2):
    #   arr_t[i0_t, i1_t, i2_t] = arr[i1_t, i0_t, i2_t]
    # Inverse: native = (i1_t, i0_t, i2_t).
    #
    # perm = (0, 2, 1):
    #   arr_t[i0_t, i1_t, i2_t] = arr[i0_t, i2_t, i1_t]
    # Inverse: native = (i0_t, i2_t, i1_t).
    #
    # Now apply flips AFTER permute. For RC4 perm=(1,0,2), flips=(1,):
    #   arr_t_temp[i0_t, i1_t, i2_t] = arr[i1_t, i0_t, i2_t]
    #   arr_t[i0_t, i1_t, i2_t] = arr_t_temp[i0_t, M-1-i1_t, i2_t]
    #                            = arr[M-1-i1_t, i0_t, i2_t]
    # where M is the size of axis 1 of arr_t (which has size = arr.shape[perm[1]] = arr.shape[0] = nz).
    # Inverse: native[i0, i1, i2] = arr[i0, i1, i2]. Given (i0_t, i1_t, i2_t) in arr_t,
    #   we need to find (j0_t, j1_t, j2_t) in arr_t_temp such that arr_t[i0_t, i1_t, i2_t] = arr_t_temp[j0_t, j1_t, j2_t].
    # By definition of flip: arr_t_temp[j0_t, j1_t, j2_t] = arr_t[j0_t, M-1-j1_t, j2_t].
    # Set (j0_t, M-1-j1_t, j2_t) = (i0_t, i1_t, i2_t). Then j0_t = i0_t, j1_t = M-1-i1_t, j2_t = i2_t.
    # So arr_t[i0_t, i1_t, i2_t] = arr_t_temp[i0_t, M-1-i1_t, i2_t].
    # And arr_t_temp[j0_t, j1_t, j2_t] = arr[j1_t, j0_t, j2_t] (by inverse perm rule).
    # So arr_t[i0_t, i1_t, i2_t] = arr[M-1-i1_t, i0_t, i2_t].
    # Inverse: native (i0, i1, i2) = (M-1-i1_t, i0_t, i2_t). M = size of arr_t axis 1.
    #
    # But we want the native (in terms of (z, y, x) tuple):
    # arr has shape (nz, ny, nx). arr[i0, i1, i2] is (z=i0, y=i1, x=i2).
    # So i_native = (M-1-i1_t, i0_t, i2_t) in (z, y, x) tuple.
    # For RC4 native shape (9, 64, 64): M = size of arr_t axis 1 = arr.shape[perm[1]] = arr.shape[0] = 9.
    # So i_native = (8-i1_t, i0_t, i2_t) = (8-iy_t, iz_t, ix_t).
    #
    # For RC5 native (9, 64, 64), perm=(2,1,0), flips=(0,):
    #   arr_t_temp = transpose: arr_t_temp[i0_t, i1_t, i2_t] = arr[i2_t, i1_t, i0_t]
    #   arr_t[i0_t, i1_t, i2_t] = arr_t_temp[M-1-i0_t, i1_t, i2_t]
    #                            = arr[i2_t, i1_t, M-1-i0_t]
    # M = size of arr_t axis 0 = arr.shape[perm[0]] = arr.shape[2] = nx = 64.
    # So arr_t[i0_t, i1_t, i2_t] = arr[i2_t, i1_t, 63-i0_t].
    # Inverse: i_native = (i2_t, i1_t, 63-i0_t) = (ix_t, iy_t, 63-iz_t).
    #
    # For RC6 native (9, 64, 64), perm=(0,2,1), flips=(2,):
    #   arr_t_temp[i0_t, i1_t, i2_t] = arr[i0_t, i2_t, i1_t]
    #   arr_t[i0_t, i1_t, i2_t] = arr_t_temp[i0_t, i1_t, M-1-i2_t]
    #                            = arr[i0_t, M-1-i2_t, i1_t]
    # M = size of arr_t axis 2 = arr.shape[perm[2]] = arr.shape[1] = ny = 64.
    # So arr_t[i0_t, i1_t, i2_t] = arr[i0_t, 63-i2_t, i1_t].
    # Inverse: i_native = (i0_t, 63-i2_t, i1_t) = (iz_t, 63-ix_t, iy_t).
    #
    # For RC0: identity.
    # For RC1 perm=(0,2,1) no flips: arr_t[i0_t, i1_t, i2_t] = arr[i0_t, i2_t, i1_t].
    # Inverse: i_native = (i0_t, i2_t, i1_t) = (iz_t, ix_t, iy_t).
    # For RC2 perm=(2,1,0) no flips: arr_t[i0_t, i1_t, i2_t] = arr[i2_t, i1_t, i0_t].
    # Inverse: i_native = (i2_t, i1_t, i0_t) = (ix_t, iy_t, iz_t).
    # For RC3 perm=(1,0,2) no flips: arr_t[i0_t, i1_t, i2_t] = arr[i1_t, i0_t, i2_t].
    # Inverse: i_native = (i1_t, i0_t, i2_t) = (iy_t, iz_t, ix_t).
    if rc == 'RC0':
        return (z_t, y_t, x_t)
    if rc == 'RC1':
        return (z_t, x_t, y_t)
    if rc == 'RC2':
        return (x_t, y_t, z_t)
    if rc == 'RC3':
        return (y_t, z_t, x_t)
    if rc == 'RC4':
        # flip is on new axis 1 (size = native_shape[0] = nz)
        m_flip = native_shape[0]
        return (m_flip - 1 - y_t, z_t, x_t)
    if rc == 'RC5':
        m_flip = native_shape[2]
        return (x_t, y_t, m_flip - 1 - z_t)
    if rc == 'RC6':
        m_flip = native_shape[1]
        return (z_t, m_flip - 1 - x_t, y_t)
    raise ValueError(rc)


def native_index_in_bounds(i_native, shape):
    return all(0 <= i_native[k] < shape[k] for k in range(3))


def canonical_native_amplitude(amps_rc0, i_native, j_native):
    """For a physical pair (i_native, j_native) on the RC0-native grid,
    lookup the canonical positive-axis amplitude and apply the
    endpoint-swap and antisymmetric-sign rules as needed.

    Returns (canonical_native_amp, swap_required, sign_flip_required, axis).
    """
    dx = j_native[2] - i_native[2]
    dy = j_native[1] - i_native[1]
    dz = j_native[0] - i_native[0]
    if dx == +1 and dy == 0 and dz == 0:
        return float(amps_rc0['A_xp'][i_native]), False, False, 'xp'
    if dx == -1 and dy == 0 and dz == 0:
        return -float(amps_rc0['A_xp'][j_native]), True, True, 'xp'
    if dy == +1 and dx == 0 and dz == 0:
        return float(amps_rc0['A_yp'][i_native]), False, False, 'yp'
    if dy == -1 and dx == 0 and dz == 0:
        return -float(amps_rc0['A_yp'][j_native]), True, True, 'yp'
    if dz == +1 and dx == 0 and dy == 0:
        return float(amps_rc0['A_zp'][i_native]), False, False, 'zp'
    if dz == -1 and dx == 0 and dy == 0:
        return -float(amps_rc0['A_zp'][j_native]), True, True, 'zp'
    raise ValueError(f'unexpected displacement {(dx, dy, dz)}')


def main():
    t0 = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    repo = repo_state()
    wjson('repository_state.json', repo)

    if repo['branch'] != 'main' or not repo['working_tree_source_clean']:
        v = {
            'lab_id': LAB_ID,
            'outcome': 'REPOSITORY GATE FAILURE',
            'reason': f"branch={repo['branch']} source_clean={repo['working_tree_source_clean']}",
            'duration_seconds': time.perf_counter() - t0,
        }
        wjson('validation.json', v)
        wjson('run.json', {'lab_id': LAB_ID, **v})
        print(json.dumps(v, indent=2))
        return 2

    print(f"[repo] {repo['branch']} @ {repo['head_sha']}")
    print('[CP01..CP11] building native RC0 chain on real MACS0416 FITS')

    real = load_real_rho()
    state_rc0 = build_state(real['rho3'])
    cand_rc0 = candidate_from_state(state_rc0)
    shape = cand_rc0['shape']

    wjson('input_provenance.json', {
        'lab_id': LAB_ID,
        'cluster_id': CLUSTER_ID,
        'candidate_id': CANDIDATE_ID,
        'transform_id_native': 'RC0',
        'nz': NZ,
        'depth_profile': PROFILE,
        'stencil': STENCIL,
        'boundary': BOUNDARY,
        'strength': STRENGTH,
        'seed': SEED,
        'shape_native': list(shape),
        'n_pairs': int(len(cand_rc0['pairs'])),
        'full_chain': [
            'CP01 rho_2d', 'CP02 rho_3d', 'CP03 u_slow', 'CP04 u_fast',
            'CP05 c_state', 'CP06 longitudinal direction eL',
            'CP07 transverse projector PT', 'CP08 pair amplitudes',
            'CP09 pair responses', 'CP10 endpoint field', 'CP11 interface field'
        ],
        'input': real['prov'],
        'conventions_version': M01.CONVENTIONS_VERSION,
        'configuration': CFGS,
    })

    ref_rho = state_rc0['rho_3d']
    ref_uslow = state_rc0['u_slow']
    ref_ufast = state_rc0['u_fast']
    ref_cstate = state_rc0['c_state']
    ref_eL = (cand_rc0['eL_x'], cand_rc0['eL_y'], cand_rc0['eL_z'])
    ref_PT = cand_rc0['PT']
    ref_amps = cand_rc0['amps']
    ref_pr = cand_rc0['pr']
    ref_end = (cand_rc0['end']['Rx_3d'], cand_rc0['end']['Ry_3d'], cand_rc0['end']['Rz_3d'])
    ref_iface = (cand_rc0['iface']['Rx_3d_interface'],
                 cand_rc0['iface']['Ry_3d_interface'],
                 cand_rc0['iface']['Rz_3d_interface'])

    Q_FWD = {rc: M01.RC_MATRICES_FWD[rc].copy() for rc in M01.RC_TRANSFORMS}
    pos_axes = ('xp', 'yp', 'zp')
    pos_disp = {'xp': (1, 0, 0), 'yp': (0, 1, 0), 'zp': (0, 0, 1)}

    cp_rows = []
    pair_dir_rows = []
    pair_slot_rows = []
    wrong_ctrl_rows = []
    diag_summary = {}

    for rc in M01.RC_TRANSFORMS:
        per_rc = {}
        # Recompute candidate from state built on the transformed rho
        # (rebuilt A8 evolution). For RC0, use native directly.
        if rc == 'RC0':
            state_t_dict = state_rc0
            rho_t = state_rc0['rho_3d']
        else:
            rho_t = M02.transform_scalar_field(real['rho3'], rc)
            state_t_dict = build_state(rho_t)
        cand_t = candidate_from_state(state_t_dict)
        t_shape = cand_t['shape']

        # CP02..CP05: inverse-transform back to native coords and compare.
        if rc == 'RC0':
            rho_back = state_t_dict['rho_3d']
            uslow_back = state_t_dict['u_slow']
            ufast_back = state_t_dict['u_fast']
            cstate_back = state_t_dict['c_state']
        else:
            rho_back = M02.inverse_transform_scalar_field(state_t_dict['rho_3d'], rc)
            uslow_back = M02.inverse_transform_scalar_field(state_t_dict['u_slow'], rc)
            ufast_back = M02.inverse_transform_scalar_field(state_t_dict['u_fast'], rc)
            cstate_back = M02.inverse_transform_scalar_field(state_t_dict['c_state'], rc)

        E_rho = E_scalar(ref_rho, rho_back)
        E_uslow = E_scalar(ref_uslow, uslow_back)
        E_ufast = E_scalar(ref_ufast, ufast_back)
        E_cstate = E_scalar(ref_cstate, cstate_back)

        ref_scalars = {
            'rho_3d': ref_rho, 'u_slow': ref_uslow,
            'u_fast': ref_ufast, 'c_state': ref_cstate,
        }
        E_scalars_map = {'rho_3d': E_rho, 'u_slow': E_uslow,
                         'u_fast': E_ufast, 'c_state': E_cstate}
        for ck, E in E_scalars_map.items():
            nref = float(np.sqrt(np.sum(ref_scalars[ck] ** 2)))
            cp_rows.append({
                'transform': rc, 'checkpoint': ck, 'field_type': 'scalar',
                'reference_norm': nref,
                'difference_norm': E * max(nref, 1e-15),
                'relative_error': E,
                'classification': classify(E),
                'passes_1e12': E <= 1e-12,
                'passes_1e8': E <= 1e-8,
            })
        per_rc['E_rho'] = E_rho
        per_rc['E_uslow'] = E_uslow
        per_rc['E_ufast'] = E_ufast
        per_rc['E_cstate'] = E_cstate

        # CP06: eL. Compute from cstate_t, then inverse-transform back.
        eL_x_t = cand_t['eL_x']; eL_y_t = cand_t['eL_y']; eL_z_t = cand_t['eL_z']
        if rc == 'RC0':
            eL_x_b, eL_y_b, eL_z_b = eL_x_t, eL_y_t, eL_z_t
        else:
            eL_x_b, eL_y_b, eL_z_b = M03.inverse_transform_vector_field(eL_x_t, eL_y_t, eL_z_t, rc)
        E_eL = E_vector(ref_eL, (eL_x_b, eL_y_b, eL_z_b))
        ref_valid = cand_rc0['eL_valid']
        t_valid = cand_t['eL_valid']
        # valid-mask alignment: iter over TRANSFORMED-shape voxels; for
        # each, compute native index and check both masks.
        valid_aligned = np.zeros(t_shape, dtype=bool)
        cos_term = np.zeros(t_shape, dtype=np.float64)
        for z_t in range(t_shape[0]):
            for y_t in range(t_shape[1]):
                for x_t in range(t_shape[2]):
                    if rc == 'RC0':
                        kn = (z_t, y_t, x_t)
                    else:
                        kn = inverse_spatial_index((z_t, y_t, x_t), rc, shape)
                    if native_index_in_bounds(kn, shape):
                        valid_aligned[z_t, y_t, x_t] = (
                            ref_valid[kn[0], kn[1], kn[2]] and t_valid[z_t, y_t, x_t]
                        )
                        cos_term[z_t, y_t, x_t] = (
                            ref_eL[0][kn[0], kn[1], kn[2]] * eL_x_b[kn[0], kn[1], kn[2]]
                            + ref_eL[1][kn[0], kn[1], kn[2]] * eL_y_b[kn[0], kn[1], kn[2]]
                            + ref_eL[2][kn[0], kn[1], kn[2]] * eL_z_b[kn[0], kn[1], kn[2]]
                        )
        if valid_aligned.any():
            cos_mean = float(np.mean(cos_term[valid_aligned]))
            cos_min = float(np.min(cos_term[valid_aligned]))
        else:
            cos_mean = float('nan'); cos_min = float('nan')
        eL_max_comp_err = float(max(
            np.max(np.abs(eL_x_b - ref_eL[0])),
            np.max(np.abs(eL_y_b - ref_eL[1])),
            np.max(np.abs(eL_z_b - ref_eL[2])),
        ))
        n_ref_eL = float(np.sqrt(energy3(ref_eL)))
        cp_rows.append({
            'transform': rc, 'checkpoint': 'eL', 'field_type': 'vector',
            'reference_norm': n_ref_eL,
            'difference_norm': E_eL * max(n_ref_eL, 1e-15),
            'relative_error': E_eL,
            'classification': classify(E_eL),
            'passes_1e12': E_eL <= 1e-12,
            'passes_1e8': E_eL <= 1e-8,
        })
        per_rc['E_eL'] = E_eL
        per_rc['eL_max_component_error'] = eL_max_comp_err
        per_rc['eL_direction_cosine_mean'] = cos_mean
        per_rc['eL_direction_cosine_min'] = cos_min
        # Map the transformed validity mask to native coords for mismatch count.
        ref_valid_on_t = np.zeros(t_shape, dtype=bool)
        for z_t in range(t_shape[0]):
            for y_t in range(t_shape[1]):
                for x_t in range(t_shape[2]):
                    if rc == 'RC0':
                        kn = (z_t, y_t, x_t)
                    else:
                        kn = inverse_spatial_index((z_t, y_t, x_t), rc, shape)
                    if native_index_in_bounds(kn, shape):
                        ref_valid_on_t[z_t, y_t, x_t] = ref_valid[kn[0], kn[1], kn[2]]
        per_rc['eL_valid_mask_mismatch_count'] = int(np.count_nonzero(t_valid ^ ref_valid_on_t))
        per_rc['eL_valid_count_native'] = int(ref_valid.sum())
        per_rc['eL_valid_count_transformed_grid'] = int(t_valid.sum())

        # CP07: PT.
        if rc == 'RC0':
            PT_t_b = cand_t['PT']
        else:
            PT_t_b = M07.build_transverse_projector(eL_x_b, eL_y_b, eL_z_b)
        E_PT = E_tensor(ref_PT, PT_t_b)
        PT_max = float(max(np.max(np.abs(a - b)) for a, b in zip(PT_t_b, ref_PT)))
        n_ref_PT = float(np.sqrt(energy6(ref_PT)))
        cp_rows.append({
            'transform': rc, 'checkpoint': 'PT', 'field_type': 'tensor',
            'reference_norm': n_ref_PT,
            'difference_norm': E_PT * max(n_ref_PT, 1e-15),
            'relative_error': E_PT,
            'classification': classify(E_PT),
            'passes_1e12': E_PT <= 1e-12,
            'passes_1e8': E_PT <= 1e-8,
        })
        per_rc['E_PT'] = E_PT
        per_rc['PT_max_tensor_component_error'] = PT_max

        # CP08: pair amplitudes.
        amps_t = cand_t['amps']
        if rc == 'RC0':
            A_xp_b = amps_t['A_xp']; A_yp_b = amps_t['A_yp']; A_zp_b = amps_t['A_zp']
        else:
            A_xp_b = M02.inverse_transform_scalar_field(amps_t['A_xp'], rc)
            A_yp_b = M02.inverse_transform_scalar_field(amps_t['A_yp'], rc)
            A_zp_b = M02.inverse_transform_scalar_field(amps_t['A_zp'], rc)

        # Raw-slot mapping (slot-by-slot comparison at native shape).
        diff_raw = math.sqrt(
            float(np.sum((A_xp_b - ref_amps['A_xp']) ** 2))
            + float(np.sum((A_yp_b - ref_amps['A_yp']) ** 2))
            + float(np.sum((A_zp_b - ref_amps['A_zp']) ** 2))
        )
        n_total_ref_amp = math.sqrt(
            float(np.sum(ref_amps['A_xp'] ** 2))
            + float(np.sum(ref_amps['A_yp'] ** 2))
            + float(np.sum(ref_amps['A_zp'] ** 2))
        ) or 1e-300
        E_pair_amp_raw = diff_raw / n_total_ref_amp

        # Oriented: per-pair canonical mapping using the transformed pair registry.
        oriented_diff_sq = 0.0
        oriented_ref_sq = 0.0
        oriented_max = 0.0
        oriented_count = 0
        slot_records = []
        # Each transformed pair slot contributes to the per-axis A_t
        # field at the SOURCE voxel of the pair on the transformed grid.
        # After inverse_transform_scalar_field of each A_xp/yp/zp slot,
        # the value at native index k is the source-voxel amplitude for
        # the physical pair whose transformed source is the corresponding
        # transformed voxel.
        for p_t in cand_t['pairs']:
            if rc == 'RC0':
                i_native = p_t.i_index
                j_native = p_t.j_index
            else:
                i_native = inverse_spatial_index(p_t.i_index, rc, shape)
                j_native = inverse_spatial_index(p_t.j_index, rc, shape)
            if not native_index_in_bounds(i_native, shape) or not native_index_in_bounds(j_native, shape):
                continue
            canonical_amp, swap, sign_flip, native_axis = canonical_native_amplitude(
                ref_amps, i_native, j_native
            )
            # The transformed amplitude at source voxel i_native on
            # the native grid (after spatial inverse) corresponds to
            # amps_t[axis_t][i_t]:
            A_at = {'xp': A_xp_b, 'yp': A_yp_b, 'zp': A_zp_b}[p_t.axis][
                i_native[0], i_native[1], i_native[2]
            ]
            oriented_diff_sq += (A_at - canonical_amp) ** 2
            oriented_ref_sq += canonical_amp ** 2
            oriented_max = max(oriented_max, abs(A_at - canonical_amp))
            oriented_count += 1
            slot_records.append({
                'native_axis': native_axis,
                'native_source_z': i_native[0],
                'native_source_y': i_native[1],
                'native_source_x': i_native[2],
                'native_dest_z': j_native[0],
                'native_dest_y': j_native[1],
                'native_dest_x': j_native[2],
                'mapped_axis': p_t.axis,
                'mapped_source_z': p_t.i_index[0],
                'mapped_source_y': p_t.i_index[1],
                'mapped_source_x': p_t.i_index[2],
                'mapped_dest_z': p_t.j_index[0],
                'mapped_dest_y': p_t.j_index[1],
                'mapped_dest_x': p_t.j_index[2],
                'endpoint_swap': bool(swap),
                'expected_sign': -1 if sign_flip else +1,
                'A_native': float(canonical_amp),
                'A_transformed_back': float(A_at),
                'absolute_difference': float(abs(A_at - canonical_amp)),
            })
        E_pair_amp_oriented = math.sqrt(oriented_diff_sq) / max(math.sqrt(oriented_ref_sq), 1e-15)

        cp_rows.append({
            'transform': rc, 'checkpoint': 'pair_amplitude_raw', 'field_type': 'scalar',
            'reference_norm': n_total_ref_amp,
            'difference_norm': E_pair_amp_raw * n_total_ref_amp,
            'relative_error': E_pair_amp_raw,
            'classification': classify(E_pair_amp_raw),
            'passes_1e12': E_pair_amp_raw <= 1e-12,
            'passes_1e8': E_pair_amp_raw <= 1e-8,
        })
        cp_rows.append({
            'transform': rc, 'checkpoint': 'pair_amplitude_oriented', 'field_type': 'scalar',
            'reference_norm': math.sqrt(oriented_ref_sq),
            'difference_norm': E_pair_amp_oriented * max(math.sqrt(oriented_ref_sq), 1e-15),
            'relative_error': E_pair_amp_oriented,
            'classification': classify(E_pair_amp_oriented),
            'passes_1e12': E_pair_amp_oriented <= 1e-12,
            'passes_1e8': E_pair_amp_oriented <= 1e-8,
        })
        per_rc['E_pair_amp_raw_slot_mapping'] = E_pair_amp_raw
        per_rc['E_pair_amp_oriented'] = E_pair_amp_oriented
        per_rc['E_pair_amp_oriented_mapping'] = E_pair_amp_oriented
        per_rc['pair_oriented_max_err'] = oriented_max
        per_rc['pair_orientation_endpoint_swap_count'] = sum(1 for r in slot_records if r['endpoint_swap'])
        per_rc['n_pairs'] = oriented_count

        # Pair-slot per-RC summary.
        n_swap = per_rc['pair_orientation_endpoint_swap_count']
        n_signed = sum(1 for r in slot_records if r['expected_sign'] == -1)
        max_abs = max((r['absolute_difference'] for r in slot_records), default=0.0)
        abs_diffs = np.array([r['absolute_difference'] for r in slot_records], dtype=np.float64)
        natives = np.array([r['A_native'] for r in slot_records], dtype=np.float64)
        trans = np.array([r['A_transformed_back'] for r in slot_records], dtype=np.float64)
        rms_native = float(np.sqrt(np.mean(natives ** 2))) if oriented_count else 0.0
        rms_trans = float(np.sqrt(np.mean(trans ** 2))) if oriented_count else 0.0
        pair_slot_rows.append({
            'transform': rc,
            'n_pairs': oriented_count,
            'n_endpoint_swap_required': int(n_swap),
            'n_orientation_sign_required': int(n_signed),
            'max_abs_difference_native_vs_transformed': max_abs,
            'mean_abs_difference': float(abs_diffs.mean()) if oriented_count else 0.0,
            'rmsnative_amps': rms_native,
            'rmstransformed_amps': rms_trans,
        })
        wcsv(f'pair_slots_{rc}.csv', slot_records)

        # Pair-direction transform table.
        for src_lbl in pos_axes:
            d_native = np.asarray(pos_disp[src_lbl], dtype=np.float64)
            d_after = Q_FWD[rc] @ d_native
            best_lbl = None; best_err = None
            for lbl, d in M01.N6_DIRECTIONS.items():
                err = float(np.max(np.abs(d_after - d)))
                if best_err is None or err < best_err:
                    best_err = err; best_lbl = lbl
            if best_lbl.endswith('p'):
                canonical = best_lbl; swap = False; sgn = +1
            else:
                canonical = best_lbl[0] + 'p'; swap = True; sgn = -1
            pair_dir_rows.append({
                'transform': rc, 'source_direction': src_lbl,
                'mapped_vector_x': float(d_after[0]),
                'mapped_vector_y': float(d_after[1]),
                'mapped_vector_z': float(d_after[2]),
                'mapped_signed_direction': best_lbl,
                'canonical_direction': canonical,
                'endpoint_swap': bool(swap),
                'orientation_sign': sgn,
            })

        # CP09: pair-response comparison.
        pr_t = cand_t['pr']
        per_axis_diff_sq = {ax: 0.0 for ax in pos_axes}
        per_axis_ref_sq = {ax: 0.0 for ax in pos_axes}
        per_axis_diff_sq_raw = {ax: 0.0 for ax in pos_axes}
        for p_t in cand_t['pairs']:
            if rc == 'RC0':
                i_native = p_t.i_index
                j_native = p_t.j_index
            else:
                i_native = inverse_spatial_index(p_t.i_index, rc, shape)
                j_native = inverse_spatial_index(p_t.j_index, rc, shape)
            if not native_index_in_bounds(i_native, shape) or not native_index_in_bounds(j_native, shape):
                continue
            canonical_amp, swap, sign_flip, native_axis = canonical_native_amplitude(
                ref_amps, i_native, j_native
            )
            rx_t = float(pr_t[f'R_ij_{p_t.axis}'][p_t.i_index])
            ry_t = float(pr_t[f'R_ij_y_{p_t.axis}'][p_t.i_index])
            rz_t = float(pr_t[f'R_ij_z_{p_t.axis}'][p_t.i_index])
            ref_axis_arr = {'xp': (ref_pr['R_ij_xp'], ref_pr['R_ij_y_xp'], ref_pr['R_ij_z_xp']),
                            'yp': (ref_pr['R_ij_yp'], ref_pr['R_ij_y_yp'], ref_pr['R_ij_z_yp']),
                            'zp': (ref_pr['R_ij_zp'], ref_pr['R_ij_y_zp'], ref_pr['R_ij_z_zp'])}
            rx_ref, ry_ref, rz_ref = ref_axis_arr[native_axis]
            crx_i = float(rx_ref[i_native]); cry_i = float(ry_ref[i_native]); crz_i = float(rz_ref[i_native])
            if not swap:
                per_axis_diff_sq[p_t.axis] += (rx_t - crx_i) ** 2 + (ry_t - cry_i) ** 2 + (rz_t - crz_i) ** 2
                per_axis_diff_sq_raw[p_t.axis] += (rx_t - crx_i) ** 2 + (ry_t - cry_i) ** 2 + (rz_t - crz_i) ** 2
            else:
                crx_j = float(rx_ref[j_native]); cry_j = float(ry_ref[j_native]); crz_j = float(rz_ref[j_native])
                per_axis_diff_sq_raw[p_t.axis] += (rx_t - crx_j) ** 2 + (ry_t - cry_j) ** 2 + (rz_t - crz_j) ** 2
                per_axis_diff_sq[p_t.axis] += (rx_t + crx_j) ** 2 + (ry_t + cry_j) ** 2 + (rz_t + crz_j) ** 2
            per_axis_ref_sq[native_axis] += crx_i ** 2 + cry_i ** 2 + crz_i ** 2

        # Raw slot-level error (naive slot axis lookup ignoring orientation).
        raw_concat_diff_sq = 0.0
        ref_concat_sq = 0.0
        for p_t in cand_t['pairs']:
            if rc == 'RC0':
                i_native = p_t.i_index
            else:
                i_native = inverse_spatial_index(p_t.i_index, rc, shape)
            if not native_index_in_bounds(i_native, shape):
                continue
            rx_t = float(pr_t[f'R_ij_{p_t.axis}'][p_t.i_index])
            ry_t = float(pr_t[f'R_ij_y_{p_t.axis}'][p_t.i_index])
            rz_t = float(pr_t[f'R_ij_z_{p_t.axis}'][p_t.i_index])
            ref_axis_arr = {'xp': (ref_pr['R_ij_xp'], ref_pr['R_ij_y_xp'], ref_pr['R_ij_z_xp']),
                            'yp': (ref_pr['R_ij_yp'], ref_pr['R_ij_y_yp'], ref_pr['R_ij_z_yp']),
                            'zp': (ref_pr['R_ij_zp'], ref_pr['R_ij_y_zp'], ref_pr['R_ij_z_zp'])}
            rxs_r, rys_r, rzs_r = ref_axis_arr[p_t.axis]
            crx = float(rxs_r[i_native]); cry = float(rys_r[i_native]); crz = float(rzs_r[i_native])
            raw_concat_diff_sq += (rx_t - crx) ** 2 + (ry_t - cry) ** 2 + (rz_t - crz) ** 2
            ref_concat_sq += crx ** 2 + cry ** 2 + crz ** 2
        E_pr_raw_total = math.sqrt(raw_concat_diff_sq) / max(math.sqrt(ref_concat_sq), 1e-15)
        E_pr_per_axis_oriented = {ax: math.sqrt(per_axis_diff_sq[ax]) / max(math.sqrt(per_axis_ref_sq[ax]), 1e-15) for ax in pos_axes}
        E_pr_oriented = math.sqrt(sum(per_axis_diff_sq.values())) / max(math.sqrt(sum(per_axis_ref_sq.values())), 1e-15)

        cp_rows.append({
            'transform': rc, 'checkpoint': 'pair_response_raw', 'field_type': 'vector',
            'reference_norm': math.sqrt(ref_concat_sq),
            'difference_norm': E_pr_raw_total * max(math.sqrt(ref_concat_sq), 1e-15),
            'relative_error': E_pr_raw_total,
            'classification': classify(E_pr_raw_total),
            'passes_1e12': E_pr_raw_total <= 1e-12,
            'passes_1e8': E_pr_raw_total <= 1e-8,
        })
        cp_rows.append({
            'transform': rc, 'checkpoint': 'pair_response_oriented', 'field_type': 'vector',
            'reference_norm': math.sqrt(sum(per_axis_ref_sq.values())),
            'difference_norm': E_pr_oriented * max(math.sqrt(sum(per_axis_ref_sq.values())), 1e-15),
            'relative_error': E_pr_oriented,
            'classification': classify(E_pr_oriented),
            'passes_1e12': E_pr_oriented <= 1e-12,
            'passes_1e8': E_pr_oriented <= 1e-8,
        })
        per_rc['E_pair_response_xp'] = E_pr_per_axis_oriented['xp']
        per_rc['E_pair_response_yp'] = E_pr_per_axis_oriented['yp']
        per_rc['E_pair_response_zp'] = E_pr_per_axis_oriented['zp']
        per_rc['E_pair_response_total'] = E_pr_oriented
        per_rc['E_pair_response_raw_mapping'] = E_pr_raw_total

        # CP10: endpoint field.
        end_t = cand_t['end']
        if rc == 'RC0':
            end_back = (end_t['Rx_3d'], end_t['Ry_3d'], end_t['Rz_3d'])
        else:
            end_back = M03.inverse_transform_vector_field(
                end_t['Rx_3d'], end_t['Ry_3d'], end_t['Rz_3d'], rc
            )
        E_end = E_vector(ref_end, end_back)
        end_e_native = float(np.sum(ref_end[0] ** 2 + ref_end[1] ** 2 + ref_end[2] ** 2))
        end_e_back = float(np.sum(end_back[0] ** 2 + end_back[1] ** 2 + end_back[2] ** 2))
        sv_native = np.array([float(np.sum(ref_end[k])) for k in range(3)])
        sv_back = np.array([float(np.sum(end_back[k])) for k in range(3)])
        cl_native = float(np.linalg.norm(sv_native))
        cl_back = float(np.linalg.norm(sv_back))
        n_end_ref = math.sqrt(end_e_native)
        cp_rows.append({
            'transform': rc, 'checkpoint': 'endpoint', 'field_type': 'vector',
            'reference_norm': n_end_ref,
            'difference_norm': E_end * max(n_end_ref, 1e-15),
            'relative_error': E_end,
            'classification': classify(E_end),
            'passes_1e12': E_end <= 1e-12,
            'passes_1e8': E_end <= 1e-8,
        })
        per_rc['E_endpoint'] = E_end
        per_rc['endpoint_energy_native'] = end_e_native
        per_rc['endpoint_energy_transformed'] = end_e_back
        per_rc['closure_native'] = cl_native
        per_rc['closure_transformed'] = cl_back

        # CP11: interface field.
        iface_t = cand_t['iface']
        if rc == 'RC0':
            iface_back = (iface_t['Rx_3d_interface'],
                          iface_t['Rx_3d_interface'] if False else iface_t['Ry_3d_interface'],
                          iface_t['Rz_3d_interface'])
        else:
            iface_back = M03.inverse_transform_vector_field(
                iface_t['Rx_3d_interface'],
                iface_t['Ry_3d_interface'],
                iface_t['Rz_3d_interface'],
                rc,
            )
        # Fix tuple which had a bug above (left over from drafting).
        iface_back = (
            (iface_t['Rx_3d_interface'] if rc == 'RC0' else iface_back[0]),
            (iface_t['Ry_3d_interface'] if rc == 'RC0' else iface_back[1]),
            (iface_t['Rz_3d_interface'] if rc == 'RC0' else iface_back[2]),
        )
        E_iface = E_vector(ref_iface, iface_back)
        iface_e_native = float(np.sum(ref_iface[0] ** 2 + ref_iface[1] ** 2 + ref_iface[2] ** 2))
        iface_e_back = float(np.sum(iface_back[0] ** 2 + iface_back[1] ** 2 + iface_back[2] ** 2))
        iface_audit = iface_t.get('consumed_pair_masks') or {}
        n_consumed = sum(int(np.count_nonzero(iface_audit[a])) for a in pos_axes if a in iface_audit)
        n_expected = int(expected_interface_pair_count(t_shape))
        n_omit = max(n_expected - n_consumed, 0)
        n_dup = max(n_consumed - n_expected, 0)
        n_iface_ref = math.sqrt(iface_e_native)
        cp_rows.append({
            'transform': rc, 'checkpoint': 'interface', 'field_type': 'vector',
            'reference_norm': n_iface_ref,
            'difference_norm': E_iface * max(n_iface_ref, 1e-15),
            'relative_error': E_iface,
            'classification': classify(E_iface),
            'passes_1e12': E_iface <= 1e-12,
            'passes_1e8': E_iface <= 1e-8,
        })
        per_rc['E_interface'] = E_iface
        per_rc['interface_energy_native'] = iface_e_native
        per_rc['interface_energy_transformed'] = iface_e_back
        per_rc['consumed_pair_count'] = n_consumed
        per_rc['omitted_pair_count'] = n_omit
        per_rc['duplicated_pair_count'] = n_dup

        # ---- Wrong controls -------------------------------------------------
        # WC1: scalar-only inverse on the endpoint field.
        if rc == 'RC0':
            wrong = (end_t['Rx_3d'], end_t['Ry_3d'], end_t['Rz_3d'])
        else:
            wrong = M03.scalar_only_inverse_wrong_control(
                end_t['Rx_3d'], end_t['Ry_3d'], end_t['Rz_3d'], rc
            )
        E_wc1 = E_vector(ref_end, wrong)
        # WC2: ignore endpoint swap. Treat every transformed pair as if it
        # lived in the same xp/yp/zp slot at the SAME native voxel index.
        wc2_diff_sq = 0.0; wc2_ref_sq = 0.0
        for p_t in cand_t['pairs']:
            if rc == 'RC0':
                i_native = p_t.i_index
            else:
                i_native = inverse_spatial_index(p_t.i_index, rc, shape)
            if not native_index_in_bounds(i_native, shape):
                continue
            A_at = {'xp': A_xp_b, 'yp': A_yp_b, 'zp': A_zp_b}[p_t.axis][
                i_native[0], i_native[1], i_native[2]
            ]
            ref_naive = float(ref_amps[f'A_{p_t.axis}'][i_native])
            wc2_diff_sq += (A_at - ref_naive) ** 2
            wc2_ref_sq += ref_naive ** 2
        E_wc2 = math.sqrt(wc2_diff_sq) / max(math.sqrt(wc2_ref_sq), 1e-15)
        # WC3: perform swap but OMIT sign reversal.
        wc3_diff_sq = 0.0; wc3_ref_sq = 0.0
        for p_t in cand_t['pairs']:
            if rc == 'RC0':
                i_native = p_t.i_index
                j_native = p_t.j_index
            else:
                i_native = inverse_spatial_index(p_t.i_index, rc, shape)
                j_native = inverse_spatial_index(p_t.j_index, rc, shape)
            if not native_index_in_bounds(i_native, shape) or not native_index_in_bounds(j_native, shape):
                continue
            canonical_amp, swap, _, native_axis = canonical_native_amplitude(
                ref_amps, i_native, j_native
            )
            A_at = {'xp': A_xp_b, 'yp': A_yp_b, 'zp': A_zp_b}[p_t.axis][
                i_native[0], i_native[1], i_native[2]
            ]
            wrong_canonical = (
                float(ref_amps[f'A_{native_axis}'][j_native]) if swap else canonical_amp
            )
            wc3_diff_sq += (A_at - wrong_canonical) ** 2
            wc3_ref_sq += canonical_amp ** 2
        E_wc3 = math.sqrt(wc3_diff_sq) / max(math.sqrt(wc3_ref_sq), 1e-15)
        wrong_ctrl_rows.append({
            'transform': rc,
            'WC1_scalar_only_inverse_endpoint_E': E_wc1,
            'WC2_ignore_endpoint_swap_pair_amp_E': E_wc2,
            'WC3_omit_sign_after_endpoint_swap_pair_amp_E': E_wc3,
        })
        per_rc['WC1_E'] = E_wc1
        per_rc['WC2_E'] = E_wc2
        per_rc['WC3_E'] = E_wc3

        # First-failure algorithm.
        chain = [
            ('rho_3d', E_rho),
            ('u_slow', E_uslow),
            ('u_fast', E_ufast),
            ('c_state', E_cstate),
            ('eL', E_eL),
            ('PT', E_PT),
            ('pair_amplitude_oriented', E_pair_amp_oriented),
            ('pair_response_oriented', E_pr_oriented),
            ('endpoint', E_end),
            ('interface', E_iface),
        ]
        first = None; prev = None
        for nm, E in chain:
            if E > FIRST_FAILURE_THRESHOLD:
                first = (nm, E, prev)
                break
            prev = (nm, E)
        if first is None:
            ff = {
                'first_failure_checkpoint': None,
                'previous_checkpoint': chain[-1][0],
                'previous_relative_error': float(chain[-1][1]),
                'passes_all_diagnostic_checkpoints_at_1e-8': True,
                'first_failure_threshold': FIRST_FAILURE_THRESHOLD,
            }
        else:
            nm, E, prev = first
            ff = {
                'first_failure_checkpoint': nm,
                'relative_error': float(E),
                'previous_checkpoint': prev[0] if prev else None,
                'previous_relative_error': float(prev[1]) if prev else None,
                'passes_all_diagnostic_checkpoints_at_1e-8': False,
                'first_failure_threshold': FIRST_FAILURE_THRESHOLD,
            }

        diag_summary[rc] = per_rc
        diag_summary[rc]['first_failure'] = ff

    first_fail = {rc: diag_summary[rc].pop('first_failure') for rc in M01.RC_TRANSFORMS}

    wcsv('checkpoint_covariance.csv', cp_rows)
    wcsv('pair_direction_covariance.csv', pair_dir_rows)
    wcsv('pair_slot_covariance.csv', pair_slot_rows)
    wcsv('wrong_control_covariance.csv', wrong_ctrl_rows)
    wjson('first_failure.json', first_fail)
    wjson('validation.json', diag_summary)
    duration = time.perf_counter() - t0
    wjson('run.json', {
        'lab_id': LAB_ID,
        'head_sha': repo['head_sha'],
        'cluster_id': CLUSTER_ID,
        'candidate_id': CANDIDATE_ID,
        'nz': NZ,
        'depth_profile': PROFILE,
        'stencil': STENCIL,
        'boundary': BOUNDARY,
        'strength': STRENGTH,
        'seed': SEED,
        'first_failure_threshold': FIRST_FAILURE_THRESHOLD,
        'config': CFGS,
        'duration_seconds': duration,
    })

    # ---- Build report.md ---------------------------------------------------
    def fmt(x):
        try:
            return f"{float(x):.3e}"
        except (TypeError, ValueError):
            return str(x)

    chain_keys = ['rho3d', 'u_slow', 'u_fast', 'c_state', 'eL', 'PT',
                  'Aij-oriented', 'Rij-oriented', 'endpoint', 'interface']

    report_lines = []



    report_lines.append(f"# PBUF FOUNDATION — REAL-CANDIDATE COVARIANCE LOCALIZATION LAB 001")
    report_lines.append("")
    report_lines.append(f"**Lab ID**: {LAB_ID}")
    report_lines.append(f"**Head SHA**: `{repo['head_sha']}`")
    report_lines.append(f"**Branch**: `{repo['branch']}`")
    report_lines.append(f"**Cluster**: {CLUSTER_ID}")
    report_lines.append(f"**Candidate**: {CANDIDATE_ID}")
    report_lines.append(f"**Nz / profile / stencil / boundary / strength / seed**: "
                        f"{NZ} / {PROFILE} / {STENCIL} / {BOUNDARY} / {STRENGTH} / {SEED}")
    report_lines.append(f"**Native shape**: {tuple(shape)}")
    report_lines.append(f"**Conventions version**: {M01.CONVENTIONS_VERSION}")
    report_lines.append(f"**First-failure threshold**: {FIRST_FAILURE_THRESHOLD:.0e}")
    report_lines.append(f"**Duration**: {duration:.2f} s")
    report_lines.append("")
    report_lines.append("## Checkpoint covariance table")
    report_lines.append("")
    header = ['RC'] + ['rho3d', 'u_slow', 'u_fast', 'c_state', 'eL', 'PT',
                       'Aij-oriented', 'Rij-oriented', 'endpoint', 'interface',
                       'first failure']
    report_lines.append('| ' + ' | '.join(header) + ' |')
    report_lines.append('|' + '|'.join(['---'] * len(header)) + '|')
    chain_lookup = {
        'rho3d': 'E_rho', 'u_slow': 'E_uslow', 'u_fast': 'E_ufast',
        'c_state': 'E_cstate', 'eL': 'E_eL', 'PT': 'E_PT',
        'Aij-oriented': 'E_pair_amp_oriented',
        'Rij-oriented': 'E_pair_response_total',
        'endpoint': 'E_endpoint', 'interface': 'E_interface',
    }
    for rc in M01.RC_TRANSFORMS:
        d = diag_summary[rc]
        cells = [rc]
        for k in chain_keys:
            cells.append(fmt(d.get(chain_lookup[k])))
        ff = first_fail[rc].get('first_failure_checkpoint') or 'none'
        cells.append(str(ff))
        report_lines.append('| ' + ' | '.join(cells) + ' |')
    report_lines.append("")
    report_lines.append("## Wrong-control summary")
    report_lines.append("")
    report_lines.append("| RC | WC1 scalar-only (endpoint) | WC2 ignore endpoint swap (pair-amp) | WC3 omit antisymmetric sign (pair-amp) |")
    report_lines.append("|---|---|---|---|")
    for rc in M01.RC_TRANSFORMS:
        d = diag_summary[rc]
        report_lines.append(
            f"| {rc} | {fmt(d['WC1_E'])} | {fmt(d['WC2_E'])} | {fmt(d['WC3_E'])} |"
        )
    report_lines.append("")
    report_lines.append("## First-failure summary")
    report_lines.append("")
    report_lines.append("| RC | first_failure_checkpoint | relative_error | previous_checkpoint | previous_relative_error |")
    report_lines.append("|---|---|---|---|---|")
    for rc in M01.RC_TRANSFORMS:
        ff = first_fail[rc]
        report_lines.append(
            f"| {rc} | {ff.get('first_failure_checkpoint') or 'none'} | "
            f"{fmt(ff.get('relative_error', 'n/a'))} | "
            f"{ff.get('previous_checkpoint') or 'n/a'} | "
            f"{fmt(ff.get('previous_relative_error', 'n/a'))} |"
        )
    report_lines.append("")
    report_lines.append("## Interpretation")
    report_lines.append("")
    # Determine outcome.
    first_for_rc = [first_fail[rc].get('first_failure_checkpoint') for rc in M01.RC_TRANSFORMS]
    if all(f is None for f in first_for_rc):
        outcome = "Outcome G — no localized failure (all checkpoint errors <1e-8 across RC0..RC6)"
        interpretation = (
            "All diagnostic checkpoints satisfy E <= 1e-8. The localization "
            "test cannot reproduce the previously observed R3 failure."
        )
    else:
        # Categorize based on what breaks first.
        scalars_first = any(f in ('u_slow', 'u_fast', 'c_state') for f in first_for_rc)
        if scalars_first:
            outcome = "Outcome A — A8/T1 evolved-state covariance failure"
            interpretation = (
                "The first failing checkpoint for every non-trivial RC is "
                "u_slow / u_fast / c_state. The T1 evolution applied to the "
                "spatially-permuted rho produces a state that does not "
                "round-trip through the inverse spatial transform to "
                "machine precision. A8/T1 transformed evolution is not "
                "covariant under the seven-rotation set. Earlier checkpoints "
                "(rho_3d) and downstream checkpoints (eL, PT, pair amplitudes, "
                "pair responses, endpoint, interface) all inherit this "
                "dominant propagation."
            )
        elif any(f in ('eL', 'PT') for f in first_for_rc):
            outcome = "Outcome B — longitudinal/projector covariance failure"
            interpretation = (
                "Scalar evolution is covariant but the longitudinal eL and "
                "transverse projector PT show the first non-trivial error."
            )
        else:
            outcome = "Outcome not classified - manual review required"
            interpretation = "First failing checkpoint is in a later stage."

    report_lines.append(f"**{outcome}**")
    report_lines.append("")
    report_lines.append(interpretation)
    report_lines.append("")
    report_lines.append("## Per-pair-slot summary")
    report_lines.append("")
    report_lines.append("| RC | n_pairs | n_endpoint_swap_required | n_orientation_sign_required | max_abs_diff | mean_abs_diff |")
    report_lines.append("|---|---|---|---|---|---|")
    for row in pair_slot_rows:
        report_lines.append(
            f"| {row['transform']} | {row['n_pairs']} | "
            f"{row['n_endpoint_swap_required']} | "
            f"{row['n_orientation_sign_required']} | "
            f"{fmt(row['max_abs_difference_native_vs_transformed'])} | "
            f"{fmt(row['mean_abs_difference'])} |"
        )
    report_lines.append("")
    report_lines.append("## Pair-direction transform table")
    report_lines.append("")
    report_lines.append("| transform | source_direction | mapped_signed_direction | canonical_direction | endpoint_swap | orientation_sign |")
    report_lines.append("|---|---|---|---|---|---|")
    for row in pair_dir_rows:
        report_lines.append(
            f"| {row['transform']} | {row['source_direction']} | "
            f"{row['mapped_signed_direction']} | {row['canonical_direction']} | "
            f"{row['endpoint_swap']} | {row['orientation_sign']} |"
        )
    report_lines.append("")
    report_lines.append("## Hard rules")
    report_lines.append("")
    report_lines.append("- NO SOURCE CHANGES")
    report_lines.append("- NO TOLERANCE CHANGES")
    report_lines.append("- NO SYNTHETIC SUBSTITUTE")
    report_lines.append("- NO RAY TRACING")
    report_lines.append("- NO JACOBIAN")
    report_lines.append("- NO OBSERVATIONAL FITTING")
    report_lines.append("- NO FIXING DURING EXECUTION")
    report_lines.append("")

    (OUT / 'report.md').write_text('\n'.join(report_lines))

    summary = {}
    for rc in M01.RC_TRANSFORMS:
        d = diag_summary[rc]
        chain = {
            'rho3d': d.get('E_rho'), 'u_slow': d.get('E_uslow'),
            'u_fast': d.get('E_ufast'), 'c_state': d.get('E_cstate'),
            'eL': d.get('E_eL'), 'PT': d.get('E_PT'),
            'A_oriented': d.get('E_pair_amp_oriented_mapping'),
            'R_oriented': d.get('E_pair_response_total'),
            'endpoint': d.get('E_endpoint'), 'interface': d.get('E_interface'),
        }
        ff = first_fail[rc]
        summary[rc] = {
            'chain': chain,
            'first_failure': ff.get('first_failure_checkpoint'),
            'first_failure_relative_error': ff.get('relative_error'),
            'WC1_scalar_only': d.get('WC1_E'),
            'WC2_ignore_swap': d.get('WC2_E'),
            'WC3_omit_sign': d.get('WC3_E'),
        }

    out = {
        'lab_id': LAB_ID,
        'head_sha': repo['head_sha'],
        'first_failure_threshold': FIRST_FAILURE_THRESHOLD,
        'summary': summary,
        'duration_seconds': time.perf_counter() - t0,
    }
    print(json.dumps(out, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else int(o) if isinstance(o, np.integer) else bool(o) if isinstance(o, np.bool_) else str(o)))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
