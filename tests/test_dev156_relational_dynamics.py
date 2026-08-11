import numpy as np
from pbuf.excitation.native_relational_state import perturbation
from pbuf.excitation.native_bond_state import antisymmetry_error, axis_antisymmetry
from pbuf.excitation.native_relational_dynamics import f01_step, f02_step, f02_inverse, f02_invariant, f03_step, f03_inverse, f03_invariant


def test_bond_orientation_is_antisymmetric_and_direction_is_derived():
    q=perturbation("P02",(9,9,9)); b=np.zeros(q.shape+(3,))
    assert antisymmetry_error(b)==0
    assert np.max(np.abs(axis_antisymmetry(q)[...,0]))>0


def test_f01_is_a_dissipative_null_control():
    q=perturbation("P01",(9,9,9)); q1=f01_step(q)
    assert np.sum(q1*q1)<np.sum(q*q)


def test_f02_and_f03_have_executable_exact_inverse_maps():
    q=perturbation("P02",(9,9,9)); b=np.zeros(q.shape+(3,)); r=np.zeros_like(q)
    q1,b1=f02_step(q,b); q0,b0=f02_inverse(q1,b1)
    assert np.allclose(q0,q) and np.allclose(b0,b)
    q1,r1=f03_step(q,r); q0,r0=f03_inverse(q1,r1)
    assert np.allclose(q0,q) and np.allclose(r0,r)


def test_f03_law_derived_quadratic_invariant():
    q=perturbation("P04",(11,11,11)); r=np.zeros_like(q); e=f03_invariant(q,r)
    for _ in range(20): q,r=f03_step(q,r)
    assert np.isclose(f03_invariant(q,r),e,rtol=1e-11,atol=1e-11)


def test_f02_node_plus_bond_invariant():
    q=perturbation("P03",(11,11,11)); b=np.zeros(q.shape+(3,)); e=f02_invariant(q,b)
    for _ in range(20): q,b=f02_step(q,b)
    assert np.isclose(f02_invariant(q,b),e,rtol=1e-11,atol=1e-11)
