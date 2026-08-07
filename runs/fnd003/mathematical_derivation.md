# FND-003 mathematical derivation and boundary

## Conditional minimality theorem

Let the microscopic state at each point be a finite-dimensional real vector space `V`, and require spatial rotations to act through a continuous, linear, faithful representation `rho: SO(3) -> GL(V)`. A one-dimensional real representation is trivial because connected `SO(3)` has no nontrivial continuous homomorphism to `R*`. A two-dimensional faithful representation would embed the three-dimensional Lie group `SO(3)` into `GL(2,R)`; its maximal compact subgroup is conjugate to `O(2)`, whose connected part is only one-dimensional, so no such embedding exists. The defining action on `R^3` is faithful. Therefore `dim(V)>=3`, and `V=R^3` realizes the minimum.

This derives **three only from the added faithful-linear-vector premises**. Three-dimensional space by itself does not require a microscopic state to transform faithfully: a scalar state is one-dimensional, a symmetric tensor has six components, and multiple fields give arbitrarily many. The PBUF ontology supplied here does not establish those added premises, so the overall milestone cannot promote the conditional theorem to an unconditional derivation.

## Scalarization obstruction

For a linear scalar map `l:R^3->R` to be rotation invariant, `l(Rq)=l(q)` for every `R`. Writing `l(q)=e.q` implies `R^T e=e` for every rotation, hence `e=0`. Thus no nonzero direction-free linear scalarization exists. CORE-001's `u=C_L[e.q]` is covariant only if `e` is additional physical data that transforms with `q`. The divergence and norm alternatives avoid a fixed direction but define different effective theories.

## Coupling identifiability

Corrected CORE-001 contains `-epsilon_* g_dev eta e.q`. The normalized microscopic source therefore depends directly on `g_dev`, and the former inverse-rescaling degeneracy is withdrawn. This makes `g_dev` identifiable inside a fully specified microscopic response calculation, but it does not derive the value `1/137`: no supplied symmetry or dynamical principle selects that number, and downstream closure and access maps remain incomplete.
