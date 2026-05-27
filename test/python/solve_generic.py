"""
Test for triqs_cthyb.solve_generic.

Structural coverage:
  * MeshImFreq input with postprocess=None, 'dyson', TailFitParams, LegendreParams
  * MeshDLRImFreq input with postprocess='dyson' and CRMParams
  * h_loc0 accepted as both an Operator and as a list of block matrices

Plus an h5diff regression case at the end that catches silent numerical drift.
"""

import numpy as np

import triqs.utility.mpi as mpi
from triqs.gfs import BlockGf, Gf, MeshImFreq, MeshDLRImFreq, iOmega_n, inverse
from triqs.operators import n, c, c_dag
from triqs.operators.util import block_matrix_from_op
from triqs.solver_utils import SolverResults
from h5 import HDFArchive
from triqs.utility.h5diff import h5diff

from triqs_cthyb import solve_generic, TailFitParams, LegendreParams, CRMParams


# ---------- impurity setup (mirrors setup_Delta_tau_and_h_loc.py) ----------
beta = 10.0
n_iw = 48
target_shape = (2, 2)

Ek    = np.array([[ 1.00, 0.75], [0.75, -1.20]])
V_hyb = np.array([[ 1.00, 0.25], [0.25, -1.00]])
E_loc = np.array([[ 0.20, 0.30], [0.30,  0.40]])

wmesh = MeshImFreq(beta=beta, statistic='Fermion', n_iw=n_iw)
g_block = Gf(mesh=wmesh, target_shape=target_shape)
g_block << inverse(iOmega_n - Ek) + inverse(iOmega_n + Ek)
g_block.from_L_G_R(V_hyb, g_block, V_hyb)

Delta_iw = BlockGf(name_list=['0'], block_list=[g_block], make_copies=True)

gf_struct = [('0', 2)]
h_loc0_op = sum(c_dag('0', i) * E_loc[i, j] * c('0', j)
                for i in range(2) for j in range(2))
# exercises both input forms: list-of-matrices (bl) and Operator (op)
h_loc0_bl = block_matrix_from_op(h_loc0_op, gf_struct)
h_int     = n('0', 0) * n('0', 1)

mc_params = dict(
    length_cycle    = 10,
    n_warmup_cycles = 50,
    n_cycles        = 50,
    random_seed     = 123 * mpi.rank + 567,
)


def assert_finite_blockgf(G, label):
    for bl, g in G:
        assert np.all(np.isfinite(g.data)), f"{label}[{bl}] contains NaN/Inf"


# ---------- case 1: postprocess=None (raw output), h_loc0 as Operator ----------
mpi.report("=== solve_generic: postprocess=None (h_loc0 = Operator) ===")
res = solve_generic(Delta_iw, h_loc0_op, h_int, postprocess=None, **mc_params)

assert isinstance(res, SolverResults)
assert res.Solver is not None
assert res.Sigma_iw is None, "raw output should not carry a Sigma_iw"

S = res.Solver
assert_finite_blockgf(S.G_iw,  "S.G_iw")
assert_finite_blockgf(S.G_tau, "S.G_tau")
for bl, g in S.G_iw:
    assert g.target_shape == target_shape
    assert len(g.mesh) == 2 * n_iw


# ---------- case 2: postprocess='dyson', h_loc0 as list of block matrices ----------
mpi.report("=== solve_generic: postprocess='dyson' (h_loc0 = list of matrices) ===")
res = solve_generic(Delta_iw, h_loc0_bl, h_int, postprocess='dyson', **mc_params)

assert isinstance(res, SolverResults)
for field in ('G_iw', 'G_tau', 'Sigma_iw', 'Sigma_dynamic', 'Sigma_HartreeFock'):
    assert getattr(res, field) is not None, f"{field} should be populated"
assert_finite_blockgf(res.G_iw,          "G_iw")
assert_finite_blockgf(res.G_tau,         "G_tau")
assert_finite_blockgf(res.Sigma_iw,      "Sigma_iw")
assert_finite_blockgf(res.Sigma_dynamic, "Sigma_dynamic")

# wiring: Sigma_iw = Sigma_dynamic + Sigma_HartreeFock per block
for i, (bl, s_iw) in enumerate(res.Sigma_iw):
    hf = res.Sigma_HartreeFock[i]
    s_dyn = res.Sigma_dynamic[bl]
    diff = np.max(np.abs(s_iw.data - s_dyn.data - hf[None, :, :]))
    assert diff < 1e-10, f"Sigma_iw - Sigma_dynamic != Sigma_HF for block {bl}: {diff}"


# ---------- case 3: postprocess=TailFitParams(...) ----------
mpi.report("=== solve_generic: postprocess=TailFitParams ===")
res = solve_generic(
    Delta_iw, h_loc0_bl, h_int,
    postprocess=TailFitParams(fit_min_n=20, fit_max_n=40, fit_max_moment=3),
    **mc_params,
)
assert isinstance(res, SolverResults)
for field in ('G_iw', 'G_tau', 'Sigma_iw', 'Sigma_dynamic', 'Sigma_HartreeFock'):
    assert getattr(res, field) is not None, f"{field} should be populated"
assert_finite_blockgf(res.Sigma_iw, "Sigma_iw (tail-fitted)")


# ---------- case 4: postprocess=LegendreParams (post-fit) ----------
mpi.report("=== solve_generic: postprocess=LegendreParams ===")
res = solve_generic(
    Delta_iw, h_loc0_bl, h_int,
    postprocess=LegendreParams(n_l=20, use_measured=False),
    **mc_params,
)
assert isinstance(res, SolverResults)
for field in ('G_iw', 'G_tau', 'G_l', 'Sigma_iw', 'Sigma_dynamic', 'Sigma_HartreeFock'):
    assert getattr(res, field) is not None, f"{field} should be populated"
assert_finite_blockgf(res.G_iw,     "G_iw (legendre)")
assert_finite_blockgf(res.Sigma_iw, "Sigma_iw (legendre)")


# ---------- DLR input path ----------
# Build Delta on a MeshDLRImFreq using the same physical parameters.
w_max   = 8.0
dlr_eps = 1e-10
dlr_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=dlr_eps)

g_block_dlr = Gf(mesh=dlr_mesh, target_shape=target_shape)
g_block_dlr << inverse(iOmega_n - Ek) + inverse(iOmega_n + Ek)
g_block_dlr.from_L_G_R(V_hyb, g_block_dlr, V_hyb)

Delta_iw_dlr = BlockGf(name_list=['0'], block_list=[g_block_dlr], make_copies=True)


# ---------- case 5: MeshDLRImFreq + postprocess='dyson' ----------
mpi.report("=== solve_generic: MeshDLRImFreq + 'dyson' ===")
res = solve_generic(Delta_iw_dlr, h_loc0_bl, h_int, postprocess='dyson', **mc_params)

assert isinstance(res, SolverResults)
for field in ('G_iw', 'G_tau', 'Sigma_iw', 'Sigma_dynamic', 'Sigma_HartreeFock'):
    assert getattr(res, field) is not None, f"{field} should be populated"
assert_finite_blockgf(res.G_iw,     "G_iw (DLR input, dyson)")
assert_finite_blockgf(res.Sigma_iw, "Sigma_iw (DLR input, dyson)")

# wiring check: Sigma_iw = Sigma_dynamic + Sigma_HF per block
for i, (bl, s_iw) in enumerate(res.Sigma_iw):
    hf = res.Sigma_HartreeFock[i]
    s_dyn = res.Sigma_dynamic[bl]
    diff = np.max(np.abs(s_iw.data - s_dyn.data - hf[None, :, :]))
    assert diff < 1e-10, f"Sigma_iw - Sigma_dynamic != Sigma_HF for block {bl}: {diff}"


# ---------- case 6: MeshDLRImFreq + CRMParams ----------
mpi.report("=== solve_generic: MeshDLRImFreq + CRMParams ===")
res = solve_generic(
    Delta_iw_dlr, h_loc0_bl, h_int,
    postprocess=CRMParams(w_max=w_max, eps=dlr_eps),
    **mc_params,
)
assert isinstance(res, SolverResults)
for field in ('G_iw', 'G_tau', 'G_tau_dlr', 'Sigma_iw', 'Sigma_dynamic',
              'Sigma_dlr', 'Sigma_HartreeFock'):
    assert getattr(res, field) is not None, f"{field} should be populated"
assert_finite_blockgf(res.G_iw,     "G_iw (CRM)")
assert_finite_blockgf(res.Sigma_iw, "Sigma_iw (CRM)")


# ---------- case 7: h5diff regression (MeshImFreq + 'dyson') ----------
# Catches silent numerical drift in the Dyson post-processing pipeline.
# Uses a deterministic seed and modest MC statistics.
mpi.report("=== solve_generic: regression (h5diff) ===")
regression_params = dict(
    length_cycle    = 10,
    n_warmup_cycles = 100,
    n_cycles        = 500,
    random_seed     = 567,
)
res = solve_generic(Delta_iw, h_loc0_bl, h_int, postprocess='dyson', **regression_params)

if mpi.is_master_node():
    with HDFArchive('solve_generic.out.h5', 'w') as A:
        A['G_tau']    = res.G_tau
        A['G_iw']     = res.G_iw
        A['Sigma_iw'] = res.Sigma_iw

    h5diff('solve_generic.out.h5', 'solve_generic.ref.h5')

mpi.report("solve_generic tests passed.")
