"""
Test an initial setup where the Delta_tau and Delta_infty are provided instead of G0_iw.
Make sure the Delta_infty and Delta_tau of the solver aren't changed. 
 """

import numpy as np
from triqs_cthyb import Solver
from triqs.operators import *
from triqs.gf import *
from triqs.utility.comparison_tests import *

beta = 10.0

gf_struct = [['0', [0, 1]]]
target_shape = [2, 2]

nw = 48
nt = 3 * nw

S = Solver(beta=beta, gf_struct=gf_struct, n_iw=nw, n_tau=nt, Delta_interface = True)

h_int = n('0', 0) * n('0', 1)

wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)
tmesh = MeshImTime(beta=beta, S='Fermion', n_max=nt)

Delta_iw = Gf(mesh=wmesh, target_shape=target_shape)

Ek = np.array([
    [1.00,  0.75],
    [0.75, -1.20],
    ])

E_loc = np.array([
    [0.2, 0.3],
    [0.3, 0.4],
    ])

V = np.array([
    [1.0, 0.25],
    [0.25, -1.0],
    ])

Delta_iw << inverse( iOmega_n - Ek ) + inverse( iOmega_n + Ek )
Delta_iw.from_L_G_R(V, Delta_iw, V)

Delta_tau = Gf(mesh=tmesh, target_shape=target_shape)
known_moments = make_zero_tail(Delta_iw, 1)
Delta_tail, Delta_tail_err = Delta_iw.fit_hermitian_tail(known_moments)
Delta_tau << Fourier(Delta_iw, Delta_tail)

S.Delta_tau['0'] << Delta_tau

S.solve(
    h_int = h_int,
    length_cycle = 10,
    n_warmup_cycles = 1,
    n_cycles = 1,
    Delta_infty = [E_loc]
    )

assert_gfs_are_close(Delta_tau, S.Delta_tau['0'])
np.testing.assert_array_almost_equal(S.Delta_infty[0], E_loc)
