"""
Test the initial setup that computes the hybridization function Delta_tau
and the local quadratic term of the Hamiltonain form G0_iw.

Author: Hugo U.R. Strand """

import numpy as np

from triqs_cthyb import SolverCore

from triqs.operators import n, c, c_dag, Operator
import triqs.utility.mpi as mpi
from triqs.gf import Gf, MeshImFreq, MeshImTime, iOmega_n, inverse, Fourier

beta = 10.0

gf_struct = [['0', 2]]
target_shape = [2, 2]

nw = 124
nt = 6 * nw + 1

S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw=nw, n_tau=nt)

h_int = n('0', 0) * n('0', 1)

wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=nw)
tmesh = MeshImTime(beta=beta, S='Fermion', n_max=nt)

Delta_iw = Gf(mesh=wmesh, target_shape=target_shape)

Ek = np.array([
    [1.00,  0.75 + 0.25j],
    [0.75 - 0.25j, -1.20],
    ])

E_loc = np.array([
    [0.2, 0.3 - 0.1j],
    [0.3 + 0.1j, 0.4],
    ])

V = np.array([
    [1.0, 0.25 - 0.3j],
    [0.25 + 0.3j, -1.0],
    ])

h_loc = c_dag('0', 0) * E_loc[0, 0] * c('0', 0) + \
        c_dag('0', 0) * E_loc[0, 1] * c('0', 1) + \
        c_dag('0', 1) * E_loc[1, 0] * c('0', 0) + \
        c_dag('0', 1) * E_loc[1, 1] * c('0', 1)

Delta_iw << inverse( iOmega_n - Ek ) + inverse( iOmega_n + Ek )
Delta_iw.from_L_G_R(V, Delta_iw, V)

Delta_tau = Gf(mesh=tmesh, target_shape=target_shape)
Delta_tau << Fourier(Delta_iw)

G0_iw = Gf(mesh=wmesh, target_shape=target_shape)
G0_iw << inverse( iOmega_n - Delta_iw - E_loc )

S.G0_iw << G0_iw

S.solve(
    h_int = h_int,
    length_cycle = 10,
    n_warmup_cycles = 1,
    n_cycles = 1,
    )

h_loc_ref = S.h_loc - h_int

print('h_loc =\n', h_loc)
print('h_loc_ref =\n', h_loc_ref)

Delta_tau_ref = S.Delta_tau['0']
Delta_iw_ref = Delta_iw.copy()
Delta_iw_ref << Fourier(Delta_tau_ref)

diff = h_loc - h_loc_ref
print('h_loc diff =', diff)
for ops, prefactor in diff:
    print(prefactor, ops)
    assert( np.abs(prefactor) < 1e-10 )

diff = np.max(np.abs(Delta_tau.data - Delta_tau_ref.data))
print('Delta_tau diff =', diff)
np.testing.assert_array_almost_equal(Delta_tau.data, Delta_tau_ref.data)
assert( diff < 1e-8 )

diff = np.max(np.abs(Delta_iw.data - Delta_iw_ref.data))
print('Delta_iw diff =', diff)
np.testing.assert_array_almost_equal(Delta_iw.data, Delta_iw_ref.data)
assert( diff < 1e-7 )

if False:
    from triqs.plot.mpl_interface import oplot, oplotr, oploti, plt
    oplotr(Delta_tau - Delta_tau_ref)
    oploti(Delta_tau - Delta_tau_ref)
    plt.show(); exit()

