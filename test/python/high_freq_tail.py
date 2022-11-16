import numpy as np

from triqs.gf import Gf, BlockGf, iOmega_n, inverse, Fourier, Wilson
from triqs_cthyb import Solver
from triqs.atom_diag import trace_rho_op
from triqs.operators import n, c, c_dag
import triqs.utility.mpi as mpi

D,V,U = 1., 0.2, 4.
ef,beta = -U/2., 50

# paramagnetic one-orbital case
H = U*n('up',0)*n('down',0)

S = Solver(beta=beta, gf_struct=[('up',1), ('down',1)])

for name, g0 in S.G0_iw: g0 << inverse(iOmega_n - ef - V**2 * Wilson(D))


S.solve(h_int = H,
        n_cycles = 100000,
        length_cycle = 200,
        n_warmup_cycles = 10000,
        measure_density_matrix=True,
        use_norm_as_weight =True,
        perform_tail_fit=True,
        fit_max_moment = 4,
        fit_min_w = 5,
        fit_max_w = 9
        )

# Sigma_Hartree = G1 - E
E = {term[0][1][0] : coef for term, coef in S.h_loc0}
for spin in ['up','down']: np.testing.assert_almost_equal(S.Sigma_moments[spin][0,0,0], S.G_moments[spin][2,0,0]-E[spin])

# exact moments
n_up = trace_rho_op(S.density_matrix, n('up',0), S.h_loc_diagonalization)
n_dn = trace_rho_op(S.density_matrix, n('down',0), S.h_loc_diagonalization)

for spin, nimp in zip(['up', 'down'], [n_dn, n_up]):
    np.testing.assert_almost_equal(S.Sigma_Hartree[spin][0,0], U*nimp)
    np.testing.assert_almost_equal(S.Sigma_moments[spin][0,0,0], U*nimp)
    np.testing.assert_almost_equal(S.Sigma_moments[spin][1,0,0], U*U*nimp*(1-nimp))


# magnetic one-orbital case
S = Solver(beta=beta, gf_struct=[('up',1), ('down',1)])

for name, g0 in S.G0_iw: 
    if 'up' in name:
        g0 << inverse(iOmega_n - ef+1 - V**2 * Wilson(D))
    else:
        g0 << inverse(iOmega_n - ef-1 - V**2 * Wilson(D))

S.solve(h_int = H,
        n_cycles = 100000,
        length_cycle = 200,
        n_warmup_cycles = 10000,
        measure_density_matrix=True,
        use_norm_as_weight =True
        )

# Sigma_Hartree = G1 - E
E = {term[0][1][0] : coef for term, coef in S.h_loc0}
for spin in ['up','down']: np.testing.assert_almost_equal(S.Sigma_moments[spin][0,0,0], S.G_moments[spin][2,0,0]-E[spin])

# exact moments
n_up = trace_rho_op(S.density_matrix, n('up',0), S.h_loc_diagonalization)
n_dn = trace_rho_op(S.density_matrix, n('down',0), S.h_loc_diagonalization)

for spin, nimp in zip(['up', 'down'], [n_dn, n_up]):
    np.testing.assert_almost_equal(S.Sigma_Hartree[spin][0,0], U*nimp)
    np.testing.assert_almost_equal(S.Sigma_moments[spin][0,0,0], U*nimp)
    np.testing.assert_almost_equal(S.Sigma_moments[spin][1,0,0], U*U*nimp*(1-nimp))
