import numpy as np

from triqs.gf import Gf, BlockGf, iOmega_n, inverse, Fourier, Wilson
from triqs_cthyb import Solver
from triqs.atom_diag import trace_rho_op
from triqs.operators import n, c, c_dag
import triqs.utility.mpi as mpi

# paramagnetic one-orbital case

D,V,U = 1., 0.2, 4.
ef,beta = -U/2., 50

H = U*n('up_0',0)*n('down_0',0)

S = Solver(beta=beta, gf_struct=[('up_0',1), ('down_0',1)])

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

# exact moments
n_up = trace_rho_op(S.density_matrix, n('up_0',0), S.h_loc_diagonalization)
n_dn = trace_rho_op(S.density_matrix, n('down_0',0), S.h_loc_diagonalization)

np.testing.assert_almost_equal(S.Sigma_Hartree['up_0'][0,0], U*n_dn)
np.testing.assert_almost_equal(S.Sigma_Hartree['down_0'][0,0], U*n_up)

np.testing.assert_almost_equal(S.Sigma_moments[0]['up_0'][0,0], U*n_dn)
np.testing.assert_almost_equal(S.Sigma_moments[0]['down_0'][0,0], U*n_up)

np.testing.assert_almost_equal(S.Sigma_moments[1]['up_0'][0,0], U*U*n_dn*(1-n_dn))
np.testing.assert_almost_equal(S.Sigma_moments[1]['down_0'][0,0], U*U*n_up*(1-n_up))


# magnetic case
S = Solver(beta=beta, gf_struct=[('up_0',1), ('down_0',1)])

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

# exact moments
n_up = trace_rho_op(S.density_matrix, n('up_0',0), S.h_loc_diagonalization)
n_dn = trace_rho_op(S.density_matrix, n('down_0',0), S.h_loc_diagonalization)

np.testing.assert_almost_equal(S.Sigma_Hartree['up_0'][0,0], U*n_dn)
np.testing.assert_almost_equal(S.Sigma_Hartree['down_0'][0,0], U*n_up)

np.testing.assert_almost_equal(S.Sigma_moments[0]['up_0'][0,0], U*n_dn)
np.testing.assert_almost_equal(S.Sigma_moments[0]['down_0'][0,0], U*n_up)

np.testing.assert_almost_equal(S.Sigma_moments[1]['up_0'][0,0], U*U*n_dn*(1-n_dn))
np.testing.assert_almost_equal(S.Sigma_moments[1]['down_0'][0,0], U*U*n_up*(1-n_up))
