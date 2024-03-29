
""" Test for bug occuring in complex G0_iw when the resulting h_loc 
has zero real-part. In this case h_loc is errorenously set to zero.

GitHub issue 81, by Robert Triebl and Gernot Kraberger

Author: Hugo Strand """

# ----------------------------------------------------------------------

import triqs.utility.mpi as mpi

from triqs.gf import *
from triqs.operators import *
from h5 import HDFArchive
from triqs.utility.comparison_tests import *

# ----------------------------------------------------------------------

from triqs_cthyb import *

# ----------------------------------------------------------------------
def get_h_loc(eps0, eps1, V, orb):
    
    h_loc_op = eps0 * n(orb, 0) + eps1 * n(orb, 1) \
                + V * c_dag(orb, 0) * c(orb, 1) \
                + V.conjugate() * c_dag(orb, 1) * c(orb, 0)

    h_loc_mat = np.array([
        [eps0, V],
        [V.conjugate(), eps1]
        ])

    return h_loc_op, h_loc_mat

# ----------------------------------------------------------------------
def check_h_loc(h_loc_ref, h_loc_mat, orb):
    
    S = Solver(beta=40, gf_struct={orb : 2})

    S.G0_iw << inverse(Omega - h_loc_mat)

    h_int = Operator()
    S.solve(h_int=h_int, n_cycles=1)

    print('h_loc_ref =', h_loc_ref)
    print('S.h_loc =', S.h_loc)

    diff = h_loc_ref - S.h_loc

    print('diff =', diff)

    for op, prefactor in diff:
        assert( np.abs(prefactor) < 1e-10 )

# ----------------------------------------------------------------------

orb = 'ud'

# -- Test quadratic G0_iw with imaginary offidagonal terms

h_loc_ref, h_loc_mat = get_h_loc(V=0.1j, eps0=1., eps1=3., orb=orb)
check_h_loc(h_loc_ref, h_loc_mat, orb=orb)

# -- Test quadratic G0_iw with imaginary offidagonal terms
# -- but with no diagonal terms

h_loc_ref, h_loc_mat = get_h_loc(V=0.1j, eps0=0., eps1=0., orb=orb)
check_h_loc(h_loc_ref, h_loc_mat, orb=orb)
