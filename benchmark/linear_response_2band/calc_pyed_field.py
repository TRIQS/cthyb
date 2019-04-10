
""" Test calculation for two-band Hubbard atom with two bath sites.

Use pyed to compute the linear response in the single-particle density 
matrix of the system in a general applied quadratic field, giving
the rank 4 tensor response of the system.

Compare with previous pyed results using one-time dynamic response functions.

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import copy
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.utility import mpi
from pytriqs.archive import HDFArchive
from pytriqs.gf import GfImTime, GfImFreq
from pytriqs.gf import MeshImTime, MeshProduct, Gf
from pytriqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# ----------------------------------------------------------------------
if __name__ == '__main__':

    if mpi.is_master_node():
        with HDFArchive('data_model.h5','r') as A: m = A["p"]
        with HDFArchive('data_pyed.h5','r') as A: p = A["p"]
    else:
        m, p = None, None
    m, p = mpi.bcast(m), mpi.bcast(p)

    p.chi_field = np.zeros((4, 4, 4, 4), dtype=np.complex)
    p.g_tau_field = {}

    g_tau = GfImTime(name=r'$g$', beta=m.beta,
                     statistic='Fermion', n_points=50,
                     target_shape=(4, 4))

    # -- The field is symmetric in (i1, i2)
    # -- only calculate upper triangle
    
    index_list = []
    for i1 in xrange(4):
        for i2 in xrange(i1, 4):
            index_list.append((i1, i2))

    #F = 0.0001 / m.beta
    F = 0.1 / m.beta
    F_vec = np.array([-F, F])
    
    work_list = np.array(index_list)
    work_list = mpi.slice_array(work_list)

    for i1, i2 in work_list:

        o1, o2 = m.op_imp[i1], m.op_imp[i2]
        O1 = dagger(o1) * o2
        O1 = O1 + dagger(O1)

        ed_p = TriqsExactDiagonalization(m.H + F * O1, m.op_full, m.beta)
        ed_m = TriqsExactDiagonalization(m.H - F * O1, m.op_full, m.beta)
        
        p.g_tau_field[(i1, i2)] = g_tau.copy()
        ed_p.set_g2_tau_matrix(p.g_tau_field[(i1, i2)], m.op_imp)

        for i3, i4 in itertools.product(range(4), repeat=2):

            o3, o4 = m.op_imp[i3], m.op_imp[i4]
            O2 = dagger(o3) * o4

            O2_vec = np.zeros(2, dtype=np.complex)
            O2_vec[0] = ed_m.get_expectation_value(O2)
            O2_vec[1] = ed_p.get_expectation_value(O2)

            p.chi_field[i1, i2, i3, i4] = -(O2_vec[1] - O2_vec[0])/(2*F)
            p.chi_field[i2, i1, i3, i4] = p.chi_field[i1, i2, i3, i4] 

            chi_tau = p.chi[i1, i2, i3, i4] + p.chi[i2, i1, i3, i4]
            chi_field = p.chi_field[i1, i2, i3, i4]

            np.testing.assert_almost_equal(chi_tau, chi_field, decimal=3)

            if(np.abs(chi_tau.real) > 1e-6):
                print i1, i2, i3, i4
                print 'chi_tau   = %+2.6f' % chi_tau.real
                print 'chi_field = %+2.6f' % chi_field.real
                print 'diff      = %+2.6E' % (chi_field.real - chi_tau.real)

    p.chi_field = mpi.all_reduce(mpi.world, p.chi_field, lambda x, y : x + y)

    chi = np.copy(p.chi).real
    chi = chi + chi.swapaxes(0, 1)

    np.testing.assert_almost_equal(chi, p.chi_field, decimal=3)
    
    # ------------------------------------------------------------------
    # -- Store to hdf5

    if mpi.is_master_node():
        with HDFArchive('data_pyed_field.h5','w') as res:
            res['p'] = p
            
# ----------------------------------------------------------------------
