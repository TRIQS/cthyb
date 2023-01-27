
  
""" Test calculation for two-band Hubbard atom with two bath sites.

Use pyed to compute static and dynamic two-particle responses. 

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import copy
import itertools
import numpy as np

# ----------------------------------------------------------------------

from triqs.utility import mpi
from h5 import HDFArchive
from triqs.gf import GfImTime, GfImFreq
from triqs.gf import MeshImTime, MeshProduct, Gf
from triqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    if mpi.is_master_node():
        with HDFArchive('data_model.h5','r') as A: m = A["p"]
    else: m = None
    m = mpi.bcast(m)

    p = ParameterCollection()
    ed = TriqsExactDiagonalization(m.H, m.op_full, m.beta)

    p.O1_exp = np.zeros((4, 4), dtype=complex)
    p.O2_exp = np.zeros((4, 4), dtype=complex)
    p.chi_dissconn = np.zeros((4, 4, 4, 4), dtype=complex)
    
    p.chi_static = np.zeros((4, 4, 4, 4), dtype=complex)
    
    p.chi_tau = GfImTime(name=r'$g$', beta=m.beta,
                         statistic='Boson', n_points=50,
                         target_shape=(4, 4, 4, 4))

    p.chi = np.zeros_like(p.chi_static)
    
    tau = np.array([float(tau) for tau in p.chi_tau.mesh])

    for i1, i2, i3, i4, in itertools.product(range(4), repeat=4):

        print(i1, i2, i3, i4)
        
        o1, o2, o3, o4 = m.op_imp[i1], m.op_imp[i2], m.op_imp[i3], m.op_imp[i4]

        O1 = dagger(o1) * o2
        O2 = dagger(o3) * o4

        p.O1_exp[i1, i2] = ed.get_expectation_value(O1)
        p.O2_exp[i3, i4] = ed.get_expectation_value(O2)

        p.chi_dissconn[i1, i2, i3, i4] = p.O1_exp[i1, i2] * p.O2_exp[i3, i4]
        
        p.chi_static[i1, i2, i3, i4] = ed.get_expectation_value(O1 * O2)

        ed.set_g2_tau(p.chi_tau[i1, i2, i3, i4], O1, O2)

        chi_tau = p.chi_tau[i1, i2, i3, i4]
        chi_tau *= -1. # cancel gf -1 prefactor

        chi_dissconn = p.chi_dissconn[i1, i2, i3, i4]
        chi_tau.data[:] -= chi_dissconn        
        
        p.chi[i1, i2, i3, i4] = np.trapz(chi_tau.data, x=tau) / m.beta

    # ------------------------------------------------------------------
    # -- Store to hdf5

    if mpi.is_master_node():
        with HDFArchive('data_pyed.h5','w') as res:
            res['p'] = p
        
# ----------------------------------------------------------------------
