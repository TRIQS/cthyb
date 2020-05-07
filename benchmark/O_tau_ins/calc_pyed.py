# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation 
    for a Hubbard atom with two bath sites. 

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

"""

# ----------------------------------------------------------------------

import os
import numpy as np

from triqs.operators import c, c_dag, n
from h5 import HDFArchive
from triqs.gf import GfImTime, GfImFreq, BlockGf, Gf, MeshImTime

# ----------------------------------------------------------------------

from triqs.utility import mpi

# ----------------------------------------------------------------------

from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization
from pyed.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
def make_calc(beta=2.0, h_field=0.0):
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    p = ParameterCollection(
        beta = beta,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.10,
        epsilon2 = 3.00,
        h_field = h_field,
        mu = 0.0,
        U = 5.0,
        ntau = 800,
        niw = 15,
        )

    # ------------------------------------------------------------------

    print('--> Solving SIAM with parameters')
    print(p)
    
    # ------------------------------------------------------------------

    up, do = 'up', 'dn'
    docc = c_dag(up,0) * c(up,0) * c_dag(do,0) * c(do,0)
    mA = c_dag(up,0) * c(up,0) - c_dag(do,0) * c(do,0)

    nA = c_dag(up,0) * c(up,0) + c_dag(do,0) * c(do,0)
    nB = c_dag(up,1) * c(up,1) + c_dag(do,1) * c(do,1)
    nC = c_dag(up,2) * c(up,2) + c_dag(do,2) * c(do,2)

    p.H = -p.mu * nA + p.U * docc + p.h_field * mA + \
        p.epsilon1 * nB + p.epsilon2 * nC + \
        p.V1 * (c_dag(up,0)*c(up,1) + c_dag(up,1)*c(up,0) + \
              c_dag(do,0)*c(do,1) + c_dag(do,1)*c(do,0) ) + \
        p.V2 * (c_dag(up,0)*c(up,2) + c_dag(up,2)*c(up,0) + \
              c_dag(do,0)*c(do,2) + c_dag(do,2)*c(do,0) )

    # ------------------------------------------------------------------

    fundamental_operators = [
        c(up,0), c(do,0), c(up,1), c(do,1), c(up,2), c(do,2)]
    
    ed = TriqsExactDiagonalization(p.H, fundamental_operators, p.beta)

    g_tau = GfImTime(beta=beta, statistic='Fermion', n_points=400, indices=[0])
    g_iw = GfImFreq(beta=beta, statistic='Fermion', n_points=10, indices=[0])

    p.G_tau = BlockGf(name_list=[up,do], block_list=[g_tau]*2, make_copies=True)
    p.G_iw = BlockGf(name_list=[up,do], block_list=[g_iw]*2, make_copies=True)
    
    ed.set_g2_tau(p.G_tau[up][0, 0], c(up,0), c_dag(up,0))
    ed.set_g2_tau(p.G_tau[do][0, 0], c(do,0), c_dag(do,0))

    ed.set_g2_iwn(p.G_iw[up][0, 0], c(up,0), c_dag(up,0))
    ed.set_g2_iwn(p.G_iw[do][0, 0], c(do,0), c_dag(do,0))

    p.magnetization = ed.get_expectation_value(0.5 * mA)

    p.O_tau = Gf(mesh=MeshImTime(beta, 'Fermion', 400), target_shape=[])
    ed.set_g2_tau(p.O_tau, n(up,0), n(do,0))
    p.O_tau.data[:] *= -1.

    p.exp_val = ed.get_expectation_value(n(up,0) * n(do,0))
    
    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    filename = 'data_pyed_h_field_%4.4f.h5' % h_field
    with HDFArchive(filename,'w') as res:
        res['p'] = p
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    make_calc(beta=10., h_field=0.0)
