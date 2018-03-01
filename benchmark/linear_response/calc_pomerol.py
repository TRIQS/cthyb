# ----------------------------------------------------------------------

""" Pomerol exact diagonalization test calculation 
    for a Hubbard atom with two bath sites. 

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

"""

# ----------------------------------------------------------------------

import os
import itertools
import numpy as np

from pytriqs.operators import c, c_dag
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------

from pytriqs.utility import mpi

# ----------------------------------------------------------------------

from pomerol2triqs import PomerolED
from pyed.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
def make_calc(nw=10, beta=2.0, h_field=0.0):
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    p = ParameterCollection(
        beta = beta,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.00,
        epsilon2 = 4.00,
        h_field = h_field,
        mu = 2.0,
        U = 5.0,
        ntau = 40,
        niw = 15,
        n_inu = nw,
        )

    # ------------------------------------------------------------------

    print '--> Solving SIAM with parameters'
    print p
    
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
    # -- Exact diagonalization

    # Conversion from TRIQS to Pomerol notation for operator indices
    # TRIQS:   block_name, inner_index
    # Pomerol: site_label, orbital_index, spin_name
    index_converter = {
        (up, 0) : ('loc', 0, 'up'),
        (do, 0) : ('loc', 0, 'down'),
        (up, 1) : ('loc', 1, 'up'),
        (do, 1) : ('loc', 1, 'down'),
        (up, 2) : ('loc', 2, 'up'),
        (do, 2) : ('loc', 2, 'down'),
        }

    # -- Create Exact Diagonalization instance
    ed = PomerolED(index_converter, verbose=True)
    ed.diagonalize(p.H) # -- Diagonalize H

    gf_struct = {up : [0], do : [0]}

    # -- Single-particle Green's functions
    p.G_iw = ed.G_iw(gf_struct, beta, n_iw=100)
    p.G_tau = ed.G_tau(gf_struct, beta, n_tau=400)

    # -- Particle-particle two-particle Matsubara frequency Green's function
    opt = dict(
        block_order='AABB',
        beta=beta, gf_struct=gf_struct,
        blocks=set([(up, up), (up, do)]),
        n_iw=1, n_inu=p.n_inu)
    
    p.G2_iw_ph = ed.G2_iw_inu_inup(channel='PH', **opt)

    # ------------------------------------------------------------------
    # -- Store to hdf5
    
    mpi.barrier()
    if mpi.is_master_node():
        filename = 'data_pomerol_h_field_%4.4f.h5' % h_field
        with HDFArchive(filename,'w') as res:
            res['p'] = p
            
# ----------------------------------------------------------------------
if __name__ == '__main__':

    nw_vec = np.array([2, 4, 8])
    beta_vec = 1./np.array([10., 20., 40., 60., 80., 100.])

    for nw, beta in itertools.product(nw_vec, beta_vec):
        path = 'pomerol_nw%i_beta%6.6f' % (nw, beta)
        print '-->path: ', path
        if mpi.is_master_node():
            os.mkdir(path)

        mpi.barrier()
        os.chdir(path)
        
        make_calc(nw=nw, beta=beta, h_field=0.0)

        os.chdir('../')
    
