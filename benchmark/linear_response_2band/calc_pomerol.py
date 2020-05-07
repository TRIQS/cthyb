
""" Test calculation for two-band Hubbard atom with two bath sites.

Use pomerol2triqs and pomerol to compute the one- and two-particle Green's 
function in the PH channel in Matsubara frequency space and zero 
bosonic transfer.

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

import triqs.utility.mpi as mpi
from h5 import HDFArchive

from pyed.ParameterCollection import ParameterCollection

from pomerol2triqs import PomerolED

# ----------------------------------------------------------------------
if __name__ == '__main__':

    if mpi.is_master_node():
        with HDFArchive('data_model.h5','r') as A: p = A["p"]
    else: p = None
    p = mpi.bcast(p)

    p.convert_keys_from_string_to_python('index_converter')
    
    pom = PomerolED(p.index_converter, verbose=True)
    pom.diagonalize(p.H)

    p.g_iw = pom.G_iw(p.gf_struct, p.beta, n_iw=400)
    p.g_tau = pom.G_tau(p.gf_struct, p.beta, n_tau=200)['0']
    p.tau = np.array([float(tau) for tau in p.g_tau.mesh])
    
    opt = dict(
        block_order='AABB',
        beta=p.beta, gf_struct=p.gf_struct,
        blocks=set([('0', '0')]),
        n_iw=1, n_inu=10)

    p.g4_ph = pom.G2_iw_inu_inup(channel='PH', **opt)['0', '0']
        
    if mpi.is_master_node():
        with HDFArchive('data_pomerol.h5','w') as A:
            A['p'] = p 
