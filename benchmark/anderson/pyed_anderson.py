import numpy as np

from pytriqs.archive import HDFArchive
from pytriqs.operators import c, c_dag, n

from pytriqs.gf import Gf
from pytriqs.gf import MeshImTime, MeshImFreq

from pyed.OperatorUtils import fundamental_operators_from_gf_struct
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# ----------------------------------------------------------------------
if __name__ == '__main__':

    orb_idxs = list(range(3))
    gf_struct = [ ['up', orb_idxs], ['do', orb_idxs] ]
    
    fundamental_operators = fundamental_operators_from_gf_struct(gf_struct)

    beta = 10.0
    U = 2.0
    mu = 1.0
    h = 0.1
    #V = 0.5
    V = 1.0
    epsilon = 2.3
    
    H = U * n('up', 0) * n('do', 0) + \
        - mu * (n('up', 0) + n('do', 0)) + \
        h * n('up', 0) - h * n('do', 0) + \
        epsilon * (n('up', 1) + n('do', 1)) - epsilon * (n('up', 2) + n('do', 2)) + \
        V * ( c_dag('up', 0) * c('up', 1) + c_dag('up', 1) * c('up', 0) ) + \
        V * ( c_dag('do', 0) * c('do', 1) + c_dag('do', 1) * c('do', 0) ) + \
        V * ( c_dag('up', 0) * c('up', 2) + c_dag('up', 2) * c('up', 0) ) + \
        V * ( c_dag('do', 0) * c('do', 2) + c_dag('do', 2) * c('do', 0) )
    
    ed = TriqsExactDiagonalization(H, fundamental_operators, beta)

    n_tau = 101    
    G_tau_up = Gf(mesh=MeshImTime(beta, 'Fermion', n_tau), target_shape=[])
    G_tau_do = Gf(mesh=MeshImTime(beta, 'Fermion', n_tau), target_shape=[])

    ed.set_g2_tau(G_tau_up, c('up',0), c_dag('up',0))
    ed.set_g2_tau(G_tau_do, c('do',0), c_dag('do',0))

    with HDFArchive('anderson.pyed.h5', 'w') as Results:
        Results['up'] = G_tau_up
        Results['dn'] = G_tau_do
