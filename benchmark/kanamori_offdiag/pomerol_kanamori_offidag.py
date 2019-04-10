
import itertools
import numpy as np
from scipy.linalg import block_diag

from pytriqs.archive import HDFArchive
from pytriqs.operators import c, c_dag, n

from pytriqs.gf import Gf
from pytriqs.gf import MeshImTime, MeshImFreq
from pytriqs.operators.util.hamiltonians import h_int_kanamori

from pyed.OperatorUtils import get_quadratic_operator
from pyed.OperatorUtils import fundamental_operators_from_gf_struct

from pytriqs.utility import mpi
from pomerol2triqs import PomerolED

# ----------------------------------------------------------------------
if __name__ == '__main__':

    norb = 2
    spin_names = ['up', 'do']
    imp_idxs = list(range(0, norb))
    bath_idxs = list(range(norb, 3*norb))
    
    fop_imp = fundamental_operators_from_gf_struct(
        [ ['up', imp_idxs], ['do', imp_idxs] ])

    fop_bath = fundamental_operators_from_gf_struct(
        [ ['up', bath_idxs], ['do', bath_idxs] ])

    print spin_names
    print imp_idxs
    print bath_idxs
    print fop_imp
    print fop_bath
    
    fop = fop_imp + fop_bath

    beta = 10.0
    mu = 1.0
    U = 2.0
    J = 0.2

    V = 1.0*np.eye(norb) + 0.1*(np.ones((norb, norb)) - np.eye(norb))
    V = np.hstack([V, V])
    
    T_imp = -mu * np.eye(norb)
    
    T_bath = np.diag([-2.3, -2.3, 2.3, 2.3])
 
    print 'V =\n', V
    print 'T_imp =\n', T_imp
    print 'T_bath =\n', T_bath

    V = block_diag(V, V)
    T_imp = block_diag(T_imp, T_imp)
    T_bath = block_diag(T_bath, T_bath)

    print 'V =\n', V
    print 'T_imp =\n', T_imp
    print 'T_bath =\n', T_bath
    
    T = np.block([
        [ T_imp, V ],
        [ V.T, T_bath]
        ])
    
    print 'T =\n', T

    H_int = h_int_kanamori(
        spin_names, imp_idxs,
        np.array([[0,U-3*J],[U-3*J,0]]),
        np.array([[U,U-2*J],[U-2*J,U]]),
        J,True)
    
    H = H_int + get_quadratic_operator(T, fop)

    up, do = spin_names
    index_converter = {
        (up, 0) : ('loc', 0, 'up'),
        (do, 0) : ('loc', 0, 'down'),
        (up, 1) : ('loc', 1, 'up'),
        (do, 1) : ('loc', 1, 'down'),
        (up, 2) : ('loc', 2, 'up'),
        (do, 2) : ('loc', 2, 'down'),
        (up, 3) : ('loc', 3, 'up'),
        (do, 3) : ('loc', 3, 'down'),
        (up, 4) : ('loc', 4, 'up'),
        (do, 4) : ('loc', 4, 'down'),
        (up, 5) : ('loc', 5, 'up'),
        (do, 5) : ('loc', 5, 'down'),
        }
    
    ed = PomerolED(index_converter, verbose=True)
    ed.diagonalize(H)

    gf_struct_up = [[up, [0, 1]]]
    gf_struct_do = [[do, [0, 1]]]

    G_tau_up = ed.G_tau(gf_struct_up, beta, n_tau=100)
    G_tau_do = ed.G_tau(gf_struct_do, beta, n_tau=100)

    with HDFArchive('kanamori_offdiag.pomerol.h5', 'w') as Results:
        Results['up'] = G_tau_up
        Results['dn'] = G_tau_do
