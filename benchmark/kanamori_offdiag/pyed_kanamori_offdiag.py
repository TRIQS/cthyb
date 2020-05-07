
import itertools
import numpy as np
from scipy.linalg import block_diag

from h5 import HDFArchive
from triqs.operators import c, c_dag, n

from triqs.gf import Gf
from triqs.gf import MeshImTime, MeshImFreq
from triqs.operators.util.hamiltonians import h_int_kanamori

from pyed.OperatorUtils import get_quadratic_operator
from pyed.OperatorUtils import fundamental_operators_from_gf_struct
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

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

    print(spin_names)
    print(imp_idxs)
    print(bath_idxs)
    print(fop_imp)
    print(fop_bath)
    
    fop = fop_imp + fop_bath

    beta = 10.0
    mu = 1.0
    U = 2.0
    J = 0.2

    V = 1.0*np.eye(norb) + 0.1*(np.ones((norb, norb)) - np.eye(norb))
    V = np.hstack([V, V])
    
    T_imp = -mu * np.eye(norb)
    
    T_bath = np.diag([-2.3, -2.3, 2.3, 2.3])
 
    print('V =\n', V)
    print('T_imp =\n', T_imp)
    print('T_bath =\n', T_bath)

    V = block_diag(V, V)
    T_imp = block_diag(T_imp, T_imp)
    T_bath = block_diag(T_bath, T_bath)

    print('V =\n', V)
    print('T_imp =\n', T_imp)
    print('T_bath =\n', T_bath)
    
    T = np.block([
        [ T_imp, V ],
        [ V.T, T_bath]
        ])
    
    print('T =\n', T)

    H_int = h_int_kanamori(
        spin_names, imp_idxs,
        np.array([[0,U-3*J],[U-3*J,0]]),
        np.array([[U,U-2*J],[U-2*J,U]]),
        J,True)
    
    H = H_int + get_quadratic_operator(T, fop)
            
    ed = TriqsExactDiagonalization(H, fop, beta)

    n_tau = 101    
    G_tau_up = Gf(mesh=MeshImTime(beta, 'Fermion', n_tau), target_shape=[2, 2])
    G_tau_do = Gf(mesh=MeshImTime(beta, 'Fermion', n_tau), target_shape=[2, 2])

    for i1, i2 in itertools.product(range(norb), repeat=2):
        ed.set_g2_tau(G_tau_up[i1, i2], c('up', i1), c_dag('up', i2))
        ed.set_g2_tau(G_tau_do[i1, i2], c('do', i1), c_dag('do', i2))

    with HDFArchive('kanamori_offdiag.pyed.h5', 'w') as Results:
        Results['up'] = G_tau_up
        Results['dn'] = G_tau_do
