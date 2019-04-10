  
""" Test calculation for two-band Hubbard atom with two bath sites.

A unitary transformed impurity model is also constructed with
off diagonal hybridization terms, inorder to be able to sample
all components of the two-particle Green's function using cthyb

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """ 

# ----------------------------------------------------------------------

import copy
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.utility import mpi
from pytriqs.archive import HDFArchive

from pytriqs.gf import Gf, inverse, iOmega_n, InverseFourier
from pytriqs.gf import MeshImTime, MeshProduct
from pytriqs.gf import GfImTime, GfImFreq

from pytriqs.operators import n, c, c_dag, Operator, dagger

from pytriqs.operators.util.op_struct import set_operator_structure
from pytriqs.operators.util.U_matrix import U_matrix_kanamori, U_matrix
from pytriqs.operators.util.hamiltonians import h_int_kanamori

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection
from pyed.OperatorUtils import relabel_operators
from pyed.OperatorUtils import op_is_fundamental, get_quadratic_operator, \
    quadratic_matrix_from_operator, operator_single_particle_transform

from pyed.GfUtils import g2_single_particle_transform

# ----------------------------------------------------------------------
def unitary_transf_4x4(t_vec):

    """ Not completely general generation of real valued 4x4
    unitary transformation. """

    assert( len(t_vec) == 6 )
    
    I = np.eye(2)
    a1 = np.array([[0, 1], [-1, 0]])
    b1 = np.array([[0, 1], [1, 0]])
    sz = np.array([[1, 0], [0, -1]])

    # Even modes
    A1 = np.kron(I, a1)
    A2 = np.kron(a1, I)
    A3 = np.kron(a1, b1)
    # Odd modes
    A4 = np.kron(b1, a1) 
    A5 = np.kron(sz, a1)
    A6 = np.kron(a1, sz)
        
    As = [A1, A2, A3, A4, A5, A6]
    
    M = sum([t * A for t, A in zip(t_vec, As)])
    from scipy.linalg import expm as scipy_expm
    T = np.mat(scipy_expm(M))

    np.testing.assert_array_almost_equal(T * T.H, np.eye(4))

    return T

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    # ------------------------------------------------------------------
    # -- Hubbard atom with two bath sites, Hamiltonian

    p = ParameterCollection(
        U = 1.0,
        J = 0.2,
        beta = 1.0,
        crystal_field = 1.0,
        num_orbitals = 2,
        num_spins = 2,
        nw = 400,
        ntau = 801,
        )
    
    up, do = '0', '1'
    p.spin_names = [up, do]
    p.orb_names = range(p.num_orbitals)

    # -- Unitary transform

    if False:
        from scipy.stats import unitary_group, ortho_group
        np.random.seed(seed=233423) # -- Reproducible "randomness"
        p.T = unitary_group.rvs(4) # General complex unitary transf
        #p.T = np.kron(ortho_group.rvs(2), np.eye(2)) # orbital only rotation

    # use parametrized "weak" realvalued 4x4 transform
    p.T = unitary_transf_4x4(t_vec=[0.3, 0.3, 0.3, 0.0, 0.0, 0.0])
    #p.T = unitary_transf_4x4(t_vec=[0.25*np.pi, 0.0, 0.0, 0.0, 0.0, 0.0])
    #p.T = np.eye(4)
    p.T = np.matrix(p.T)
    print p.T
    
    # -- Different types of operator sets
    
    p.op_imp = [c('0',0), c('0',1), c('0',2), c('0',3)]
    p.op_full = p.op_imp + [c('1',0), c('1',1), c('1',2), c('1',3)]

    p.spin_block_ops = [c('0', 0), c('1', 0), c('0', 1), c('1', 1)]
    p.spin_block_ops += [ dagger(op) for op in p.spin_block_ops ]

    p.org_ops = copy.copy(p.op_imp)
    p.org_ops += [ dagger(op) for op in p.org_ops ]
    
    p.diag_ops = [c('0_0',0), c('0_1',0), c('0_2',0), c('0_3',0)]
    p.diag_ops += [ dagger(op) for op in p.diag_ops ]
    
    p.gf_struct = [['0', [0, 1, 2, 3]]]
    p.gf_struct_diag = [['0_0',[0]], ['0_1',[0]], ['0_2',[0]], ['0_3',[0]]]

    p.index_converter = {
        ('0', 0) : ('loc', 0, 'up'),
        ('0', 1) : ('loc', 0, 'down'),
        ('0', 2) : ('loc', 1, 'up'),
        ('0', 3) : ('loc', 1, 'down'),
        ('1', 0) : ('bath', 0, 'up'),
        ('1', 1) : ('bath', 0, 'down'),
        ('1', 2) : ('bath', 1, 'up'),
        ('1', 3) : ('bath', 1, 'down'),
        }
    
    # ------------------------------------------------------------------
    # -- Setup Hamiltonian

    # -- Local crystal field
    
    h_loc = np.array([
        [1.0, 0.0],
        [0.0, -1.0],
        ])

    I = np.eye(p.num_spins)
    h_loc = np.kron(h_loc, I)

    p.H_loc = p.crystal_field * get_quadratic_operator(h_loc, p.op_imp)

    # -- Bath sites and hybridizations
    
    V2 = 1.0
    V3 = 3.0
    E2 = -1.0
    E3 = -4.0
    
    h_bath = np.array([
        [0,  0, V2, 0 ],
        [0,  0, 0,  V3],
        [V2, 0, E2, 0 ],
        [0,  V3, 0,  E3],
        ])

    h_bath = np.kron(h_bath, I)
    p.H_bath = get_quadratic_operator(h_bath, p.op_full)

    # -- Weiss field of the bath
    
    h_tot = quadratic_matrix_from_operator(p.H_bath + p.H_loc, p.op_full)

    g0_iw = GfImFreq(beta=p.beta, statistic='Fermion',
                      n_points=p.nw, target_shape=(8, 8))

    g0_iw << inverse(iOmega_n - h_tot)

    p.g0_iw = g0_iw[:4, :4] # -- Cut out impurity Gf
    p.g0t_iw = g2_single_particle_transform(p.g0_iw, p.T.H)

    p.g0_tau = GfImTime(beta=p.beta, statistic='Fermion',
                        n_points=p.ntau, target_shape=(4, 4))

    p.g0_tau << InverseFourier(p.g0_iw)
    p.g0t_tau = g2_single_particle_transform(p.g0_tau, p.T.H)

    p.g0_tau_ref = g2_single_particle_transform(p.g0t_tau, p.T)
    np.testing.assert_array_almost_equal(p.g0_tau_ref.data, p.g0_tau.data)
    
    # -- Interaction Hamiltonian: Kanamori interaction
    
    U_ab, UPrime_ab = U_matrix_kanamori(n_orb=p.num_orbitals, U_int=p.U, J_hund=p.J)
    
    H_int = h_int_kanamori(
        p.spin_names, p.orb_names, U_ab, UPrime_ab, J_hund=p.J,
        off_diag=True, map_operator_structure=None, H_dump=None)

    p.H_int = relabel_operators(H_int, p.spin_block_ops, p.org_ops)

    # -- Complete 2site + 2bath Hamiltonian
    
    p.H = p.H_loc + p.H_bath + p.H_int

    # -- Single particle transformed Hamiltonian
    
    p.Ht = operator_single_particle_transform(p.H, p.T, p.op_imp)
    p.Ht_int = operator_single_particle_transform(p.H_int, p.T, p.op_imp)

    # ------------------------------------------------------------------
    # -- Store to hdf5

    if mpi.is_master_node():        
        with HDFArchive('data_model.h5','w') as res:
            res['p'] = p
        
# ----------------------------------------------------------------------
