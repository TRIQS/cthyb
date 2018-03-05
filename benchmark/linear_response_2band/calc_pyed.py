  
""" Test calculation for Hubbard atom with two bath sites.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

 """ 

# ----------------------------------------------------------------------

import copy
import itertools
import numpy as np

# ----------------------------------------------------------------------

from pytriqs.utility import mpi
from pytriqs.archive import HDFArchive
from pytriqs.gf import GfImTime, GfImFreq
from pytriqs.gf import MeshImTime, MeshProduct, Gf

# ----------------------------------------------------------------------

from pyed.ParameterCollection import ParameterCollection
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization
from pyed.GfUtils import g2_single_particle_transform
from pyed.GfUtils import g4_single_particle_transform

# ----------------------------------------------------------------------
if __name__ == '__main__':
    
    if mpi.is_master_node():
        with HDFArchive('data_model.h5','r') as A: p = A["p"]
    else: p = None
    p = mpi.bcast(p)

    p.T = np.matrix(p.T)
    
    ed = TriqsExactDiagonalization(p.H, p.op_full, p.beta)
    edt = TriqsExactDiagonalization(p.Ht, p.op_full, p.beta)

    if mpi.is_master_node():
        print p.T
        print p.T * p.T.H
        print 'Ground state energies:',
        print ed.get_ground_state_energy(),
        print edt.get_ground_state_energy()
        print 'max(abs(dE)) =', np.max(np.abs(ed.ed.E - edt.ed.E))

    N = ed.get_expectation_value(p.N_tot)
    Nt = edt.get_expectation_value(p.N_tot)

    print 'Total density:', N, Nt
        
    np.testing.assert_array_almost_equal(ed.ed.E, edt.ed.E, decimal=4)

    # ------------------------------------------------------------------
    # -- Single-particle Green's functions

    g_tau = GfImTime(name=r'$g$', beta=p.beta,
                     statistic='Fermion', n_points=15,
                     target_shape=(4, 4))

    gt_tau = g_tau.copy()

    ed.set_g2_tau_matrix(g_tau, p.op_imp)
    edt.set_g2_tau_matrix(gt_tau, p.op_imp)

    g_iwn = GfImFreq(name=r'$g$', beta=p.beta,
                     statistic='Fermion', n_points=15,
                     target_shape=(4, 4))
    
    gt_iwn = g_iwn.copy()

    ed.set_g2_iwn_matrix(g_iwn, p.op_imp)
    edt.set_g2_iwn_matrix(gt_iwn, p.op_imp)
    
    # ------------------------------------------------------------------
    # -- Transform single particle Green's function

    g_tau_ref = g2_single_particle_transform(gt_tau, p.T)
    np.testing.assert_array_almost_equal(g_tau_ref.data, g_tau.data)
    if mpi.is_master_node(): print '--> pyed g_tau and gt_tau agree'

    calc_g4 = False
    
    if calc_g4:
        # ------------------------------------------------------------------
        # -- Two particle Green's functions

        ntau = 4
        imtime = MeshImTime(p.beta, 'Fermion', ntau)
        prodmesh = MeshProduct(imtime, imtime, imtime)

        g4_tau = Gf(name='g4_tau', mesh=prodmesh, target_shape=[4, 4, 4, 4])
        g40_tau = g4_tau.copy()

        g4t_tau = g4_tau.copy()
        g40t_tau = g4_tau.copy()

        if mpi.is_master_node():
            print '--> g4 calc'

        ed.set_g40_tau_matrix(g40_tau, g_tau)
        edt.set_g40_tau_matrix(g40t_tau, gt_tau)

        ed.set_g4_tau_matrix(g4_tau, p.op_imp)
        edt.set_g4_tau_matrix(g4t_tau, p.op_imp)

        # ------------------------------------------------------------------
        # -- Transform two particle Green's function

        g40_tau_ref = g4_single_particle_transform(g40t_tau, p.T)
        np.testing.assert_array_almost_equal(g40_tau.data, g40_tau_ref.data)
        if mpi.is_master_node(): print '--> pyed g40_tau and g40t_tau agree'    

        g4_tau_ref = g4_single_particle_transform(g4t_tau, p.T)
        np.testing.assert_array_almost_equal(g4_tau.data, g4_tau_ref.data)
        if mpi.is_master_node(): print '--> pyed g4_tau and g4t_tau agree'    

    # ------------------------------------------------------------------
    # -- Store to hdf5

    if mpi.is_master_node():

        p = ParameterCollection()
        
        p.g_iwn = g_iwn
        p.g_tau = g_tau

        p.gt_iwn = gt_iwn
        p.gt_tau = gt_tau

        if calc_g4:
            p.g4_tau = g4_tau
            p.g40_tau = g40_tau

            p.g4t_tau = g4t_tau
            p.g40t_tau = g40t_tau
        
        with HDFArchive('data_pyed.h5','w') as res:
            res['p'] = p
        
# ----------------------------------------------------------------------
