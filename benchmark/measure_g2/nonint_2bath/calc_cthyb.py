# ----------------------------------------------------------------------    

""" CTHYB test calculation 
    for Hubbard atom with two bath sites.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

"""

# ----------------------------------------------------------------------    

import time
import numpy as np

# ----------------------------------------------------------------------    

import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.operators import *
from h5 import HDFArchive
from triqs_cthyb import Solver

# ----------------------------------------------------------------------    
if __name__ == '__main__':

    params = dict(
        beta = 2.0,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.0,
        epsilon2 = 4.0,
        mu = 2.0,
        U = 0.0,
        )
    
    # ------------------------------------------------------------------

    class Dummy():
        def __init__(self):
            pass
        
    d = Dummy() # storage space
    d.params = params

    if mpi.is_master_node():
        print('--> Solving SIAM with parameters')
        for key, value in list(params.items()):
            print('%10s = %-10s' % (key, str(value)))
            globals()[key] = value # populate global namespace

    # ------------------------------------------------------------------

    solv = Solver(
        beta=beta,
        gf_struct=[['up', [0]], ['do', [0]]],
        n_iw=15,
        n_tau=4*128+1,
        n_l=20,
        )

    # Local Hamiltonian
    H = U*n("up",0)*n("do",0)

    # -- and update the Weiss field of the impurity
    for name, g0 in solv.G0_iw:
        g0 << inverse(iOmega_n + mu
                      - V1**2*inverse(iOmega_n - epsilon1)
                      - V2**2*inverse(iOmega_n - epsilon2)
                     )

    starttime = time.time()
    
    # ------------------------------------------------------------------
    # -- Solve the impurity model
    solv.solve(
        h_int=H,
        max_time=-1,
        random_name="",
        #random_seed=123 * mpi.rank + 567, # default uses mpi.rank
        measure_G_l=True,
        move_double=True,
        measure_pert_order=True,
        use_norm_as_weight=True, # needed for density matrix
        measure_density_matrix=True,
        # -- measurements
        length_cycle=1000,
        n_warmup_cycles=int(1e3),
        n_cycles=int(4e3),
        # -- g2 measurements
        measure_G2_tau=True,
        measure_G2_iw=True,
        measure_G2_iw_pp=True,
        measure_G2_iw_ph=True,
        measure_G2_n_tau=40,
        measure_G2_n_bosonic=15,
        measure_G2_n_fermionic=15,
        #nfft_buf_sizes=dict(up=64, do=64),
        )

    runtime = time.time() - starttime
    if mpi.is_master_node():
        print('--> runtime:', runtime)

    G_iw = solv.G0_iw.copy()
    G_tau = solv.G_tau.copy()
    for name, g0 in solv.G0_iw:
        G_iw[name] << LegendreToMatsubara(solv.G_l[name])
        G_tau[name] << Fourier(G_iw[name])

    # ------------------------------------------------------------------
    # -- Collect results
    
    d.Sigma_iw = solv.Sigma_iw['up']
    d.G_tau = solv.G_tau['up']
    d.G_l = solv.G_l['up']
    d.G0_iw = solv.G0_iw['up']

    d.G_iw = G_iw['up']
    d.Gl_tau = G_tau['up']

    d.runtime = runtime
    d.G2_tau = solv.G2_tau[('up', 'do')]
    d.G2_iw = solv.G2_iw[('up', 'do')]
    d.G2_iw_pp = solv.G2_iw_pp[('up', 'do')]
    d.G2_iw_ph = solv.G2_iw_ph[('up', 'do')]

    d.mpi_size = mpi.size

    d.perturbation_order = solv.perturbation_order
    d.perturbation_order_total = solv.perturbation_order_total
        
    # ------------------------------------------------------------------
    # -- Store results
    if mpi.is_master_node():
        filename = 'data_cthyb.h5'
        with HDFArchive(filename,'w') as res:
            for key, value in d.__dict__.items():
                res[key] = value



