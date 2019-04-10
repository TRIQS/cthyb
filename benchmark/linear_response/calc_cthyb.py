# ----------------------------------------------------------------------    

""" CTHYB test calculation 
    for Hubbard atom with two bath sites.

Author: Hugo U.R. Strand (2017) hugo.strand@gmail.com

"""

# ----------------------------------------------------------------------    

import os
import itertools
import numpy as np

# ----------------------------------------------------------------------    

import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive

# ----------------------------------------------------------------------    

from triqs_cthyb import Solver

# ----------------------------------------------------------------------    

from pyed.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------    
def make_calc(nw=2, nc=1e5, beta=2.0, h_field=0.0, rand=1):

    p = ParameterCollection(
        beta = beta,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.0,
        epsilon2 = 4.0,
        mu = 2.0,
        U = 5.0,
        h_field = h_field,
        n_cycles = nc,
        )

    p.init = ParameterCollection(
        beta = p.beta,
        gf_struct = [['up',[0]],['do',[0]]],
        n_iw = 1000,
        n_tau = 6*1000+1,
        n_l = 20,
        )

    p.solve = ParameterCollection(
        h_int = p.U*n('up',0)*n('do',0) - p.h_field*(n('up',0) - n('do',0)),
        max_time = -1,
        random_name = "",
        random_seed = 123 * mpi.rank + rand, # default uses mpi.rank
        measure_G_l = True,
        measure_G_tau = True,
        move_double = True,
        # -- measurements
        length_cycle = 50,
        n_warmup_cycles = int(5e4),
        n_cycles = int(p.n_cycles / mpi.size),
        # -- G2 measurements
        measure_G2_iw_ph = True,
        measure_G2_n_bosonic = 1,
        measure_G2_n_fermionic = nw,
        )

    print 'p.solve.random_seed =', p.solve.random_seed
    
    # ------------------------------------------------------------------
    
    if mpi.is_master_node():
        print '--> Solving SIAM with parameters'
        print p

    # ------------------------------------------------------------------

    solv = Solver(**p.init.dict())

    # -- and update the Weiss field of the impurity
    for name, g0 in solv.G0_iw:
        g0 << inverse(iOmega_n + p.mu
                      - p.V1**2*inverse(iOmega_n - p.epsilon1)
                      - p.V2**2*inverse(iOmega_n - p.epsilon2)
                     )

    # ------------------------------------------------------------------
    # -- Solve the impurity model
    solv.solve(**p.solve.dict())

    p.Gl_iw = solv.G0_iw.copy()
    p.Gl_tau = solv.G_tau.copy()
    for name, g0 in solv.G0_iw:
        p.Gl_iw[name] << LegendreToMatsubara(solv.G_l[name])
        p.Gl_tau[name] << InverseFourier(p.Gl_iw[name])

    # ------------------------------------------------------------------
    # -- Collect results

    attribs = [
        'G0_iw',
        'G_tau', 'G_l', 'G_iw',
        'G2_iw_ph'
        ]
    
    for key in attribs:
        val = getattr(solv, key)
        setattr(p, key, val)

    # ------------------------------------------------------------------
    # -- Store results
    if mpi.is_master_node():
        filename = 'data_cthyb.h5'
        with HDFArchive(filename, 'w') as res:
            res['p'] = p

# ----------------------------------------------------------------------
if __name__ == '__main__':

    #nw_vec = np.array([2, 4, 8])
    #beta_vec = 1./np.array([10., 20., 40., 60., 80., 100.])

    repeats = 10
    nw_vec = np.array([8])
    nc_vec = np.array([5e7])
    beta_vec = 1./np.array([10.])

    for nw, nc, beta in itertools.product(nw_vec, nc_vec, beta_vec):

        for repeat in xrange(repeats):

            if mpi.is_master_node():
                max = 10000
                rand = np.random.randint(10000, high=10*max - 1)
            else:
                rand = None

            rand = mpi.bcast(rand)

            path = 'cthyb_nw%i_nc%i_beta%6.6f_rand%i' % (nw, int(np.log10(nc)), beta, rand)
            print '-->path: ', path
            if mpi.is_master_node():
                os.mkdir(path)

            mpi.barrier()
            os.chdir(path)

            make_calc(nw=nw, nc=nc, beta=beta, h_field=0.0, rand=rand)

            os.chdir('../')
