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
def make_calc():

    p = ParameterCollection(
        beta = 1.0,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.0,
        epsilon2 = 4.0,
        mu = 2.0,
        U = 5.0,
        n_cycles = 1e6,
        )

    p.init = ParameterCollection(
        beta = p.beta,
        gf_struct = [['up',[0]],['do',[0]]],
        n_iw = 40,
        n_tau = 2*40+1,
        n_l = 20,
        )

    p.solve = ParameterCollection(
        h_int = p.U*n('up',0)*n('do',0),
        max_time = -1,
        random_name = "",
        random_seed = 123 * mpi.rank + 22222, # default uses mpi.rank
        measure_G_l = True,
        measure_G_tau = True,
        move_double = True,
        # -- measurements
        length_cycle = 20,
        n_warmup_cycles = int(1e3),
        n_cycles = int(p.n_cycles / mpi.size),
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
        'O_tau',
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

    make_calc()
