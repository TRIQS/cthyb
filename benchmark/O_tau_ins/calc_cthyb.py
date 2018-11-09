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
        beta = 10.,
        V1 = 2.0,
        V2 = 5.0,
        epsilon1 = 0.1,
        epsilon2 = 3.0,
        mu = 0.0,
        U = 5.0,
        #n_cycles = 1e8,
        n_cycles = 1e6,
        #n_cycles = 1e5,
        )

    p.init = ParameterCollection(
        beta = p.beta,
        gf_struct = [['up',[0]],['do',[0]]],
        n_iw = 1000,
        n_tau = 2*1000+1,
        n_l = 20,
        )

    ops = (n('up',0), n('do',0))
    
    p.solve = ParameterCollection(
        h_int = p.U*n('up',0)*n('do',0),
        max_time = -1,
        random_name = "",
        #random_seed = 123 * mpi.rank + 22222, # default uses mpi.rank
        random_seed = 123 * mpi.rank + 42211, # default uses mpi.rank
        #measure_G_l = True,
        measure_G_tau = True,
        move_double = True,
        # -- measurements
        length_cycle = 200,
        n_warmup_cycles = int(1e3),
        n_cycles = int(p.n_cycles / mpi.size),
        measure_O_tau = ops,
        # --
        measure_density_matrix = True,
        use_norm_as_weight = True,
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

    from pytriqs.atom_diag import trace_rho_op, AtomDiagReal

    o1, o2 = ops
    p.exp_val = trace_rho_op(solv.density_matrix, o1 * o2, solv.h_loc_diagonalization)
    p.n_exp = trace_rho_op(solv.density_matrix, n('up', 0) + n('do', 0), solv.h_loc_diagonalization)

    print 'p.exp_val =', p.exp_val
    
    # ------------------------------------------------------------------
    # -- Collect results

    attribs = [
        'G0_iw',
        'G_tau',
        #'G_l', 'G_iw',
        'O_tau',
        'density_matrix',
        'h_loc_diagonalization',
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
