""" 
Sampling of the density density correlator by operator insertion
regression test derived from the benchmark ./benchmark/O_tau_ins/

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------    

import numpy as np

# ----------------------------------------------------------------------    

from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive

from pytriqs.utility.h5diff import h5diff
import pytriqs.utility.mpi as mpi

# ----------------------------------------------------------------------    

from triqs_cthyb import Solver

# ----------------------------------------------------------------------
if __name__ == '__main__':

    solv = Solver(
        beta = 2.1,
        gf_struct = [['up',[0]],['do',[0]]],
        n_iw = 30,
        n_tau = 2*30+1,
        )

    # -- Weiss field of the impurity
    
    V1 = 2.0
    V2 = 5.0
    epsilon1 = 0.0
    epsilon2 = 4.0
    mu = 2.0
    
    for name, g0 in solv.G0_iw:
        g0 << inverse(iOmega_n + mu
                      - V1**2*inverse(iOmega_n - epsilon1)
                      - V2**2*inverse(iOmega_n - epsilon2)
                     )

    # -- Solve the impurity model
    
    solv.solve(
        h_int = 5.0*n('up',0)*n('do',0),
        measure_G_tau = True,
        move_double = True,
        # -- measurements
        length_cycle = 20,
        n_warmup_cycles = int(1e4),
        n_cycles = int(1e4),
        # -- measure density-density correlator
        measure_O_tau = (n('up',0), n('do',0)),
        )

    # -- Store results
    
    filename = 'O_tau_ins.out.h5'
    with HDFArchive(filename, 'w') as res:
        res['O_tau'] = solv.O_tau

    h5diff(filename, 'O_tau_ins.ref.h5')
