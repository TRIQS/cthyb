
""" Test calculation for two-band Hubbard atom with two bath sites.

Use cthyb to compute one- and two-particle Green's function in 
unitary transformed basis with off-diagonal hybridization.

Author: Hugo U.R. Strand (2018) hugo.strand@gmail.com """

# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------

import pytriqs.utility.mpi as mpi
from h5 import HDFArchive
from pytriqs.operators import n, c, c_dag, Operator, dagger

# ----------------------------------------------------------------------

from triqs_cthyb import SolverCore

# ----------------------------------------------------------------------

from pyed.OperatorUtils import relabel_operators
from pyed.ParameterCollection import ParameterCollection

# ----------------------------------------------------------------------
if __name__ == '__main__':

    if mpi.is_master_node():
        with HDFArchive('data_model.h5','r') as A: p = A["p"]
    else: p = None
    p = mpi.bcast(p)
    
    S = SolverCore(beta=p.beta, gf_struct=p.gf_struct, n_tau=p.ntau, n_iw=p.nw)

    S.G0_iw['0'] << p.g0t_iw

    solve_parameters = dict(
        h_int = p.Ht_int,
        max_time = -1,
        random_name = "",
        random_seed = 123 * mpi.rank + 567,
        length_cycle = 100,
        n_warmup_cycles = int(1e4),
        #n_cycles = int(1e9) / mpi.size,
        n_cycles = int(1e7) / mpi.size,
        move_double = True,
        measure_G2_iw_ph = True,
        measure_G2_n_fermionic = 10,
        measure_G2_n_bosonic = 1,
        )

    S.solve(**solve_parameters)

    if mpi.is_master_node():
        p.g_tau = S.G_tau
        p.g0_iw = S.G0_iw
        p.g4_ph = S.G2_iw_ph[('0', '0')]
        with HDFArchive('data_cthyb.h5','w') as A:
            A['p'] = p 
