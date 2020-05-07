import numpy as np
import triqs.utility.mpi as mpi
from triqs.gf import GfImFreq, BlockGf, SemiCircular, inverse, iOmega_n
from triqs.operators import n, c, c_dag
from h5 import HDFArchive
from triqs_cthyb import SolverCore

cp = dict(
    beta=10.0,
    gf_struct=[['up',[0]], ['do',[0]]],
    n_iw=1025, n_tau=2500, n_l=20
    )

solver = SolverCore(**cp)

# Set hybridization function
mu = 0.5
half_bandwidth = 1.0
delta_w = GfImFreq(indices = [0], beta=cp['beta'])
delta_w << (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)
for name, g0 in solver.G0_iw:
    g0 << inverse(iOmega_n + mu - delta_w)

sp = dict(
    h_int = n('up',0)*n('do',0) + c_dag('up',0)*c('do',0) + c_dag('do',0)*c('up',0),
    max_time = -1,
    length_cycle = 50,
    n_warmup_cycles = 50,
    n_cycles = 500,
    move_double = False,
    )
    
solver.solve(**sp)

filename = 'h5_read_write.h5'

with HDFArchive(filename, 'w') as A:
    A['solver'] = solver

with HDFArchive(filename, 'r') as A:
    solver_ref = A['solver']

assert( solver.last_constr_parameters == solver_ref.last_constr_parameters )
assert( solver.last_solve_parameters == solver.last_solve_parameters )
#assert( solver.last_container_set == solver_ref.last_container_set ) # want to write this

# -- Poor mans version of comparison of the container sets
for key in dir(solver_ref):
    if 'G' in key or 'Delta' in key:
        print('comparing', key)
        
        val = getattr(solver, key)
        val_ref = getattr(solver_ref, key)

        if isinstance(val, BlockGf):
            for (n1, g1), (n1, g2) in zip(val, val_ref):
                np.testing.assert_array_almost_equal(g1.data, g2.data)
        elif val == None or isinstance(val, list):
            assert(val == val_ref)
        else:
            raise Exception("Invalid type in comparison")
