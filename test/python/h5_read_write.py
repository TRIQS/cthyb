import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf import GfImFreq, SemiCircular, inverse, iOmega_n
from pytriqs.operators import n, c, c_dag
from pytriqs.archive import HDFArchive
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
    
print type(sp['h_int'])
solver.solve(**sp)

sp = solver.last_solve_parameters
cp = solver.last_constr_parameters

filename = 'h5_read_write.h5'

with HDFArchive(filename, 'w') as A:
    A['solver'] = solver

with HDFArchive(filename, 'r') as A:
    solver_ref = A['solver']

cp_h5_ref = solver_ref.last_constr_parameters
sp_h5_ref = solver_ref.last_solve_parameters

assert( cp_h5_ref == cp )
assert( sp_h5_ref == sp )
