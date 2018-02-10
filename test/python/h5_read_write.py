import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.operators import n, Operator
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from pytriqs.operators.util.op_struct import set_operator_structure
from pytriqs.archive import HDFArchive
from cthyb import SolverCore
from pytriqs.utility.comparison_tests import *

def cf_elements_in_A_with_B(a, b, verbose=False):
    """ Compare keys in a with b """
    if verbose:
        print 'a.keys() =', a.keys()
        print 'b.keys() =', b.keys()
        print 'a =', a
        print 'b =', b
    assert( set(a.keys()) <= set(b.keys()) )
    for key in a.keys():
        if verbose: print key,
        va, vb = a[str(key)], b[str(key)]
        if verbose: print va, vb, type(va), type(vb)
        if type(va) == Operator:
            if verbose: print 'va - vb = ', va - vb
            if not (va - vb).is_zero():
                print key, va, vb
                return False
        else:
            if va != vb:
                print key, va, vb
                return False
    return True

cp_in = dict(
    beta=10.0,
    gf_struct={'up' : [0], 'do' : [0]},
    n_iw=1025, n_tau=2500, n_l=20
    )

solver = SolverCore(**cp_in)

cp = solver.last_constr_parameters

assert( cf_elements_in_A_with_B(cp_in, cp) )

# Set hybridization function
mu = 0.5
half_bandwidth = 1.0
delta_w = GfImFreq(indices = [0], beta=cp['beta'])
delta_w << (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)
for name, g0 in solver.G0_iw:
    g0 << inverse(iOmega_n + mu - delta_w)

sp_in = dict(
    h_int = n('up',0)*n('do',0),
    max_time = -1,
    length_cycle = 50,
    n_warmup_cycles = 50,
    n_cycles = 500,
    move_double = False,
    )
    
solver.solve(**sp_in)

sp = solver.last_solve_parameters

assert( cf_elements_in_A_with_B(sp_in, sp) )

filename = 'h5_read_write.h5'

with HDFArchive(filename, 'w') as A:
    A['solver'] = solver

with HDFArchive(filename, 'r') as A:
    solver_ref = A['solver']

cp_h5_ref = solver_ref.constr_parameters
sp_h5_ref = solver_ref.solve_parameters

assert( cf_elements_in_A_with_B(cp_h5_ref, cp) )
assert( cf_elements_in_A_with_B(sp_h5_ref, sp) )
