import numpy as np
import triqs.utility.mpi as mpi
from triqs.gf import *
from triqs.operators import n
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs.operators.util.op_struct import set_operator_structure
from h5 import HDFArchive
from triqs_cthyb import *
from triqs.utility.comparison_tests import *
from itertools import product

# H_loc parameters
beta = 10.0
norb = 2
mu = 1.0
U = 2.0
J = 0.2
J_h = 2.0

# Poles of delta
epsilon = 2.3 * np.eye(norb)
offdiag = 0.1 * (np.ones(norb) - np.eye(norb))

# Hybridization matrices
V = 1.0 * np.eye(norb) + 2.0 * (np.ones(norb) - np.eye(norb))

# Define Hybridization Function
delta_w = GfImFreq(target_shape=(norb, norb), beta=beta)
delta_w << inverse(iOmega_n - epsilon - offdiag) + inverse(iOmega_n + epsilon - offdiag)
delta_w.from_L_G_R(V, delta_w, V)

# Block structure of GF
spin_names = ('up','down')
gf_struct = set_operator_structure(spin_names,norb,True)

# Hamiltonian
orb_names = list(range(norb))
H = h_int_kanamori(spin_names,orb_names,
                   np.array([[0,U-3*J],[U-3*J,0]]),
                   np.array([[U,U-2*J],[U-2*J,U]]),
                   J_h,off_diag=True)

# Solve Parameters
sp = {}
sp["max_time"] = -1
sp["random_name"] = ""
sp["random_seed"] = 123 * mpi.rank + 567
sp["length_cycle"] = 50
sp["n_warmup_cycles"] = 5000
sp["n_cycles"] = 5000
sp["measure_G_l"] = True
sp["move_double"] = False
sp["move_shift"] = False
sp["pauli_prob"] = 0.5

# Helper Function to write to hdf5 and compare results
def write_and_compare(S, fname):
    if mpi.is_master_node():
        with HDFArchive("{}.out.h5".format(fname),'w') as Results:
            Results["G_tau"] = S.G_tau
            Results["G_leg"] = S.G_l

    if mpi.is_master_node():
        with HDFArchive("{}.ref.h5".format(fname),'r') as Results:
            assert_block_gfs_are_close(Results["G_tau"], S.G_tau)
            assert_block_gfs_are_close(Results["G_leg"], S.G_l)

    from triqs.utility.h5diff import h5diff
    h5diff("{}.out.h5".format(fname),"{}.ref.h5".format(fname))

# ==== G0_iw Interface ====

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=2500, n_l=50)

# check if solver properties for complex support can be accessed
mpi.report('Is cthyb compiled with support for complex hybridization? '+str(S.hybridisation_is_complex))
mpi.report('Is cthyb compiled with support for complex local Hamiltonian? '+str(S.local_hamiltonian_is_complex))

# Set Weiss Field
S.G0_iw << inverse(iOmega_n + mu - 2. * offdiag - delta_w)

S.solve(h_int=H, **sp)

write_and_compare(S, "pauli")

