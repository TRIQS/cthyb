import numpy as np
import triqs.utility.mpi as mpi
from triqs.gf import *
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs.operators.util.op_struct import set_operator_structure
from h5 import HDFArchive
from triqs_cthyb import *
from triqs.utility.comparison_tests import *

# Test for Hilbert space truncation feature
# Uses a 2-orbital Kanamori model with energy truncation

# H_loc parameters
beta = 10.0
n_orb = 2
mu = 1.0
U = 2.0
J = 0.2

# Poles of delta
epsilon = 2.3

# Hybridization matrices
V = 1.0 * np.eye(n_orb) + 0.1 * (np.ones(n_orb) - np.eye(n_orb))

# Define Hybridization Function
delta_w = GfImFreq(target_shape=(n_orb, n_orb), beta=beta)
delta_w << inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon)
delta_w.from_L_G_R(V, delta_w, V)

# Block structure of GF
spin_names = ('up','down')
gf_struct = set_operator_structure(spin_names, n_orb, True)

# Hamiltonian
H = h_int_kanamori(spin_names, n_orb,
                   np.array([[0, U-3*J], [U-3*J, 0]]),
                   np.array([[U, U-2*J], [U-2*J, U]]),
                   J, off_diag=True)

# Solve Parameters with truncation
sp = {}
sp["max_time"] = -1
sp["random_name"] = ""
sp["random_seed"] = 123 * mpi.rank + 567
sp["length_cycle"] = 50
sp["n_warmup_cycles"] = 50
sp["n_cycles"] = 5000
sp["measure_G_l"] = True
sp["move_double"] = False
# Truncation: keep only states within 2.0 of the ground state energy
sp["truncate_energy_cutoff"] = 2.0

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=2500, n_l=50)

# Set Weiss Field
S.G0_iw << inverse(iOmega_n + mu - delta_w)

S.solve(h_int=H, **sp)

if mpi.is_master_node():
    with HDFArchive("truncation.out.h5", 'w') as Results:
        Results["G_tau"] = S.G_tau
        Results["G_leg"] = S.G_l

if mpi.is_master_node():
    with HDFArchive("truncation.ref.h5", 'r') as Results:
        assert_block_gfs_are_close(Results["G_tau"], S.G_tau)
        assert_block_gfs_are_close(Results["G_leg"], S.G_l)

from triqs.utility.h5diff import h5diff
h5diff("truncation.out.h5", "truncation.ref.h5")
