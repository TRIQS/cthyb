import pytriqs.utility.mpi as mpi
from pytriqs.operators import *
from pytriqs.operators.util.op_struct import set_operator_structure
from pytriqs.operators.util.U_matrix import U_matrix
from pytriqs.operators.util.hamiltonians import h_int_slater
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
from pytriqs.utility.comparison_tests import *

beta = 100.0
# H_loc parameters
L = 2 # angular momentum
U = 5.0
J = 0.1
F0 = U
F2 = J*(14.0/(1.0 + 0.63))
F4 = F2*0.63
half_bandwidth = 1.0
mu = 32.5  # 3 electrons in 5 bands

spin_names = ("up","down")
cubic_names = map(str,range(2*L+1))
U_mat = U_matrix(L, radial_integrals=[F0,F2,F4], basis="cubic")

# Parameters
p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50
p["n_cycles"] = 5000
p["measure_g_l"] = True
p["move_double"] = False

# Block structure of GF
gf_struct = set_operator_structure(spin_names,cubic_names,False)

# Local Hamiltonian
H = h_int_slater(spin_names,cubic_names,U_mat,False)

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=100000)

# Set hybridization function
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w << (half_bandwidth/2.0)**2 * SemiCircular(half_bandwidth)
for name, g0 in S.G0_iw:
    g0 << inverse(iOmega_n + mu - delta_w)


S.solve(h_int=H, **p)

if mpi.is_master_node():
    with HDFArchive("slater.out.h5",'w') as Results:
        Results["G_tau"] = S.G_tau
        Results["G_leg"] = S.G_l

if mpi.is_master_node():
    with HDFArchive("slater.ref.h5",'r') as Results:
        assert_block_gfs_are_close(Results["G_tau"], S.G_tau)
        assert_block_gfs_are_close(Results["G_leg"], S.G_l)
