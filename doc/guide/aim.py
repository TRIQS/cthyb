from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import Solver
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi

# Parameters
D, V, U = 1.0, 0.2, 4.0
e_f, beta = -U/2.0, 50

# Construct the impurity solver with the inverse temperature
# and the structure of the Green's functions
S = Solver(beta = beta, gf_struct = [ ['up',[0]], ['down',[0]] ], n_l = 100)

# Initialize the non-interacting Green's function S.G0_iw
for name, g0 in S.G0_iw: g0 << inverse(iOmega_n - e_f - V**2 * Wilson(D))

# Run the solver. The results will be in S.G_tau, S.G_iw and S.G_l
S.solve(h_int = U * n('up',0) * n('down',0),     # Local Hamiltonian
        n_cycles  = 500000,                      # Number of QMC cycles
        length_cycle = 200,                      # Length of one cycle
        n_warmup_cycles = 10000,                 # Warmup cycles
        measure_g_l = True)                      # Measure G_l

# Save the results in an HDF5 file (only on the master node)
if mpi.is_master_node():
    with HDFArchive("aim_solution.h5",'w') as Results:
        Results["G_tau"] = S.G_tau
        Results["G_iw"] = S.G_iw
        Results["G_l"] = S.G_l
