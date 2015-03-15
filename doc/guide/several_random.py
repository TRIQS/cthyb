from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import Solver
import pytriqs.utility.mpi as mpi

D, V, U = 1.0, 0.2, 4.0
e_f, beta = -U/2.0, 50

# Construct the impurity solver
S = Solver(beta = beta, gf_struct = {'up':[0],'down':[0]})

# Loop for two random generators
for random_name in ['mt11213b','lagged_fibonacci607']:

    for spin, g0 in S.G0_iw:
        g0 << inverse( iOmega_n - e_f - V**2 * Wilson(D) )
  
    # Solve using random_name as a generator
    S.solve(h_loc = U * n('up',0) * n('down',0),   # Local Hamiltonian
            n_cycles = 100000,                     # Number of QMC cycles
            length_cycle = 200,                    # Length of one cycle
            n_warmup_cycles = 10000,               # Warmup cycles
            random_name = random_name)             # Name of the random generator
  
    # Save the results in an hdf5 file (only on the master node)
    if mpi.is_master_node():
        Results = HDFArchive("random.h5")
        Results["G_iw%s"%(random_name)] = S.G_iw
        del Results
