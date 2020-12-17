from triqs.gf import *
from triqs.operators import *
from h5 import *
import triqs.utility.mpi as mpi
from triqs.applications.impurity_solvers.cthyb import Solver

# Set up a few parameters
U = 2.5
half_bandwidth = 1.0
chemical_potential = U/2.0
beta = 100
n_loops = 5

# Construct the CTQMC solver
S = Solver(beta = beta, gf_struct = [ ('up',1), ('down',1) ])

# Initalize the Green's function to a semi circular DOS
S.G_iw << SemiCircular(half_bandwidth)

# Now do the DMFT loop
for i_loop in range(n_loops):

    # Compute new S.G0_iw with the self-consistency condition while imposing paramagnetism
    g = 0.5 * (S.G_iw['up'] + S.G_iw['down'])
    for name, g0 in S.G0_iw:
        g0 << inverse(iOmega_n + chemical_potential - (half_bandwidth/2.0)**2 * g)

    # Run the solver
    S.solve(h_int = U * n('up',0) * n('down',0),    # Local Hamiltonian
            n_cycles = 5000,                        # Number of QMC cycles
            length_cycle = 200,                     # Length of a cycle
            n_warmup_cycles = 1000)                 # How many warmup cycles

    # Some intermediate saves
    if mpi.is_master_node():
        with HDFArchive("dmft_solution.h5") as ar:
            ar["G_tau-%s"%i_loop] = S.G_tau
            ar["G_iw-%s"%i_loop] = S.G_iw
            ar["Sigma_iw-%s"%i_loop] = S.Sigma_iw

    # Here we could write some convergence test
    # if converged : break
