@preamble@
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.utility.mpi import *
import numpy

# Set up a few parameters
beta = 20.0
U = 2.0
V = 1.0
h = 0.0
mu = 1.0
epsilon = 2.3

# Construct a CTQMC solver
from pytriqs.applications.impurity_solvers.ctqmc_hyb import Solver
S = Solver(beta = beta, gf_struct = [ ('up',[0]), ('down',[0]) ])

# Compute G0
S.G0['up'] <<= inverse( iOmega_n - h + mu - V*V*inverse(iOmega_n - epsilon) - V*V*inverse(iOmega_n + epsilon) )
S.G0['down'] <<= inverse( iOmega_n + h + mu - V*V*inverse(iOmega_n - epsilon) - V*V*inverse(iOmega_n + epsilon) )

# Solve
from pytriqs.applications.impurity_solvers.operators import *
S.solve(H_local = U * N('up',0) * N('down',0),
        quantum_numbers = { 'Nup' : N('up',0), 'Ndown' : N('down',0) },
        n_cycles  = 400000000,
        length_cycle = 50,
        n_warmup_cycles = 20000000,
        n_time_slices_delta = 1000,
        n_time_slices_gtau = 1000,
        legendre_accumulation = False,
        time_accumulation = True,
        random_name = "",
	random_seed = 123*rank + 567,
        use_segment_picture = True)

# Calculation is done. Now save a few things
if is_master_node():
	Results = HDFArchive("anderson.hyb.h5",'w')
	Results["G"] = S.G
	Results["G_tau"] = S.G_tau

