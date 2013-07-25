@preamble@
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.utility.mpi import rank, is_master_node
import numpy

# Set up a few parameters
beta = 10.0
U = 2.0
mu = 1.0
t = 0.1
epsilon = 2.3

# Construct a CTQMC solver
from pytriqs.applications.impurity_solvers.ctqmc_hyb import Solver
S = Solver(beta = beta, gf_struct = [ ('tot', [0,1]) ])

# Compute Delta
Delta = S.G0.copy()
Delta['tot'] <<= inverse(iOmega_n - numpy.matrix([[epsilon,-t],[-t,epsilon]])) + inverse(iOmega_n - numpy.matrix([[-epsilon,-t],[-t,-epsilon]]))

# Compute G0
tmat = numpy.array([[mu,t],[t,mu]])
S.G0['tot'] <<= inverse(iOmega_n + tmat - Delta['tot'])

# Solve
from pytriqs.applications.impurity_solvers.operators import *
S.solve(H_local = U * N('tot',0) * N('tot',1),
        quantum_numbers = {'N' : N('tot',0)+N('tot',1)},
        n_cycles  = 30000000,
        length_cycle = 50,
        n_warmup_cycles = 50000,
        n_time_slices_delta = 1000,
        n_time_slices_gtau = 1000,
        legendre_accumulation = False,
        time_accumulation = True,
        random_name = "",
	random_seed = 123*rank + 567,
        use_segment_picture = False)

# Calculation is done. Now save a few things
if is_master_node():
	Results = HDFArchive("spinless.hyb.h5",'w')
	Results["G"] = S.G
	Results["Delta_tau"] = S.Delta_tau
	Results["G_tau"] = S.G_tau

