import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive
from pytriqs.utility.comparison_tests import *

from cthyb import *

#  Example of DMFT single site solution with CTQMC

# set up a few parameters
half_bandwidth = 1.0
U = 2.5
mu = (U/2.0)+0.2
beta = 100.0

gf_struct = {'up':[0],'down':[0]}

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=3000, n_l=30)

# Local Hamiltonian
H = U*n("up",0)*n("down",0)

# init the Green function
S.G_iw << SemiCircular(half_bandwidth)

# Compute G0
for name, g0 in S.G0_iw:
  g0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2  * S.G_iw[name])

# Parameters
p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 500
p["n_warmup_cycles"] = 5000
p["n_cycles"] = 5000
p["measure_g_l"] = True
p["move_double"] = False
p["perform_tail_fit"] = True
p["fit_max_moment"] = 2

S.solve(h_int=H, **p)

# Calculation is done. Now save a few things
if mpi.is_master_node():
    with HDFArchive("single_site_bethe.out.h5",'w') as Results:
        Results["Sigma_iw"] = S.Sigma_iw
        Results["G_tau"] = S.G_tau
        Results["G_iw"] = S.G_iw
        Results["G_l"] = S.G_l

# Check against reference
if mpi.is_master_node():
    with HDFArchive("single_site_bethe.ref.h5",'r') as Results:
        assert_block_gfs_are_close(Results["Sigma_iw"], S.Sigma_iw)
        assert_block_gfs_are_close(Results["G_tau"], S.G_tau)
        assert_block_gfs_are_close(Results["G_iw"], S.G_iw)
        assert_block_gfs_are_close(Results["G_l"], S.G_l)

# Redondant with previous check
from pytriqs.utility.h5diff import h5diff
h5diff("single_site_bethe.out.h5","single_site_bethe.ref.h5")

