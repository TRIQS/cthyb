import triqs.utility.mpi as mpi
from triqs.gf import *
from triqs.operators import *
from h5 import HDFArchive
from triqs.utility.comparison_tests import *

from triqs_cthyb import *

#  Example of DMFT single site solution with CTQMC

# set up a few parameters
half_bandwidth = 1.0
U = 2.5
mu = (U/2.0)+0.2
beta = 80.0

gf_struct = [['down',1], ['up',1]]

# Parameters
p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 500
p["n_warmup_cycles"] = 5000
p["n_cycles"] = 5000
p["measure_G_l"] = True
p["move_double"] = False
p["perform_tail_fit"] = True
p["fit_max_moment"] = 3
p["fit_min_w"] = 1.2
p["fit_max_w"] = 3.0

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=8001, n_l=30)

# Local Hamiltonian
H = U*n("up",0)*n("down",0)

# init the Green function
S.G_iw << SemiCircular(half_bandwidth)


for i in range(2):

    g = 0.5 * ( S.G_iw['up'] + S.G_iw['down'] )
    # Compute G0
    for name, g0 in S.G0_iw:
      g0 << inverse(iOmega_n + mu - (half_bandwidth/2.0)**2  * g)

    S.solve(h_int=H, **p)

# Calculation is done. Now save a few things
if mpi.is_master_node():
    with HDFArchive("single_site_bethe.out.h5",'w') as Results:

        Results["Delta_infty"] = S.Delta_infty
        Results["Delta_tau"] = S.Delta_tau

        Results["G0_iw"] = S.G0_iw

        Results["G_tau"] = S.G_tau
        Results["G_l"] = S.G_l

        Results["G_iw"] = S.G_iw
        Results["G_iw_raw"] = S.G_iw_raw

        # we store Sigma_iw_arw here, but comparing noisy Sigma from Dyson after 2 DMFT
        # iteration is pointless. The tail fitted Sigma_iw can be compared to some degree
        Results["Sigma_iw"] = S.Sigma_iw
        Results["Sigma_iw_raw"] = S.Sigma_iw_raw

# Check against reference
if mpi.is_master_node():
    with HDFArchive("single_site_bethe.ref.h5",'r') as Results:

        assert_block_gfs_are_close(Results["G_tau"], S.G_tau)
        assert_block_gfs_are_close(Results["G_l"], S.G_l)

        assert_block_gfs_are_close(Results["G_iw"], S.G_iw)
        assert_block_gfs_are_close(Results["G_iw_raw"], S.G_iw_raw)

        assert_block_gfs_are_close(Results["Sigma_iw"], S.Sigma_iw, precision=1e-4)

        assert_block_gfs_are_close(Results["Delta_tau"], S.Delta_tau)

