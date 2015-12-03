import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
import numpy as np

# Input parameters
beta = 10.0
U = 2.0
mu = 1.0
h = 0.1
V = 1.0
t = 0.1
epsilon = 2.3

n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 1000
p["n_cycles"] = 10000
p["measure_g_tau"] = False
p["move_double"] = False
p["measure_density_matrix"] = True
p["use_norm_as_weight"] = True

H = U*n("up",1)*n("dn",1) + U*n("up",2)*n("dn",2)
H = H + 0.5*h*(n("up",1) - n("dn",1)) + 0.5*h*(n("up",2) - n("dn",2))

# Construct the solver
S = Solver(beta=beta, gf_struct={"up":[1,2], "dn":[1,2]}, n_tau=n_tau, n_iw=n_iw)

# Set hybridization function
delta_w = GfImFreq(indices = [1,2], beta=beta)
delta_w << (V**2)*(inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon))
for bn, g in S.G0_iw: g << inverse(iOmega_n - np.matrix([[-mu,t],[t,-mu]]) - delta_w)

# Solve!
S.solve(h_int=H, **p)

if mpi.is_master_node():
    # Measure expectation values
    dm = S.density_matrix
    static_observables = {"N1_up" : n("up",1), "N1_dn" : n("dn",1),
                          "N2_up" : n("up",2), "N2_dn" : n("dn",2)}
    with HDFArchive('measure_static.out.h5','w') as ar:
        for name,op in static_observables.iteritems():
            ave = trace_rho_op(dm,op,S.h_loc_diagonalization)
            ar[name] = ave

from pytriqs.utility.h5diff import h5diff
h5diff("measure_static.out.h5","measure_static.ref.h5")
