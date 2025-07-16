import triqs.utility.mpi as mpi
from h5 import HDFArchive
from triqs.operators import *
#from atom_diag import trace_rho_op
from triqs.atom_diag import trace_rho_op
from triqs_cthyb import *
from triqs.gf import *
import numpy as np

# Input parameters
beta = 10.0
U = 2.0
mu = 1.0
h = 0.1
norb = 2
V = 1.0 * np.eye(norb) + 0.1 * (np.ones(norb) - np.eye(norb))
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
p["measure_G_tau"] = False
p["use_norm_as_weight"] = True
p["measure_density_matrix"] = True
p["time_invariance"] = True

gm = {}
gm['flip_spins'] = {("up",0) : ("dn",0), ("dn",0) : ("up",0), ("up",1) : ("dn",1), ("dn",1) : ("up",1)}
gm['swap_orbs']  = {("up",0) : ("up",1), ("up",1) : ("up",0), ("dn",0) : ("dn",1), ("dn",1) : ("dn",0)}
p["move_global"] = gm
p["move_global_prob"] = 0.06

qn = [n("up",0) + n("up",1),n("dn",0) + n("dn",1)]
p["quantum_numbers"] = qn
p["partition_method"] = "quantum_numbers"

H = U*n("up",0)*n("dn",0) + U*n("up",1)*n("dn",1)
H = H + 0.5*h*(n("up",0) - n("dn",0)) + 0.5*h*(n("up",1) - n("dn",1))

# Construct the solver
S = Solver(beta=beta, gf_struct=[["dn",2], ["up",2]], n_tau=n_tau, n_iw=n_iw)

# Set hybridization function
delta_w = GfImFreq(beta=beta, target_shape=(2,2))
delta_w << inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon)
delta_w.from_L_G_R(V, delta_w, V)

S.G0_iw << inverse(iOmega_n + mu - delta_w)

# Solve!
S.solve(h_int=H, **p)

if mpi.is_master_node():
    # Measure expectation values
    dm = S.density_matrix
    static_observables = {"N1_up" : n("up",0), "N1_dn" : n("dn",0),
                          "N2_up" : n("up",1), "N2_dn" : n("dn",1),
                          "N12_up": c_dag("up",0) * c("up",1),
                          "N12_dn": c_dag("dn",0) * c("dn",1)}
    with HDFArchive('measure_static_time_invariance.out.h5','w') as ar:
        for name,op in static_observables.items():
            ave = trace_rho_op(dm,op,S.h_loc_diagonalization)
            assert( np.abs(ave.imag) < 1e-10 )
            ar[name] = ave.real

from triqs.utility.h5diff import h5diff
h5diff("measure_static_time_invariance.out.h5","measure_static_time_invariance.ref.h5")
