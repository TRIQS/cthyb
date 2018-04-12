import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf import *
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from pytriqs.operators.util.op_struct import set_operator_structure
from pytriqs.operators.util.observables import S_op
from pytriqs.archive import HDFArchive
from cthyb import *
from pytriqs.utility.comparison_tests import *

# H_loc parameters
beta = 60.0
num_orbitals = 2
mu = 1.0
U = 2.0
J = 0.2
h = 0.1

# Poles of delta
epsilon = 2.3

# Hybridization matrices
V = 2.0 * np.eye(num_orbitals) + 0.2 * (np.ones(num_orbitals) - np.eye(num_orbitals))

# Block structure of GF
spin_names = ('up','dn')
orb_names = range(num_orbitals)
gf_struct = set_operator_structure(spin_names,orb_names,True)
gf_struct.reverse() # the reference data was computed with reversed block order

# Construct solver
S = Solver(beta=beta, gf_struct=gf_struct, n_iw=1025, n_tau=2500, n_l=50)

# Hamiltonian
H = h_int_kanamori(spin_names,orb_names,
                   np.array([[0,U-3*J],[U-3*J,0]]),
                   np.array([[U,U-2*J],[U-2*J,U]]),
                   J,True)
H += h*S_op('z',spin_names,orb_names,True)

# Set hybridization function
delta_w = GfImFreq(indices = orb_names, beta=beta)
delta_w << inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon)
delta_w.from_L_G_R(V, delta_w, V)
S.G0_iw << inverse(iOmega_n + mu - delta_w)

# Parameters
p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 1000
p["n_cycles"] = 5000
p["measure_g_tau"] = True
p["measure_g_l"] = False
p["move_double"] = False

# Global moves
gm = {}
gm['flip_spins'] = {("up",0) : ("dn",0), ("dn",0) : ("up",0), ("up",1) : ("dn",1), ("dn",1) : ("up",1)}
gm['swap_orbs']  = {("up",0) : ("up",1), ("up",1) : ("up",0), ("dn",0) : ("dn",1), ("dn",1) : ("dn",0)}
p["move_global"] = gm
p["move_global_prob"] = 0.06

S.solve(h_int=H, **p)

if mpi.is_master_node():
    with HDFArchive("move_global.out.h5",'w') as Results:
        Results["G_tau"] = S.G_tau

if mpi.is_master_node():
    with HDFArchive("move_global.ref.h5",'r') as Results:
        assert_block_gfs_are_close(Results["G_tau"], S.G_tau)

from pytriqs.utility.h5diff import h5diff
h5diff("move_global.out.h5","move_global.ref.h5")
