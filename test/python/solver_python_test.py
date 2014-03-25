#!/bin/python

import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *
from pytriqs.parameters.parameters import Parameters
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb_matrix.cthyb_matrix import *
import itertools
from collections import OrderedDict

# H_loc parameters
beta = 10.0
num_orbitals = 2
mu = 1.0
U = 2.0
J = 0.2
n_n_only = False

# Poles of delta
epsilon = 2.3

# Hybridization matrices
V = 1.0 * np.eye(num_orbitals) + 0.1 * (np.ones(num_orbitals) - np.eye(num_orbitals))

# Parameters
p = Parameters()
p["beta"] = beta
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50
p["n_cycles"] = 500
p["n_tau_delta"] = 1000
p["n_tau_g"] = 1000
p["krylov_bs_use_cutoff"] = True
p["krylov_bs_prob_cutoff"] = -1.0
p["krylov_small_matrix_size"] = 25

# Block structure of GF
gf_struct = OrderedDict()
gf_struct['up'] = range(0,num_orbitals)
gf_struct['down'] = range(0,num_orbitals)

# Hamiltonian
H = Operator()
for o in range(0,num_orbitals):
    H += -mu*(N("up",o) + N("down",o))

for o in range(0,num_orbitals):
    H += U*N("up",o)*N("down",o)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o1==o2: continue
    H += (U-2*J)*N("up",o1)*N("down",o2)

for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
    if o2>=o1: continue;
    H += (U-3*J)*N("up",o1)*N("up",o2)
    H += (U-3*J)*N("down",o1)*N("down",o2)
      
if not n_n_only: # spin flips and pair hopping
    for o1,o2 in itertools.product(range(0,num_orbitals),range(0,num_orbitals)):
        if o1==o2: continue
        H += -J*C_dag("up",o1)*C_dag("down",o1)*C("up",o2)*C("down",o2)
        H += -J*C_dag("up",o1)*C_dag("down",o2)*C("up",o2)*C("down",o1)

# Quantum numbers
qn = [Operator(),Operator()]
for o in range(0,num_orbitals):
    qn[0] += N("up",o)
    qn[1] += N("down",o)
    
S = Solver(parameters=p, H_local=H, quantum_numbers=qn, gf_struct=gf_struct)

# Set hybridization function
delta_w = GfImFreq(indices = range(0,num_orbitals), beta=beta)
delta_w <<= inverse(iOmega_n - epsilon) + inverse(iOmega_n + epsilon)
delta_w.from_L_G_R(V, delta_w, V)

S.Delta_tau["up"] <<= InverseFourier(delta_w)
S.Delta_tau["down"] <<= InverseFourier(delta_w)

S.solve(parameters=p)
  
if mpi.rank==0:
    Results = HDFArchive("solver_python_test.output.h5",'w')
    Results["G_up"] = S.G_tau["up"]
    Results["G_down"] = S.G_tau["down"]
