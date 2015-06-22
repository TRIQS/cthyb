#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *

spin_names = ("up","dn")
mkind = lambda spin: (spin,0) if use_blocks else ("tot",spin)

# Input parameters
beta = 10.0
U = 2.0
mu = 1.0
h = 0.0
V = 1.0
epsilon = 2.3

# Use block structure of G
use_blocks = True
# Use quantum numbers
use_qn = True

n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50000
p["n_cycles"] = 5000000

results_file_name = "anderson"
if use_blocks: results_file_name += ".block"
if use_qn: results_file_name += ".qn"
results_file_name += ".h5"

mpi.report("Welcome to Anderson (1 correlated site + symmetric bath) test.")

H = U*n(*mkind("up"))*n(*mkind("dn"))

QN = []
if use_qn:
    for spin in spin_names: QN.append(n(*mkind(spin)))
    p["quantum_numbers"] = QN
    p["partition_method"] = "quantum_numbers"

gf_struct = {}
for spin in spin_names:
    bn, i = mkind(spin)
    gf_struct.setdefault(bn,[]).append(i)

mpi.report("Constructing the solver...")

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

mpi.report("Preparing the hybridization function...")

# Set hybridization function    
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w << (V**2) * inverse(iOmega_n - epsilon) + (V**2) * inverse(iOmega_n + epsilon)
for spin in spin_names:
    bn, i = mkind(spin)
    S.G0_iw[bn][i,i] << inverse(iOmega_n + mu - {'up':h,'dn':-h}[spin] - delta_w)

mpi.report("Running the simulation...")

# Solve the problem
S.solve(h_int=H, **p)

# Save the results  
if mpi.is_master_node():
    with HDFArchive(results_file_name,'w') as Results:
        Results['G_tau'] = S.G_tau
