#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi

spin_names = ("up","dn")
def mkind(spin,n): return spin + "-" + str(n),0

# Input parameters
beta = 10.0
num_orbitals = 2
mu = 1.0
U = 2.0
J = 0.2
V = 1.0
epsilon = 2.3

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
p["n_cycles"] = 3000000

results_file_name = "kanamori." + ("qn." if use_qn else "") + "h5"
matrix_results_file_name = "kanamori.matrix.h5"
ref_file_name = "kanamori.ed.h5"
