#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi

spin_names = ("up","dn")
def mkind(spin): return (spin,0)

# Input parameters
beta = 10.0
U = 2.0
mu = 1.0
h = 0.0
V = 1.0
epsilon1 = 2.1
epsilon2 = -2.4

n_iw = 1025
n_tau = 10001
n_l = 50

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50000
p["n_cycles"] = 5000000
p["measure_g_tau"] = False
p["measure_g_l"] = True

results_file_name = "legendre.h5"
matrix_results_file_name = "legendre.matrix.h5"
ref_file_name = "legendre.ed.h5"
