#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi

spin_names = ("up","dn")
def mkind(spin):
    return (spin,0) if use_blocks else ("tot",spin)

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
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50000
p["n_cycles"] = 5000000

def results_file_name(use_blocks,use_qn):
    name = "anderson."
    if use_blocks: name += "block."
    if use_qn: name += "qn."
    name += "h5"
    return name

matrix_results_file_name = "anderson.matrix.h5"
ref_file_name = "anderson.ed.h5"
