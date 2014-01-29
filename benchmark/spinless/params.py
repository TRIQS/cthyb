#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi

# Input parameters
beta = 10.0
U = 2.0
mu = 1.0
epsilon = 2.3
t = 0.1

# Use quantum numbers
use_qn = True

p = {}
p["beta"] = beta
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 20000
p["n_cycles"] = 1000000
p["n_tau_delta"] = 1000
p["n_tau_g"] = 200
p["use_truncation"] = True
p["use_old_trace"] = False
p["use_quick_trace_estimator"] = False
p["trace_estimator_n_blocks_guess"] = -1
p["krylov_gs_energy_convergence"] = 1e-8
p["krylov_small_matrix_size"] = 100

def results_file_name(use_qn):
    name = "spinless."
    if use_qn: name += "qn."
    name += "h5"
    return name

matrix_results_file_name = "spinless.matrix.h5"
ref_file_name = "spinless.ed.h5"