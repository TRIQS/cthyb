#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.gf import Gf, MeshImFreq, iOmega_n, inverse
from pytriqs.operators import n
from h5 import HDFArchive
from triqs_cthyb import SolverCore

mpi.report("Welcome to asymm_bath test (1 band with a small asymmetric hybridization function).")
mpi.report("This test helps to detect sampling problems.")

# H_loc parameters
beta = 40.0
ed = -1.0
U = 2.0

epsilon = [0.0,0.1,0.2,0.3,0.4,0.5]
V = 0.2

# Parameters
n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 20000
p["n_cycles"] = 1000000
p["performance_analysis"] = True
p["measure_pert_order"] = True

# Block structure of GF
gf_struct = [['up', [0]], ['dn', [0]]]

# Hamiltonian
H = U*n("up",0)*n("dn",0)

# Quantum numbers
qn = [n("up",0),n("dn",0)]
p["partition_method"] = "quantum_numbers"
p["quantum_numbers"] = qn

# Construct solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

def read_histo(f,type_of_col_1):
    histo = []
    for line in f:
        cols = [s for s in line.split(' ') if s]
        histo.append((type_of_col_1(cols[0]),float(cols[1]),float(cols[2])))
    return histo

if mpi.is_master_node():
    arch = HDFArchive('asymm_bath.h5','w')

# Set hybridization function
for e in epsilon:
    delta_w = Gf(mesh=S.G0_iw.mesh, target_shape=[1, 1]) 
    delta_w << (V**2) * inverse(iOmega_n - e)

    S.G0_iw["up"] << inverse(iOmega_n - ed - delta_w)
    S.G0_iw["dn"] << inverse(iOmega_n - ed - delta_w)

    S.solve(h_int=H, **p)

    if mpi.is_master_node():
        arch.create_group('epsilon_' + str(e))
        gr = arch['epsilon_' + str(e)]
        gr['G_tau'] = S.G_tau
        gr['beta'] = beta
        gr['U'] = U
        gr['ed'] = ed
        gr['V'] = V
        gr['e'] = e
        gr['perturbation_order'] = S.perturbation_order
        gr['perturbation_order_total'] = S.perturbation_order_total
        gr['performance_analysis'] = S.performance_analysis
