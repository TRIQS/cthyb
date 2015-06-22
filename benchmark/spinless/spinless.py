#!/bin/env pytriqs

import numpy as np

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.operators import *
from pytriqs.gf.local import *

# Input parameters
beta = 10.0
U = 2.0
mu = 1.0
epsilon = 2.3
t = 0.1

# Use quantum numbers
use_qn = True

n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 20000
p["n_cycles"] = 1000000

results_file_name = "spinless"
if use_qn: results_file_name += ".qn"
results_file_name += ".h5"

mpi.report("Welcome to spinless (spinless electrons on a correlated dimer) test.")

H = U*n("tot","A")*n("tot","B")

QN = []
if use_qn:
    QN.append(n("tot","A")+n("tot","B"))
    p["partition_method"] = "quantum_numbers"
    p["quantum_numbers"] = QN

gf_struct = {"tot":["A","B"]}

mpi.report("Constructing the solver...")

## Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw = n_iw, n_tau = n_tau)

mpi.report("Preparing the hybridization function...")

## Set hybridization function
delta_w = GfImFreq(indices = ["A","B"], beta=beta)
delta_w << inverse(iOmega_n - np.array([[epsilon,-t],[-t,epsilon]])) + inverse(iOmega_n - np.array([[-epsilon,-t],[-t,-epsilon]]))
S.G0_iw["tot"] << inverse(iOmega_n - np.array([[-mu,-t],[-t,-mu]]) - delta_w)

mpi.report("Running the simulation...")

## Solve the problem
S.solve(h_int=H, **p)

## Save the results  
if mpi.is_master_node():
    with HDFArchive(results_file_name,'w') as Results:
        Results["tot"] = S.G_tau["tot"]
