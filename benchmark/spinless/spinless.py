#!/bin/env pytriqs

import numpy as np

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.operators import *
from pytriqs.gf.local import *

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

print_master("Welcome to spinless (spinless electrons on a correlated dimer) test.")

H = U*n("tot","A")*n("tot","B")

QN = []
if use_qn:
    QN.append(n("tot","A")+n("tot","B"))
    p["partition_method"] = "quantum_numbers"
    p["quantum_numbers"] = QN

gf_struct = {"tot":["A","B"]}

print_master("Constructing the solver...")

## Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw = n_iw, n_tau = n_tau)

print_master("Preparing the hybridization function...")

## Set hybridization function
delta_w = GfImFreq(indices = ["A","B"], beta=beta)
delta_w << inverse(iOmega_n - np.array([[epsilon,-t],[-t,epsilon]])) + inverse(iOmega_n - np.array([[-epsilon,-t],[-t,-epsilon]]))
S.G0_iw["tot"] << inverse(iOmega_n - np.array([[-mu,-t],[-t,-mu]]) - delta_w)

print_master("Running the simulation...")

## Solve the problem
S.solve(h_loc=H, **p)

## Save the results  
if mpi.rank==0:
    Results = HDFArchive(results_file_name(use_qn),'w')
    Results["tot"] = S.G_tau["tot"]
