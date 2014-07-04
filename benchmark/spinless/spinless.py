#!/bin/env pytriqs

import numpy as np

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.operators.operators2 import *
from pytriqs.gf.local import *

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

pp = SolverCore.solve_parameters()
for k in p: pp[k] = p[k]

print_master("Welcome to spinless (spinless electrons on a correlated dimer) test.")

H = U*n("tot",0)*n("tot",1)
  
Qn = []
if use_qn: Qn.append(n("tot",0)+n("tot",1))

gf_struct = {"tot":[0,1]}

print_master("Constructing the solver...")

## Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_iw = n_iw, n_tau = n_tau)

print_master("Preparing the hybridization function...")

## Set hybridization function    
delta_w = GfImFreq(indices = [0,1], beta=beta)
delta_w <<= inverse(iOmega_n - np.array([[epsilon,-t],[-t,epsilon]])) + inverse(iOmega_n - np.array([[-epsilon,-t],[-t,-epsilon]]))
S.G0_iw["tot"] <<= inverse(iOmega_n - np.array([[-mu,-t],[-t,-mu]]) - delta_w)

print_master("Running the simulation...")

## Solve the problem
S.solve(h_loc=H, params=pp, quantum_numbers=Qn, use_quantum_numbers=True)

## Save the results  
if mpi.rank==0:
    Results = HDFArchive(results_file_name(use_qn),'w')
    Results["tot"] = S.G_tau["tot"]
