#!/bin/env pytriqs

from datetime import datetime
import numpy

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb_matrix import *
from pytriqs.gf.local import *

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

# Remove Krylov-specific parameters 
krylov_only_params = [
    'beta',
    'verbosity',
    'use_old_trace',
    'trace_estimator_n_blocks_guess',
    'use_truncation',
    'use_quick_trace_estimator',
    'krylov_gs_energy_convergence',
    'krylov_small_matrix_size']
p = {k:v for k, v in p.items() if not k in krylov_only_params}


print_master("Welcome to spinless (spinless electrons on a correlated dimer) test, matrix version.")

H = U*N("tot",0)*N("tot",1) -mu*(N("tot",0) + N("tot",1)) -t*(Cdag("tot",0)*C("tot",1) + Cdag("tot",1)*C("tot",0))
  
QN = {}
if use_qn: QN["N"] = N("tot",0)+N("tot",1)

gf_struct = {"tot":[0,1]}

print_master("Constructing the solver...")

## Construct the solver
S = Solver(beta=beta, gf_struct=gf_struct.items())

print_master("Preparing the hybridization function...")

## Set hybridization function    
delta_w = GfImFreq(indices = [0,1], beta=beta)
delta_w <<= inverse(iOmega_n - numpy.matrix([[epsilon,-t],[-t,epsilon]])) + inverse(iOmega_n - numpy.matrix([[-epsilon,-t],[-t,-epsilon]]))
S.G0["tot"] <<= inverse(iOmega_n - delta_w)

print_master("Running the simulation...")

# Solve the problem
p['H_local'] = H
p['quantum_numbers'] = QN

p['legendre_accumulation'] = False
p['time_accumulation'] = True

start_time = datetime.now()
S.solve(**p)
print_master("Simulation lasted: " + str((datetime.now() - start_time).total_seconds()) + " seconds")

# Save the results  
if mpi.rank==0:
    Results = HDFArchive(matrix_results_file_name,'w')
    Results["tot"] = S.G_tau["tot"]