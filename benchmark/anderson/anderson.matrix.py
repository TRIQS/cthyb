#!/bin/env pytriqs

from datetime import datetime
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

# Remove cthyb-specific parameters 
cthyb_only_params = [
    'beta',
    'verbosity']
p = {k:v for k, v in p.items() if not k in cthyb_only_params}

print_master("Welcome to Anderson (1 correlated site + symmetric bath) test, matrix version.")

H = U*N("up",0)*N("dn",0)
  
QN = {}
for spin in spin_names: QN['N_'+spin] = N(spin,0)

gf_struct = {}
for spin in spin_names: gf_struct[spin] = [0]

print_master("Constructing the solver...")

# Construct the solver
S = Solver(beta=beta, gf_struct=gf_struct.items())

print_master("Preparing the hybridization function...")

# Set hybridization function    
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w <<= (V**2) * inverse(iOmega_n - epsilon) + (V**2) * inverse(iOmega_n + epsilon)
for spin in spin_names:
    S.G0[spin][0,0] <<= inverse(iOmega_n + mu -{'up':h,'dn':-h}[spin] - delta_w)

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
    for b in gf_struct: Results[b] = S.G_tau[b]
