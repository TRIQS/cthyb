#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb_krylov import *
from pytriqs.gf.local import *

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

pp = Parameters()
for k in p: pp[k] = p[k]

print_master("Welcome to Anderson (1 correlated site + symmetric bath) test.")

H = U*N(*mkind("up"))*N(*mkind("dn")) + (-mu+h)*N(*mkind("up")) + (-mu-h)*N(*mkind("dn"))
  
QN = []
if use_qn:
    for spin in spin_names: QN.append(N(*mkind(spin)))

gf_struct = {}
for spin in spin_names:
    bn, i = mkind(spin)
    gf_struct.setdefault(bn,[]).append(i)

print_master("Constructing the solver...")

# Construct the solver
S = Solver(parameters=pp, H_local=H, quantum_numbers=QN, gf_struct=gf_struct)

print_master("Preparing the hybridization function...")

# Set hybridization function    
delta_w = GfImFreq(indices = [0], beta=beta)
delta_w <<= (V**2) * inverse(iOmega_n - epsilon) + (V**2) * inverse(iOmega_n + epsilon)
for spin in spin_names:
    bn, i = mkind(spin)
    S.Delta_tau[bn][i,i] <<= InverseFourier(delta_w)

print_master("Running the simulation...")

# Solve the problem
S.solve(parameters=pp)

# Save the results  
if mpi.rank==0:
    Results = HDFArchive(results_file_name(use_blocks,use_qn),'w')
    for b in gf_struct: Results[b] = S.G_tau[b]
