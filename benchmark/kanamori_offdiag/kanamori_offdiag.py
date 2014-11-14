#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
from itertools import product

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

print_master("Welcome to Kanamori benchmark.")

gf_struct = {}
for spin in spin_names:
    for o in range(num_orbitals):
        bn, i = mkind(spin,o)
        gf_struct.setdefault(bn,[]).append(i)

# Hamiltonian
H = Operator()

for o in range(num_orbitals):
    H += U *n(*mkind("up",o))*n(*mkind("dn",o))

for o1, o2 in product(range(num_orbitals),range(num_orbitals)):
    if o1==o2: continue
    H += (U-2*J)*n(*mkind("up",o1))*n(*mkind("dn",o2))

for o1, o2 in product(range(num_orbitals),range(num_orbitals)):
    if o2>=o1: continue
    H += (U-3*J)*n(*mkind("up",o1))*n(*mkind("up",o2))
    H += (U-3*J)*n(*mkind("dn",o1))*n(*mkind("dn",o2))

for o1, o2 in product(range(num_orbitals),range(num_orbitals)):
    if(o1==o2): continue
    H += -J*c_dag(*mkind("up",o1))*c_dag(*mkind("dn",o1))*c(*mkind("up",o2))*c(*mkind("dn",o2))
    H += -J*c_dag(*mkind("up",o1))*c_dag(*mkind("dn",o2))*c(*mkind("up",o2))*c(*mkind("dn",o1))

if use_qn:
    QN=[]
    QN.append(sum([n(*mkind("up",o)) for o in range(num_orbitals)],Operator()))
    QN.append(sum([n(*mkind("dn",o)) for o in range(num_orbitals)],Operator()))
    for o in range(num_orbitals):
        dn = n(*mkind("up",o)) - n(*mkind("dn",o))
        QN.append(dn*dn)
    p["partition_method"] = "quantum_numbers"
    p["quantum_numbers"] = QN

print_master("Constructing the solver...")

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

print_master("Preparing the hybridization function...")

# Set hybridization function    
delta_w = GfImFreq(indices = range(num_orbitals), beta=beta, n_points=n_iw)
delta_w_part = delta_w.copy()
for e, v in zip(epsilon,V):
    delta_w_part <<= inverse(iOmega_n - e)
    delta_w_part.from_L_G_R(np.transpose(v),delta_w_part,v)
    delta_w += delta_w_part

for spin in spin_names:
    S.G0_iw[spin] <<= inverse(iOmega_n + mu - delta_w)

print_master("Running the simulation...")

# Solve the problem
S.solve(h_loc=H, **p)

# Save the results
if mpi.rank==0:
    Results = HDFArchive(results_file_name,'w')
    for b in gf_struct: Results[b] = S.G_tau[b]
