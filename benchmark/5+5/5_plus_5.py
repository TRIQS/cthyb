#!/bin/env pytriqs

from multiorbital import *

from itertools import product
import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.operators.operators2 import *
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
from collections import OrderedDict

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

print_master("Welcome to 5+5 (5 orbitals + 5 bath sites) test.")

cubic_names, W = spherical2cubic(L)
N_comp = len(cubic_names)

gf_struct = OrderedDict()
for sn, cn in product(spin_names,cubic_names):
    bn, i = mkind(sn,cn)
    gf_struct[bn] = [i]

# Local Hamiltonian
H = Operator()
H_term = Operator()

if H_dump: H_dump_file = open(H_dump,'w')

if use_interaction:
    print_master("Preparing the interaction matrix...")

    U_matrix = transform_U_matrix(U_matrix_spherical([F0,F2,F4]),W)

    # Interaction terms of the Hamiltonian
    a_range=range(N_comp)
    for s1, s2 in product(spin_names,spin_names):
        for ap1, ap2, a1, a2 in product(a_range,a_range,a_range,a_range):
            U_val = U_matrix[ap1,ap2,a1,a2]
            if abs(U_val.imag) > 1e-10:
                raise RuntimeError("Cubic harmonics are real, so should be the matrix elements of U.")

            H_term = 0.5*U_val.real*c_dag(*mkind(s1,cubic_names[ap1]))*c_dag(*mkind(s2,cubic_names[ap2]))*c(*mkind(s2,cubic_names[a1]))*c(*mkind(s1,cubic_names[a2]))
            H += H_term

            # Dump quartic terms of H
            if H_dump and not H_term.is_zero():
                H_dump_file.write(mkind(s1,cubic_names[ap1])[0] + '\t')
                H_dump_file.write(mkind(s2,cubic_names[ap2])[0] + '\t')
                H_dump_file.write(mkind(s2,cubic_names[a1])[0] + '\t')
                H_dump_file.write(mkind(s1,cubic_names[a2])[0] + '\t')
                H_dump_file.write(str(U_val.real) + '\n')

# Quantum numbers (N_up and N_down)
QN=[Operator(),Operator()]
for cn in cubic_names:
    for i, sn in enumerate(spin_names):
        QN[i] += n(*mkind(sn,cn))
        
print_master("Constructing the solver...")

# Construct the solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

print_master("Preparing the hybridization function...")

# Set hybridization function
if Delta_dump: Delta_dump_file = open(Delta_dump,'w')
for sn, cn in product(spin_names,cubic_names):
    bn, i = mkind(sn,cn)
    V = delta_params[cn]['V']
    e = delta_params[cn]['e']
    
    delta_w = GfImFreq(indices = [i], beta=beta)
    delta_w <<= (V**2) * inverse(iOmega_n - e)

    S.G0_iw[bn][i,i] <<= inverse(iOmega_n +mu - atomic_levels[(bn,i)] - delta_w)

    # Dump Delta parameters
    if Delta_dump:
        Delta_dump_file.write(bn + '\t')
        Delta_dump_file.write(str(V) + '\t')
        Delta_dump_file.write(str(e) + '\n')

print_master("Running the simulation...")

# Solve the problem
S.solve(h_loc=H, params=pp, quantum_numbers=QN, use_quantum_numbers=use_quantum_numbers)

# Save the results  
if mpi.rank==0:
    Results = HDFArchive(results_file_name,'w')
    for b in gf_struct: Results[b] = S.G_tau[b]
