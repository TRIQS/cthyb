#!/bin/env pytriqs

from multiorbital import *

from itertools import product
import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb_krylov import *
from pytriqs.gf.local import *

from params import *

def print_master(msg):
	if mpi.rank==0: print msg

pp = Parameters()
for k in p: pp[k] = p[k]

print_master("Welcome to 5+5 (5 orbitals + 5 bath sites) test.")

cubic_names, W = spherical2cubic(L)
N_comp = len(cubic_names)

gf_struct = {}
for sn, cn in product(spin_names,cubic_names):
    bn, i = mkind(sn,cn)
    gf_struct[bn] = [i]

# Local Hamiltonian
H = Operator()

# Chemical potential
for b in gf_struct:
    for c in gf_struct[b]:
        H += -mu*N(b,c)

# Atomic levels
for i in atomic_levels:
    H += atomic_levels[i] * N(*i)

# Dump quadratic terms of H
if H_dump:
    H_dump_file = open(H_dump,'w')
    for b in gf_struct:
        H_dump_file.write(b + '\t')
        H_dump_file.write(b + '\t')
        H_dump_file.write(str(atomic_levels[(b,0)]-mu) + '\n')

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

            H_term = 0.5*U_val.real*C_dag(*mkind(s1,cubic_names[ap1]))*C_dag(*mkind(s2,cubic_names[ap2]))*C(*mkind(s2,cubic_names[a1]))*C(*mkind(s1,cubic_names[a2]))
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
    for n, sn in enumerate(spin_names):
        QN[n] += N(*mkind(sn,cn))

# Use PS quantum numbers (see arXiv:1209.0915)
if use_PS_quantum_numbers:
    for cn in cubic_names:
        QN += [Operator()]
        dN = N(*mkind(spin_names[0],cn)) - N(*mkind(spin_names[1],cn))
        QN[-1] = dN*dN
        
print_master("Constructing the solver...")

# Construct the solver
S = Solver(parameters=pp, H_local=H, quantum_numbers=QN, gf_struct=gf_struct)

print_master("Preparing the hybridization function...")

# Set hybridization function
if Delta_dump: Delta_dump_file = open(Delta_dump,'w')
for sn, cn in product(spin_names,cubic_names):
    bn, i = mkind(sn,cn)
    V = delta_params[cn]['V']
    e = delta_params[cn]['e']
    
    delta_w = GfImFreq(indices = [i], beta=beta)
    delta_w <<= (V**2) * inverse(iOmega_n - e)
    S.Delta_tau[bn][i,i] <<= InverseFourier(delta_w)

    # Dump Delta parameters
    if Delta_dump:
        Delta_dump_file.write(bn + '\t')
        Delta_dump_file.write(str(V) + '\t')
        Delta_dump_file.write(str(e) + '\n')

print_master("Running the simulation...")

# Solve the problem
S.solve(parameters=pp)

# Save the results  
if mpi.rank==0:
    Results = HDFArchive(results_file_name,'w')
    for b in gf_struct: Results[b] = S.G_tau[b]
