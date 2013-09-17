#!/bin/env pytriqs
#$ -N five_plus_five.matrix
#$ -pe mpi 128
#$ -q amd
#$ -cwd

from multiorbital import *

from itertools import product
import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.operators import *
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb_matrix import *
from pytriqs.gf.local import *

from params import *

# Remove Krylov-specific parameters 
p = {k:v for k, v in p.items() if not k.startswith('krylov')}

print "Welcome to 5+5 (5 orbitals + 5 bath sites) test, matrix version."

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

if use_interaction:
    print "Preparing the interaction matrix..."
    
    U_matrix = transform_U_matrix(U_matrix_spherical([F0,F2,F4]),W)

    # Interaction terms of the Hamiltonian
    a_range=range(N_comp)
    for s1, s2 in product(spin_names,spin_names):
        for ap1, ap2, a1, a2 in product(a_range,a_range,a_range,a_range):
            U_val = U_matrix[ap1,ap2,a1,a2]
            if abs(U_val.imag) > 1e-10:
                raise RuntimeError("Cubic harmonics are real, so should be the matrix elements of U.")

            H += 0.5*U_val.real*Cdag(*mkind(s1,cubic_names[ap1]))*Cdag(*mkind(s2,cubic_names[ap2]))*C(*mkind(s2,cubic_names[a1]))*C(*mkind(s1,cubic_names[a2]))

# Quantum numbers (N_up and N_dn)
QN={'N_up':Operator(),'N_dn':Operator()}
for cn in cubic_names:
    for sn in spin_names:
        QN['N_'+sn] += N(*mkind(sn,cn))

print "Constructing the solver..."

# Construct the solver
S = Solver(beta=beta, gf_struct=gf_struct.items())

print "Preparing the hybridization function..."

# Set hybridization function
for sn, cn in product(spin_names,cubic_names):
    bn, i = mkind(sn,cn)
    V = delta_params[cn]['V']
    e = delta_params[cn]['e']
    
    delta_w = GfImFreq(indices = [i], beta=beta)
    delta_w <<= (V**2) * inverse(iOmega_n - e)
    S.G0[bn][i,i] <<= inverse(iOmega_n - delta_w)

print "Running the simulation..."

# Solve the problem
p['H_local'] = H
p['quantum_numbers'] = QN

S.solve(**p)

# Save the results  
if mpi.rank==0:
    Results = HDFArchive(matrix_results_file_name,'w')
    for b in gf_struct: Results[b] = S.G_tau[b]